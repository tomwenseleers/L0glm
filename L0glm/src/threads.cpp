#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#include <R_ext/BLAS.h>
#endif

#include <Rcpp.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

#ifdef _WIN32
using Library = HMODULE;
void* lookup(Library library, const char* name) {
  return reinterpret_cast<void*>(GetProcAddress(library, name));
}
#else
using Library = void*;
void* lookup(Library library, const char* name) {
  return dlsym(library, name);
}
#endif

// Read and write a count only when the loaded BLAS provides both operations.
template<class Int>
struct ThreadApi {
  using Get = Int (*)();
  using Set = void (*)(Int);
  Get get = nullptr;
  Set set = nullptr;
  std::array<Get, 5> blis_ways{};

  bool available() const { return get && set; }

  int count() const {
    if (!available()) return NA_INTEGER;
    // A single count cannot restore a manual BLIS loop partition.
    for (const auto way : blis_ways) {
      if (way && way() > 1) return NA_INTEGER;
    }
    const Int value = get();
    if (value < 1 || value > std::numeric_limits<int>::max()) return NA_INTEGER;
    return static_cast<int>(value);
  }
};

template<class Int>
ThreadApi<Int> thread_api(Library library, const char* get, const char* set) {
  return {reinterpret_cast<typename ThreadApi<Int>::Get>(lookup(library, get)),
          reinterpret_cast<typename ThreadApi<Int>::Set>(lookup(library, set))};
}

template<class Int>
ThreadApi<Int> blis_api(Library library) {
  auto api = thread_api<Int>(library, "bli_thread_get_num_threads", "bli_thread_set_num_threads");
  const char* names[] = {"bli_thread_get_jc_nt", "bli_thread_get_pc_nt", "bli_thread_get_ic_nt",
                         "bli_thread_get_jr_nt", "bli_thread_get_ir_nt"};
  for (std::size_t i = 0; i < api.blis_ways.size(); ++i) {
    api.blis_ways[i] = reinterpret_cast<typename ThreadApi<Int>::Get>(lookup(library, names[i]));
    if (!api.blis_ways[i]) return {};
  }
  return api;
}

struct BlasThreads {
  ThreadApi<int> standard;
  ThreadApi<std::int64_t> wide;

  bool available() const { return standard.available() || wide.available(); }
  int count() const { return wide.available() ? wide.count() : standard.count(); }

  void change(int previous, int requested) const {
    if (previous == NA_INTEGER || requested <= 0 || previous == requested) return;
    if (wide.available()) wide.set(requested);
    else standard.set(requested);
  }
};

BlasThreads find_thread_api(Library library) {
  const char* symbols[][2] = {
    {"flexiblas_get_num_threads", "flexiblas_set_num_threads"},
    {"openblas_get_num_threads", "openblas_set_num_threads"},
    {"MKL_Get_Max_Threads", "MKL_Set_Num_Threads"},
    {"mkl_get_max_threads", "mkl_set_num_threads"}
  };
  for (const auto& names : symbols) {
    auto api = thread_api<int>(library, names[0], names[1]);
    if (api.available()) return {api, {}};
  }
  // BLIS can be configured with either 32- or 64-bit dimension arguments.
  const auto width = reinterpret_cast<const char* (*)()>(
    lookup(library, "bli_info_get_int_type_size_str"));
  if (width) {
    const char* bits = width();
    if (std::strcmp(bits, "64") == 0) return {{}, blis_api<std::int64_t>(library)};
    if (std::strcmp(bits, "32") == 0) return {blis_api<int>(library), {}};
  }
  return {};
}

BlasThreads r_blas_threads() {
  // Follow R's BLAS, which may differ from another library in this process.
#ifdef _WIN32
  const Library library = GetModuleHandleW(L"Rblas.dll");
  if (!library) return {};
  const auto api = find_thread_api(library);
  if (api.available()) return api;
  const auto dgemm = GetProcAddress(library, "dgemm_");
  Library owner = nullptr;
  if (dgemm && GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                                 GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                                 reinterpret_cast<LPCSTR>(dgemm), &owner) && owner != library) {
    return find_thread_api(owner);
  }
#elif defined(RTLD_NOLOAD)
  Dl_info info{};
  if (dladdr(reinterpret_cast<void*>(&F77_CALL(dgemm)), &info) && info.dli_fname) {
    const Library library = dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD);
    if (library) {
      const auto api = find_thread_api(library);
      // The module remains loaded by R after releasing our temporary reference.
      dlclose(library);
      return api;
    }
  }
#endif
  return {};
}

} // namespace

// [[Rcpp::export]]
Rcpp::List l0glm_threads_cpp(int blas, int openmp) {
  static const auto api = r_blas_threads();
  const int previous_blas = api.count();
  int previous_omp = NA_INTEGER;
#ifdef _OPENMP
  previous_omp = omp_get_max_threads();
#endif
  auto previous = Rcpp::List::create(Rcpp::Named("blas") = previous_blas,
                                     Rcpp::Named("openmp") = previous_omp);
  api.change(previous_blas, blas);
#ifdef _OPENMP
  if (openmp > 0 && openmp != omp_get_max_threads()) omp_set_num_threads(openmp);
#endif
  return previous;
}
