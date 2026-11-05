#ifndef L0GLM_EIGEN_DIAGNOSTICS_H
#define L0GLM_EIGEN_DIAGNOSTICS_H

// GCC can report uninitialized-buffer false positives in Eigen's products.
// Suppress them in dependency headers and restore warnings for package code.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <RcppEigen.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#endif
