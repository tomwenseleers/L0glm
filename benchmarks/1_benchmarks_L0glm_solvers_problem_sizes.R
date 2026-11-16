# One fit per solver, problem and constraint strategy; no untimed fitting warmup.
# Run from the repository root, or from benchmarks/. Data generation, process
# startup, namespace loading and metrics are outside the fitting timer.
# Dense/sparse deconvolution labels describe input storage. Native solvers may
# compress exact zeros. The Gaussian normal-design example is genuinely dense.

suppressPackageStartupMessages(library(L0glm))
suppressPackageStartupMessages(library(ggplot2))
if(!"nonnegative.method"%in%names(formals(L0glm::controlSolverGen)))
  stop("Install a current L0glm version supporting controlSolverGen(nonnegative.method=...).")

repository_dir <- if(file.exists("L0glm/DESCRIPTION"))"."else".."
helper <- normalizePath(file.path(repository_dir,"benchmarks/solver_problem_sizes_helpers.R"),winslash="/")
source(helper)
output_dir <- file.path(repository_dir,"benchmarks/1_benchmarks_solvers_problem_sizes/nonnegative_methods")

# Main settings: both strategies for nonnegative Gaussian/Poisson deconvolution.
nonnegative_methods <- c("clamp","constrained")
deconvolution_sizes <- c(100L,500L,1000L,10000L)
# Two dense-normal sweeps, sharing the n=p=1000 condition.
dense_sizes <- c(100L,1000L,10000L)
dense_fixed_n <- dense_fixed_p <- 1000L
config <- list(criterion="gic",maxit=500L,tol=1e-10,ncores=1L,blas_threads=1L,
  max_memory=NULL,timeout=900,deconvolution_seed=1L,dense_seed=123L,
  dense_noise_seed=1L,dense_effects=60L,dense_screen="none")   # or "ISIS"

# Use .libPaths() before sourcing to select a private package installation.
# A fresh R worker per fit isolates caches and permits a hard per-fit timeout.
arguments <- commandArgs(TRUE)
plots_only <- "--plots-only" %in% arguments
resume <- "--resume" %in% arguments
smoke <- "--smoke" %in% arguments
output_argument <- grep("^--output-dir=",arguments,value=TRUE)
if(length(output_argument)>1L)stop("Supply --output-dir only once.")
if(length(output_argument))output_dir <- sub("^--output-dir=","",output_argument)
if(!nzchar(output_dir))stop("The output directory must not be empty.")
if(smoke) {
  deconvolution_sizes <- 100L;dense_sizes <- c(100L,200L)
  dense_fixed_n <- dense_fixed_p <- 100L;config$dense_effects <- 10L
  output_dir <- file.path(output_dir,"smoke")
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
output_dir <- normalizePath(output_dir,winslash="/")
screening_label <- if(config$dense_screen=="none")"screen=none"else paste(config$dense_screen,"when p>n")
screening_caption <- if(config$dense_screen=="none")"Screening disabled for every problem."else
  sprintf("The complete %s workflow is timed when p>n.",config$dense_screen)
solver_grid <- data.frame(Solver=c("OSQP","bigGlm","Eigen LLT","Eigen LDLT","Eigen CG-IC","Eigen LSCG-IC"),
  Solver_backend=c("osqp","bigGlm",rep("eigen",4L)),Solver_method=c(1L,1L,1L,2L,3L,4L))
solver_colours <- setNames(c("blue","red2","green3","orange","purple","cyan3"),solver_grid$Solver)
solver_shapes <- setNames(c(1,17,0,15,18,4),solver_grid$Solver)
solver_linetypes <- setNames(c("dashed","solid","solid","dotted","longdash","dotdash"),solver_grid$Solver)
metric_columns <- c("Time","GIC","adj_R2_weighted","R2","RMSE_coefficients","Precision","Recall","F1","deviance")
metric_labels <- c(Time="Fitting time (s)",GIC="GIC",adj_R2_weighted="Adjusted weighted prediction R-squared",
  R2="Coefficient R-squared",RMSE_coefficients="Coefficient RMSE",Precision="Precision",Recall="Recall",
  F1="Support F1",deviance="True-signal deviance / n")

makeCases <- function() {
  d <- expand.grid(Noise_family=c("poisson","gaussian"),Sparse=c(TRUE,FALSE),N=deconvolution_sizes,
    Nonnegative_method=nonnegative_methods,stringsAsFactors=FALSE)
  d$Experiment <- "deconvolution";d$P <- d$N
  d$Matrix_density <- ifelse(d$Sparse,"sparse","dense")
  dense <- unique(rbind(data.frame(N=dense_fixed_n,P=dense_sizes),data.frame(N=dense_sizes,P=dense_fixed_p)))
  dense$Experiment <- "dense_normal";dense$Noise_family <- "gaussian";dense$Sparse <- FALSE
  dense$Nonnegative_method <- "unrestricted";dense$Matrix_density <- "dense"
  cases <- rbind(d,dense[,names(d)])
  cases$Case <- paste(cases$Experiment,cases$Noise_family,cases$Matrix_density,
    paste0("n",cases$N),paste0("p",cases$P),cases$Nonnegative_method,sep="-")
  cases
}

results_file <- file.path(output_dir,"metrics_all_conditions.csv")
dir.create(file.path(output_dir,"fits"),showWarnings=FALSE)
dir.create(file.path(output_dir,"logs"),showWarnings=FALSE)
if(!plots_only) {
  cases <- makeCases()
  signature <- list(cases=cases,solvers=solver_grid,config=config,
    L0glm=as.character(packageVersion("L0glm")),R=R.version.string,
    DLL=unname(tools::md5sum(list.files(file.path(find.package("L0glm"),"libs"),
      pattern="[.](dll|so)$",full.names=TRUE,recursive=TRUE))),
    Rcode=unname(tools::md5sum(list.files(file.path(find.package("L0glm"),"R"),
      pattern="[.](rdb|rdx)$",full.names=TRUE))),
    helper=unname(tools::md5sum(helper)))
  manifest <- file.path(output_dir,"run_manifest.rds")
  if(resume)stopifnot(file.exists(manifest),identical(signature,readRDS(manifest)))
  saveRDS(signature,manifest)
  jobs <- expand.grid(Case_index=seq_len(nrow(cases)),Solver_index=seq_len(nrow(solver_grid)))
  jobs <- jobs[order(cases$Experiment[jobs$Case_index],cases$Noise_family[jobs$Case_index],
    cases$Sparse[jobs$Case_index],cases$N[jobs$Case_index],cases$P[jobs$Case_index],
    jobs$Solver_index,match(cases$Nonnegative_method[jobs$Case_index],nonnegative_methods)),]
  rows <- list()
  for(k in seq_len(nrow(jobs))) {
    i <- jobs$Case_index[k];j <- jobs$Solver_index[k]
    case <- cases[i,,drop=FALSE];solver <- solver_grid[j,,drop=FALSE]
    key <- paste(case$Case,solver$Solver_backend,solver$Solver_method,sep="-")
    file <- file.path(output_dir,"fits",paste0(key,".rds"))
    message(sprintf("[%d/%d] %s; %s",k,nrow(jobs),case$Case,solver$Solver))
    if(resume&&file.exists(file))out <- readRDS(file)else {
      out <- tryCatch(callr::r(solverBenchmarkWorker,args=list(helper=helper,case=case,solver=solver,config=config),
        libpath=.libPaths(),timeout=config$timeout,
        stdout=file.path(output_dir,"logs",paste0(key,".out")),
        stderr=file.path(output_dir,"logs",paste0(key,".err"))),error=function(e) {
          row <- cbind(case,solver,data.frame(Nonnegative=case$Experiment=="deconvolution",
            Start=if(case$Experiment=="deconvolution")"1"else"NULL",
            Screen=if(case$Experiment=="dense_normal"&&case$P>case$N)config$dense_screen else "none",Local_search=FALSE,
            Ncores=config$ncores,BLAS_threads=config$blas_threads,Timing_repetitions=1L,Warmup_fits=0L,
            Fit_calls=NA_integer_,Timing_clock="microbenchmark",Maxit=config$maxit,Tolerance=config$tol,
            Iterations=NA_integer_,Converged=NA,Dispersion=NA_real_,
            Status=if(inherits(e,"callr_timeout_error"))"timeout"else"worker_error",Warnings=NA_character_,Error=conditionMessage(e)))
          for(name in c(metric_columns,"MCC","TP","FP","FN","TN"))row[[name]] <- NA_real_
          list(row=row,coefficients=NULL,truth=NULL,convergence=NULL)
        })
      saveRDS(out,file)
    }
    rows[[length(rows)+1L]] <- out$row
    results <- do.call(rbind,rows)
    write.csv(results,results_file,row.names=FALSE)
    saveRDS(results,file.path(output_dir,"metrics_all_conditions.rds"))
    message(sprintf("  %s; %.4f s; GIC %.3f; coefficient R2 %.4f; F1 %.4f; converged %s",
      out$row$Status,out$row$Time,out$row$GIC,out$row$R2,out$row$F1,out$row$Converged))
  }
}else results <- read.csv(results_file,stringsAsFactors=FALSE)

results$Solver <- factor(results$Solver,levels=solver_grid$Solver)
results$Constraint_handling <- ifelse(!results$Nonnegative,"unrestricted",results$Nonnegative_method)
write.csv(results,results_file,row.names=FALSE)
saveRDS(results,file.path(output_dir,"metrics_all_conditions.rds"))
results_long <- tidyr::pivot_longer(results,cols=tidyselect::all_of(metric_columns),names_to="Metric",values_to="Value")
results_long$Metric <- factor(results_long$Metric,levels=metric_columns)
write.csv(results_long,file.path(output_dir,"metrics_all_conditions_long.csv"),row.names=FALSE)

reportPlots <- list()
savePlot <- function(plot,name,width=11,height=8) {
  reportPlots[[name]] <<- plot
  writePlot <- function(file)withCallingHandlers(
    ggplot2::ggsave(file,plot,width=width,height=height,dpi=180),warning=function(w) {
      if(grepl("could not write|cannot open file",conditionMessage(w),ignore.case=TRUE))
        stop(conditionMessage(w),call.=FALSE)
    })
  for(extension in c("png","pdf")) {
    file <- file.path(output_dir,paste0(name,".",extension))
    tryCatch(writePlot(file),error=function(e) {
      alternate <- file.path(output_dir,paste0(name,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),".",extension))
      warning("Could not replace ",basename(file),"; saving ",basename(alternate),call.=FALSE)
      writePlot(alternate)
    })
  }
}
solverLayers <- function(data) {
  has_path <- length(unique(data$N))>1L||length(unique(data$P))>1L
  unlist(lapply(list(!data$Solver_backend%in%c("osqp","bigGlm"),data$Solver_backend=="bigGlm",data$Solver_backend=="osqp"),
    function(use)list(if(has_path)geom_line(data=data[use,],linewidth=.8,na.rm=TRUE),
      geom_point(data=data[use,],size=2,stroke=.65,na.rm=TRUE))),recursive=FALSE)
}
solverScales <- function()list(scale_colour_manual(values=solver_colours,drop=FALSE),
  scale_shape_manual(values=solver_shapes,drop=FALSE),scale_linetype_manual(values=solver_linetypes,drop=FALSE))
basePlot <- function(data,x,label) {
  ggplot(data,aes(x=.data[[x]],y=Value,colour=Solver,shape=Solver,linetype=Solver,group=Solver))+
    solverLayers(data)+solverScales()+theme_minimal(base_size=11)+
    theme(legend.position="bottom")+guides(colour=guide_legend(nrow=2),shape=guide_legend(nrow=2),linetype=guide_legend(nrow=2))+
    scale_x_log10(breaks=sort(unique(data[[x]])),labels=scales::label_number(big.mark=","))+
    labs(x=label,y=NULL,colour="Solver",shape="Solver",linetype="Solver")
}
caption <- paste("One timed fit; no fitting warmup. Native dispatch; local search disabled.",
  "OSQP dashed overlays bigGlm solid. Constrained LSCG uses checked bounded polishing.")
caption <- paste(strwrap(caption,width=130),collapse="\n")
for(mode in nonnegative_methods) {
  data <- subset(results_long,Experiment=="deconvolution"&Nonnegative_method==mode)
  data$Noise_family <- factor(data$Noise_family,levels=c("poisson","gaussian"))
  for(noise in c("poisson","gaussian"))for(storage in c("sparse","dense")) {
    condition <- subset(data,Noise_family==noise&Matrix_density==storage)
    plot <- basePlot(condition,"N","Problem size (n = p)")+
      facet_wrap(~Metric,scales="free_y",ncol=3,labeller=as_labeller(metric_labels))+
      labs(title=paste(tools::toTitleCase(noise),"deconvolution;",storage,"input;",mode),
        subtitle=sprintf("GIC preset; start=1; screen=none; maxit=%d; ncores=%d; BLAS threads=%d",config$maxit,config$ncores,config$blas_threads),caption=caption)
    savePlot(plot,paste("timings_metrics",storage,noise,mode,sep="_"))
  }
  timing <- subset(data,Metric=="Time"&Value>0)
  plot <- basePlot(timing,"N","Problem size (n = p)")+facet_grid(Matrix_density~Noise_family)+scale_y_log10()+
    labs(title=paste("Deconvolution solver timings:",mode),
      subtitle=sprintf("GIC preset; start=1; screen=none; maxit=%d; one core and one BLAS thread",config$maxit),
      y="Fitting time (s)",caption=caption)
  savePlot(plot,paste0("timings_all_conditions_",mode),width=11,height=7)
}

for(sweep in c("vary_p","vary_n")) {
  data <- subset(results_long,Experiment=="dense_normal")
  if(sweep=="vary_p") {data <- subset(data,N==dense_fixed_n);axis <- "P";fixed <- paste("n =",dense_fixed_n)} else {
    data <- subset(data,P==dense_fixed_p);axis <- "N";fixed <- paste("p =",dense_fixed_p)}
  plot <- basePlot(data,axis,if(axis=="P")"Predictors (p)"else"Observations (n)")+
    facet_wrap(~Metric,scales="free_y",ncol=3,labeller=as_labeller(metric_labels))+
    labs(title=paste("Dense Gaussian regression:",fixed),
      subtitle=sprintf("%d signed effects; unit Gaussian noise; GIC preset; start=NULL; %s; maxit=%d",config$dense_effects,screening_label,config$maxit),
      caption=paste(strwrap(paste("One fresh-worker fit; no fitting warmup; one core and one BLAS thread.",screening_caption,
        "Local search disabled. OSQP dashed overlays bigGlm solid; coefficients are unrestricted."),width=130),collapse="\n"))
  savePlot(plot,paste0("timings_metrics_dense_normal_",sweep))
}

comparison <- subset(results,Experiment=="deconvolution")
clamp <- subset(comparison,Nonnegative_method=="clamp")
constrained <- subset(comparison,Nonnegative_method=="constrained")
paired <- merge(clamp,constrained,by=c("Experiment","Noise_family","Matrix_density","Sparse","N","P","Solver","Solver_backend","Solver_method"),
  suffixes=c("_clamp","_constrained"))
paired$Time_ratio_constrained_clamp <- paired$Time_constrained/paired$Time_clamp
for(metric in c("GIC","R2","RMSE_coefficients","F1","TP","FP"))
  paired[[paste0("Change_",metric)]] <- paired[[paste0(metric,"_constrained")]]-paired[[paste0(metric,"_clamp")]]
write.csv(paired,file.path(output_dir,"clamp_vs_constrained.csv"),row.names=FALSE)
write.csv(subset(results,Experiment=="dense_normal"),file.path(output_dir,"dense_normal_results.csv"),row.names=FALSE)
writeLines(c(paste("L0glm:",packageVersion("L0glm")),paste("R:",R.version.string),
  paste("Fits:",nrow(results)),paste("Successful:",sum(results$Status=="ok")),
  paste("Converged:",sum(results$Converged,na.rm=TRUE)),paste("Completed:",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")),
  caption,"Deconvolution: Gaussian/Poisson, sparse/dense input, unit starts, both nonnegative strategies.",
  screening_caption,
  sprintf("Dense normal design: %d signed effects drawn uniformly on [-0.8,0.8]; errors N(0,1); NULL starts.",config$dense_effects),
  "Gaussian deconvolution retains the original simulation's zero clipping of negative observations.",
  "Each input is regenerated from fixed seeds in an isolated worker; only l0glmFit is timed.",
  "Failed or timed-out fits remain explicit NA rows; maxit and convergence are recorded.",
  "GIC uses fitted dispersion and coefficient-count degrees of freedom.",
  "The dense signed example has no bounds, so clamp/constrained modes are not duplicated."),
  file.path(output_dir,"run_summary.txt"))
grDevices::pdf(file.path(output_dir,"solver_size_comparison.pdf"),width=11,height=8)
tryCatch(for(plot in reportPlots)print(plot),finally=grDevices::dev.off())
message("Saved benchmark outputs in ",output_dir)
if(any(results$Status!="ok"))warning("Some fits failed or timed out; inspect Status and Error.",call.=FALSE)
if(any(!results$Converged,na.rm=TRUE))warning("Some fits reached the iteration ceiling; inspect Converged and Iterations.",call.=FALSE)
