# Render the measured large examples without fitting models again.
website_plot_examples <- function(cases=c("million", "single", "multi", "selection")) {
  pkgload::load_all("L0glm", recompile=FALSE, quiet=TRUE)
  cache <- "codex_sandbox/website_large_examples"
  output <- "benchmarks/plots/website"
  draw <- function(name, width, height, code) {
    grDevices::png(file.path(output, paste0(name, ".png")), width=width,
      height=height, res=160, type="cairo")
    on.exit(grDevices::dev.off())
    graphics::par(mar=c(4,4.5,3,1),mgp=c(2.5,.7,0),family="sans")
    force(code)
    message("Rendered ", name)
  }
  labels <- c(L0glm="L0glm",lasso="LASSO",relaxed_lasso="Relaxed LASSO")
  for (case in cases) {
    if (case == "million") {
      s <- readRDS(file.path(cache,"million-summary.rds"))
      z <- structure(list(fits=list(L0glm=structure(list(coefficients=s$coefficients),
        class="L0glm_fit")),beta_true=s$beta_true,
        timings=data.frame(seconds=s$metrics$seconds),data=list(intercept=FALSE)),class="sparseglm")
      draw("million-coefficients",1400,1000,plot(z))
      path <- readRDS(file.path(cache,"million-path.rds"))
      draw("million-path",1700,1000,L0glm:::L0glm_plot_tuning(path,s$beta_true))
    } else if (case %in% c("single","multi")) {
      z <- readRDS(file.path(cache,paste0(case,"-comparison.rds")))
      scores <- read.csv(file.path(output,paste0(case,"-metrics.csv")))
      best_other <- which.max(replace(scores$MCC,scores$method=="L0glm",-Inf))
      panel_labels <- labels
      if (case == "multi") panel_labels[c("lasso","relaxed_lasso")] <- c("Group LASSO","Relaxed group LASSO")
      names(z$fits) <- unname(panel_labels[names(z$fits)])
      if (case == "single") {
        draw("deconvolution-single-large",1700,1300,{
          graphics::par(mfrow=c(2,1),mar=c(4.5,4.5,3,1))
          observed <- as.numeric(z$data$y)
          truth <- as.numeric(z$beta_true)
          true_mean <- as.numeric(z$data$X %*% z$beta_true)
          for (i in c(1L,best_other)) {
            estimate <- as.numeric(coef(z$fits[[i]]))
            mean <- as.numeric(z$data$X %*% estimate)
            limit <- max(observed,truth,estimate,mean,true_mean)
            graphics::plot(seq_along(observed),observed,type="l",col="grey65",
              ylim=c(-limit,limit),xlab="Time index",ylab="Signal and spikes",
              main=sprintf("%s (%.3g s)",names(z$fits)[i],z$timings$seconds[i]))
            graphics::lines(true_mean,col="#173e65",lwd=.8)
            graphics::lines(truth,type="h",col="navy",lwd=.7)
            graphics::lines(-observed,col="grey85",lwd=.6)
            graphics::lines(-estimate,type="h",col="red3",lwd=.7)
            graphics::lines(-mean,col="red",lwd=1)
            graphics::abline(h=0,col="grey70",lty=3)
            graphics::legend("topleft",c("Observed signal","True mean / spikes","Fitted mean / spikes (mirrored)"),
              col=c("grey65","navy","red"),lty=1,bty="n",cex=.7,horiz=TRUE)
          }
        })
        draw("deconvolution-single-coefficients-large",1700,850,
          plot(z,which=c(1L,best_other)))
      } else {
        draw("deconvolution-multi-large",1700,1800,
          plot(z,type="heatmap",which=c(1L,best_other)))
        draw("deconvolution-multi-coefficients-large",1100,1800,
          plot(z,which=1:3,channels="all"))
        draw("deconvolution-multi-home-coefficients",1200,1300,
          plot(z,which=c(1L,best_other),channels="all"))
      }
    } else if (case == "selection") {
      s <- readRDS(file.path(cache,"selection-paths.rds"))
      fields <- c("gic","TP","FP","FN","TN","F1_score","MCC")
      preset <- read.csv(file.path(output,"selection-metrics.csv"))[fields]
      fits <- list(s$gic$fit_best,s$cv$fit_best,s$cv$fit_1se)
      measures <- do.call(rbind,c(list(preset),lapply(fits,function(f)
        L0glm.metrics(f,beta_true=s$beta_true)[fields])))
      measures$method <- c("GIC preset","Tuned GIC","CV minimum","CV one-SE")
      write.csv(measures,file.path(output,"selection-comparison.csv"),row.names=FALSE)
      draw("selection-gic",1700,1000,L0glm:::L0glm_plot_tuning(s$gic,s$beta_true))
      draw("selection-cv",1700,1000,L0glm:::L0glm_plot_tuning(s$cv,s$beta_true))
    }
  }
  invisible(NULL)
}
if (sys.nframe()==0L) website_plot_examples()
