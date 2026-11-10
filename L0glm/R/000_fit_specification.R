# Preserve argument laziness while retaining effective options for resampling.
L0glm_recorded_fitter <- function(kernel,engine) {
  wrapper <- function()NULL
  formals(wrapper) <- formals(kernel)
  body(wrapper) <- substitute({
    supplied <- names(match.call())[-1L]
    forwarded <- stats::setNames(lapply(supplied,as.name),supplied)
    call <- as.call(c(list(quote(KERNEL)),forwarded))
    source <- inherits(X,'L0glm_data')
    if(source)X <- L0glm_data_training(X,weights)
    if(ENGINE=='l0glmFit') {
      if(!inherits(family,'family'))family <- family()
      if(family$family %in% c("gaussian","Gamma","inverse.gaussian")) {
        family$consistent.scale <- TRUE
        call[['family']] <- quote(family)
      }
      if('control.family'%in%supplied && is.list(control.family) && !is.null(control.family$dispersion.calibration))
        L0glm_calibration_mode(control.family$dispersion.calibration)
    }
    result <- eval(call,environment())
    if(source)result <- L0glm_data_restore(result,X)
    if(ENGINE=='l0glmFit' && 'control.family'%in%supplied)force(control.family)
    recorded <- setdiff(supplied,c('X','y'))
    recorded <- recorded[!rlang::env_binding_are_lazy(environment(),recorded)]
    specification <- mget(recorded,envir=environment(),inherits=FALSE)
    result$engine <- ENGINE
    result$specification <- specification
    result
  },list(KERNEL=substitute(kernel),ENGINE=engine))
  environment(wrapper) <- environment(kernel)
  wrapper
}
