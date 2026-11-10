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
    result <- eval(call,environment())
    if(source)result <- L0glm_data_restore(result,X)
    if(ENGINE=='l0glmFit' && 'control.family'%in%supplied)force(control.family)
    recorded <- setdiff(supplied,c('X','y'))
    recorded <- recorded[!rlang::env_binding_are_lazy(environment(),recorded)]
    specification <- mget(recorded,envir=environment(),inherits=FALSE)
    if(ENGINE=='l0glmFit')result <- L0glm_calibrated_fit(result,specification)
    result$engine <- ENGINE
    result$specification <- specification
    result
  },list(KERNEL=substitute(kernel),ENGINE=engine))
  environment(wrapper) <- environment(kernel)
  wrapper
}
