# Clone the mutable NB shape cell, retaining immutable callbacks and their data.
# Serializing an entire family can recursively duplicate its caller's workspace.
L0glm_copy_family <- function(family) {
  if(!is.function(family$getTheta)||!is.function(family$putTheta))return(family)
  source <- environment(family$getTheta)
  if(!identical(environment(family$putTheta),source)||
     !identical(environment(family$variance),source))
    return(unserialize(serialize(family,NULL)))
  target <- new.env(parent=source)
  target$.Theta <- family$getTheta()
  answer <- family
  for(name in names(answer))if(is.function(answer[[name]])&&identical(environment(answer[[name]]),source))
    environment(answer[[name]]) <- target
  answer
}
