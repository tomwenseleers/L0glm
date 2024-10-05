# EXTENSION: TEST L0glm FIT WITH OPTION FOR RANDOM FACTORS ADDED ####

devtools::load_all(".")
library(scales)

# IDEA: I think we can incorporate something akin to random effects
# by allowing an option random.factor.groups which would be a factor
# with the IDs of each group of one-hot encoded variables corresponding to each
# individual factor level (see https://stackoverflow.com/questions/62116652/one-hot-encoding-using-model-matrix)
# and NAs for variables not included in any random factor
# and in the IRLS loop then updating the prior.mean vector to the means on the coefficients
# estimated for all the variables within each random.factor.groups
# in each iteration
# so that it will shrink the coefficients to the center of each re-estimated group
# average coefficient in each iteration
# (something like this is done in mgcv gam fits with s(bs="re") terms)

# If we would allow random.factor.groups to be a list with nested hierarchical
# levels, e.g. individuals nested within colonies and colonies nested within populations,
# then this averaging & fitting could also be done sequentially, e.g. 
# in each IRLS iteration first do lmridge fit with random.factor.groups given
# individual IDs & shrink towards calculated per-individual mean,
# then do another lmridge fit with random.factor.groups given colony IDs & shrink
# towards calculated per-colony means given by average coefs for individuals within
# each group & finally do a last lmridge fit with random.factor.groups given
# population IDs & shrink towards calculated per-population means given by average
# coefs for different colonies within each population
# This would then require 3 lmridge fits within each IRLS iteration, which would
# still be OK.
# If we would like these nested random effects to be independent we could just
# calculate the averages independently for each nested level, i.e. not use
# averages of averages.

# An alternative would be to use fused L2/ridge type penalties for this & put in
# all penalties together in one single lmridge fit (even allowing for e.g.
# fused L1 penalties for the random effects, which would then correspond to
# Laplace random effects instead of Gaussian ones), but this would then
# introduce additional lambda penalties, which I would rather like to avoid.

# EXTRA OPTIONS NEEDED:
# random.factor.groups = named list with random factors with random factor group IDs
# random.factor.type = "nested" or "crossed"

# Predict method would probably need to allow to predict either for 
# original factor levels or allow original factor levels to be randomly
# resampled, so that predictions would be for the original individuals/colonies/populations
# or for an average individual/colony/population
# (the latter would probably make most sense)


