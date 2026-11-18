# L0glm 0.2.0.9004

* Statistical fitting, tuning, screening and generic simulation remain in L0glm.
* Sparse matrix factorization, PSF estimation, blind deconvolution and the experimental chromatogram windows are provided by logconNMF.
* GCMSprime provides chromatogram files, deconvolution workflows, spectral libraries and reports.
* Shared native likelihood and conditional background interfaces preserve the fitting engines across packages.
* Coupled Poisson backgrounds recognize stationary starts before another
  Newton solve; ordinary conditional initialization retains its existing updates.
  Shared scalar kernels provide density-free derivatives without changing
  absolute likelihoods or information criteria.

