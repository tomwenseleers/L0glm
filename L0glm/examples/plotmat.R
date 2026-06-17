## Dense matrix example: plot a GC/MS chromatogram
M = simulate_GCMS()$D # simulated GC/MS chromatogram
plotmat(M = M, 
        gamma = 0.13, 
        xlab = "Elution time (mins)", 
        ylab = "Mass (m/z)",
        main = "GC/MS chromatogram")
        
