## Dense matrix example: plot a simulated multichannel signal
M <- simul(n=400L, p=400L, n_nonzero=30L, m=24L)$y
plotmat(M = M, 
        gamma = 0.13, 
        xlab = "Time index", 
        ylab = "Channel",
        main = "Simulated multichannel signal")
        
