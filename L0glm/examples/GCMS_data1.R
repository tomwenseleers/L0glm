data(GCMS_data1)
plotmat(GCMS_data1, 
      xlab = "Elution time (min)",
      ylab = "Mass (m/z)",
      main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
matplot(x = as.numeric(rownames(GCMS_data1)), 
        y = GCMS_data1^0.2, type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
        xlab = "Elution time (min)",
        ylab = "Intensity ^ 0.2",
        main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
plot(x=as.numeric(names(attr(GCMS_data1, "alkanes_TIC"))), 
     y=attr(GCMS_data1, "alkanes_TIC"), 
     xlab="Retention time (min)", 
     ylab="Total ion intensity", 
     main="C9 to C40 alkane ladder standard", 
     type="l")
