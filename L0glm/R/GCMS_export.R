GCMS_xlsx_path <- function(path) {
  if(!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)||
      !grepl("[.]xlsx$",path,ignore.case=TRUE)) stop("xlsx_file must end in .xlsx.",call.=FALSE)
  if(!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx for Excel export.",call.=FALSE)
  invisible(path)
}

GCMS_write_xlsx <- function(model,time,mass,compounds,path,cdf_file,alkanes_file,seconds) {
  GCMS_xlsx_path(path)
  wb<-openxlsx::createWorkbook(creator="L0glm")
  header<-openxlsx::createStyle(fgFill="#102d4d",fontColour="#FFFFFF",textDecoration="bold",
    fontSize=11,wrapText=TRUE,valign="center")
  S<-model$S[compounds$compound,,drop=FALSE]
  if(nrow(S)) S<-S/rowSums(S)
  spectra<-cbind(compounds[,c("sample","compound","name","retention_time"),drop=FALSE],
    as.data.frame(S,optional=TRUE))
  names(spectra)[4L+seq_along(mass)]<-paste0("mz_",mass)
  metadata<-data.frame(field=c("sample","source CDF","calibration CDF","noise family",
    "dispersion","theta","scans","mass channels","compounds","retention time unit",
    "area unit","area_min unit","spectral normalization","fitting seconds"),
    value=c(if(nrow(compounds)) compounds$sample[1L] else tools::file_path_sans_ext(basename(cdf_file)),
      cdf_file,if(is.null(alkanes_file)) "" else alkanes_file,model$spec$name,
      format(model$spec$dispersion),format(model$spec$theta),length(time),length(mass),nrow(compounds),
      "minutes","intensity-scans","intensity-minutes","unit L1 sum per compound",
      if("total"%in%names(seconds)) seconds[["total"]] else sum(seconds)),
    stringsAsFactors=FALSE)
  tables<-list(Compounds=compounds,Spectra=spectra,Run=metadata)
  for(name in names(tables)) {
    openxlsx::addWorksheet(wb,name,gridLines=FALSE)
    openxlsx::writeData(wb,name,tables[[name]],headerStyle=header,withFilter=name!="Run")
    openxlsx::freezePane(wb,name,firstActiveRow=2L,firstActiveCol=if(name=="Spectra") 5L else 1L)
    openxlsx::setRowHeights(wb,name,1L,32)
    openxlsx::setColWidths(wb,name,seq_len(ncol(tables[[name]])),if(name=="Run") c(28,90) else 16)
  }
  openxlsx::setColWidths(wb,"Compounds",3L,38)
  openxlsx::setColWidths(wb,"Spectra",3L,38)
  if(nrow(compounds)) {
    rows<-seq_len(nrow(compounds))+1L
    openxlsx::addStyle(wb,"Compounds",openxlsx::createStyle(numFmt="0.00000"),rows=rows,cols=4L,gridExpand=TRUE)
    openxlsx::addStyle(wb,"Compounds",openxlsx::createStyle(numFmt="0.000E+00"),rows=rows,cols=6:7,gridExpand=TRUE)
    openxlsx::addStyle(wb,"Spectra",openxlsx::createStyle(numFmt="0.00000"),rows=rows,
      cols=4L+seq_along(mass),gridExpand=TRUE)
  }
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  openxlsx::saveWorkbook(wb,path,overwrite=TRUE)
  invisible(path)
}
