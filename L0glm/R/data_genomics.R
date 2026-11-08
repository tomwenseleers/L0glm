L0glm_data_pgen <- function(x,rows,columns) {
  selected <- sort(unique(rows));key <- digest::digest(selected,algo='xxhash64')
  if(is.null(x$cache$pgen))x$cache$pgen <- list()
  handle <- x$cache$pgen[[key]]
  valid <- if(is.null(handle))FALSE else tryCatch(pgenlibr::GetVariantCt(handle)>0,error=function(e)FALSE)
  if(!valid) {
    handle <- pgenlibr::NewPgen(x$object,sample_subset=selected)
    x$cache$pgen[[key]] <- handle
  }
  pgenlibr::ReadList(handle,columns,meanimpute=FALSE)[match(rows,selected),,drop=FALSE]
}

L0glm_data_seurat <- function(path,group) {
  if(!requireNamespace('rhdf5',quietly=TRUE))stop('Install rhdf5.',call.=FALSE)
  table <- rhdf5::h5ls(path,recursive=TRUE)
  full <- paste0(ifelse(table$group=='/','',table$group),'/',table$name)
  entry <- match(paste0('/',sub('^/','',group)),full)
  if(is.na(entry))stop('Assay dataset is absent from the h5Seurat file.',call.=FALSE)
  if(table$otype[entry]=='H5I_DATASET')return(HDF5Array::HDF5Array(path,group))
  dimensions <- rhdf5::h5readAttributes(path,group)$dims
  if(is.null(dimensions))stop('Sparse h5Seurat assay needs its dims attribute.',call.=FALSE)
  HDF5Array::H5SparseMatrix(HDF5Array::H5SparseMatrixSeed(path,group,
    dim=as.integer(dimensions),sparse.layout='CSC'))
}

L0glm_data_convert_vcf <- function(path,control) {
  if(length(path)!=1L||!file.exists(path)||is.null(control$cache))
    stop('VCF input needs a file and control$cache output prefix for a PGEN conversion.',call.=FALSE)
  exe <- if(is.null(control$plink))Sys.which('plink2')else control$plink
  if(length(exe)!=1L||!nzchar(exe)||!file.exists(exe))stop('Supply control$plink pointing to PLINK 2.',call.=FALSE)
  prefix <- normalizePath(control$cache,winslash='/',mustWork=FALSE)
  if(!dir.exists(dirname(prefix)))stop('Create the VCF cache directory first.',call.=FALSE)
  path <- normalizePath(path,winslash='/',mustWork=TRUE)
  information <- file.info(path)
  provenance <- list(path=path,bytes=information$size,modified=information$mtime,dosage=control$dosage)
  manifest <- paste0(prefix,'.l0glm.rds');output <- paste0(prefix,'.pgen')
  if(file.exists(output)) {
    if(!file.exists(manifest)||!identical(readRDS(manifest),provenance))
      stop('VCF cache already exists with different provenance; choose a new prefix.',call.=FALSE)
    return(output)
  }
  if(any(file.exists(paste0(prefix,c('.pvar','.pvar.zst','.psam','.log','.conversion.log','.l0glm.rds')))))
    stop('The VCF cache prefix already contains files; choose a new prefix.',call.=FALSE)
  if(!is.null(control$dosage) && !identical(control$dosage,'DS'))
    stop('dosage must be NULL for GT hardcalls or "DS" for expected ALT dosage.',call.=FALSE)
  args <- c('--vcf',shQuote(path),if(!is.null(control$dosage))'dosage=DS',
    '--make-pgen','--out',shQuote(prefix))
  status <- system2(exe,args,stdout=paste0(prefix,'.conversion.log'),stderr=paste0(prefix,'.conversion.log'))
  if(!identical(status,0L)||!file.exists(output))stop('PLINK VCF conversion failed; inspect the conversion.log.',call.=FALSE)
  saveRDS(provenance,manifest);output
}
