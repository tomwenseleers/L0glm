# MSP IMPORT USING OWN FAST READER ####


# helper function to import msp as (Erah style) list ####
library(stringr)
library(stringi)
library(readr)

import_msp <- function(filepath, 
                       rds_path = paste0(filepath, ".rds"), 
                       default_column = "SemiStdNP") { # assumed RI if column unspecified, SemiStdNP / StdNP / StdPolar / Any
  # Check if cached RDS exists
  if (file.exists(rds_path)) {
    message("Reading from cached RDS file...")
    return(readRDS(rds_path))
  }
  
  msp_text <- readr::read_file(filepath)
  msp_text <- str_replace_all(msp_text, "\r\n|\r", "\n")
  compound_blocks <- str_split(msp_text, regex("(?im)^\\s*NAME:", multiline = TRUE))[[1]]
  compound_blocks <- compound_blocks[nzchar(str_trim(compound_blocks))]
  
  compounds <- vector("list", length(compound_blocks))
  pb <- txtProgressBar(min = 0, max = length(compound_blocks), style = 3)
  
  for (i in seq_along(compound_blocks)) {
    # print(i)
    setTxtProgressBar(pb, i)
    column <- default_column
    block <- compound_blocks[i]
    lines <- stri_split_lines(block)[[1]]
    lines <- stri_trim_both(lines)
    lines <- lines[lines != ""]
    
    compound <- list()
    synon_vec <- c()
    
    compound[["NAME"]] <- lines[1]  # Use first line as NAME
    line_idx <- 2
    while (line_idx <= length(lines)) {
      line <- lines[line_idx]
      if (stri_detect_fixed(line, ":")) {
        parts <- stri_split_fixed(line, ":", n = 2)[[1]]
        field <- stri_trans_toupper(stri_trim_both(parts[1]))
        field <- stri_replace_last_fixed(field, "#", "") 
        if (field=="SYNONYMS") field <- "SYNON"
        if (field=="CAS") field <- "CASNO"
        if (field=="COMMENT") field <- "COMMENTS"
        if (field=="RETENTION_INDEX") field <- "RI" # for MoNa MSP files
        value <- stri_trim_both(parts[2])
        if (field %in% c("RI", "RT", "RRT", "MW", "EXACTMASS"))  {
          if (value=="") { value = NA } else {
            if (stri_detect_fixed(value, "polar")) { # deal with RI: 1777 (polar) etc in MACE
              column <- "StdPolar"
              value <- stri_extract_first_regex(value, "^[^\\s]+")
            }
            if (stri_detect_fixed(value, ".")&stri_detect_fixed(value, ","))  { # deal with European or US style nrs
              value <- stri_replace_all_fixed(value, ".", "") }
            if (is.na(as.numeric(stri_replace_last_fixed(value, ",", ".")))) { print(i);print(value) }
            value <- as.numeric(stri_replace_last_fixed(value, ",", ".")) }
            if (field=="RI") field <- paste0(field, "_", column)
        }
        
        if (stri_cmp_equiv(field, "SYNON")) {
          synon_vec <- c(synon_vec, value)
        } else if (stri_cmp_equiv(field, "NUM PEAKS")) {
          num_peaks <- as.integer(value)
          peak_lines <- lines[(line_idx + 1):(line_idx + num_peaks)]
          peak_lines <- peak_lines[!is.na(peak_lines)]
          if (stri_count_fixed(peak_lines[[1]], ";")>1) { # MACE formatted MSP
            peak_lines <- unlist(stri_split_regex(peak_lines, "(?<=;)\\s+"))
          } 
          if (stri_detect_fixed(peak_lines[[1]], ";", max_count = 1L)) { # remove trailing semicolon
            peak_lines <- stri_replace_last_fixed(peak_lines, ";", "") 
          }
          if (stri_detect_fixed(peak_lines[[1]], ":", max_count = 1L)) { # GOLM formatted MSP
            peak_lines <- unlist(stri_split_fixed(peak_lines, " ", simplify = FALSE))
            peak_lines <- stri_replace_first_fixed(peak_lines, ":", " ")
          }
          peak_matrix <- stri_split_fixed(peak_lines, " ", simplify = TRUE) # stri_split_regex(peak_lines, "\\s+", simplify = TRUE)
          spectrum_data <- matrix(as.numeric(peak_matrix[,1:2]), ncol=2)
          colnames(spectrum_data) <- c("mz", "int")
          compound[["SPECTRUM"]] <- spectrum_data
          line_idx <- line_idx + num_peaks  # Skip peak lines
        } else {
          compound[[field]] <- value
        }
      }
      line_idx <- line_idx + 1
    }
    
    if (length(synon_vec) > 0) compound[["SYNON"]] <- stri_join(synon_vec, collapse = "|")
    compounds[[i]] <- compound
  }
  close(pb)
  
  saveRDS(compounds, rds_path)
  message("Parsed MSP file saved as RDS.")
  return(compounds)
}



library(stringi)
library(stringr)
library(tibble)
library(data.table)
library(fst)

# ------------------------------------------------------------------
# helper: aggregate a numeric 2-column matrix to unit-mass
# ------------------------------------------------------------------
collapse_spectrum <- function(mat, round_fun = round) {
  if (is.null(mat) || nrow(mat) == 0) return(NA_character_)
  mass  <- round_fun(mat[, 1])
  inten <- mat[, 2]
  agg   <- rowsum(inten, mass, reorder = FALSE)           # fast
  paste(paste(rownames(agg), agg[, 1], sep = " "),
        collapse = "; ")
}

# ------------------------------------------------------------------
# main import function
# ------------------------------------------------------------------
import_msp_tbl <- function(filepath,
                           fst_path       = paste0(filepath, ".fst"),
                           default_column = "SemiStdNP", # assumed RI if column unspecified, options SemiStdNP / StdNP / StdPolar / Any
                           fields         = c("NAME", "MW", "INCHIKEY",
                                              "FORMULA", "CASNO",
                                              "RI_SemiStdNP", "RI")) {
  
  ## ---------- fast return if cache exists ----------
  if (file.exists(fst_path)) {
    message("Reading cached .fst …")
    return(read_fst(fst_path, as_tibble = TRUE))
  }
  
  ## ---------- read file and split into blocks ----------
  msp_text <- readr::read_file(filepath)
  msp_text <- str_replace_all(msp_text, "\r\n|\r", "\n")
  blocks   <- str_split(msp_text, regex("(?im)^\\s*NAME:", multiline = TRUE))[[1]]
  blocks   <- blocks[nzchar(str_trim(blocks))]
  n_blocks <- length(blocks)
  
  # data.table list to fill
  rows <- vector("list", n_blocks)
  pb   <- txtProgressBar(min = 0, max = n_blocks, style = 3)
  
  for (idx in seq_len(n_blocks)) {
    
    setTxtProgressBar(pb, idx)
    column <- default_column
    
    lines <- stri_split_lines(blocks[idx])[[1]] |> stri_trim_both()
    lines <- lines[lines != ""]
    
    rec   <- list()
    rec$NAME <- lines[1]
    syn   <- character()
    
    line_idx <- 2
    while (line_idx <= length(lines)) {
      
      ln <- lines[line_idx]
      
      if (stri_detect_fixed(ln, ":")) {
        
        parts <- stri_split_fixed(ln, ":", n = 2)[[1]]
        field <- stri_trans_toupper(str_trim(parts[1]))
        field <- stri_replace_last_fixed(field, "#", "")
        if (field == "SYNONYMS")        field <- "SYNON"
        if (field == "CAS")             field <- "CASNO"
        if (field == "COMMENT")         field <- "COMMENTS"
        if (field == "RETENTION_INDEX") field <- "RI"
        
        val <- str_trim(parts[2])
        
        ## numeric normalisation ------------------------------------------
        if (field %in% c("RI", "RT", "RRT", "MW", "EXACTMASS")) {
          if (val == "")  val <- NA
          if (!is.na(val)) {
            if (stri_detect_fixed(val, "polar")) {
              column <- "StdPolar"
              val <- stri_extract_first_regex(val, "^[^\\s]+")
            }
            if (stri_detect_fixed(val, ".") && stri_detect_fixed(val, ",")) {
              val <- stri_replace_all_fixed(val, ".", "")
            }
            val <- as.numeric(stri_replace_last_fixed(val, ",", "."))
          }
          if (field == "RI") field <- paste0(field, "_", column)
        }
        
        if (field == "SYNON") {
          syn <- c(syn, val)
        } else if (field == "NUM PEAKS") {
          
          n_peak   <- as.integer(val)
          pk_lines <- lines[(line_idx + 1):(line_idx + n_peak)]
          
          ##  ---------- reformat peak lines (handles all flavours) ----------
          pk_lines <- pk_lines[!is.na(pk_lines)]
          
          # MACE : many peaks per line with semicolons
          if (stri_count_fixed(pk_lines[1], ";") > 1)
            pk_lines <- unlist(stri_split_regex(pk_lines, "(?<=;)\\s+"))
          
          # strip trailing semicolons
          if (stri_detect_fixed(pk_lines[1], ";", max_count = 1L))
            pk_lines <- stri_replace_last_fixed(pk_lines, ";", "")
          
          # Golm : mz:int pairs separated by spaces
          if (stri_detect_fixed(pk_lines[1], ":", max_count = 1L)) {
            pk_lines <- unlist(stri_split_fixed(pk_lines, " "))
            pk_lines <- stri_replace_first_fixed(pk_lines, ":", " ")
          }
          
          pk_mat  <- stri_split_fixed(pk_lines, " ", simplify = TRUE)
          pk_mat  <- matrix(as.numeric(pk_mat[, 1:2]), ncol = 2)
          colnames(pk_mat) <- c("mz", "int")
          
          rec$SPECTRUM <- collapse_spectrum(pk_mat)
          line_idx <- line_idx + n_peak
        } else {
          rec[[field]] <- val
        }
      }
      line_idx <- line_idx + 1
    }
    
    if (length(syn)) rec$SYNON <- paste(syn, collapse = "|")
    rows[[idx]] <- rec
  }
  close(pb)
  
  tbl <- rbindlist(rows, fill = TRUE) |> as_tibble()
  
  ## keep only requested fields + SPECTRUM -------------------------------
  needed <- unique(c(fields, "SPECTRUM"))
  missing <- setdiff(needed, names(tbl))
  tbl[, missing] <- NA
  tbl <- tbl[, needed]
  
  ## ---------- write cache (.fst) & return ------------------------------
  write_fst(tbl, fst_path)
  message("Parsed MSP file cached as .fst")
  tbl
}

# note: Ierah::importMSP & mspReader::readMsp did not work on MSP files below






# helper function to stack spectra in sparse matrix, with annotations stored in tidy table annot attribute ####

library(Matrix)
library(data.table)
library(tibble)

spectra_to_matrix <- function(lib,
                              round_fun = round,
                              fields = c("NAME", "MW", "INCHIKEY"))
{
  n_comp <- length(lib)
  
  ## ---------- 1.  collect unique integer masses ----------
  unique_masses <- sort(unique(unlist(
    lapply(lib, function(z)
      if (!is.null(z$SPECTRUM))
        round_fun(z$SPECTRUM[, 1])
      else integer())
  )))
  col_index <- setNames(seq_along(unique_masses),
                        as.character(unique_masses))
  
  ## ---------- 2.  build triplet list (with progress bar) ----------
  triplets <- vector("list", n_comp)
  pb <- txtProgressBar(min = 0, max = n_comp, style = 3)
  
  for (idx in seq_len(n_comp)) {
    
    setTxtProgressBar(pb, idx)
    
    spec <- lib[[idx]]$SPECTRUM
    if (is.null(spec)) next
    
    unit_mz <- round_fun(spec[, 1])
    ints    <- spec[, 2]
    
    agg  <- rowsum(ints, unit_mz, reorder = FALSE)     # fast aggregation
    mass <- as.numeric(rownames(agg))
    
    triplets[[idx]] <- data.table(
      i = idx,
      j = col_index[as.character(mass)],
      x = as.numeric(agg[, 1])
    )
  }
  close(pb)
  
  triplets_dt <- rbindlist(triplets, use.names = FALSE)
  
  ## ---------- 3.  build sparse matrix ----------
  spmat <- sparseMatrix(
    i    = triplets_dt$i,
    j    = triplets_dt$j,
    x    = triplets_dt$x,
    dims = c(n_comp, length(unique_masses)),
    dimnames = list(
      as.character(seq_len(n_comp)),           # row names: 1…n
      as.character(unique_masses)              # col names: unit masses
    ),
    giveCsparse = TRUE
  )
  
  ## ---------- 4.  annotation tibble (atomic) ----------
  annot <- as_tibble(setNames(
    lapply(fields, function(f) {
      vapply(lib, function(rec) {
        v <- rec[[f]]
        if (is.null(v) || length(v) == 0) {              # missing
          if (f %in% c("MW", "EXACTMASS")) NA_real_ else NA_character_
        } else if (length(v) == 1) {
          v
        } else {                                         # multiple -> collapse
          paste(v, collapse = "|")
        }
      },
      if (f %in% c("MW", "EXACTMASS")) numeric(1) else character(1))
    }),
    nm = fields))
  
  attr(spmat, "annot") <- annot
  spmat
}


# Example usage:
# CASNO -> insert hyphens & convert "0" to NA

# TODO import CSV files with in-silico predicted spectra & RIs from GCMSid

# NIST23 lib imported from Shimadzu format
system.time(lib1a <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/NIST_2023_converted/from_shimadzu/NIST23M1.msp")) # main1, 111s # TODO still put hyphens in CAS nrs like 287230 for Cyclobutane
system.time(lib1b <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/NIST_2023_converted/from_shimadzu/NIST23M2.msp")) # main2, 33s
system.time(lib1c <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/NIST_2023_converted/from_shimadzu/NIST23R.msp")) # repl, 25s
# NIST23 lib imported from Agilent format
system.time(lib1d <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/NIST_2023_converted/from_agilent/NIST23.msp")) # 111s # TODO still put hyphens in CAS nrs like 287230 for Cyclobutane

# GOLM db
system.time(lib2 <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/GOLM/GMD_20111121_VAR5_ALK_MSP.txt")) # 2.1s # TODO parse RI from Synon: part "RI: 973,23"

# MoNA db
system.time(lib3 <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/MONA/MoNA-export-GC-MS_Spectra.msp")) # 7s # TODO still parse INCHI & SMILES & CAS in COMMENTS

# MACE db
system.time(lib4 <- import_msp(filepath="C:/Users/bherr/Dropbox/GCMS/GCMS_libraries/MACE/MACE_r05.mace")) # 0.19s # TODO parse Smiles from comments " Smiles=CCCCCCCCCCCCC[C@@H](C[C@@H](C(OCCCCCCCCCCCCCC)=O)C)C "

# check available fields
names(lib1a[[1]])
# [1] "NAME"         "COMMENTS"     "RT"           "RRT"          "RI_SemiStdNP" "FORMULA"      "MW"          
# [8] "CASNO"        "SMILES"       "DB"           "REFID"        "SPECTRUM"     "SYNON"  
names(lib1b[[1]])
# [1] "NAME"         "COMMENTS"     "RT"           "RRT"          "RI_SemiStdNP" "FORMULA"      "MW"          
# [8] "CASNO"        "SMILES"       "DB"           "REFID"        "SPECTRUM" 
names(lib1c[[1]])
# [1] "NAME"         "COMMENTS"     "RT"           "RRT"          "RI_SemiStdNP" "FORMULA"      "MW"          
# [8] "CASNO"        "SMILES"       "DB"           "REFID"        "SPECTRUM"     "SYNON" 
names(lib1d[[9]]) # 262s
# [1] "NAME"         "COMMENTS"     "RT"           "RRT"          "RI_SemiStdNP" "FORMULA"      "MW"          
# [8] "CASNO"        "SMILES"       "DB"           "REFID"        "SPECTRUM" 
# but length(which(unlist(sapply(1:length(lib3), function (i) lib3[[i]]$RI_SemiStdNP[[1]]))!=0)) only 3255
# TODO Parse COMMENTS e.g. COMMENTS:  SemiStdNP=317 StdNP=323/4/2 StdPolar=1053/14/2
# or COMMENTS:  SemiStdNP=200 StdNP=200 StdPolar=200 ; NIST MS# 18824, Seq# R83

names(lib2[[1]])
# "NAME"     "CASNO"    "COMMENTS" "DB"       "SPECTRUM" "SYNON"  
# PARSE RI from Synon: RI: 988 (replace comma with .)
# PARSE CAS from Synon: METB CAS: 298-12-4 (but CAS of precursor / underivatised)
# PARSE INCHI from Synon: METB InChI: InChI=1S/C2H2O3/c3-1-2(4)5/h1H,(H,4,5)
# PARSE INCHIKEY from Synon: METB InChIKey: HHLFWLYXYJOTON-UHFFFAOYSA-N
names(lib3[[1]]) 
# [1] "NAME"            "DB"              "INCHIKEY"        "SPECTRUM_TYPE"   "INSTRUMENT_TYPE"
# [6] "INSTRUMENT"      "ION_MODE"        "FORMULA"         "MW"              "EXACTMASS"      
# [11] "COMMENTS"        "SPECTRUM"        "SYNON"
# -> PARSE SMILES FROM COMMENT (use computed SMILES?)
# -> PARSE CAS FROM COMMENT
# change FIELD Retention_index OCCURING IN SOME RECS to RI_SemiStdNP (and make sure RI is numeric)
# length(which(sapply(1:length(lib3), function (id) grepl("retention index", lib3[[id]]$COMMENTS)))) # 1104, but mix of kovats and fame
# 
# parse from comments "kovats retention index=978" "column type=ZB-1MS non-polar" "column=5%-phenyl-95%-dimethylpolysiloxane capillary column"
names(lib4[[1]])
# [1] "NAME"         "CASNO"        "INCHIKEY"     "FORMULA"      "MW"           "EXACTMASS"    "RI_SemiStdNP"
# [8] "COMMENTS"     "SPECTRUM"     "SYNON" 
# parse from comments
# Smiles=O=C1CC(C)CC(C)O1
# RI: 1793 (polar)

# TODO add chemical structures if not available, maybe using Chemical Identifier Resolver??
# library(webchem)
# cir_query(identifier=lib1d[[1001]]$CASNO, representation="smiles", resolver="cas_number", match="first")
# 
# TODO add RIs from Pubchem or GCMSid
