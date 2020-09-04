###################################################################################
# Author: Carson Callahan
#
# Date: 2020-09-03
#
# Purpose: Functions for integrating incucyte data with GRmetrics
#
# Notes:
# 1) GRmetrics: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4887336/
# 2) Web Application: http://www.grcalculator.org/grtutorial/Home.html
# 3) R package: https://www.bioconductor.org/packages/release/bioc/html/GRmetrics.html
#################################################################################


#' @title parseIncucyte
#' 
#' @description Converts raw incucyte text output into a GRmetrics-friendly format, then writes them to disk with
#' the "*for_GR.tsv" suffix/extension. MAY NEED TO BE TWEAKED BY THE USER. NOT BROADLY TESTED.
#' 
#' @param dir Directory where incucyte text outputs are stored
#' @param skip Number of lines to skip before reading in data. Default = 2.
#' Incucyte usually has some unneeded rows that make reading in files difficult.
#' @param pattern Pattern to use for listing input files. Example: "*_data.txt"
#' @param treatment Name of treatment. For now, just a single treatment is accepted.
#' @param max_conc Integer; maximum concentration of the drug used in uM
#' @param dilution Integer; fold dilution (e.g. if well 1 is 50uM and well 2 is 25uM, etc., this would be a 2)
#' @param series Vector; if not using a consistent dilution, manually enter concentrations (e.g. c(100, 50, 10, 5, 2)); Default = NULL
#' @param control_well_num Integer; the column number used as control (e.g. DMSO in column 2)
#' @param control2_well_num Integer; column number for the second control (if used); Default = NULL
#' @param rm_artefact Integer or vector; columns with obvious technical artefacts to be excluded
#' @param time_0 Integer; the time recorded by incucyte where you have set your experimental time_0;
#' Default = minimum value from "Elapsed" column.
#' @param time_F Integer; the time recorded by incucyte where you have set your experimental time_final;
#' Default = maximum value from "Elapsed" column
#' @param time_E Integer; the amount of time your experiment took; Default = time_F;
#' Change if incucyte times are not correct
#'
#' @return Nothing. Writes files to disk.
#' @export
#'
parseIncucyte <- function(dir, skip = 2, pattern, treatment,
                          max_conc = NULL, dilution = NULL, series = NULL,
                          control_well_num, control2_well_num = NULL, rm_artefact = NULL,
                          time_0 = NULL, time_F = NULL, time_E = NULL){
  
  inc_files <- list.files(path = dir, pattern = pattern)
  
  for (file in inc_files){
    
    # ------------------------------ #
    #### read in file, grab times ####
    # ------------------------------ #
    dat <- read.table(file = paste0(dir, file), sep = '\t', header = TRUE, skip = skip)
    
    # run checks for the time values
    if (is.null(time_0)){
      time_start <- min(dat$Elapsed)
    } else{
      time_start <- time_0
    }
    
    if (is.null(time_F)){
      time_final <- max(dat$Elapsed)
    } else{
      time_final <- time_F
    }
    
    if (is.null(time_E)){
      time_E <- max(dat$Elapsed)
    } else{
      time_E <- time_E
    }
    # time_final <- time_F
    # time_start <- time_0
    
    # we generally don't need intermediate times for this analysis
    # filter for time0 and timeFinal
    dat_timepoints <- dat[dat$Elapsed == time_final | dat$Elapsed == time_start, ]
    
    # -------------------------------------------- #
    #### remove std dev columns and dates/times ####
    # -------------------------------------------- #
    dat_prune <- dat_timepoints[, (!grepl("Std", colnames(dat_timepoints)))]
    dat_prune <- dat_prune[, -c(1:2)]
    
    
    # -------------------------------------- #
    #### remove any technically bad wells ####
    # -------------------------------------- #
    
    if (!is.null(rm_artefact)){
      dat_prune_final <- dat_prune[, !grepl(paste(rm_artefact, collapse = "|"), colnames(dat_prune))]
    } else {
      dat_prune_final <- dat_prune
    }
    
    # ------------------------------- #
    #### calculate dilution series ####
    # ------------------------------- #
    # I assume the follwing plate setup:
    # Outer ring is PBS for edge effect
    # 2 columns are controls, giving 8 total concentrations
    # Can be manually set
    if(is.null(series)){
      if(is.null(max_conc) || is.null(dilution)){
        stop("Must provide max concentation and dilution factor if no dilution series provided")
      } else {
        message("Series not provided, using provided max and dilution arguments")
        #### change dilution number/strategy as needed ####
        conc <- c(max_conc, max_conc/dilution, max_conc/(dilution^2),
                  max_conc/(dilution^3), max_conc/(dilution^4), max_conc/(dilution^5), max_conc/(dilution^6), max_conc/(dilution^7))
      }
    }
    
    if(is.null(max_conc) || is.null(dilution)){
      if(is.null(series)){
        stop("Must provide dilution series if no max and dilution are provided")
      } else {
        message("Using provided dilution series")
        conc <- series
      }
    }
    
    # --------------------------------------- #
    #### control well NUMBER (i.e. column) ####
    # --------------------------------------- #
    control <- control_well_num
    
    # isolate just the drug dilution data
    dat_drug <- dat_prune_final[, !grepl(control, colnames(dat_prune_final))]
    
    # allow for dual controls
    if(!is.null(control2_well_num)){
      control2 <- control2_well_num
      dat_drug <- dat_drug[, !grepl(control2, colnames(dat_drug))]
    }
    
    
    # ---------------------------- #
    #### get cell counts and CV ####
    # ---------------------------- #
    cell_count__time0 <- apply(dat_prune_final, 1, mean)[1]
    time0_stdev <- apply(dat_prune_final, 1, sd)[1]
    time0_cv <- time0_stdev/cell_count__time0*100
    message("Coefficient of variation for initial ", file, " cell counts is ", round(time0_cv, digits = 3), "%")
    
    # average control counts at timeFinal
    cell_count__ctrl <- apply(dat_prune_final[2, grepl(control, colnames(dat_prune_final))], 1, mean)
    
    # take timeFinal counts for drugs
    dat_drug_final <- dat_drug[2,]
    
    # -------------------------------------- #
    #### match counts with concentrations ####
    # -------------------------------------- #
    tmp <- as.data.frame(t(dat_drug_final))
    
    # I assume 3 replicates
    #### change reps if needed ####
    tmp$conc <- as.data.frame(rep(conc, times = 3))
    colnames(tmp) <- c("cell_count", "concentration")
    rownames(tmp) <- NULL
    
    # ----------------------------------- #
    #### make final dataframe and save ####
    # ----------------------------------- #
    
    ## extract cell line
    cell_line <- str_split(file, pattern = "_", simplify = TRUE)[1]
    
    ## make output df
    final <- data.frame(cell_line = rep(cell_line, times = nrow(tmp)),
                        treatment = rep(treatment, times = nrow(tmp)),
                        perturbation = rep(0, times = nrow(tmp)),
                        #### change replicates here too if needed ####
                        replicate = c(rep(1, times = nrow(tmp)/3), rep(2, times = nrow(tmp)/3), rep(3, times = nrow(tmp)/3)),
                        time = rep(time_E, times = nrow(tmp)),
                        concentration = tmp$concentration,
                        cell_count = tmp$cell_count,
                        cell_count__ctrl = rep(cell_count__ctrl, times = nrow(tmp)),
                        cell_count__time0 = rep(cell_count__time0, times = nrow(tmp)))
    colnames(final)[6] <- "concentration"
    
    write.table(x = final, file = paste0(dir, file_path_sans_ext(file), "_for_GR.tsv"), sep = '\t',
                quote = FALSE, row.names = FALSE)
  }
}



#' @title calcGRmetrics
#' 
#' @description Combines the outputs from `parseIncucyte` into a single dataframe, and runs a few simple
#' functions from the `GRmetrics` package. Writes GRvalues and GRmetrics files to disk. Writes dose-response curve plot
#' to disk if argument set to "TRUE".
#'
#' @param dir Directory where `parseIncucyte` text outputs are stored
#' @param pattern Pattern to use for listing input files. Example: ".tsv"
#' @param doseResponse Boolean indicating whether or not to produce the dose response curve plot. Default = FALSE
#'
#' @return A SummarizedExperiment object containing GR metrics. The same as returned by `GRfit`.
#' @export
#'
calcGRmetrics <- function(dir, pattern,
                          doseResponse = FALSE){
  
  # ------------------------------------------ #
  #### read in files and merge for GR input ####
  # ------------------------------------------ #
  
  # get files
  # full.names = TRUE for full path
  gr_files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
  
  # need to cat all files together for convenience and plots
  list_files <- lapply(gr_files, read.table, header = TRUE)
  frame <- do.call(rbind, list_files)
  
  # save this table
  write.table(x = frame, file = paste0(dir, "merged_incucyte_data.tsv"), sep = '\t',
              quote = FALSE, row.names = FALSE)
  
  
  # ------------------------------------------- #
  #### calculate GR metrics and make DR plot ####
  # ------------------------------------------- #
  drc_output <- GRfit(inputData = frame, groupingVariables = c("cell_line", "treatment"))
  
  # save GR values
  write.table(GRgetValues(drc_output), file = paste0(dir, "GR_values.tsv"), sep = '\t', quote = FALSE,
              row.names = FALSE)
  
  # save GR metrics
  write.table(GRgetMetrics(drc_output), file = paste0(dir, "GR_metrics.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE)
  
  # make dose response curve plot and save file
  if (doseResponse == TRUE){
    dr_plot <- GRdrawDRC(drc_output, plotly = FALSE)
    ggsave(dr_plot, filename = paste0(dir, "dose_response_curve.png"), dpi = 600,
           width = 12, height = 10)
  }
  
  # return the summarizedExperiment object
  return(drc_output)
  
}



#' @title incucytePipe
#' 
#' @description Simply pipes the workflow of `parseIncucyte` and `calcGRmetrics`.
#'
#' @param dir Directory where incucyte text outputs are stored
#' @param skip Number of lines to skip before reading in data. Default = 2.
#' Incucyte usually has some unneeded rows that make reading in files difficult.
#' @param pattern Pattern to use for listing input files. Example: "*_data.txt"
#' @param treatment Name of treatment. For now, just a single treatment is accepted.
#' @param max_conc Integer; maximum concentration of the drug used in uM
#' @param dilution Integer; fold dilution (e.g. if well 1 is 50uM and well 2 is 25uM, etc., this would be a 2)
#' @param series Vector; if not using a consistent dilution, manually enter concentrations (e.g. c(100, 50, 10, 5, 2)); Default = NULL
#' @param control_well_num Integer; the column number used as control (e.g. DMSO in column 2)
#' @param control2_well_num Integer; column number for the second control (if used); Default = NULL
#' @param rm_artefact Integer or vector; columns with obvious technical artefacts to be excluded
#' @param time_0 Integer; the time recorded by incucyte where you have set your experimental time_0;
#' Default = minimum value from "Elapsed" column.
#' @param time_F Integer; the time recorded by incucyte where you have set your experimental time_final;
#' Default = maximum value from "Elapsed" column
#' @param time_E Integer; the amount of time your experiment took; Default = time_F;
#' Change if incucyte times are not correct 
#' @param GR_pattern Directory where `parseIncucyte` text outputs are stored
#' @param doseResponse Boolean indicating whether or not to produce the dose response curve plot. Default = FALSE
#'
#' @return
#' @export
#'
#' @examples
incucytePipe <- function(dir, skip = 2, pattern,
                         treatment, max_conc = NULL, dilution = NULL, series = NULL,
                         control_well_num, control2_well_num = NULL, rm_artefact = NULL,
                         time_0 = NULL, time_F = NULL, time_E = NULL,
                         GR_pattern, doseResponse = FALSE){
  
  parseIncucyte(dir = dir, skip = skip, pattern = pattern, treatment = treatment,
                max_conc = max_conc, dilution = dilution, series = series,
                control_well_num = control_well_num, control2_well_num = control2_well_num, rm_artefact = rm_artefact,
                time_0 = time_0, time_F = time_F, time_E = time_E)
  
  calcGRmetrics(dir = dir, pattern = GR_pattern, doseResponse = doseResponse)
  
}













