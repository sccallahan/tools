###################################################################################
# Author: Carson Callahan
#
# Date: 2020-08-27
#
# Purpose: Parse Incucyte confluency data and format for input into
#         GRmetrics package/web application
#
# Notes:
# 1) GRmetrics: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4887336/
# 2) Web Application: http://www.grcalculator.org/grtutorial/Home.html
# 3) R package: https://www.bioconductor.org/packages/release/bioc/html/GRmetrics.html
#################################################################################

#### Brief overview ####
# In general, this takes the .txt outputs from an Incucyte and formats them for GRmetrics.
# The function will read in each file of a certain pattern (using `list.files`), format it, then save to disk.


#### NB: this code is not fully tested ####
##    several sections of this code may need to be changed
##    in order to fit individual needs; this is just
##    the setup that works for my usual layout --
##    areas that may need to be changed are usually marked
##    with 4 hashes


#### Explanation of arguments ####
# 1. dir = directory where .txt outputs are stored.
# 2. skip = how many lines to skip before reading. Incucyte seems to have some blank spaces and metadata
#           that aren't strictly needed for this.
# 3. cell_line = cell line name
# 4. treatment = treatment name
# 5. max_conc = integer; maximum concentration used (in uM)
# 6. dilution = integer; the fold dilution used between wells
# 7. series = concentrations used; only if not using max_conc + dilution;
#             e.g. c(100, 80, 60, 10)
# 8. control_well_num = integer; the column number used as control (e.g. DMSO in column 2)
# 9. control2_well_num = integer; column number for the second control (if used)
# 10. time = integer; the duration of the treatment; By default, this is pulled from the maximum incucyte "elapsed" value


#### Usage ####
# The easiest way to use this is probably as follows:
# 1) Save parseIncucyte.R to whatever directory you choose.
# 2) Source the file in your code by running:
#       source("path/to/script/directory/parseIncucyte.R")


parseIncucyte <- function(dir, skip = 2, cell_line, treatment,
                          max_conc = NULL, dilution = NULL, series = NULL,
                          control_well_num, control2_well_num = NULL,
                          time = max(dat$Elapsed)){
  #### change pattern if needed ####
  files <- list.files(path = dir, pattern = "*_well.txt")
  
  for (file in files){
    
    # ------------------------------ #
    #### read in file, grab times ####
    # ------------------------------ #
    # time can also be manually set if desired
    dat <- read.table(file = paste0(dir, file), sep = '\t', header = TRUE, skip = skip)
    time_final <- max(dat$Elapsed)
    time_start <- min(dat$Elapsed)
    
    # we generally don't need intermediate times for this analysis
    # filter for time0 and timeFinal
    dat_timepoints <- dat[dat$Elapsed == time_final | dat$Elapsed == time_start, ]
    
    # -------------------------------------------- #
    #### remove std dev columns and dates/times ####
    # -------------------------------------------- #
    dat_prune <- dat_timepoints[, (!grepl("Std", colnames(dat_timepoints)))]
    dat_prune <- dat_prune[, -c(1:2)]
    
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
    dat_drug <- dat_prune[, !grepl(control, colnames(dat_prune))]
    
    # allow for dual controls
    if(!is.null(control2_well_num)){
      control2 <- control2_well_num
      dat_drug <- dat_drug[, !grepl(control2, colnames(dat_drug))]
    }
    
    # ---------------------------- #
    #### get cell counts and CV ####
    # ---------------------------- #
    cell_count__time0 <- apply(dat_prune, 1, mean)[1]
    time0_stdev <- apply(dat_prune, 1, sd)[1]
    time0_cv <- time0_stdev/cell_count__time0*100
    message("Coefficient of variation for initial cell counts is ", round(time0_cv, digits = 3), "%")
    
    # average control counts at timeFinal
    cell_count__ctrl <- apply(dat_prune[2, grepl(control, colnames(dat_prune))], 1, mean)
    
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
    final <- data.frame(cell_line = rep(cell_line, times = nrow(tmp)),
                        treatment = rep(treatment, times = nrow(tmp)),
                        perturbation = rep(0, times = nrow(tmp)),
                        #### change replicates here too if needed ####
                        replicate = c(rep(1, times = nrow(tmp)/3), rep(2, times = nrow(tmp)/3), rep(3, times = nrow(tmp)/3)),
                        time = rep(time_final, times = nrow(tmp)),
                        concentration = tmp$concentration,
                        cell_count = tmp$cell_count,
                        cell_count__ctrl = rep(cell_count__ctrl, times = nrow(tmp)),
                        cell_count__time0 = rep(cell_count__time0, times = nrow(tmp)))
    colnames(final)[6] <- "concentration"
    
    write.table(x = final, file = paste0(dir, file_path_sans_ext(file), "_for_GR.txt"), sep = '\t',
                quote = FALSE, row.names = FALSE)
  }
}