# Function converts a dataframe into bed interval data
# Requires valr package:
#   https://cran.r-project.org/web/packages/valr/index.html
#   https://github.com/rnabioco/valr

convertToBed <- function(object, columns = c(1:3)){
  object_field_select <- object[, columns]
  tmp_file <- file.path(tempdir(), "temp_bed.tsv")
  write.table(x = object_field_select,
              file = tmp_file,
              sep = '\t', col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  bed_format <- read_bed(tmp_file,
                         n_fields = length(columns))
  unlink(tmp_file)
  return(bed_format)
}
