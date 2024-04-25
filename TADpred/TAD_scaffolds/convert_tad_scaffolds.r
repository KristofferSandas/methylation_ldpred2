# Script for converting files from UCSC LiftOver to the
# format they need to be in TADpred:
# chromosome  start end 
# For example:
# chr1	735380	1285380

# Open the input file for reading
input_file <- file('TAD_scaffolds/tad_scaffold_conversion.txt', 'r')

# Open a new file for writing
output_file <- file('TAD_scaffolds/tad_scaffold_hg19.txt', 'w')

# Read each line from the input file
while (length(line <- readLines(input_file, n = 1)) > 0) {
  # Split the line by tabs
  parts <- unlist(strsplit(line, '\t'))
  # Keep only the first three columns
  modified_line <- paste(parts[1:3], collapse = '\t')
  # Write the modified line to the output file
  cat(modified_line, '\n', file = output_file)
}

# Close the files
close(input_file)
close(output_file)

# Remove variables
rm(input_file, output_file, line, modified_line, parts)
