# code to take [project_name]_picked_mips.txt file from MIPGEN output
# and output a .csv file for an initial panel compatible with MIPGEN_scoring

# USER INPUTS
## filepath to initial picked_mips.txt file from MIPGEN
in_filepath <- "mipgen_practice/practice_design.picked_mips.txt"
## output path
out_filepath <- "practice_design_mipgen.csv"

library(data.table)
library(tidyverse)

initial_panel <- fread(in_filepath)

the_chosen %>% 
  # length of ligation and extension arms
  mutate(lig_len = str_length(lig_probe_sequence), 
         ext_len = str_length(ext_probe_sequence)) %>% 
  select(mip_start=mip_scan_start_position, 
         mip_end=mip_scan_stop_position, 
         ext_len, 
         lig_len, 
         probe_strand) %>%
  write.csv(out_filepath, row.names=FALSE)