



# 1.0 Load required libraries----
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(viridis)
library(gridExtra)



# 2.0 This loads a .txt file that consists of one column populated with primate genera. Depending on the study, this could be exchanged, or added to, with genera from other taxa i.e. birds, fish, etc.
primate_genera<- read.table("working_data/primate_genera.txt", quote="\"", comment.char="") 

# 2.1 This creates a character string that makes a pattern out of the column of genera. 
pattern<- paste0("\\b", primate_genera$V1, "\\b", collapse = "|")



# 3.0 Define a function that processes files in the Mascot output format. This will process all files found in the input directory.
process_file <- function(file_path) {
  
  # 3.1 Extract the sample ID from the file name. All of the files in this study have a six-character alphanumeric code at the start of their files name: two letters followed by four numbers.
  # This extraction steps allows the data to be "tagged" with it's relevant meta-data. In an effort to make this the standard in ancient proteomics, I have hard-coded this extraction criterion.
  # Therefore, when the pipeline is developed, anyone who wishes to use it will have to follow this naming protocol. This should also facilitate different labs being able to reproduce each others'
  # results more easily.
  sample_id <- str_extract(basename(file_path), "(?i)^[A-Z]{2}\\d{4}")
  
  # 3.2 Read the data for a specific sample file----
  data <- fread(file_path, header = TRUE)
  
  # 3.3 Check if the required columns exist. If they don't, the code prints the warning message "File XX does not contain the required columns. Skipping. Check your data!"
  if (!all(c("prot_desc", "pep_seq") %in% colnames(data))) {
    warning(paste("File", file_path, "does not contain the required columns. Skipping. Check your data!"))
    return(NULL)
  }
  
  # 3.4 The following code performs the actual filtering and transformations.
  data %>%
    # 3.5 Select the protein description and peptide sequence columns from the data file.
    select("prot_desc", "pep_seq") %>%
    # 3.6 Filter the protein description column so as to select only those genera defined in the pattern character string in section 2.1
    filter(grepl(pattern, prot_desc)| grepl("\\|[^|]*_HUMAN", prot_desc)) %>% 
    # 3.7 Remove records that are fragments, trypsin catalyst, and contaminated.
    filter(!grepl("TRYP_PIG", prot_desc)) %>%
    mutate(combinations = paste(prot_desc, pep_seq, sep = "-")) %>%
    
    # 3.8 Normalize whitespace in prot_desc
    mutate(prot_desc = str_squish(prot_desc)) %>%
    
    # 3.9 Group records by protein description. This allows for different peptides that code for the same protein to be identified.
    group_by(prot_desc) %>%
    
    # 3.10 Filter the data so that only proteins that are recovered by at least 2 unique peptide sequences are kept in the dataset.
    filter(n_distinct(pep_seq) >= 2) %>%
    ungroup() %>%
    group_by(prot_desc) %>%
    
    # 3.11 Calculate PSM.
    summarise(value = n()) %>%
    mutate(prot_desc = ifelse(
      str_detect(prot_desc, "\\|[^|]*_HUMAN"),
      str_extract(prot_desc, "(?<=\\|)[^|]*_HUMAN"),
      prot_desc)) %>% 
    
    # 3.12 Remove unnecessary characters from the protein description.
    mutate(across(prot_desc, str_remove, pattern = c("\\OS.*"))) %>%
    mutate(sample_ID = sample_id) %>%
    
    # 3.13 Rename variables so they are more intuitive.
    rename(protein_name = "prot_desc",
           psm = "value") %>%
    
    # 3.14 Select the variables required for further downstream analyses.
    select("protein_name", "sample_ID", "psm") %>% 
    
    # 3.15 Change some protein names to a more descriptive nomenclature.
    mutate(protein_name = if_else(protein_name %in% "ALBU_HUMAN", "Albumin", protein_name),
           
           protein_name = if_else(protein_name %in% "CATG_HUMAN", "Cathepsin G", protein_name),
           
           protein_name = if_else(protein_name %in% "CRP_HUMAN", "C-reactive protein", protein_name))
}



# 4.0 Process all files in the input directory and store results in a list.
# 4.1 Set input directory where the data files are stored.
input_directory<- "working_data/Seli_Data/"

# 4.2 Create a list of only those files that contain the Mascot output suffix. This ensures that only the required data files (and not other files such as a meta-data files) are processed.
csvs<- list.files(input_directory, pattern = "\\.csv.0.01.out.csv", full.names = TRUE)

# 4.3 Iterativily run the process_file() function on each data file listed in the "csvs" list.
data_Seili<- map(csvs, process_file) %>% 
  compact() %>% # Removes NULL elements (files that were skipped)
  bind_rows() %>% # Combine all processed data into a single dataframe
  
  # 4.4 Select only those records with a PSM equal or greater to two
  filter(psm >= 2) %>% 
  mutate(study = "Seili") 


# 5.0 Combine each separate study's records into a single dataframe for downstream analyses.
data_all<- rbind.data.frame(data_Mickleburgh, data_Ntasi, data_Sawafuji, data_Seili)
setwd("output/")
write.csv(data_all, "all_studies_and_samples.csv", row.names = FALSE)
