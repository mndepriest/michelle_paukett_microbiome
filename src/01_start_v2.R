

# Load packages ####

  # Define packages needed
  list.of.packages <- c("phyloseq",
                        "ggplot2",
                        "plyr",
                        "dplyr",
                        # "microbiome",
                        "hrbrthemes",
                        "gcookbook",
                        "tidyverse",
                        "BiocManager",
                        "proto",
                        "vegan",
                        "RColorBrewer",
                        "reshape2",
                        "knitr",
                        # "DECIPHER",
                        "phangorn",
                        "readr",
                        "here",
                        "stringr",
                        "stringi",
                        "microbiome",
                        "readxl", 
                        "writexl",
                        "openxlsx")
  
  # Install new packages if necessary
  new.packages <-
    list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages))
    install.packages(new.packages)
  
  # Load all packages
  packages_load <-
    lapply(list.of.packages, require, character.only = TRUE)
  
  # Check that all packages were loaded successfully
  if (any(as.numeric(packages_load) == 0)) {
    warning(paste("Package/s: ", paste(list.of.packages[packages_load != TRUE], sep = ", "), "not loaded!"))
  } else {
    print("All packages were successfully loaded.")
  }
  
  # Clear global environment
  rm(list=ls())

# Define directories and import data #### 

  # Directories  
    base.dir <- here(here())
    data.dir <- paste0(base.dir, "/data")
    out.dir <- paste0(base.dir, "/out")
    
  # Import metadata, seqtabs, and tax objects
    metadata <- read.csv(paste0(data.dir, "/metadata_MD.csv"), row.names = 1)
    
    b.seqtab <- as.matrix(readRDS(paste0(data.dir, "/16S_seqtab_nochim.rds")))
    f.seqtab <- as.matrix(readRDS(paste0(data.dir, "/forwards_seqtab_250bp.rds")))
    
    b.tax <- as.matrix(readRDS(paste0(data.dir, "/16S_tax.rds")))
    f.tax <- as.matrix(readRDS(paste0(data.dir, "/forwards_taxa_250bp.rds")))


# Bacterial phyloseq object #### 

  # Remove controls and outliers from b.seqtab
    b.seqtab <- b.seqtab[7:382,]
    b.seqtab <- b.seqtab[!(row.names(b.seqtab) %in% "45D8_S134"),]
  
  # Remove controls and outliers from b.metadata
    b.metadata <- metadata[row.names(metadata) %in% row.names(b.seqtab),]
  
  # Create the phyloseq object 
    b.ps <- phyloseq(otu_table(b.seqtab, taxa_are_rows = FALSE),
                     sample_data(b.metadata), tax_table(b.tax))
    
  # Filter anything not from Kingdom Bacteria
    b.ps <- subset_taxa(b.ps, Kingdom == "Bacteria")
    
  # Save the bacterial phyloseq object
    saveRDS(b.ps, paste0(base.dir, "/data/bac_ps.RDS"))


# Fungal phyloseq object #### 
  # Remove controls from f.seqtab
    f.seqtab <- f.seqtab[1:370,]
  
  # Remove controls from f.metadata
    f.metadata <- metadata[row.names(metadata) %in% row.names(f.seqtab),]
  
  # Create the phyloseq object
    f.ps <- phyloseq(otu_table(f.seqtab, taxa_are_rows = FALSE),
                     sample_data(f.metadata), tax_table(f.tax))
    
  # Filter out anything not from Kingdom Fungi
    f.ps <- subset_taxa(f.ps, Kingdom == "Fungi")

  # Save the fungal phyloseq object
    saveRDS(f.ps, paste0(base.dir, "/data/fun_ps.RDS"))
    