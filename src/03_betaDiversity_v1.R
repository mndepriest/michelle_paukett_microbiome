# Notes/references #### 

# Set up environment, import data, create phyloseq objects
    source("C:/Users/miran/OneDrive - The Pennsylvania State University/Esker lab/Michelle_ITSand16S/source/01_start_v2.R")


# Prepare for ordination and modeling #### 
  # Define the ps you want to work with
    # For bacterial data:
    ps <- b.ps
    
    # Make extra sure to remove anything that is not bacteria
    ps <- subset_taxa(ps, Kingdom == "Bacteria")
    
    # For fungal data:
    # ps <- f.ps
    
  # # Subset the data (optional)
    ps <- subset_samples(ps, Experiment == "Greenhouse")
  #   # Be sure to modify the model below if you do subset. You can't model by 
  #   # Experiment if you're only looking at Greenhouse, or by Treatment if you're
  #   # only looking at the control. 
    
    # Subset by phyla
    # ps <- subset_taxa(ps, Phylum == "Bacteroidetes")
    
  # Rarefy
    # Set random seed
    set.seed(123456)
    
    # Rarefy
    rarefied <- rarefy_even_depth(ps, sample.size = 500, replace = TRUE)
    
  # Create the Bray-Curtis distance matrix
    bray <- phyloseq::distance(rarefied, method = "bray")
    
    # Format a sample data frame
    sampledf <- data.frame(sample_data(rarefied))
    
  # Run the full PERMANOVA. We include Plant_no to estimate the natural variation.
    # We remove Experiment if working within a single Experiment
    adonis2(bray ~ Plant_no*Day*Treatment*shoot_dis*root_dis*soil, data = sampledf)
    
    
  # Is there a way to directly write the PERMANOVA results to a new sheet
    # in the already existing full_bacterial, full_fungal, etc. docs?
    