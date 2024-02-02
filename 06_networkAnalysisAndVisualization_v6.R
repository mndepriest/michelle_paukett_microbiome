# Notes/resources #### 
  
  # This is the current guide to SpiecEasi:
  # https://github.com/zdk123/SpiecEasi#spieceasi 
  
  # For the strongest relationships - I've ordered the taxa by Phylum so that
  # you can potentially draw brackets on the graphs to denote which Phyla each
  # taxa belongs to. 


  
# Packages #### 
  
  # Define all necessary packages 
  list.of.packages <- c("devtools",
                        "BiocManager",
                        "Biostrings",
                        "phyloseq",
                        "dplyr",
                        "ggplot2",
                        "Rcpp",
                        "here",
                        "R.utils",
                        "SpiecEasi",
                        "conflicted",
                        "Matrix",
                        "corrr", 
                        "corrplot")
  
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
    warning(paste("Package/s: ", paste(list.of.packages[packages_load != TRUE], 
                                       sep = ", "), "not loaded!"))
  } else {
    print("All packages were successfully loaded.")
  }
  
  # Command conflicts 
  conflicts_prefer(SpiecEasi::tril)
  conflicts_prefer(SpiecEasi::triu)
  conflicts_prefer(igraph::make_graph)
  
  # Clear global environment
  rm(list.of.packages, new.packages, packages_load)

# Import data and get files ####

  # Define the base directory
    base.dir <- here(here())
  
  # Get files
    
    # Read taxonomy file, filter, make into matrix
    tax <- data.frame(readRDS(paste0(base.dir, "/data/16S_tax.rds")))
    tax <- tax[which(tax$Kingdom == "Bacteria"),]
    
    # Remove any Phylum name containing 'Unknown'
    tax$Phylum[grepl("nknown", tax$Phylum)] <- NA
    tax <- as.matrix(tax)
    
    # Replace any hyphens (-) with underscores (_)
    tax <- apply(tax, 2, function(x) gsub("-", "_", x))
    
    # Read metadata file, remove outlier
    meta <- read.csv(paste0(base.dir, "/data/metadata_MD.csv"), row.names = 1)
    meta <- meta[!(rownames(meta) %in% "45D8_S134"),]
    
    # Read seqtab_nochim file, remove controls & outliers, filter by read number
    seqtab <- data.frame(readRDS(paste0(base.dir, "/data/16S_seqtab_nochim.rds")))
    seqtab <- seqtab[(rownames(seqtab) %in% rownames(meta)),]
    seqtab <- dplyr::filter(seqtab, rowSums(seqtab) > 500)
    

# Create ps, get objects needed for model #### 
    
  # # Create, filter, glom, and pull seqtab and tax back out
  #   ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
  #                  sample_data(meta), tax_table(tax)) %>%
  #     subset_taxa(Kingdom == "Bacteria") %>%
  #     tax_glom(taxrank = "Phylum", NArm = TRUE)
  #   saveRDS(ps, paste0(base.dir, "/out/network_analysis/ps_Phylum.RDS"))

    # If ps has already been created:
    ps <- readRDS(paste0(base.dir, "/out/network_analysis/ps_Phylum.RDS"))
  
    # Get seqtab and tax out of ps object
    seqtab1 <- data.frame(otu_table(ps))
    tax1 <- data.frame(tax_table(ps))
    
    # Format taxa with identical names
      # Look for duplicates and see how often & where they occur
      n_occur <- data.frame(table(tax1$Phylum))
      n_occur <- n_occur[n_occur$Freq > 1,]
      # There are duplicate Phylum names. These are Families with the same name
      # but from different Phylums
      
      # If there are duplicates, append the duplicate names with their Phyla 
      # (if possible). If not, append their duplicate names with the first
      # unique taxonomic level.
      if (nrow(n_occur) > 0){
        for (i in as.numeric(1:nrow(n_occur))){
          for (j in as.numeric(1:nrow(tax1))){
            if (tax1$Phylum[j] == n_occur$Var[i]){
              
              ord <- substr(tax1$Class[j], 1, 3)
              tax1$Phylum[j] <- paste0(ord, "-", tax1$Phylum[j])
            }
          }
        }
      }
    
   
# Run correlation matrix model ####

  # Run sparcc using SpiecEasi package and save it
    sp <- sparcc(seqtab1)
    saveRDS(sp, paste0(base.dir, "/out/network_analysis/sparcc_Phylum.RDS"))
    
    # If the sparcc model has already been created, load it in
    sp <- readRDS(paste0(base.dir, "/out/network_analysis/sparcc_Phylum.RDS"))

    # Get correlation matrix and save for each taxonomic level
    cor <- sp[["Cor"]]

    # Set row and column names from seqs to Phylum name
    rownames(cor) <- tax1$Phylum
    colnames(cor) <- tax1$Phylum

    # Remove self:self relationships
    diag(cor) <- 0

    # Save the correlation matrix for each taxonomic level
    write.csv(cor, paste0(base.dir, "/out/network_analysis/sparcc_correlation_Phylum.csv"))

      
# Visualizing: Determine which relationships to show ####   
  
    
# Keep taxa with the strongest relationships ####
    
    # Create a copy of cor that is absolute values
    abs.cor <- abs(cor)

    # Get max relationship strength for each taxon and look at distribution
    max.rel <- data.frame(apply(abs.cor, 2, max))
    colnames(max.rel) <- "max"
    rownames(max.rel) <- tax1$Phylum

    # Filtering for the strongest 50 relationships
    max.rel <- arrange(max.rel, desc(max))
    
    # # Choose the strongest 50 relationships. Not available for Phyla
    # max.rel <- dplyr::filter(max.rel, max.rel$max >= max.rel$max[50])

    # Removing taxa with weaker relationships from correlation matrix
    cor1 <- cor[(rownames(cor) %in% rownames(max.rel)),]
    cor1 <- cor1[,(colnames(cor1) %in% rownames(max.rel))]


    # # Display taxa in order. Comment out in phyla
    #   max.rel <- tax1[(tax1$Phylum %in% rownames(cor1)),]
    #   max.rel <- arrange(max.rel, Phylum)
    #   tax.Phylum <- as.vector(max.rel$Phylum)
    #   
    # # Arrange columns and rows according to higher taxonomic level
    #   cor1 <- cor1[, tax.Phylum]
    #   cor1 <- cor1[tax.Phylum, ]
    #   
    # # Make and save a .csv showing which lower taxa belong to which higher taxa
    #   # so that you can draw brackets grouping them. Not available for phy.lum
    #   max.rel <- dplyr::select(max.rel, Phylum, Phylum)
    #   # write.csv(max.rel, paste0(base.dir, 
    #   #                           "/out/network_analysis/correlation_Phylum_by_Phylum.csv"), 
    #   #           row.names = FALSE)
    
  # Visualize and save #### 
    
    # # Create relationship strength heatmap and save it 
    # png(paste0(base.dir, 
    #            "/out/network_analysis/sparcc_correlation_Phylum_strongest_relationships.png"),
    #     width = 1600, height = 1600, units = 'px')
    corrplot(cor1, type = 'lower', tl.cex = 0.5, tl.col = "black")
    # adjust tl.cex to change font size
    dev.off()
  
    
# Alternative: Keep the most abundant taxa #### 
    
    # Arrange seqtab by number of reads per taxa
    seq.abund <- data.frame(t(seqtab1))
    seq.abund <- arrange(seq.abund, desc(rowSums(seq.abund)))
    
    # Keep only the top 50 most abundant taxa
    seq50 <- seq.abund[1:50,]
    tax50 <- tax1[rownames(tax1) %in% rownames(seq50),]
    
    # Remove less abundant taxa from the correlation matrix
    cor1 <- cor[(rownames(cor) %in% tax50$Phylum),]
    cor1 <- cor1[,(colnames(cor1) %in% tax50$Phylum)]
    
    # Remember to change the image name! 
   
  # Visualize and save #### 
    
  # # As heatmap 
  #   png(paste0(base.dir, "/out/network_analysis/sparcc_correlation_Phylum_50_most_abundant.png"),
  #       width = 1600, height = 1600, units = 'px')
    corrplot(cor1, type = 'lower', tl.cex = 0.5, tl.col = "black")
    # dev.off()
  
