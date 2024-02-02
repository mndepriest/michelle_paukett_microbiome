
# Set up environment, import data, create phyloseq objects ####
  source("C:/Users/miran/OneDrive - The Pennsylvania State University/Esker lab/Michelle_ITSand16S/source/01_start_v2.R")

# Choose data to work with
  # For bacteria:
  ps <- b.ps

  # Make sure to only get bacteria 
  ps <- subset_taxa(ps, Kingdom == "Bacteria")

# # For fungi:
  # ps <- f.ps
  # 
  # # Make sure to only get fungi
  # ps <- subset_taxa(ps, Kingdom == "Fungi")

# Inspect and manipulate the phyloseq data frame (optional)
  ps0 <- tax_glom(ps, taxrank = "Phylum") 
  df0 <- psmelt(ps0)

# Get the unique phyla 
  phy <- data.frame(unique(df0$Phylum))

  phy.counts <- c()

  for (i in 1:nrow(phy)){
    i <- as.numeric(i)
    result <- df0[which(df0$Phylum == phy[i,]),]
    result1 <- sum(result$Abundance)
    
    phy.counts <- append(phy.counts, result1)
  }
  
  phy$abundance <- phy.counts
  phy <- arrange(phy, by_group = desc(abundance))
  phy.top20 <- phy[1:20,]
  phyla <- phy.top20$unique.df0.Phylum.
  
  write.csv(phy, paste0(base.dir, "/out/bacterial_phylum_counts.csv"))