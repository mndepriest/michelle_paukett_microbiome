# Notes/resources #### 

# Some of the phyla that show up in phyla object are only present above
# 500 reads in a few samples. We can't really run analysis on that. 

# Set up environment, import data, create phyloseq objects ####
    source("C:/Users/miran/OneDrive - The Pennsylvania State University/Esker lab/Michelle_ITSand16S/source/01_start_v2.R")
    

# List relevant variables #### 
variables <- colnames(b.metadata)

# We want to remove Plant_no from the modelling because the results are not
# meaningful due to a 1 vs. 1 comparison. 
variables <- variables[!variables == 'Plant_no']


# Prepare the phyloseq object ####
       
    # If you want to work on bacterial data:
      ps <- b.ps

      # Make extra sure that there are no eukaryotes or metazoans
      ps <- subset_taxa(ps, Kingdom == "Bacteria")
      
    # # If you want to work on fungal data:
    # ps <- f.ps
    # 
    # ps <- subset_taxa(ps, Kingdom == "Fungi")

# If you want to subset: #### 
        
    # # Subset by variables
    #   ps1 <- subset_samples(ps, Experiment == "Greenhouse")
    # 
    #   # If you are working within a single group in a variable (i.e., "Field"),
    #   # remove that variable from the variables list
    #   variables <- variables[ !variables == 'Experiment']
      
      
    # # Subset by taxonomic group
      
        # Get most abundant phyla
          phyla <- read.csv(paste0(base.dir, "/out/bacterial_phylum_counts.csv"),
                            row.names = 1)
          phyla$unique.df0.Phylum.

          ps1 <- subset_taxa(ps, Phylum == "Latescibacteria")
      
    # # If you DON'T want to subset:
    #   ps1 <- ps

        
# Set seed, rarefy, and calculate Shannon alpha diversity metric
       
    # Prepare the data   
      # Set a random seed for reproducible results
      set.seed(123456)
      
      # Randomly take 500 reads (with replacement) from samples to correct for
      # uneven sampling depth
      rarefied <- rarefy_even_depth(ps1, sample.size = 500, replace = TRUE)
      
    # Calculate Shannon diversity indices for each sample
      Shannon <- as.data.frame(vegan::diversity((otu_table(rarefied, taxa_are_rows = FALSE)),
                                                  index = "shannon"))
      colnames(Shannon) <- c("Shannon")
      
    
    # Pair Shannon values with metadata
      meta <- data.frame(sample_data(ps1))
      meta <- meta[(rownames(meta) %in% rownames(Shannon)),]
      meta <- cbind(rownames(meta), meta, Shannon)
      
    # Run ANOVA 
      
      # Create an empty dataframe to save ANOVA results to
      ANOVA <- data.frame(matrix(ncol = 8, nrow = 0))
      
      # Run ANOVA, looking for significant differences in alpha diversity by
      # each variable. 
        for (j in variables){
          meta.sub <- dplyr::select(meta, as.character(j), Shannon)
          colnames(meta.sub) <- c("Value", "Shannon")
          result <- aov(Shannon ~ Value, data = meta.sub)
          result.summ <- summary(result)
          result1 <- data.frame(result.summ[[1]])
          part1 <- result1[1,]
          part2 <- result1[2,]
          result2 <- cbind(part1, part2)
          result2 <- result2[,1:8]
          result2$variable <- as.character(j)
          rownames(result2) <- as.character(j)
          ANOVA <- rbind(ANOVA, result2)
        }
      colnames(ANOVA) <- c("Df", "Sum sq", "Mean sq", "F", "p", "Residual df", 
                           "Residual sum sq", "Residual Mean Sq", "variable")
      
      # Adjust the p-value to correct for multiple hypothesis testing
      ANOVA$adj.p <- p.adjust(ANOVA$p, method = "fdr", n = nrow(ANOVA))
      
      # Order the ANOVA object by smallest adjust p value
      ANOVA <- arrange(ANOVA, by_group = adj.p)
      
    # Run Tukey post-hoc 
      
      # Collect variables with significant adjusted 
      ANOVA.sig <- ANOVA[(ANOVA$adj.p < 0.05),]
      vars <- rownames(ANOVA.sig)
      
      Tukey.df <- data.frame(matrix(ncol = 7, nrow = 0))
      colnames(Tukey.df) <- c("difference", "lower", "upper", "p", "comparison", "test", "adj.p")
      
      # Tukey post-hoc for variables with significant ANOVA results
        for (k in vars){
          meta.sub <- dplyr::select(meta, as.character(k), Shannon)
          colnames(meta.sub) <- c("Value", "Shannon")
          meta.sub$Value <- as.character(meta.sub$Value)
          result <- aov(Shannon ~ Value, data = meta.sub)  
          tukey.result <- TukeyHSD(result)
          tukey <- data.frame(tukey.result[["Value"]])
          tukey$comparison <- rownames(tukey)
          tukey$test <- as.character(k)
          colnames(tukey) <- c("difference", "lower", "upper", "p", "comparison", "test")
          tukey$adj.p <- p.adjust(tukey$p, method = "fdr", n = nrow(tukey))
          Tukey.df <- rbind(Tukey.df, tukey)
        }
      
      # Order Tukey.df by significance
      Tukey.df <- arrange(Tukey.df, by_group = adj.p)
      
# Save results to an appropriately-named Excel sheet ####
      
      # Direct Excel to save the results to different sheets within the same
      # Excel document
      sheets <- list("Reads and Shannon values" = meta, "ANOVA" = ANOVA, "Tukey" = Tukey.df)
      openxlsx::write.xlsx(sheets, file = paste0(base.dir, "/out/Latescibacteria.xlsx"))
      

      
