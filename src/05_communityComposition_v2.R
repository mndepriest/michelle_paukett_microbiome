# Notes/references ####

# Used this to make palette: https://mokole.com/palette.html

# Set up the environment with start_v2.R
  source("G:/My Drive/Work/Esker lab/Michelle_ITSand16S/source/01_start_v2.R")

# 16 visually distinct colors
palette <- c('#2f4f4f', '#8b4513', '#008000', '#00008b', '#ff0000', '#ffd700', 
             '#00ff00', '#00ffff', '#ff00ff', '#6495ed', '#f5deb3', '#ff69b4',
             'yellow', 'lightblue', 'red3', "grey1")

# This means you can plot as many as sixteen different groups, including 'Other/Unclassified'

# Visualize community composition patterns to identify possible patterns ####

      # Use for bacteria
      # ps <- b.ps

      # Make sure to only get bacteria
      # ps <- subset_taxa(ps, Kingdom == "Bacteria")

      # Use for fungi
      ps <- f.ps
      ps <- subset_taxa(ps, Kingdom == "Fungi")
      
      # And then change any instances of 'Bacterial' or 'Fungi' in the script 
      # with Ctrl + F
      

    # And set the taxonomic level you want to work at
      ps <- tax_glom(ps, taxrank = "Phylum")
      

  # By Experiment and Treatment #### 
      
      # Create a copy of ps
      ps1 <- ps
      
      # Create a new column in sample_data(ps1) to show Experiment and Treatment
      sample_data(ps1)$group <- paste0(sample_data(ps1)$Experiment, "_", 
                                      sample_data(ps1)$Treatment)
      
      # Group ps1 by groups
      ps.merged <- merge_samples(ps1, group = sample_data(ps1)$group, fun = mean)
      # You'll get feedback that says 'NAs introduced by coercion'. That's okay.
      # We will fill the relevant variables back in.

        
      # Transform counts to 100% relative abundance
      ps.rel <- transform_sample_counts(ps.merged, function(x) x/sum(x))
        
      # Turn the ps object into a dataframe
      df.rel <- psmelt(ps.rel)
      
      # Remove any possible NaNs (only for fungi, but it doesn't hurt the bacteria)
      # Remove rows with NaN in them. These are samples that were too poor to analyze.
      df.rel <- subset(df.rel, df.rel$Abundance != 'NaN')
        
      # If counts are below 1%, turn the Phylum label into "Other/Unclassified"
      for (i in 1:nrow(df.rel)){
        i <- as.numeric(i)
      
      # If you don't have enough colors to finish your plot, adjust 0.01 
      if(df.rel$Abundance[i] < 0.02){
         df.rel$Phylum[i] <- "Other/Unclassified"
          }
        }
        
      # Remove the columns you're going to replace 
      df.rel <- dplyr::select(df.rel, -Experiment, -Treatment)
        
      # Separate groups back into Experiment and Treatment
      df.rel <- separate_wider_delim(df.rel, cols = Sample,
                                       names = c("Experiment", "Treatment"),
                                       delim = "_", 
                                       too_few = "align_start")
        
    # Visualize ####
        # Grouping by Experiment
        ggplot(df.rel, aes(x = Treatment, y = Abundance, fill = Phylum)) +
          geom_bar(stat="identity") +
          scale_fill_manual(values = palette) +
          facet_wrap(~Experiment, scales = 'free_x', nrow = 1) +
          labs(title = "Fungal relative abundance by Experiment and Treatment",
               subtitle = "Phyla above 2% relative abundance shown \nTreatment:",
               ylab = "Abundance", xlab = "Experiment")
        # You can view this in different ways by changing x = Treatment and
        # facet_wrap(~Experiment)

        
  # By Experiment, Treatment, and Days 1, 7, and 14 ####
        
        # Create a copy of ps
        ps1 <- ps
        
        # Remove Days that aren't 1, 7, and 14
        ps.day1 <- subset_samples(ps1, Day == 1)
        ps.day7 <- subset_samples(ps1, Day == 7)
        ps.day14 <- subset_samples(ps1, Day == 14)
        
        # Recombine Days 1, 7, and 14
        ps1 <- merge_phyloseq(ps.day1, ps.day7, ps.day14)
        
        # Create a new column in sample_data(ps1) to show Experiment, Treatment, and Day
        sample_data(ps1)$group <- paste(sample_data(ps1)$Experiment, 
                                         sample_data(ps1)$Treatment,
                                        sample_data(ps1)$Day,
                                        sep = "_")
        
        # Group ps1 by groups
        ps.merged <- merge_samples(ps1, group = sample_data(ps1)$group, fun = mean)
        
        
        # Transform counts to 100% relative abundance
        ps.rel <- transform_sample_counts(ps.merged, function(x) x/sum(x))
        
        # Turn the ps object into a dataframe
        df.rel <- psmelt(ps.rel)
        
        # Remove rows with NaN in them. These are samples that were too poor to analyze.
        df.rel <- subset(df.rel, df.rel$Abundance != 'NaN')
        
        # If counts are below 2%, turn the Phylum label into "Other/Unclassified"
        for (i in 1:nrow(df.rel)){
          i <- as.numeric(i)
          
          # If you don't have enough colors to finish your plot, adjust 0.02 
          if(df.rel$Abundance[i] < 0.02){
            df.rel$Phylum[i] <- "Other/Unclassified"
          }
        }
        
        # Remove the columns you're going to replace 
        df.rel <- dplyr::select(df.rel, -Experiment, -Treatment, -Day)
        
        # Separate groups back into Experiment and Treatment
        df.rel <- separate_wider_delim(df.rel, cols = Sample,
                                       names = c("Experiment", "Treatment", "Day"),
                                       delim = "_", 
                                       too_few = "align_start")
        
        # Turn Day into a factor
        df.rel$Day <- factor(df.rel$Day, levels = c(1, 7, 14))
        
      
    # Visualize ####
        # These show the same data in several different conformations
        ggplot(df.rel, aes(x = Treatment, y = Abundance, fill = Phylum)) +
          geom_bar(stat="identity") +
          scale_fill_manual(values = palette) +
          facet_wrap(~Experiment*Day, scales = 'free_x', nrow = 3, ncol = 3) +
          labs(title = "Fungal relative abundance by Experiment, Treatment, and Day",
               subtitle = "Phyla above 2% relative abundance shown \nExperiment:",
               ylab = "Abundance", xlab = "Treatment")
        # You can change which variables are shown where by changing
        # x = Treatment and facet_wrap(~Experiment*Day). Use the following
        # line if the x-axis labels are close to each other.
          # + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
  # By root disease, treatment, and experiment #### 
        
        # Make a copy of ps
        ps1 <- ps
        
        #  Group ps1 by root_dis, Treatment, and Experiment
        sample_data(ps1)$group <- paste(sample_data(ps1)$root_dis, 
                                        sample_data(ps1)$Experiment,
                                        sample_data(ps1)$Treatment,
                                        sep = "_")
        
        # Merge by group
        ps.merged <- merge_samples(ps1, group = sample_data(ps1)$group, fun = mean)
        
        # Transform into relative abundance
        ps.rel <- transform_sample_counts(ps.merged, function(x) x/sum(x))
        
        # Turn into a dataframe
        df.rel <- psmelt(ps.rel)
        
        # Remove rows with NaN in them. These are samples that were too poor to analyze.
        df.rel <- subset(df.rel, df.rel$Abundance != 'NaN')
        
        # If counts are below 2%, turn the Phylum label in 'Other/Unclassified'
        for (i in 1:nrow(df.rel)){
          i <- as.numeric(i)
          
          # If you don't have enough colors to finish your plot, adjust 0.02 
          if(df.rel$Abundance[i] < 0.02){
            df.rel$Phylum[i] <- "Other/Unclassified"
          }
        }
        
        # Remove the columns you're going to replace 
        df.rel <- dplyr::select(df.rel,  -root_dis, -Experiment, -Treatment)
        
        # Separate groups back into Experiment and Treatment
        df.rel <- separate_wider_delim(df.rel, cols = Sample,
                                       names = c("root_dis", "Experiment", "Treatment"),
                                       delim = "_", 
                                       too_few = "align_start")
        
        # Order root_dis
        df.rel$root_dis <- factor(df.rel$root_dis, levels = c("Low", "High"))
        
      
      # Visualize #### 
        # These show the same data in several different conformations
        ggplot(df.rel, aes(x = root_dis, y = Abundance, fill = Phylum)) +
          geom_bar(stat="identity") +
          scale_fill_manual(values = palette) +
          facet_wrap(~Experiment*Treatment, scales = 'free_x') +
          labs(title = "Fungal relative abundance by root_dis, Experiment, and Treatment",
               subtitle = "Phyla above 2% relative abundance shown \nExperiment:",
               ylab = "Abundance", xlab = "root_dis")
        
        # A problem with the root_dis and shoot_dis categorical ratings is that
        # since they're binary, you're not guaranteed to have high and low ratings. 
        # An alternative would be to base low and high ratings on the total
        # min and maximum ratings or to include more classifications (ex. low,
        # mid-low, mid, mid-high, high). 
        
    # Root_dis for each Day 14 sample #### 
        # Make a copy of the ps object
        ps1 <- ps
        
        # Remove samples that are not from Day 14
        ps1 <- subset_samples(ps1, Day == 14)
        
        # Aggregate at the phylum level
        ps1 <- 
        
        # Transform the counts 
        ps.rel <-  transform_sample_counts(ps1, function(x) x/sum(x))
        
        # Turn ps.rel into a dataframe
        df.rel <- psmelt(ps.rel)
        
        # Remove rows with NaN in them. These are samples that were too poor to analyze.
        df.rel <- subset(df.rel, df.rel$Abundance != 'NaN')
        
        # If counts are below 2%, turn the Phylum label in 'Other/Unclassified'
        for (i in 1:nrow(df.rel)){
          i <- as.numeric(i)
          
          # If you don't have enough colors to finish your plot, adjust 0.02 
          if(df.rel$Abundance[i] < 0.02){
            df.rel$Phylum[i] <- "Other/Unclassified"
          }
        }
        
        # Set factor levels
        df.rel$root_dis <- factor(df.rel$root_dis, levels = c("Low", "High"))
        
        # Visualize #### 
        ggplot(df.rel, aes(x = Sample, y = Abundance, fill = Phylum)) +
          geom_bar(stat="identity") +
          scale_fill_manual(values = palette) +
          facet_wrap(~root_dis*Experiment, scales = 'free_x') +
          labs(title = "Fungal relative abundance at day 14 by root_dis, Experiment, and Treatment",
               subtitle = "Phyla above 2% relative abundance shown \nExperiment:",
               ylab = "Abundance", xlab = "Treatment") + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        # It's possible that community evenness is driving alpha diversity more so 
        # than community richness in healthy/low disease plants. Maybe we should
        # try calculating the Chao scores for each of these samples and running
        # ANOVA on that. 
      
        
  # By shoot disease, treatment, and experiment #### 
        
        # Make a copy of ps
        ps1 <- ps
        
        #  Group ps1 by shoot_dis, Treatment, and Experiment
        sample_data(ps1)$group <- paste(sample_data(ps1)$shoot_dis, 
                                        sample_data(ps1)$Experiment,
                                        sample_data(ps1)$Treatment,
                                        sep = "_")
        
        # Merge by group
        ps.merged <- merge_samples(ps1, group = sample_data(ps1)$group, fun = mean)
        
        # Transform into relative abundance
        ps.rel <- transform_sample_counts(ps.merged, function(x) x/sum(x))
        
        # Turn into a dataframe
        df.rel <- psmelt(ps.rel)
        
        # Remove rows with NaN in them. These are samples that were too poor to analyze.
        df.rel <- subset(df.rel, df.rel$Abundance != 'NaN')
        
        # If counts are below 2%, turn the Phylum label in 'Other/Unclassified'
        for (i in 1:nrow(df.rel)){
          i <- as.numeric(i)
          
          # If you don't have enough colors to finish your plot, adjust 0.02 
          if(df.rel$Abundance[i] < 0.02){
            df.rel$Phylum[i] <- "Other/Unclassified"
          }
        }
        
        # Remove the columns you're going to replace 
        df.rel <- dplyr::select(df.rel,  -shoot_dis, -Experiment, -Treatment)
        
        # Separate groups back into Experiment and Treatment
        df.rel <- separate_wider_delim(df.rel, cols = Sample,
                                       names = c("shoot_dis", "Experiment", "Treatment"),
                                       delim = "_", 
                                       too_few = "align_start")
        
        # Order shoot_dis
        df.rel$shoot_dis <- factor(df.rel$shoot_dis, levels = c("Low", "High"))
        
        
      # Visualize #### 
        # These show the same data in several different conformations
        ggplot(df.rel, aes(x = shoot_dis, y = Abundance, fill = Phylum)) +
          geom_bar(stat="identity") +
          scale_fill_manual(values = palette) +
          facet_wrap(~Experiment*Treatment, scales = 'free_x') +
          labs(title = "Fungal relative abundance by shoot_dis, Experiment, and Treatment",
               subtitle = "Phyla above 2% relative abundance shown \nExperiment:",
               ylab = "Abundance", xlab = "shoot_dis")

        
    # shoot_dis for each Day 14 sample #### 
        # Make a copy of the ps object
        ps1 <- ps
        
        # Remove samples that are not from Day 14
        ps1 <- subset_samples(ps1, Day == 14)
        
        # Aggregate at the phylum level
        ps1 <- 
          
          # Transform the counts 
          ps.rel <-  transform_sample_counts(ps1, function(x) x/sum(x))
        
        # Turn ps.rel into a dataframe
        df.rel <- psmelt(ps.rel)
        
        # Remove rows with NaN in them. These are samples that were too poor to analyze.
        df.rel <- subset(df.rel, df.rel$Abundance != 'NaN')
        
        # If counts are below 2%, turn the Phylum label in 'Other/Unclassified'
        for (i in 1:nrow(df.rel)){
          i <- as.numeric(i)
          
          # If you don't have enough colors to finish your plot, adjust 0.02 
          if(df.rel$Abundance[i] < 0.02){
            df.rel$Phylum[i] <- "Other/Unclassified"
          }
        }
        
        # Set factor levels
        df.rel$shoot_dis <- factor(df.rel$shoot_dis, levels = c("Low", "High"))
        
      # Visualize #### 
        ggplot(df.rel, aes(x = Sample, y = Abundance, fill = Phylum)) +
          geom_bar(stat="identity") +
          scale_fill_manual(values = palette) +
          facet_wrap(~shoot_dis*Experiment, scales = 'free_x') +
          labs(title = "Fungal relative abundance at day 14 by shoot_dis, Experiment, and Treatment",
               subtitle = "Phyla above 2% relative abundance shown \nExperiment:",
               ylab = "Abundance", xlab = "Treatment") + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
