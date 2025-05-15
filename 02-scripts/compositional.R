library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling  

ps1 <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
print(ps1)
ps1.com <- ps1
taxa_names(ps1.com) <- paste0("ASV_", rownames(tax_table(ps1.com)))
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.
taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"
guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))

ps1.com@phy_tree <- NULL

# Transform to relative abundance first
ps1.com <- microbiome::transform(ps1.com, "compositional")


# Filter out taxa with mean relative abundance < 0.1% (0.001)
ps1.com <- filter_taxa(ps1.com, function(x) mean(x) > 0.001, prune = TRUE)


ps_fam <- aggregate_taxa(ps1.com, level = "Family")
top_families <- names(sort(taxa_sums(ps_fam), decreasing = TRUE))[1:10]
ps1.com.fam <- prune_taxa(top_families, ps_fam)

ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")

plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
                                             sample.sort = "scientific_name",
                                             x.label = "env_material") 

plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Family", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))
print(plot.composition.relAbun)




# Assume ps1.com@phy_tree is NULL and you've done the aggregation
ps_fam <- aggregate_taxa(ps1.com, level = "Family")
top_families <- names(sort(taxa_sums(ps_fam), decreasing = TRUE))[1:10]
ps1.com.fam <- prune_taxa(top_families, ps_fam)

# Transform to relative abundance
ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")
custom_order <- c("K1",  "K2",  "K3", "K4", "K5",   "K6",  "K7",  "K8",   "K9", "K10",
                  "DM1", "DM2",  "DM3", "DM4", "DM5",  "DM6", "DM7", "DM8", "DM9", "DM10", 
                  "DM11", "DM12", "DM13", "DM14", "DM15", "DM16", "DM17", "DM18", "DM19", "DM20", "DM21",
                  "PDM1", "PDM2", "PDM4", "PDM5", "PDM6", "PDM7", "PDM8", "PDM9", "PDM10", 
                  "PDM11", "PDM12", "PDM13", "PDM14", "PDM15", "PDM16", "PDM17", "PDM18")

# Assign ordered factor to sample_information in the ps1.com.fam.rel object
sample_data(ps1.com.fam.rel)$sample_information <- factor(
  sample_data(ps1.com)$sample_information,
  levels = custom_order
)



# Plot
plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
                                             sample.sort = "sample_information",  # use this for x-axis ordering
                                             x.label = "sample_information")      # label x-axis accordingly

# Add aesthetics
plot.composition.relAbun <- plot.composition.relAbun +
  theme(legend.position = "bottom") +
  scale_fill_brewer("Family", palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +

  labs(y = "Relative abundance") +  # <- sets y-axis label
  guide_italics +
  theme(legend.title = element_text(size = 18))

# Print
print(plot.composition.relAbun)

#ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/compositional_plots_family.svg", height = 10, width = 10)
#ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/compositional_plots_family.png", height = 10, width = 10)
##################### heatmap 

data.com <- plot.composition.relAbun$data
# Add sample_information to data.com
#data.com$sample_information <- sample_data(ps1)$sample_information[as.character(data.com$Sample)]



data.com$sample_information <- factor(data.com$xlabel, levels = custom_order)

colnames(data.com)

taxon_means <- data.com %>%
  group_by(Tax) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE))

# Define threshold (e.g., keep taxa with mean rel. abundance > 0.5%)
threshold <- 0.001

# Filter taxa
keep_taxa <- taxon_means %>% filter(mean_abundance > threshold) %>% pull(Tax)
data.com.filtered <- data.com %>% filter(Tax %in% keep_taxa)
# base plot
p.heat <- ggplot(data.com, aes(x = sample_information, y = Tax)) + geom_tile(aes(fill = Abundance)) 

# Change color
p.heat <- p.heat + scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw() 

# Make bacterial names italics
p.heat <- p.heat + theme(axis.text.y = element_text(colour = 'black', 
                                                    size = 10, 
                                                    face = 'italic')) 
p.heat <- p.heat + ylab("Family") 



# Clean the facet label box
p.heat <- p.heat + theme(legend.key = element_blank(), 
                         strip.background = element_rect(colour="black", fill="white"))

p.heat <- p.heat + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



print(p.heat)

#ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap_relab.svg", height = 10, width = 10)
#ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap_relab.png", height = 10, width = 10)
#############################################################################


