library("devtools")
library(microbiome)
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

head(prevalence(ps1.com.fam.rel, detection = 1/100, sort = TRUE, count = TRUE))


head(prevalence(ps1.com.fam.rel, detection = 1/100, sort = TRUE))

core.taxa.standard <- core_members(ps1.com.fam.rel, detection = 0, prevalence = 50/100)

pseq.core <- core(ps1.com.fam.rel, detection = 0, prevalence = .5)

pseq.core2 <- aggregate_rare(ps1.com.fam.rel, "Family", detection = 0, prevalence = .5)

core.taxa <- taxa(pseq.core)

core.abundance <- sample_sums(core(ps1.com.fam.rel, detection = .01, prevalence = .95))

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(ps1.com.fam.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

################################################################################
library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
library(RColorBrewer)
p <- plot_core(ps1.com.fam.rel,
               plot.type = "heatmap", 
               #colours = gray,
               colours = rev(brewer.pal(5, "RdBu")),
               prevalences = prevalences, 
               detections = detections, 
              )

p +  theme_bw()


library(viridis)
print(p + scale_fill_viridis())

################################################################################
detections <- seq(from = 50, to = round(max(abundances(ps1))/10, -1), by = 100)

p <- plot_core(ps1, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE) +
  theme(axis.text.x= element_text(size=8, face="italic", hjust=1),
        axis.text.y= element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)
################################################################################

# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(ps1.com.fam.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

