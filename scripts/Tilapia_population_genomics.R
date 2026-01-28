
#Tilapia Population Genomics Analysis
#SNP-based population structure, PCA, DAPC and FST analyses

# Load libraries
library(adegenet)
library(ape)
library(StAMPP)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(jcolors)
library(gaston)
library(reshape2)
library(RColorBrewer)
library(vcfR)
library(cowplot)
library(gdata)
library(vegan)
library(geosphere)
library(dplyr)
library(readxl)
library(grid)


# ------------------------------
# Set working directory

setwd("/home/redempta/work_2025/tilapia/new_analysis")

# Load data

til_vcf <- read.vcfR("populations.snps.vcf")
populations_til <- read.table("populations_info_lines_.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Preprocess populations data

colnames(populations_til)[colnames(populations_til) == "Lines"] <- "Populations"
populations_til$Populations <- recode(populations_til$Populations,
                                      exotic_Nile_tilapia = "introduced_Nile_tilapia",
                                      local_Nile_tilapia = "native_Nile_tilapia",
                                      local_Rufiji_tilapia = "native_Rufiji_tilapia",
                                      local_Rufiji = "local_Rufiji_tilapia")

# Count individuals per population
df_counts <- populations_til %>%
  group_by(Pop_Id) %>%
  summarise(n_individuals = n())
write.csv(df_counts, "population_counts.csv", row.names = FALSE)


# Convert VCF to genlight

til <- vcfR2genlight(til_vcf)
pop(til) <- populations_til$Pop_Id

# rufiji tilapia
pop(til) <- populations_til$Species
til_rufiji <- til[til$pop=="Rufiji_tilapia",]
pop_rufiji <- populations_til[populations_til$Species=="Rufiji_tilapia",]

# niloticus
pop(til) <- populations_til$Species
til_nile <- til[til$pop=="Nile_tilapia",]
pop(til_nile) <- populations_til[populations_til$Species=="Nile_tilapia",]$Origin


# allele freq 
myFreq <- glMean(til)
myFreq
myFreq <- c(myFreq, 1-myFreq)
windows()
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)



#PCA species and origin
pca <- glPca(til)
#6
var_exp <- pca$eig / sum(pca$eig)
pc1_var <- round(var_exp[1]*100, 1)
pc2_var <- round(var_exp[2]*100, 1)
pc3_var <- round(var_exp[3]*100, 1)

pca_df <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  PC3 = pca$scores[,3],
  Populations = populations_til$Populations,
  Pop_Id = populations_til$Pop_Id,
  Origin = populations_til$Origin
)

group_colors <- c(
  "introduced_Nile_tilapia" = "#E69F00", 
  "native_Nile_tilapia" = "#0072B2", 
  "native_Rufiji_tilapia" = "#009E73"
)
group_shapes <- c(
  "introduced_Nile_tilapia" = 16, 
  "native_Nile_tilapia" = 17, 
  "native_Rufiji_tilapia" = 15
)

# PCA nile
pcan <- glPca(til_nile)
#5
# Variance explained
var_expn <- pcan$eig / sum(pcan$eig)  
pc1_varn <- round(var_expn[1]*100, 1)
pc2_varn <- round(var_expn[2]*100, 1)
pc3_varn <- round(var_expn[3]*100, 1)

pc1_var; pc2_var; pc3_var

# dataframe for ggplot 
pca_n <- data.frame(
  PC1 = pcan$scores[,1],
  PC2 = pcan$scores[,2],
  PC3 = pcan$scores[,3],
  Populations = Nile_introduced_local$Populations,
  Pop_Id = Nile_introduced_local$Pop_Id,
  Origin = Nile_introduced_local$Origin  
)


# Colors for Origin ----
origin_colors <- c(
  "Netherland" = "#CC79A7",
  "Tanzania" = "#56B4E9",
  "Thailand" = "#D55E00",
  "Uganda" = "#882E72"
)


origin_shapes <- c(
  "Netherland" = 16,  
  "Tanzania"  = 17,   
  "Thailand"  = 15,   
  "Uganda"    = 18    
)


p1<- ggplot(pca_df, aes(x=PC1, y=PC2, color=Populations, shape=Populations)) +
  geom_point(size=3) +
  geom_text(aes(label=Pop_Id), vjust=-0.5, hjust=0.5, size=3) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  scale_color_manual(values=group_colors) +
  scale_shape_manual(values=group_shapes) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
p1

# Panel 2 by Origin 
p2<- ggplot(pca_n, aes(x=PC1, y=PC2, color=Origin, shape=Origin)) +
  geom_point(size=3) +
  #geom_text(aes(label=Pop_Id), vjust=-0.5, hjust=0.5, size=3) +
  xlab(paste0("PC1 (", pc1_varn, "%)")) +
  ylab(paste0("PC2 (", pc2_varn, "%)")) +
  scale_color_manual(values=origin_colors) +
  scale_shape_manual(values=origin_shapes) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
p2

# Add A and B labels 
combined <- plot_grid(p1, p2,
                      labels = c("A", "B"),
                      label_size = 14,  # adjust font size
                      nrow = 1)

ggsave("tilapia_pca_combined.tiff", combined, width=14, height=6, dpi=300)


# PCA rufiji
pcaj <- glPca(til_rufiji)
#5
# Variance explained
var_expj <- pcaj$eig / sum(pcaj$eig)  
pc1_varj <- round(var_expj[1]*100, 1)
pc2_varj <- round(var_expj[2]*100, 1)
pc3_varj <- round(var_expj[3]*100, 1)

pca_j <- data.frame(
  PC1 = pcaj$scores[,1],
  PC2 = pcaj$scores[,2],
  PC3 = pcaj$scores[,3],
  Populations = pop_rufiji$Populations,
  Pop_Id = pop_rufiji$Pop_Id,
  Origin = pop_rufiji$Wild  # <- add this
)

pop_rufiji$Origin <- ifelse(pop_rufiji$Wild == 1, "Wild", "Farmed")

pca_j$Origin <-  ifelse(pca_j$Origin == 1, "Wild", "Farmed")


group_colors <- c(
  "Farmed" = "#E69F00", 
  "Wild" = "#009E73"
)

# ---- Shapes ----
group_shapes <- c(
  "Farmed" = 16, 
  "Wild" = 15
)

#PCA

p3<- ggplot(pca_j, aes(x=PC1, y=PC2, color=Origin, shape=Origin)) +
  geom_point(size=3) +
  geom_text(aes(label=Pop_Id), vjust=-0.5, hjust=0.5, size=3) +
  xlab(paste0("PC1 (", pc1_varj, "%)")) +
  ylab(paste0("PC2 (", pc2_varj, "%)")) +
  scale_color_manual(values=group_colors) +
  scale_shape_manual(values=group_shapes) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
p3

ggsave("rufiji_tilapia.png", p3, width=12, height=6, dpi=300)


# DAPC Analysis

# Impute missing data
mat <- tab(til, NA.method="mean")

# Define groups
grp <- populations_til$Species

# Cross-validation to choose optimal number of PCs
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.75,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)

optimal_pcs <- xval$`Number of PCs Achieving Highest Success`

# Run DAPC with optimal PCs
dapc_result <- dapc(mat, grp, n.pca = optimal_pcs, n.da = length(unique(grp)) - 1)

# Scatterplot
scatter(dapc_result)

#density plot
#DAPC plot
df1 <- dapc_result$ind.coord[,1]
df <- data.frame(df1=df1, species=grp)

top <- ggplot(df, aes(x=df1, fill=species)) +
  geom_density(alpha=0.5) +
  geom_rug(aes(color=species), sides="b") +
  xlab("Discriminant function 1") +
  ylab("Density") +
  theme_classic() +  
  theme(
    panel.border = element_blank(),  
    legend.background = element_blank(),  
    legend.key = element_blank()
  )

top

ggsave("tdapc_species_plot.tiff", top, width=12, height=6, dpi=300)


dapc1 <- dapc(til_nile, til_nile$pop)
100
3
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", xlab="",legend=TRUE,solid=0.5)
compoplot(dapc1, posi="bottomright",
          txt.leg=c("Netherland","Tanzania", "Thailand", "Uganda"), lab="",cleg=0.7,xlab="Tilapias")


## structure like plot
# Convert dapc results
dapc.results <- as.data.frame(dapc1$posterior)
dapc.results$pop <- til_nile$pop
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_cluster","Posterior_membership_probability")

dapc.results <- dapc.results %>% arrange(Original_Pop)

# Define colors for populations
origin_colors <- c(
  "Netherland" = "#CC79A7",
  "Tanzania"   = "#56B4E9",
  "Thailand"   = "#D55E00",
  "Uganda"     = "#882E72"
)

pop_labels <- c(
  "Netherland" = "NL (n=20)",
  "Tanzania"   = "TZ (n=163)",
  "Thailand"   = "TH (n=80)",
  "Uganda"     = "UG (n=56)"
)

# DAPC admixture plot
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_cluster)) +
  geom_bar(stat='identity', width=1) +
  scale_fill_manual(
    values = origin_colors,
    labels = c(
      "Netherland" = "Netherland (introduced)",
      "Tanzania"   = "Tanzania (native)",
      "Thailand"   = "Thailand (introduced)",
      "Uganda"     = "Uganda (introduced)"
    )
  ) +
  facet_grid(~Original_Pop, scales = "free", space="free_x",
             labeller = labeller(Original_Pop = pop_labels)) +
  labs(
    x = "Nile tilapia samples",
    y = "Posterior membership probability",
    fill = "Population"
  ) +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size=14, face="bold"),
    legend.text  = element_text(size=12)
  )

print(p)

ggsave("Nile_tilapia_admix.tiff",p, width=14, height=6, dpi=300)



# Pairwise FST

Results<-read_tsv("Results.txt")
Results[1] <- NULL
Results$Fst
Results

fst_mat <- matrix(0,nrow=27,ncol=27)
fst_mat
upperTriangle(fst_mat,byrow = TRUE ) <- Results$Fst
lowerTriangle(fst_mat) <- NA

fst_result
diag(fst_mat) <- 0
colnames(fst_mat) <- unique(populations_til$Pop_Id)
rownames(fst_mat) <- unique(populations_til$Pop_Id)
populations_til$Pop_Id
fst_mat
nrow(fst_mat)
ncol(fst_mat)
fst_mat_melt <- melt(fst_mat)
colnames(fst_mat_melt) <- c("Population_1","Population_2","Fst")
windows()
fst_mat_melt
fst_map <- ggplot(fst_mat_melt, aes(Population_1,Population_2,fill=Fst)) +
  geom_tile() +
  scale_fill_viridis_c() +
  #scale_fill_gradient(low="blue",high="red") +
  labs(x="Population",color="Fst") +
  theme(plot.title=element_text(size=18,hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=12),
        #axis.title.y=element_text(size=16)) 
        #axis.text.y=element_blank(),
        axis.title.y=element_blank())

fst_map



#labels
fst_map +
  # Group 1
  annotation_custom(
    grob = textGrob(
      "Native Rufiji tilapia",
      rot = 90,
      gp = gpar(fontsize = 10, fontface = "bold")
    ),
    xmin = -3.5, xmax = -1.5,  
    ymin = 0, ymax = 10
  ) +
  # Group 2
  annotation_custom(
    grob = textGrob(
      "Introduced Nile tilapia",
      rot = 90,
      gp = gpar(fontsize = 10, fontface = "bold")
    ),
    xmin = -4.5, xmax = -1.5,
    ymin = 11, ymax = 18
  ) +
  # Group 3
  annotation_custom(
    grob = textGrob(
      "Native Nile tilapia",
      rot = 90,
      gp = gpar(fontsize = 10, fontface = "bold")
    ),
    xmin = -4, xmax = -1.5,
    ymin = 20, ymax = 27
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(10, 10, 10, 80),  # enough left space
    axis.text.y = element_text(size = 11)
  )


ggsave("fst_grounat_newfinal.png", width = 12, height = 6, dpi = 300)


