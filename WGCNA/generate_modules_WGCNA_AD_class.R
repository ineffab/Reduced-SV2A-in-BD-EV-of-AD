library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph)
library(RColorBrewer)
library(readxl)
library(org.Hs.eg.db)
library(HGNChelper)
library(devtools)
library(GO.db)
library(AnnotationDbi)
library(GOstats)
library(stringr) 
library(dplyr)
library(tidyr)
library(sigora)
library(biomaRt)
library(curl)
library(gprofiler2)
library(rrvgo)
library(GOSemSim)
library(limma)
library(clusterProfiler)
library(rlang)
library(xlsx)


#set working directory
setwd("C:/Users/aatmi/Desktop/WGCNA")

options (stringsAsFactors = FALSE)

allowWGCNAThreads()

#input the data 
data <-  read_xlsx("C:/Users/aatmi/Desktop/WGCNA/Pick_PSP_trans_proteins.xlsx", .name_repair = "universal")

data$ribaq <- factor(data$ribaq)
data$Stage <- factor(data$Stage)
#data$type <- factor(data$type)

dim(data)

#clean the data
new_data <- data[,-c(1:2)]

rownames(new_data) = data$ribaq

#check for missing values
gsg = goodSamplesGenes(new_data, verbose = 3)
gsg$allOK



#if not TRUE for gsg$allok
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(new_data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(new_data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  new_data <- new_data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

rownames(new_data) = data$ribaq


#cluster the samples
sampleTree = hclust(dist(new_data), method = "average")

#plot the sample tree
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detcted outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#plot a line to show the cut
abline(h = 0.055, col = "red")

#Determine the cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.055)
table(clust)

#clust 1 contains samples we want to keep 
keepSamples = (clust == 0)
new_data1 = new_data[keepSamples, ] 
nGenes = ncol(new_data1)
nSamples = nrow(new_data1)
rownames(new_data1) = data$ribaq

# # Identify genes with zero variance
# zeroVarianceGenes <- apply(new_data1, 2, var) == 0
# # Remove zero variance genes
# new_data1 <- new_data1[,!zeroVarianceGenes]

symbols = names(new_data1)

# Map the gene symbols to ensure they are approved
mappedSymbols <- mapIds(org.Hs.eg.db, keys = symbols, column = "SYMBOL", keytype = "SYMBOL", multiVals = "first")

# Identify and remove non-approved gene symbols
nonApprovedSymbols <- which(is.na(mappedSymbols))
if (length(nonApprovedSymbols) > 0) {
  new_data1 <- new_data1[,-nonApprovedSymbols]
  warning("Removed non-approved gene symbols from the data.")
}

# Update the row names with approved symbols
names(new_data1) <- mappedSymbols[!is.na(mappedSymbols)]

##Chosing soft threshold power
#choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from= 12, to =20, by =2))

#call the network topology analysis fucntion
sft = pickSoftThreshold(new_data1, powerVector = powers, verbose = 5)

#plot the results
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9

#scale free topology fit index as a function of soft thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[, 3])*sft$fitIndices[,2],
     xlab = "Soft threshold (power)", 
     ylab = "Scale free Topology Model Fit, signed R^2", type = "n", 
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")

#the line corresponds to using R^2 cut-off of h
abline(h = 0.90, col = "red")

#Mean connectivity as a fucntion of the thresholding power 
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft threshold (power)", 
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

#the line corresponds to using R^2 cut-off of h
abline(h = 0.90, col = "red")

#We calculate adjancies using soft thresholding power 5 
softPower = sft$powerEstimate
adjacency = adjacency(new_data1, power = 16)

correlate = corAndPvalue(new_data1)

write.csv(correlate$cor, "All_correlate.csv")

write.csv(correlate$p, "Correlate_p_value.csv")

#topological overlap matrix
#Turn adjacency to topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

#clustering using TOM
#Call the hierarichical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

#plot the resulting clustering tree
sizeGrWindow(12,9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)


#For large modules we set min module size relatively high
minModuleSize = 30

#Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro= geneTree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize
)
table(dynamicMods)

#cut numerical labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#plot the dendogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")



#Merging modules whose expression profiles are very similar 
#Calculate eigengenes
MEList = moduleEigengenes(new_data1, colors = dynamicColors)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

#plot the results
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes", 
     xlab = "", sub = "")

MEDissThres = 0.15    #1.5 for PiCK and PSP paper

#plot the cut line into dendogram
abline(h = MEDissThres, col = "red")

#call an automatic merging function
net = mergeCloseModules(new_data1, dynamicColors, cutHeight = MEDissThres, verbose = 3)

#the merged module colors
mergedColors = net$colors

#Eigengenes of the new merged modules
mergedMEs = net$newMEs

sizeGrWindow(12,9)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


table(mergedColors)

#Rename to module Colors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;

#Calculating the kME
kMEdat = signedKME(new_data1, mergedMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
write.csv(kMEdat, "Human_kME_PSP_PICK_new_modules.csv")

write.csv(mergedMEs, "Eigengene_values_PSP_PICK_modules.csv")


##########################################################
# #Module-trait relationship
# 
# behavior <- read_xlsx("Behaviour_no.xlsx")
# rownames(behavior) = behavior$Sample
# behavior = behavior[-c(1,14,22),]
# 
# nGenes = ncol(new_data1)
# nSamples = nrow(new_data1)
# 
# #recalculate MEs with color labels
# MEs0 = moduleEigengenes(new_data1, moduleColors)$eigengenes
# #MEs0 = MEs0[MEs0!='grey']
# MEs_1 =orderMEs(MEs0)
# 
# moduleTraitCor = cor(MEs_1, behavior[, -c(1:4)], use = "p")
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 
# sizeGrWindow(10,6)
# 
# testMatrix = paste(signif(moduleTraitCor, 2), "\n(",
#                    signif(moduleTraitPvalue, 1), ")", sep = "")
# 
# dim(testMatrix) = dim(moduleTraitCor)
# par(mar = c(9, 8.5, 3, 3))
# 
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = colnames(behavior[,-c(1:4)]),
#                yLabels = names(MEs_1),
#                ySymbols = names(MEs_1),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = testMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.6,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# 

#Boxplot
new_data1$ribaq = data$ribaq
new_data1$Stage = data$Stage
#new_data1$type = data$type
# new_data1$'T' = data$'T'
# new_data1$'A' = data$'A'
# new_data1$'AT' = data$'AT'
# new_data1$'V' = data$'V'
# new_data1$ATV = data$ATV

new_data1 = new_data1 %>% relocate('ribaq', .before = 'A2M')
new_data1 = new_data1 %>% relocate('Stage', .before = 'A2M')
#new_data1 = new_data1 %>% relocate('type', .before = 'A2M')
# new_data1 = new_data1 %>% relocate('T', .before = 'A1BG')
# new_data1 = new_data1 %>% relocate('A', .before = 'A1BG')
# new_data1 = new_data1 %>% relocate('AT', .before = 'A1BG')
# new_data1 = new_data1 %>% relocate('V', .before = 'A1BG')
# new_data1 = new_data1 %>% relocate('ATV', .before = 'A1BG')

new_data1$ribaq = factor(new_data1$ribaq)
new_data1$Stage = factor(new_data1$Stage)
#new_data1$type = factor(new_data1$type)
# new_data1$'T' = factor(new_data1$'T')
# new_data1$'A' = factor(new_data1$'A')
# new_data1$'AT' = factor(new_data1$'AT')
# new_data1$'V' = factor(new_data1$'V')
# new_data1$ATV = factor(new_data1$ATV)

toplot=t(mergedMEs)
cols=substring(colnames(mergedMEs),3,20)
par(mfrow=c(3,3))
#par(mar = c(bottom, left, top, right))
par(mar=c(4,3,4,3))


for (i in 1:nrow(toplot)) {
  #boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$Stage)),c('CTL','AD12','AD34', 'AD56')),
  #boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$'T')),c('T_no','T_yes')),
  #boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$'A')),c('A_yes', 'A_no', 'A_2_yes')),
          #boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$type)),c('A_yes_T_yes', 'A_no_T_no', 'A_yes_T_no', 'A_no_T_yes')),
                  #boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$Type)),c('CAA_no', 'CAA_yes')),
                          #boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$Stage)),c('A-T-CAA-', 'A+T+CAA-', 'A-T+CAA-','A+T-CAA-', 'A+T+CAA+', 'A+T-CAA+')),
                              boxplot(toplot[i,]~factor(as.vector(as.factor(new_data1$Stage)), c("CTL", "PSP", "PiD")),
          col=cols[i],ylab="ME", 
          main=rownames(toplot)[i],xlab=NULL,las=2)
  # verboseScatterplot(x=as.numeric(targets.Ref.AD_1$Age),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  # boxplot(toplot[i,]~factor(targets.Ref.AD_1$Gender),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}


# # Standardize gene symbols
# standardized_results <- checkGeneSymbols(symbols, unmapped.as.na = TRUE, species = "rat")
# standardized_symbols <- standardized_results$Suggested.Symbol
# 
# # Convert standardized symbols to sentence case (first letter capitalized)
# standardized_symbols <- str_to_title(standardized_symbols)
# 
# # Print out the number of gene symbols after standardization
# cat("Number of gene symbols after standardization: ", length(standardized_symbols), "\n")

# Connect to the Ensembl BioMart database
ensembl <- useMart("ensembl")
ensembl_hs <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Fetch Ensembl IDs from Ensembl
results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'external_gene_name',
                 values = symbols,
                 mart = ensembl_hs)

annot = data.frame(results$external_gene_name, results$ensembl_gene_id)

colnames(annot) = c("Gene_symbols", "ENSEMBL_ID")

probes2annot <- match(symbols, annot$Gene_symbols)


allENT = annot$ENSEMBL_ID[probes2annot]

dev.off()

ME_1 <- as.data.frame(cbind(allENT, symbols, moduleColors))


#getting file with modules and ribaq values
genes_exp = t(new_data1[, -c(1:2)])

genes_exp = as.data.frame(genes_exp)

names(genes_exp) = new_data1$ribaq

genes_exp$symbols = rownames(genes_exp)


modules_genes = merge(ME_1, genes_exp, by = "symbols")
write.csv(modules_genes, "All_module_PSP_PICK_modules.csv")

data_1 = new_data1[,-c(1:2)]
rownames(data_1) = new_data1$ribaq

data_2 = t(data_1)

# Assuming you have groups defined in traitData
#group <- factor(new_data1$Stage, levels = c("CTL", "AD12", "AD34", "AD56"))
#group <- factor(new_data1$'T', levels = c('T_no','T_yes'))   #no differentially expressed proteins with this classification
#group <- factor(new_data1$'A', levels = c('A_yes', 'A_no', 'A_2_yes'))    #no differentially expressed proteins with this classification
#group <- factor(new_data1$type, levels = c('A_no_T_no', 'A_yes_T_no', 'A_no_T_yes', 'A_yes_T_yes'))
#group <- factor(new_data1$type, levels = c('CAA_no', 'CAA_yes'))
# group <- factor(new_data1$Stage, levels =  c('A-T-CAA-', 'A-T+CAA-', 'A+T-CAA-', 'A+T-CAA+', 'A+T+CAA-', 'A+T+CAA+'))

group <- factor(new_data1$Stage, levels = c("CTL", "PSP", "PiD"))

design <- model.matrix(~ 0 + group)
# 
# #colnames(design) <- c("CTL", "AD12", "AD34", "AD56")
# #colnames(design) <- c('T_no','T_yes')
# #colnames(design) <- c('A_yes', 'A_no', 'A_2_yes')
# #colnames(design) <- c('A_no_T_no', 'A_yes_T_no', 'A_no_T_yes', 'A_yes_T_yes')
# #colnames(design) <- c('CAA_no', 'CAA_yes')
# colnames(design) <- make.names(c('A-T-CAA-', 'A-T+CAA-', 'A+T-CAA-', 'A+T-CAA+', 'A+T+CAA-', 'A+T+CAA+'), unique = TRUE)

colnames(design) <- c("CTL", "PSP", "PiD")

fit <- lmFit(data_2, design)

# # Create contrast matrix for all pairwise comparisons
contrast.matrix <- makeContrasts(
   PSP_vs_CTL = PSP - CTL,
   PiD_vs_CTL = PiD - CTL,
   PSP_vs_PiD = PSP - PiD,
   PiD_vs_PSP = PiD - PSP,
  levels = design
 )


# # Create contrast matrix for all pairwise comparisons
# # contrast.matrix <- makeContrasts(
# #   AD12_vs_CTL = AD12 - CTL,
# #   AD34_vs_CTL = AD34 - CTL,
# #   AD56_vs_CTL = AD56 - CTL,
# #   AD34_vs_AD12 = AD34 - AD12,
# #   AD56_vs_AD12 = AD56 - AD12,
# #   AD56_vs_AD34 = AD56 - AD34,
# #    levels = design
# #  )
# 
# # c('A_yes', 'A_no', 'A_2_yes')
#  # contrast.matrix <- makeContrasts(
#  #   A_yes_vs_A_no = A_yes - A_no,
#  #   A_yes_vs_A_2_yes = A_yes - A_2_yes,
#  #   A_no_vs_A_2_yes = A_no - A_2_yes,
#  #   levels = design
#  # )
# 
# #c('CAA_no', 'CAA_yes')
#  # contrast.matrix <- makeContrasts(
#  #    CAA_yes_vs_CAA_no = CAA_yes - CAA_no,
#  #    levels = design
#  #  )
# 
# #('T_no','T_yes')
#  # contrast.matrix <- makeContrasts(
#  #   T_no_vs_T_yes = T_no - T_yes,
#  #   levels = design
#  # )
# 
# #('A_no_T_no', 'A_yes_T_no', 'A_no_T_yes', 'A_yes_T_yes')
#  # contrast.matrix <- makeContrasts(
#  #   A_no_T_no_vs_A_no_T_yes = A_no_T_no - A_no_T_yes,
#  #   A_no_T_no_vs_A_yes_T_no = A_no_T_no - A_yes_T_no,
#  #   A_no_T_no_vs_A_yes_T_yes = A_no_T_no - A_yes_T_yes,
#  #   A_no_T_yes_vs_A_yes_T_no = A_no_T_yes - A_yes_T_no,
#  #   A_no_T_yes_vs_A_yes_T_yes = A_no_T_yes - A_yes_T_yes,
#  #   A_yes_T_no_vs_A_yes_T_yes = A_yes_T_no - A_yes_T_yes,
#  #   levels = design
#  # )
# 
# #('A-T-CAA-', 'A-T+CAA-', 'A+T-CAA-', 'A+T-CAA+', 'A+T+CAA-', 'A+T+CAA+')
#  contrast.matrix <- makeContrasts(
#   A_no_T_no_CAA_no_vs_A_no_T_yes_CAA_no = A-T-CAA- - A-T+CAA-,
#   A_no_T_no_CAA_no_vs_A_yes_T_no_CAA_no = A-T-CAA- - A+T-CAA-,
#   A_no_T_no_CAA_no_vs_A_yes_T_no_CAA_yes = A-T-CAA- - A+T-CAA+,
#   A_no_T_no_CAA_no_vs_A_yes_T_yes_CAA_no = A-T-CAA- - A+T+CAA-,
#   A_no_T_no_CAA_no_vs_A_yes_T_yes_CAA_yes = A-T-CAA- - A+T+CAA+,
#   A_no_T_yes_CAA_no_vs_A_yes_T_yes_CAA_no = A-T+CAA- - A+T+CAA- ,
#   A_no_T_yes_CAA_no_vs_A_yes_T_no_CAA_no = A-T+CAA- - A+T-CAA- ,
#   A_no_T_yes_CAA_no_vs_A_yes_T_yes_CAA_yes = A-T+CAA- - A+T+CAA+ ,
#   A_no_T_yes_CAA_no_vs_A_yes_T_no_CAA_yes = A-T+CAA- - A+T-CAA+ ,
#   A_yes_T_no_CAA_no_vs_A_yes_T_no_CAA_yes = A+T-CAA- - A+T+CAA+ ,
#   A_yes_T_no_CAA_no_vs_A_yes_T_yes_CAA_no = A+T-CAA- - A+T+CAA- ,
#   A_yes_T_no_CAA_no_vs_A_yes_T_yes_CAA_yes = A+T-CAA- - A+T+CAA+ ,
#   A_yes_T_no_CAA_yes_vs_A_yes_T_yes_CAA_no = A+T-CAA+ - A+T+CAA- ,
#   A_yes_T_no_CAA_yes_vs_A_yes_T_yes_CAA_yes = A+T-CAA+ - A+T+CAA+ ,
#   A_yes_T_yes_CAA_no_vs_A_yes_T_yes_CAA_yes = A+T+CAA- - A+T+CAA+ ,
#     levels = design
# )


fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
 
# # Extract differentially expressed proteins for each contrast
results_list <- list()
for (contrast in colnames(contrast.matrix)) {
   topTable(fit2, coef = contrast, adjust = "fdr", number = Inf) -> diffExprResults
   diffExprResults$Gene <- rownames(diffExprResults)
   diffExprResults$Comparison <- contrast
   results_list[[contrast]] <- diffExprResults[diffExprResults$adj.P.Val < 0.05, c("Gene", "Comparison")]
 }
 
# Combine differentially expressed proteins from all contrasts
diffExprProteins <- unique(unlist(results_list))
# 
# 
# # Combine differentially expressed proteins from all contrasts
diffExprProteinsDF <- do.call(rbind, results_list)

write.csv(diffExprProteinsDF, "diffExprProt_stage.csv")

#########################################################
# Initialize a data frame to store hub proteins and their kME values
# Ensure the column names in datExpr match the moduleColors
if (length(names(data_1)) != length(mergedColors)) {
  stop("The number of columns in data_1 does not match the length of moduleColors")
}

# Initialize a data frame to store hub proteins and their kME values
hubProteinsDF <- data.frame(Gene = character(), Module = character(), kME = numeric(), stringsAsFactors = FALSE)

for (module in unique(mergedColors)) {
  moduleGenes <- colnames(data_1)[mergedColors == module]
  eigengeneName <- paste0("ME", module)
  
  if (eigengeneName %in% colnames(mergedMEs)) {
    if (length(moduleGenes) > 0) {
      kME <- cor(data_1[, moduleGenes], mergedMEs[, eigengeneName], use = "p")
      
      # Check if kME values are correctly calculated
      if (length(kME) == length(moduleGenes)) {
        hubGenes <- moduleGenes[order(kME, decreasing = TRUE)[1:min(200, length(moduleGenes))]] # Top 200 hub proteins
        hubProteinsDF <- rbind(hubProteinsDF, data.frame(Gene = hubGenes, Module = module, kME = kME[match(hubGenes, moduleGenes)]))
      } else {
        warning(paste("Mismatch in kME values for module:", module))
      }
    }
  } else {
    warning(paste("Module eigengene", eigengeneName, "not found in MEs"))
  }
}
      
# Intersect with differentially expressed proteins
prospectiveBiomarkers <- intersect(hubProteinsDF$Gene, diffExprProteins)

# Filter hubProteinsDF to include only prospective biomarkers
prospectiveBiomarkersDF <- hubProteinsDF[hubProteinsDF$Gene %in% prospectiveBiomarkers, ]

# Display the resulting data frame
write.csv(prospectiveBiomarkersDF, "differentially_expressed_hubproteins_PSP_PICK_modules.csv")

# Intersect with differentially expressed proteins
prospectiveBiomarkersDF_new <- merge(hubProteinsDF, diffExprProteinsDF, by = "Gene")

write.csv(prospectiveBiomarkersDF_new, "differentially_expressed_hubproteins_groups__PSP_PICK__modules.csv")

##########################################################################
#extracting the genes from 





