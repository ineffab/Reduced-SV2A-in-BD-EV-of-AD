source("C:/Users/jespo/Desktop/WGCNA/Functional_annotation_WGCNA_AD_class.R")

#ORA pruning for GO:BP

#extract the GO:BP for each module 

######################################

#try another ORA pruning methods

#####################################

##first on the rat data
#required libraries
library(org.Hs.eg.db)
library(HGNChelper)
library(devtools)
library(GO.db)
library(AnnotationDbi)
library(AnnotationHub)
library(AnnotationForge)
library(GOstats)
library(GOSemSim)
library(biomaRt)
library(ensembldb)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)


# Load GO semantic similarity data
godata_BP <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
godata_MF <- godata('org.Hs.eg.db', ont="MF", computeIC=TRUE)
godata_CC <- godata('org.Hs.eg.db', ont="CC", computeIC=TRUE)


# Define the pruning function
prune_terms <- function(results, threshold=0.7, ontology_type="BP") {
  godata <- switch(ontology_type,
                   BP = godata_BP,
                   MF = godata_MF,
                   CC = godata_CC,
                   stop("Invalid ontology type"))
  
  keep <- rep(TRUE, nrow(results))
  for (i in 1:(nrow(results)-1)) {
    for (j in (i+1):nrow(results)) {
      if (!is.na(keep[j]) && keep[j]) {
        sim <- tryCatch({
          mgoSim(results$term_id[i], results$term_id[j], semData=godata, measure="Wang")
        }, error = function(e) {
          NA
        })
        if (!is.na(sim) && sim > threshold) {
          if (!is.na(results$p_value[i]) && !is.na(results$p_value[j])) {
            if (results$p_value[i] < results$p_value[j]) {
              keep[j] <- FALSE
            } else {
              keep[i] <- FALSE
            }
          }
        }
      }
    }
  }
  return(results[keep, , drop = FALSE])
}

# Function to extract GO terms and prune them
prune_go_terms_in_module <- function(module_data) {
  go_types <- c("GO:BP", "GO:MF", "GO:CC")
  pruned_data <- lapply(go_types, function(go_type) {
    go_data <- module_data %>% dplyr::filter(grepl(go_type, source))
    if (nrow(go_data) > 0) {
      ontology_type <- substr(go_type, 4, 5)
      print(paste("Processing ontology:", go_type, "with type:", ontology_type))
      pruned_subset <- prune_terms(go_data, ontology_type = ontology_type)
      return(pruned_subset)
    }
    return(NULL)
  })
  pruned_data <- do.call(rbind, pruned_data)
  non_go_data <- module_data %>% dplyr::filter(!source %in% go_types)
  combined_data <- bind_rows(pruned_data, non_go_data)
  return(combined_data)
}

# Apply pruning to each module
pruned_results <- lapply(filteredResults, prune_go_terms_in_module)


combinedResults <- combineResults(pruned_results)

# Save the combined results to a CSV file
saveToCSV(combinedResults, "C:/Users/jespo/Desktop/WGCNA/HS_AD_ORA_pruned.csv")

#extract genes from the ORA pruned data




##################################

#DO for humans

library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)

# Check available keys for SYMBOL
valid_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
cat("First 10 valid keys for SYMBOL:\n")
print(head(valid_keys, 10))

# Function to safely convert gene symbols to Entrez IDs
safe_bitr <- function(genes, fromType, toType, OrgDb) {
  gene_entrez <- tryCatch(
    bitr(genes, fromType = fromType, toType = toType, OrgDb = OrgDb),
    error = function(e) {
      warning(paste("Error in mapping genes:", e))
      return(data.frame(SYMBOL = genes, ENTREZID = NA))
    }
  )
  return(gene_entrez)
}

# Initialize a list to store DO annotation results
do_results_list <- list()


# # Loop through each module and perform DO annotation
for (module in names(moduleGenes)) {
  cat(paste0("Processing ", module, "...\n"))
  #   
  #   # Check if gene symbols are valid
  invalid_genes <- setdiff(moduleGenes[[module]], valid_keys)
  if (length(invalid_genes) > 0) {
    cat("Invalid gene symbols found:\n")
    print(invalid_genes)
    next
  }
# Convert gene symbols to Entrez IDs
  gene_entrez <- safe_bitr(moduleGenes[[module]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  #   
  # Report unmapped genes
  unmapped_genes <- gene_entrez[is.na(gene_entrez$ENTREZID), "SYMBOL"]
  if (length(unmapped_genes) > 0) {
    cat("Unmapped genes:\n")
    print(unmapped_genes)
  }
 
  
  # Perform Disease Ontology (DO) annotation if there are mapped genes
  if (sum(!is.na(gene_entrez$ENTREZID)) > 0) {
    do_results <- enrichDO(gene = gene_entrez$ENTREZID[!is.na(gene_entrez$ENTREZID)], ont = "HDO", organism ="hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 5, readable = TRUE)
    
    # Convert to data frame and store the results in the list
    do_results_df <- as.data.frame(do_results)
    if (nrow(do_results_df) > 0) {
      do_results_df$Module <- module  # Add a column for the module name
      do_results_list[[module]] <- do_results_df
    }
  } else {
    cat("No valid Entrez IDs found for ", module, "\n")
  }
} 

  combineDOResults <- function(filteredResults) {
    combinedDO_df <- do.call(rbind, filteredResults)
    
    # Ensure all columns are of atomic type
    combinedDO_df <- as.data.frame(lapply(combinedDO_df, function(x) {
      if (is.list(x)) {
        sapply(x, paste, collapse = ", ")
      } else {
        x
      }
    }))
    
    return(combinedDO_df)
  }
  
  # Combine the results into a single data frame
  combinedDO_results <- combineDOResults(do_results_list)
  
  
  # Save the combined results to a CSV file
  saveToCSV(combinedDO_results, "Human_DO_annotations.csv")

