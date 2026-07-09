source("C:/Users/jespo/Desktop/WGCNA/generate_modules_WGCNA_AD_class.R")

#change the name of columns
colnames(ME_1)[1] = "ENSEMBL.Gene.ID"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
colnames(ME_1)[2] = "GeneSymbol"
colnames(ME_1)[3] = "Initially.Assigned.Module.Color"

# Number of NAs in each column
na_counts_per_column <- colSums(is.na(ME_1))
print(na_counts_per_column)

#remove nas (939 ensembl ids missing)
ME_1 = na.omit(ME_1)

#create input files for gProfiler 
# Create a list of genes for each module
moduleGenes <- split(ME_1$GeneSymbol, ME_1$Initially.Assigned.Module.Color)

# Function to perform GO enrichment analysis for a given module
performGOEnrichment <- function(genes) {
  gostres <- gost(query = genes, organism = "hsapiens", sources = c('GO:BP', 'GO:MF',
                                                                       'GO:CC', 'KEGG', 'REAC', 
                                                                       'TF', 'MIRNA', 'CORUM', 'WP'), correction_method = "fdr",
                  domain_scope = "annotated", user_threshold = 0.5, significant = FALSE)
  return(gostres$result)
}

# Apply the function to each module
goResults <- lapply(moduleGenes, performGOEnrichment)

# Combine the results for each module into a single data frame
combineResults <- function(filteredResults) {
  combined_df <- do.call(rbind, lapply(names(filteredResults), function(module) {
    if (!is.null(filteredResults[[module]]) && nrow(filteredResults[[module]]) > 0) {
      df <- filteredResults[[module]]
      df$Module <- module  # Add a column for the module name
      return(df)
    } else {
      return(NULL)
    }
  }))
  
  # Ensure all columns are of atomic type
  combined_df <- as.data.frame(lapply(combined_df, function(x) {
    if (is.list(x)) {
      sapply(x, paste, collapse = ", ")
    } else {
      x
    }
  }))
  
  return(combined_df)
}


combinedGOResults <- combineResults(goResults)

# Save the combined data frame to a CSV file
saveToCSV <- function(combinedGOResults, filename) {
  if (!is.null(combinedGOResults) && nrow(combinedGOResults) > 0) {
    write.csv(combinedGOResults, file = filename, row.names = FALSE)
    message("Results saved to ", filename)
  } else {
    message("No results to save.")
  }
}

# Save the combined results to a CSV file
saveToCSV(combinedGOResults, "C:/Users/jespo/Desktop/WGCNA/HS_AD_All_GOResults.csv")


# Function to filter GO results by a significance threshold
filterGOResults <- function(go_results, pval_threshold = 0.05) {
  if (is.null(go_results) || nrow(go_results) == 0) {
    return(NULL)
  }
  filtered_results <- go_results[go_results$p_value <= pval_threshold, ]
  return(filtered_results)
}

# Apply the filter function to each module's results
filteredResults <- lapply(goResults, filterGOResults)

combinedResults <- combineResults(filteredResults)

# Save the combined results to a CSV file
saveToCSV(combinedResults, "C:/Users/jespo/Desktop/WGCNA/HS_AD_All_filteredResults.csv")

