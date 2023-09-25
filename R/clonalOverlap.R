#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the calculation and visualizations of 
#' various overlap metrics for clonotypes.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' clonalOverlap(combined, 
#'               cloneCall = "gene", 
#'               method = "jaccard")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param method The method to calculate the "overlap", 
#  "morisita", "jaccard" indices, "cosine" similarity or "raw" 
#' for the base numbers.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @param pcoa Plot PCoA plot instead of heatmap
#' @param group.by grouping variable in PCoA plot, default to NULL
#' @param point.size point size in PCoA plot
#' @importFrom stringr str_sort str_to_title
#' @importFrom reshape2 melt
#' @importFrom stats quantile
#' @importFrom ape pcoa
#' @export
#' @return ggplot of the clonotypic overlap between elements of a list
clonalOverlap <- function(df, 
                          cloneCall = "strict", 
                          method = NULL, 
                          chain = "both", 
                          split.by = NULL,
                          exportTable = FALSE,
                          palette = "inferno",
                          pcoa = FALSE,
                          group.by = NULL,
                          point.size = 3) {
    if(method == "morisita") {
      return_type = "freq"
    } else {
      return_type = "unique"
    }
    if (!is.null(group.by)) {
        group <- vector(mode="character", length=length(df))
        for (i in seq_along(df)) {
            group[i] <- unique(df[[i]][, group.by])
        }
        names(group) <- names(df)
        group_len <- length(unique(group))
    } else {
        group <- NULL; group_len <- 0;
    }
    cloneCall <- theCall(cloneCall)
    df <- .data.wrangle(df, split.by, cloneCall, chain)
    df <- df[order(names(df))]
    values <- str_sort(as.character(unique(names(df))), numeric = TRUE)
    df <- df[quiet(dput(values))]
    num_samples <- length(df[])
    names_samples <- names(df)
    length <- seq_len(num_samples)
    
    if (chain != "both") {
      for (i in seq_along(df)) {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
    }
    
    #Selecting Index Function
    indexFunc <- switch(method,
                        "morisita" = .morisitaCalc,
                        "jaccard"  = .jaccardCalc,
                        "raw"      = .rawCalc,
                        "overlap"  = .overlapCalc,
                        "cosine"  = .cosineCalc,
                        stop("Invalid method provided"))
    
    #Calculating Index 
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    coef_matrix <- .calculateIndex(df, length, cloneCall, coef_matrix, indexFunc, return_type)
    
    #Data manipulation
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples

    if (exportTable == TRUE) { 
      return(coef_matrix) 
    }
    if (pcoa) {
        m <- as.matrix(coef_matrix)
        m[lower.tri(m)] <- t(m)[lower.tri(m)]
        m <- as.dist(1-m, upper=TRUE)
        # res_pcoa <- cmdscale(as.dist(1-m, upper=TRUE), k=2, eig=TRUE)
        # eig_sum <- sum(res_pcoa$eig)
        # expv <- ((res_pcoa$eig / eig_sum) %>% round(2)) * 100
        res_pcoa <- ape::pcoa(m, correction="lingoes")
        if ("Rel_corr_eig" %in% colnames(res_pcoa$values)) {
            expv <- round(res_pcoa$values[,"Rel_corr_eig"] * 100, 2)        
        } else {
            expv <- round(res_pcoa$values[,"Relative_eig"] * 100, 2)
        }
        plot <- res_pcoa$vectors %>%
            data.frame() %>%
            .[,c(1,2)] %>%
            `colnames<-`(c("PC1","PC2")) %>%
            mutate(sample=row.names(.)) %>%
            mutate(group=group[sample]) %>%
            ggplot(aes(x=PC1, y=PC2, fill=group))+
            geom_point(shape=21, size=point.size) +
            scale_fill_manual(values=.colorizer(palette, group_len),
                na.value = "white",
                name=group.by)+
            xlab(paste0("PC1 (",expv[1]," %)"))+
            ylab(paste0("PC2 (",expv[2]," %)"))+
            theme_classic()
        return(plot)
    }
    mat_melt <- suppressMessages(melt(as.matrix(coef_matrix)))
    
    mean_value <- mean(na.omit(mat_melt[,"value"]))
    
    plot <- ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
                geom_tile() + 
                geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], fill = NA, lwd= 0.25, color = "black") +
                labs(fill = str_to_title(method)) +
                geom_text(aes(label = round(value, digits = 3), 
                              color = ifelse(value <= mean_value,
                                             "white", "black"))) +
                scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
                scale_color_identity() +
                theme_classic() + 
                theme(axis.title = element_blank())
    return(plot) 
}

# Helper function to prepare data
.prepareDataFrame <- function(df, cloneCall, return_type = "unique") {
  if (return_type == "unique") {
    return(unique(df[, cloneCall]))
  } else if (return_type == "freq") {
    temp_df <- data.frame(table(df[, cloneCall]))
    colnames(temp_df) <- c(cloneCall, 'Count')
    temp_df[, 2] <- as.numeric(temp_df[, 2])
    return(temp_df)
  }
}

# Helper function for common loop and conditional structure
.calculateIndex <- function(df, length, cloneCall, coef_matrix, indexFunc, return_type = "unique") {
  for (i in seq_along(length)) {
    df_i <- .prepareDataFrame(df[[i]], cloneCall, return_type)
    for (j in seq_along(length)) {
      if (i >= j) { next }
      df_j <- .prepareDataFrame(df[[j]], cloneCall, return_type)
      coef_matrix[i, j] <- indexFunc(df_i, df_j)
    }
  }
  return(coef_matrix)
}

# Morisita Index calculation function
.morisitaCalc <- function(df_i, df_j) {
  merged <- merge(df_i, df_j, by = names(df_i)[1], all = TRUE)
  merged[is.na(merged)] <- 0
  
  X <- sum(merged[, 2])
  Y <- sum(merged[, 3])
  
  num <- 2 * sum(merged[, 2] * merged[, 3])
  den <- ((sum(df_i[, 2]^2) / (X^2)) + (sum(df_j[, 2]^2) / (Y^2))) * X * Y
  
  return(num / den)
}

# Jaccard Index calculation function
.jaccardCalc <- function(df_i, df_j) {
  overlap <- length(intersect(df_i, df_j))
  return(overlap / (length(df_i) + length(df_j) - overlap))
}

# Raw Index calculation function
.rawCalc <- function(df_i, df_j) {
  return(length(intersect(df_i, df_j)))
}

# Overlap Index calculation function
.overlapCalc <- function(df_i, df_j) {
  overlap <- length(intersect(df_i, df_j))
  return(overlap / min(length(df_i), length(df_j)))
}

# Overlap Index calculation function
.cosineCalc <- function(df_i, df_j) {
  all_species <- unique(c(df_i, df_j))
  vector_location1 <- as.integer(all_species %in% df_i)
  vector_location2 <- as.integer(all_species %in% df_j)
  
  return(sum(vector_location1 * vector_location2) / 
           (sqrt(sum(vector_location1^2)) * sqrt(sum(vector_location2^2))))
}
