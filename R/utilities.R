


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                  Filter
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' Filter markers for Seurat form marker table
#'
#' @param data object
#'
#' @return data frame
#' @export
#'
FilterMarkers_SeuratForm <- function(
    data = NULL,
    gene_col = 'gene',
    p_val_col = 'p_val_adj',
    log2fc_col = 'avg_log2FC',
    pct1_col = 'pct.1',
    pval = 0.0001,
    log2fc = 0.585,
    pct1 = 0.3
){
    diff_exp <- data %>% 
                    filter(!grepl('^MT-|^MTRNR', .data[[gene_col]])) %>%  
                    filter(!grepl('^RPL|^RPS|^MRP', .data[[gene_col]])) %>% 
                    filter(!grepl('^AC[0-9]|^AL[0-9]', .data[[gene_col]])) %>% 
                    filter(!grepl('^LINC', .data[[gene_col]])) %>%
                    filter(.data[[p_val_col]] < pval & .data[[log2fc_col]] > log2fc & .data[[pct1_col]] > pct1)

    return(diff_exp)
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                  Others
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' Convert the Seurat marker table format to the scType format
#'
#' @param data object
#'
#' @return data frame
#' @export
#'
SeuratMarkerForm_to_scType <- function(
    data = NULL,
    group_col = 'cluster',
    gene_col = 'gene',
    tissueType = 'Brain',
    source_name = '',
    outfile = 'sctypeForm.xlsx'
) {
    merged_data <- data %>%
        group_by(.data[[group_col]]) %>% 
        summarise( merged_symbols = paste(unique(.data[[gene_col]]), collapse = ",") )
    colnames(merged_data) <- c('cellName', 'geneSymbolmore1')

    merged_data$tissueType <- tissueType
    merged_data <- merged_data[, c('tissueType', 'cellName', 'geneSymbolmore1')]

    merged_data$geneSymbolmore2 <- ''
    merged_data$shortName <- merged_data$cellName
    merged_data$shortName2 <- ''
    merged_data$Class <- ''
    merged_data$Subclass <- ''
    merged_data$Region <- ''
    merged_data$Subregion <- ''
    merged_data$source <- source_name

    # Save to Excel with options
    openxlsx::write.xlsx(merged_data, file = outfile, sheetName = "Sheet1", rowNames = FALSE, overwrite = TRUE)
}





#' CompareMarkers
#'
#' @param data object
#' @param col_ct cell column name
#' @param col_gene gene column name
#'
#' @return data frame
#' @export
#'
CompareMarkersByPublicDB <- function(
    data=NULL, 
    fdb=NULL,
    col_ct='cellName',
    col_gene='geneSymbol'
) { 
    d_db <- readxl::read_excel(fdb)
    
    df <- NULL
    for (x in unique(data$cluster)){
        print(x)
        geneset <- data$gene[data$cluster == x]
        print(length(geneset))

        for (ct in unique(d_db[[col_ct]])){
            markers <- d_db[[col_gene]][ d_db[[col_ct]] == ct ]
            markers <- str_split(markers, ',')[[1]]
            overlap_genes <- intersect(geneset, markers)
            geneSymbol <- paste(overlap_genes, collapse = ",")
            avg_log2fc <- round(mean(data$avg_log2FC[data$cluster == x & data$gene %in% overlap_genes]), 3)
            n <- length(overlap_genes)
            if (n == 0){ next }
            new_row <- c(cluster=x, cellType=ct, geneCount=n, avg_log2fc=avg_log2fc, geneSymbol=geneSymbol)

            df <- rbind(df, new_row)
        }
    }

    df <- as.data.frame(df)
    print(dim(df))


    d_res <- df %>% filter(as.numeric(geneCount)>1)
    d_res <- d_res %>% arrange(cluster, desc(as.numeric(avg_log2fc)))
    print(dim(d_res))

    return(d_res)
}







