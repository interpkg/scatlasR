# Hua Sun
# Edited original 'ScType' codes (https://github.com/IanevskiAleksandr/sc-type)
# v3.2 2024-07-26 (Updated to support Seurat v4 & v5)


#library(dplyr)
#library(HGNChelper)




#' gene_sets_prepare
#'
#' @param path_to_db_file db path
#' @param cell_type cell type tissue
#' @export
#'
gene_sets_prepare <- function(path_to_db_file, cell_type){
    cell_markers = openxlsx::read.xlsx(path_to_db_file)
    cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    # correct gene symbols from the given DB (up-genes)
    cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
        
        markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
        markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
        markers_all = sort(markers_all)
        
        if(length(markers_all) > 0){
        markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
        paste0(markers_all, collapse=",")
        } else {
        ""
        }
    })
    
    # correct gene symbols from the given DB (down-genes)
    cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
        
        markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
        markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
        markers_all = sort(markers_all)
        
        if(length(markers_all) > 0){
        markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
        paste0(markers_all, collapse=",")
        } else {
        ""
        }
    })
    
    cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
    gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
    
    list(gs_positive = gs, gs_negative = gs2)
}





#' sctype_score
#'
#' @param scRNAseqData matrix
#'
#' @export
#'
sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
    # check input matrix
    if(!is.matrix(scRNAseqData)){
        warning("scRNAseqData doesn't seem to be a matrix")
    } else {
        if(sum(dim(scRNAseqData))==0){
        warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
        }
    }
    
    # marker sensitivity
    marker_stat = sort(table(unlist(gs)), decreasing = T); 
    marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                        gene_ = names(marker_stat), stringsAsFactors = !1)

    # convert gene names to Uppercase
    if(gene_names_to_uppercase){
        rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
    }
    
    # subselect genes only found in data
    names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
    gs = lapply(1:length(gs), function(d_){ 
        GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
    gs2 = lapply(1:length(gs2), function(d_){ 
        GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
    names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
    cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
    
    # z-scale if not
    if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
    
    # multiple by marker sensitivity
    for(jj in 1:nrow(cell_markers_genes_score)){
        Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
    }
    
    # subselect only with marker genes
    Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
    
    # combine scores
    es = do.call("rbind", lapply(names(gs), function(gss_){ 
        sapply(1:ncol(Z), function(j) {
        gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
        sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
        if(is.na(sum_t2)){
            sum_t2 = 0;
        }
        sum_t1 + sum_t2
        })
    })) 
    
    dimnames(es) = list(names(gs), colnames(Z))
    es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
    
    es.max
}




#' CellType Anno SCType
#'
#' @param obj seurat object
#'
#' @export
#'
ScTypeAnnotation <- function(obj, assay, reduction, f_marker_db, tissue, out_path)
{
    ##-- 1. Loding DB
    # load libraries
    lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
    
    # DB file
    if (f_marker_db == ''){
        # default
        f_marker_db <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    }
    
    print(f_marker_db)
    
    # prepare gene sets
    gs_list <- gene_sets_prepare(f_marker_db, tissue)

    # //-------- 2024-07-26
    # check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
    seurat_package_v5 <- isFALSE('counts' %in% names(attributes(obj[[assay]])));
    print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

    # extract scaled scRNA-seq matrix
    scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(obj[[assay]]$scale.data) else as.matrix(obj[[assay]]@scale.data)
    # --------//

    ##-- 2. Annotation cell types
    # get cell-type by cell matrix
    es.max = sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    # testing code (the results same as @scale.data, but counts will take long time due to z-score calculation)
    #count.data <- as.data.frame(obj[[assay]]@counts)
    #es.max = sctype_score(scRNAseqData = count.data, scaled = FALSE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

    # merge by cluster
    cL_resutls = do.call("rbind", lapply(unique(obj$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj$seurat_clusters==cl, ])]), decreasing = !0)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj$seurat_clusters==cl)), 10)
    }))
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

    # set low-confident (low ScType score) clusters to "unknown"
    #sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unclassified"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < 100] = "Unclassified"
    print(sctype_scores[,1:3])

    ##-- 3. UMAP
    # cell types on UMAP plot
    obj$cell_type = ""
    obj$cell_type2 = ""
    marker_db <- readxl::read_excel(f_marker_db)
    for(j in unique(sctype_scores$cluster)){
        # add cell type
        celltype = as.character(sctype_scores$type[sctype_scores$cluster==j][1])
        obj$cell_type[obj$seurat_clusters == j] = celltype
        
        # add cell type short name (shortName) -- 'cell_type2'
        shortname = as.character(marker_db$shortName[marker_db$cellName==celltype][1])
        obj$cell_type2[obj$cell_type == celltype] = shortname

        # add cell type 'shortName2' -- 'cell_type3' (for v3)
        #shortname2 = as.character(marker_db$shortName2[marker_db$cellName==celltype][1])
        #obj$cell_type3[obj$cell_type == celltype] = shortname2
    }

    obj$cell_type2[is.na(obj$cell_type2)] <- 'Unclassified'
    # v3.1
    obj$cluster_plus <- paste0('C', obj$seurat_clusters, '.', obj$cell_type2)
    obj$sctype.score <- sctype_scores$scores[match(obj$seurat_clusters, sctype_scores$cluster)]
    obj$sctype.score <- format(round(obj$sctype.score, 1), nsmall = 1)

    # 
    obj$UMAP_1 <- as.data.frame(obj@reductions[[reduction]]@cell.embeddings)[[1]]
    obj$UMAP_2 <- as.data.frame(obj@reductions[[reduction]]@cell.embeddings)[[2]]

    return(obj)
}





#' Auto detect tissue: Guess tissue type from seurat data
#'
#' @param path_to_db_file db path
#' @param seuratObject seurat object
#'
#' @export
#'
auto_detect_tissue_type <- function(path_to_db_file, seuratObject, scaled, assay = "RNA", ...){
  
  # get all tissue types in DB
    db_read = openxlsx::read.xlsx(path_to_db_file); tissues_ = unique(db_read$tissueType); result_ = c()
    
    for(tissue in tissues_){ print(paste0("Checking...", tissue));
        
    # prepare gene sets
    gs_list = gene_sets_prepare(path_to_db_file, tissue);
    
    # prepare obj
    if(scaled){
      obj = as.matrix(seuratObject[[assay]]@scale.data)
    } else {
      obj = as.matrix(seuratObject[[assay]]@counts)
    }
                            
    es.max = sctype_score(scRNAseqData = obj, scaled = scaled, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                          marker_sensitivity = gs_list$marker_sensitivity, verbose=!0);
  
    cL_resutls = do.call("rbind", lapply(unique(seuratObject@meta.data$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es.max[ ,rownames(seuratObject@meta.data[seuratObject@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
    }))
                          
    dt_out = cL_resutls %>% group_by(cluster) %>% top_n(n = 1)
    
    # return mean score for tissue
    result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
  }
  
  # order by mean score
  result_ = result_[order(-result_$score),]
  
  # plot 
  barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
          xlab="Tissue", ylab="Summary score",  main="The higher summary score, the more likely tissue type is")
  
  result_
}





