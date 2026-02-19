

#' Gene expression markers for all group
#'
#' @param obj object
#' @param assay_use assay
#' @return data frame
#' @export
#'
RunFindAllMarkers <- function(
    obj=NULL, 
    assay_use='SCT', 
    idents='seurat_clusters', 
    logfc=.25, 
    min_pct=0.1, 
    min_diff_pct=0, 
    recorrect_umi=TRUE
) { 
    Idents(obj) <- idents

    diff_exp <- Seurat::FindAllMarkers(
        object = obj,
        only.pos = TRUE,
        logfc.threshold = logfc,
        min.pct = min_pct,
        min.diff.pct = min_diff_pct,
        recorrect_umi = recorrect_umi
    )

    diff_exp$diff_pct <- diff_exp$pct.1 - diff_exp$pct.2

    # filter genes not in SCT
    all_genes <- rownames(obj@assays[[assay_use]]@data)
    # wrong gene name from duplicated gene markers
    wrong_gene_name <- setdiff(rownames(diff_exp), all_genes)

    diff_exp$gene <- rownames(diff_exp)
    diff_exp$gene[diff_exp$gene %in% wrong_gene_name] <- str_sub(diff_exp$gene[diff_exp$gene %in% wrong_gene_name], end = -2)

    return(diff_exp)
}





#' Seruat Anchors
#'
#' @param que_obj seurat object of query
#' @param ref_obj seurat object of reference
#' @export
#'
SeuratAnchors <- function(
    que_obj=NULL, 
    ref_obj=NULL, 
    nfeatures=5000,
    norm_method='LogNormalize',
    reduction='cca',
    dim1=1,
    dim2=30,
    weight_reduction='pca'
){
    ref_obj <- Seurat::UpdateSeuratObject(ref_obj)

    ref_obj <- Seurat::FindVariableFeatures(
      object = ref_obj,
      nfeatures = nfeatures
    )

    transfer.anchors <- Seurat::FindTransferAnchors(
      reference = ref_obj,
      query = que_obj,
      normalization.method = norm_method,
      reduction = reduction,
      dims = dim1:dim2
    )

    predicted.labels <- Seurat::TransferData(
      anchorset = transfer.anchors,
      refdata = ref_obj$subclass,
      weight.reduction = que_obj[[weight_reduction]],
      dims = dim1:dim2
    )

    que_obj <- Seurat::AddMetaData(object = que_obj, metadata = predicted.labels)

    return(que_obj)
}






#' Seurat Add Score
#'
#' @param obj seurat object
#' @param name score title 
#' @param features gene vectors 
#'
#' @return seurat object
#'
#' @import Seurat
#'
#' @export
#'
SeuratAddScorePlus <- function(
    obj=NULL, 
    assay='SCT',
    slot='data',
    name='', 
    features='',
    fdata='',
    form='seurat'
){
    library(tidyr)

    df <- NULL
    if (form == 'seurat'){
        df <- read.table(fdata, sep='\t', header=T)
        df <- df[,c('cluster', 'gene')]
    }
    if (form == 'sctype'){
        df <- d_db <- readxl::read_excel(fdata)
        df <- df[,c('cellName', 'geneSymbolmore1')]

        df <- df %>% separate_rows(geneSymbolmore1, sep = ",")
        df <- as.data.frame(df)
        colnames(df) <- c('cluster', 'gene')
    }


    all_gene <- rownames(obj[[assay]]@data)

    for (group in unique(df$cluster)){
        genelist <- df[df$cluster == group, 'gene']
        gene_intersect <- intersect(genelist, all_gene)

        obj <- Seurat::AddModuleScore(
            object = obj,
            assay = assay,
            slot = slot,
            features = list(gene_intersect),
            name = group
        )
    }

    return(obj)
}






#' Seurat Add Score
#'
#' @param obj seurat object
#' @param name score title 
#' @param features gene vectors 
#'
#' @return seurat object
#'
#' @import Seurat
#'
#' @export
#'
SeuratAddScore <- function(
    obj=NULL, 
    assay='SCT',
    slot='data',
    name='', 
    features=''
){
    all_gene <- rownames(obj[[assay]]@data)
    gene_intersect <- intersect(features, all_gene)

    obj <- AddModuleScore(
        object = obj,
        assay = assay,
        slot = slot,
        features = list(gene_intersect),
        name = name
    )

    return(obj)
}






