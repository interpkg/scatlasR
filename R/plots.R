


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                  Seurat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#' Seurat DotPlot
#'
#' @param obj object
#' @param group
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
Seurat_DotPlotMarkers <- function(
    obj=NULL, 
    assay='SCT', 
    d_markers=NULL,
    scale=TRUE,
    col_min = -2.5,
    n=5, 
    group=''
) {
    topX.markers <- data.frame(d_markers %>% dplyr::group_by(cluster) %>% dplyr::slice(1:n))

    meta <- obj@meta.data
    sorted_name <- unique(topX.markers$cluster)

    color.scheme <- rev(brewer.pal(9,"RdBu"))
    p <- Seurat::DotPlot(obj, assay = assay, group.by = group, features = unique(topX.markers$gene), scale = scale, col.min = col_min, dot.min = 0) +
            theme_bw(base_line_size=0.1) +
            scale_size_area(max_size = 3) +
            scale_color_gradientn(colors=color.scheme, limits = c(-2.5, 2.5)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            labs(title='', x='', y='') +
            scale_y_discrete(limits=sorted_name) +
            theme(
                text=element_text(size=7), 
                axis.text=element_text(colour="black", size=7), 
                axis.title.x=element_blank()) +
            theme(
                    legend.position="bottom", 
                    legend.key.width = unit(4, 'mm'),
                    legend.key.height = unit(2.5, 'mm'),
                    legend.title = element_text(size=7),
                    legend.text=element_text(size=7)
            )

    p

}







