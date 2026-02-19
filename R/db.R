
#' database
#'
#' @param ct string
#'
#' @return vector
#' @export
#'
MarkerDB <- function(ctype){
    genes <- ''

    if ( ctype == 'male' ) { genes <- c('Sry','Eif2s3y','Ddx3y','Uty','Kdm5d') } 
    if ( ctype == 'immune' ) { genes <- c('Ptprc', 'Cd69') }
    if ( ctype == 'cycling' ) { genes <- c('Top2a', 'Mki67', 'Foxm1') }
    if ( ctype == 'astrocyte' ) { genes <- c('Gfap', 'Sox9', 'Aqp4') }

    return(genes)
}







