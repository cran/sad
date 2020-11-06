

#' @rdname sadforecast
#' @export
plot.sadforecast <- function(x,  mfrow=NULL, col=NULL,  ...){
    if( is.null(col) ) col <- raincols
    oldpar <- graphics::par( no.readonly=TRUE )
    on.exit( graphics::par(oldpar) )
    
    nc <- ceiling( sqrt( length(x) ) )
    nr <- ceiling( length(x)/nc )
    if( is.null(mfrow) ) mfrow <- c(nr,nc)
    
    graphics::par( mfrow=mfrow, mar=c(1,1,2,1) )
    for( i in 1:length(x) ){
        graphics::image( x[[i]], col=col, xaxt="n", yaxt="n", xlab="", ylab="", 
                         useRaster=TRUE, main=names(x)[i], ... )
    }
    
    
}
