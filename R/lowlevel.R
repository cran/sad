#' class for a list of forecasts
#'
#' check that a list of forecasts fulfills all requirements to be verified by our method
#' @param x a list of 2 or more 2D matrices with equal sizes and no missing or inifinite values
#' @param mfrow vector with the number of rows and columns you would like in the plot
#' @param col color scale for the plot
#' @param ... further arguments passed to \code{image}
#' @return an object of class \code{sadforecast} 
#' @details \code{as.sadforecast} does nothing except check that everything is as it should be, add the attributes that can be changed by \code{prepare_sad} and provide a method for quick plots of the data.
#' @examples
#' data( rrain )
#' ra <- list( rrain[1,1,,], rrain[4,5,,], rrain[2,7,,] )
#' ra <- as.sadforecast(ra)
#' plot(ra)
#' @name sadforecast
NULL

#' @rdname sadforecast
#' @export
as.sadforecast <- function(x){
    # check that we got a list
    if( !is.list(x) ) stop( "please provide observation and forecasts as a list of matrices" )
    
    # check that we got at least one pair
    n <- length(x)
    if( n<2 ) stop( "need at least two fields to compare" )
    
    # check the dimensions
    dims <- unlist( lapply(x, function(x) length(dim(x))) )
    if( any( dims != 2 ) ) stop( "all inputs must be 2D matrices" )
    nx <- unlist( lapply( x, nrow ) )
    ny <- unlist( lapply( x, ncol ) )
    if( any( nx[-1] != nx[1] ) | any( ny[-1] != ny[1] )  ) stop( "all inputs must have the same dimensions" )

    # check that we got no missing values
    if( any( is.na( unlist( x ) ) ) )      stop( "no missing values allowed" ) 
    if( any( !is.finite( unlist( x ) ) ) ) stop( "no infinite values allowed" ) 

    # maybe add names
    if( is.null( names(x) ) ) names(x) <- c( "obs", paste0( "f", 1:(n-1) ) )

    # add attributes if they aren't given already
    defaults <- list( log=FALSE, px=1:nx[1], py=1:ny[1], rsm=0, xmin=-Inf, boundaries="none" )
    for( n in names( defaults ) ){
        if( is.null( attr(x, n) ) ) attr( x, n ) <- defaults[[n]]
    }
    
    return( structure( x, class="sadforecast" ) )
}

#' prepare a sad forecast for verification
#'
#' remove small values, apply log-transform, smooth borders, handle boundary conditions
#' @param x a list of 2 or more 2D matrices with equal sizes and no missing or inifinite values, as required by \code{as.sadforecast}
#' @param xmin values smaller than \code{xmin} are set to zero
#' @param log  logical, do you want to log-transfrom the data? (recommended for precipitation)
#' @param rsm  number of pixels which are linearly smoothed at the edge
#' @param Nx   size to which the data is extended in x-direction
#' @param Ny   size to which the data is extended in y-direction
#' @param boundaries how to handle the boundary conditions, either "pad", "mirror" or "periodic"
#' @return an object of class \code{sadforecast} which has been prepared in the desired way.
#' @details the positions within the extended field where the original field resides are output as attributes "px", "py" of the result. The other input parameters are saved as attributes of the result as well.
#' @examples
#' data( rrain )
#' ra <- list( rrain[2,4,,], rrain[3,9,,] )
#' ra <- prepare_sad( ra, rsm=0, Nx=256, boundaries="mirror", log=FALSE )
#' plot(ra)
#' @export
prepare_sad <- function( x, xmin=0.1, log=TRUE, rsm=0, Nx=NULL, Ny=NULL, boundaries="pad" ){
    x <- as.sadforecast( x )
    old_att <- attributes( x )
    
    if( is.null(Nx) ) Nx <- nrow( x[[1]] )
    if( is.null(Ny) ) Ny <- ncol( x[[1]] )
    
    # remove small values
    x <- lapply(x, function(x) x*1*(x>=xmin) )
    # smooth the boundaries
    x <- lapply(x, dualtrees::smooth_borders, r=rsm)
    
    # log or no?
    if( log & !old_att$log ) logfun <- function(x) log2( x + xmin ) else logfun <- identity
    # select boundary type
    bound <- switch( boundaries, 
                     mirror=dualtrees::put_in_mirror, 
                     pad=dualtrees::pad, 
                     periodic=dualtrees::period_bc,
                     stop( "unknown boundaries, you can use 'mirror', 'pad' or 'periodic'" )
                    )
    # remember where we were in the original field
    bc <- bound( logfun(x[[1]]), N=Nx, Ny=Ny )
    # apply boundaries and maybe logarithm
    x  <- lapply( x, function(x) bound( logfun(x), N=Nx, Ny=Ny )$res )
    
    # set attributes
    class( x ) <- "sadforecast"
    attr( x, "px" ) <- bc$px
    attr( x, "py" ) <- bc$py
    attr( x, "xmin" ) <- max( xmin, old_att$xmin )
    if(log) attr( x, "log" )  <- TRUE
    attr( x, "rsm" ) <- rsm
    attr( x, "boundaries" ) <- boundaries
    return( x )
}

#' rain color scale
#'
#' eight shades of blue used in \code{plot.sadforecast}
#' @export
raincols <- c( "#FCFCFC", "#CAD5EB", "#97AFD9", "#608AC8", "#3267A5", "#214671", "#112641", "#040404" )

angdif <- function(x,y){
    d <- x - y
    if( d > 90 )  d <- d - 180
    if( d < -90 ) d <- d + 180
    return(d)
}


checkJ <- function(J, Nx, Ny){
    if( J < 1 ) stop( "you have to use at least one scale (J>=1)" )
    N    <- min( Nx, Ny )
    Jmax <- floor( log2( N ) )
    if( J > Jmax ) stop( paste0( "J must be no greater than ", Jmax ) )
}
