#' correct structure errors
#'
#' use the inverse 'dtcwt' to correct errors in scale, anisotropy and direction
#' @param x a list of equally sized matrices, the first element is assumed to be the observation
#' @param xmin values smaller than \code{xmin} are set to zero
#' @param log  logical, do you want to log-transfrom the data? (recommended for precipitation)
#' @param rsm  number of pixels which are linearly smoothed at the edge
#' @param Nx   size to which the data is extended in x-direction, has to be a whole power of 2
#' @param Ny   size to which the data is extended in y-direction, has to be a whole power of 2
#' @param J    largest scale considered
#' @param boundaries how to handle the boundary conditions, either "pad", "mirror" or "periodic"
#' @param direction if \code{TRUE}, scale and direction are corrected, otherwise only scale
#' @return an object of class \code{sadforecast}
#' @details The algorithm performs the following steps:
#' \enumerate{
#'    \item remove values below \code{xmin}
#'    \item if \code{log=TRUE} log-transform all fields
#'    \item  set all fields to zero mean, unit variance
#'    \item apply \code{dtcwt} to all fields
#'    \item loop over forecasts and scales. If \code{direction=TRUE} loop over the six directions. Multiply forecast energy at each location by the ratio of total observed energy to total forecast energy at that scale (and possibly direction)
#'    \item apply \code{idtcwt} to all forecasts
#'    \item reset means and variance of the forecasts to their original values
#'    \item if \code{log=TRUE} invert the log-transform 
#'    \item return the list of corrected fields
#'}
#' @examples
#' data(rrain)
#' ra   <- as.sadforecast( list( rrain[2,1,,], rrain[3,1,,], rrain[3,2,,], rrain[3,3,,] ) )
#' ra_c <- sadcorrect( ra, rsm=10 )
#' plot(ra_c)
#' @export
sadcorrect <- function( x, xmin=0.1, log=TRUE, rsm=0, Nx=NULL, Ny=NULL, J=NULL, boundaries="pad", direction=TRUE ){
    
    # check the input
    x   <- as.sadforecast( x )
    n   <- length( x )
    nam <- names( x )
    
    # get defaults for the dimensions etc
    if( is.null( Nx ) ) Nx <- 2**ceiling( log2( max( dim( x[[1]] ) )  ) )
    if( is.null( Ny ) ) Ny <- Nx
    if( is.null( J ) ) J <- log2( min(Nx,Ny) ) - 3
    
    # check that J is valid
    checkJ( J, Nx, Ny )
    
    # remember the original observation
    obs0 <- x[[1]]
    
    # handle thresholding, log, boundaries
    x <- prepare_sad( x, xmin, log, rsm, Nx, Ny, boundaries )
    
    # standardize the margins
    xm    <- lapply( x, mean )
    xsd   <- lapply( x, stats::sd )
    for( i in 1:n ) x[[i]] <- ( x[[i]] - xm[[i]] ) / xsd[[i]] 
    
    # apply wavelet transform
    dt    <- lapply( x, dualtrees::dtcwt, dec=TRUE, 
                        fb1 = dualtrees::near_sym_b_bp, 
                        fb2 = dualtrees::qshift_b_bp )
    

    res <- list( obs0 )
    
    # correct each forecast
    dto <- dt[[1]]
    for( i in 2:n ){ # loop over forecasts
        dtf <- dt[[i]]
        for( j in 1:J ){ # loop over all scales
            if( direction ){ # correct each direction individually
                for( d in 1:6 ){
                    fac       <- sum(Mod(dto[[j]][,,d])**2) / 
                                 sum(Mod(dtf[[j]][,,d])**2)
                    dtf[[j]][,,d] <- dtf[[j]][,,d]*fac
                } 
            }else{ # correct all directions together
                fac       <- sum(Mod(dto[[j]])**2) / 
                                 sum(Mod(dtf[[j]])**2)
                    dtf[[j]] <- dtf[[j]]*fac
            }
            
        } 
        # transform back
        fbc <- idtcwt( dtf, fb1 = dualtrees::near_sym_b_bp, 
                            fb2 = dualtrees::qshift_b_bp 
                     )
        # restore margins
        fbc <- ( fbc - mean(fbc) )/ stats::sd(fbc) * xsd[[i]] + xm[[i]]
        # cut out original domain
        fbc <- with( attributes(x), fbc[px,py] )
        # invert logarithm (maybe)
        if( log ) fbc <- 2**fbc - xmin
        res[[i]] <- fbc
    }

    
    return( structure(res, class="sadforecast" ) )
}
