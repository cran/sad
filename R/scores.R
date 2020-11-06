#' spectral emd 
#'
#' earth mover's distance between dual-tree wavelet spectra
#' @param dt1,dt2 forecast and observed spectrum
#' @param a ratio between scale- and directional component
#' @param dir whether or not to include direction information
#' @return a single value, the emd. If \code{dir=FALSE}, the value is signed, indicating the direction of the scale error.
#' @export
semd <- function( dt1, dt2, a=1, dir=TRUE ){
    if( sum(dt1)==0 | sum(dt2)==0 ) return( NA )
    J <- nrow( dt1 )
    
    if( dir ){
    
    dt1 <- c( dt1 )/sum( dt1 )
    dt2 <- c( dt2 )/sum( dt2 )
    
    d   <- rep( 1:6, each=J )
    
    x   <- a*cos( (d-1)*60*pi/180 )*( J-1 )/2
    y   <- a*sin( (d-1)*60*pi/180 )*( J-1 )/2
    z   <- rep( 1:J, 6 )
    
    dt1 <- cbind( dt1, x, y, z )
    dt2 <- cbind( dt2, x, y, z )
        
    }else{
        z1  <- dualtrees::dt2cen(dt1)[,,3]
        z2  <- dualtrees::dt2cen(dt2)[,,3]
        dt1 <- rowMeans( dt1 )
        dt2 <- rowMeans( dt2 )
        dt1 <- cbind( dt1/sum(dt1), 1:J )
        dt2 <- cbind( dt2/sum(dt2), 1:J )
        
    }
    res <- emdist::emd( dt1, dt2 )
    if(!dir) res <- res*sign( z1 - z2 )
    
    return( res )
}

#' histogram emd
#'
#' Earth Mover's Distance between two histograms, given as vectors
#' @param h1,h2 vectors of non-negtaive numbers representing two histograms
#' @param mids the bin mids corresponding to the histograms. Can also be given via the names of \code{h1}.
#' @return the value of the EMD
#' @export
hemd <- function( h1, h2, mids=NULL ){
    if( is.null(mids) ) mids <- as.numeric( names(h1) )
    if( is.null(mids) ) stop( "give me the mids" )
    
    z1 <- sum( h1*mids ) / sum(h1)
    z2 <- sum( h2*mids ) / sum(h2)
    
    res <- emdist::emd( cbind(h1,mids), cbind(h2,mids) )*sign( z1-z2 )
    return(res)
}


energyscore <- function( obs, forc ){
    nx <- nrow(forc)
    nm <- ncol(forc)
    
    en1 <- mean(  sqrt( colSums( ( forc - obs )**2) ) )
    en2 <- 0
    for( i in 1:nm ) en2 <- en2 + sum(  sqrt( colSums( ( forc[,-i,drop=FALSE] - forc[,i] )**2) ) )/( nm**2 )
    res <- en1 - en2 / 2
    return(res)
}

#' Find the pareto set 
#' 
#' Determine the set of pareto optimal forecasts in a matrix of scores
#' @param scores a matrix of negatively oriented scores where the rows correspond to different forecasts and the columns denote different scores.
#' @return a vector of indices indicating all members of the pareto set.
#' @details The Pareto set contains all those forecasts for which no other forecast is better in every respect. In this function, we assume that all scores are negatively oriented, "better" therefore means lower values.
#' @note This function becomes very memory hungry if you have more than 1000 forecasts, be careful.
#' @export
getpareto <- function( scores ){
    scores <- scores[ !is.na( rowSums( scores ) ), , drop=FALSE ]
    if( nrow(scores) < 2 ){
        winners <- NULL
    }else{
        ns <- ncol(scores)
        isworse <- TRUE
        for( i in 1:ns ) isworse <- isworse & outer( scores[,i,drop=FALSE], 
                                                     scores[,i,drop=FALSE], 
                                                     function(x,y) x > y 
                                                     )
        winners <- which( rowSums( isworse )  < 1 )        
    }
    return( winners )
}
