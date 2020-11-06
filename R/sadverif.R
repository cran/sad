#' dual-tree verification
#'
#' verify the scale, anisotropy and direction of a number of forecasts
#' @param x  a list of equally sized matrices, the first element is assumed to be the observation
#' @param dec  logical, do you want to use the decimated transform
#' @param xmin values smaller than \code{xmin} are set to zero
#' @param log  logical, do you want to log-transfrom the data? (recommended for precipitation)
#' @param a    relative weight of directional errors compared to scale errors in \code{semdd}
#' @param nbr  number of breaks for the scale histograms, has no effect if \code{dec=TRUE}
#' @param rsm  number of pixels which are linearly smoothed at the edge
#' @param Nx   size to which the data is extended in x-direction
#' @param Ny   size to which the data is extended in y-direction
#' @param J    largest scale considered
#' @param boundaries how to handle the boundary conditions, either "pad", "mirror" or "periodic"
#' @param return_specs if \code{TRUE}, the spatial mean spectra are returned as well
#' @param object object of class sadverif
#' @param ... further arguments, currently ignored.
#' @details each element of x is transformed via \code{dtcwt} from the 'dualtrees' package. Scores and centres based on the mean spectra are calculated. If \code{dec=FALSE}, scale histograms and the corresponding score \code{hemd} are calcualted as well.
#' @return an object of class \code{sadverif}, containing the following elements
#' \describe{
#'    \item{settings}{ a dataframe containing the parameters that were originally passed to dtverif }
#'    \item{centres}{ a matrix cotaining the anisotropy \code{rho}, angle \code{phi} and central scale \code{z} derived from the mean spectra. Rain area and sum are included as well.}
#'    \item{detscores}{ a matrix containing the differences in centre components, the direction/anisotropy score \code{dxy}, the emd between direction-averaged spectra (\code{semd}) and the emd between the directional spectra (\code{semdd}). If \code{dec=FALSE}, the emd between the scale histograms, hemd, is included as well. }
#'    \item{time}{ the time the calculation took in seconds }
#' }
#' if there is more than one forecast, the ensemble scores SpEn and (if available), hemd are computed as well, treating all forecasts as members of the ensemble to be verified.
#' @references Selesnick, I.W., R.G. Baraniuk, and N.C. Kingsbury (2005) <doi:10.1109/MSP.2005.1550194>
#' Buschow et al. (2019) <doi:10.5194/gmd-12-3401-2019>
#' Buschow and Friederichs (2020) <doi:10.5194/ascmo-6-13-2020>
#' @examples
#' oldpar <- par(no.readonly=TRUE)
#' on.exit(par(oldpar))
#' data(rrain)
#' ra <- as.sadforecast( list( rrain[1,1,,], rrain[1,2,,], rrain[2,1,,], rrain[3,1,,] ) )
#' plot(ra)
#' verif <- sadverif( ra, log=FALSE, xmin=0 )
#' summary(verif)
#' par( mfrow=c(2,2) )
#' plot( verif )
#' @import dualtrees
#' @name sadverif
NULL

#' @rdname sadverif
#' @export
sadverif <- function( x, dec=TRUE, xmin=0.1, log=TRUE, a=1, nbr=33, rsm=0, Nx=NULL, Ny=NULL, J=NULL, boundaries="pad", return_specs=FALSE  ){
    
    # start timing
    t0 <- proc.time()
    
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
    
    # get area and sum
    xar <- lapply( x, function(x) mean( x>xmin ) )
    xsu <- lapply( x, function(x) sum( x*1*(x>=xmin) ) )
    
    # handle thresholding, log, boundaries
    x <- prepare_sad( x, xmin, log, rsm, Nx, Ny, boundaries )
    
    # do the dual-tree transform
    dt    <- lapply( x, dualtrees::dtcwt, J=J, dec=dec, 
                        fb1 = dualtrees::near_sym_b_bp, 
                        fb2 = dualtrees::qshift_b_bp )
    
    if( !dec ){
        # cut out original domain
        dt <- with( attributes(x), lapply(dt, function(y) return(y[,px,py,]) )  )
        # correct bias
        dt <- lapply(dt, FUN=get_en, correct = "b_bp", N = 2^(J + 3) )
        # remove negative energy
        dt <- lapply( dt, function(x) x*1*(x>0) )
        # prepare for the histograms
        br <- seq( 1, J, , nbr )
        mi <- ( br[-1] + br[-nbr] ) / 2
        # calculate the histograms
        hists <- array( dim=c(n, nbr-1), dimnames=list(nam, mi) )
        for( i in 1:n ){
            # get the mask of non-zero rain
            th <- if(log) log2(xmin) else xmin
            ma <- x[[i]] <= th
            # map of central scales
            zc <- with( attributes(x), dt2cen( dt[[i]], mask=ma[px,py] )[ ,,3 ] )
            hists[ i, ] <- graphics::hist( zc, breaks=br, plot=FALSE )$density
        }
    }
    
    # average over space
    dt   <- lapply( dt, dualtrees::dtmean )      
    # get the mean centre
    dtcen <- lapply( dt, dualtrees::dt2cen )
    # get the x- and y-component
    xy    <- lapply( dtcen, function(x) dualtrees::cen2xy(x)[1:2] )

    # prepare the results
    cennames   <- c( "rho", "phi", "z", "area", "sum" )
    scorenames <- c( paste0("d", cennames), "dxy", "semd", "semdd" )
    if( !dec ) scorenames <- c( scorenames, "hemd" )
    
    cen       <- array( dim=c(n, length(cennames) ), 
                        dimnames=list(nam, cennames ) ) 
    detscores <- array( dim=c(n-1, length(scorenames) ), 
                        dimnames=list(nam[-1], scorenames) )
    
    # calculate the centres and deterministic scores
    for( i in 1:n ){
        # get the centres
        cen[ i, ] <- c( dtcen[[i]], xar[[i]], xsu[[i]] )
        j <- i-1
        if( j > 0 & xar[[i]] > 0 ){
            # differences in centre
            detscores[ j, 1:5 ] <- cen[i,] - cen[1,]
            # special treatment for the angle
            detscores[ j, 2 ]   <- angdif( cen[i,2], cen[1,2] )
            
            # combined anisotropy / direction score
            detscores[ j, 6 ]   <- sqrt( sum( ( xy[[i]] - xy[[1]] )**2 ) )
            
            # direction-less and directional emd
            detscores[ j, 7 ] <- semd( dt[[i]], dt[[1]], dir=FALSE )
            detscores[ j, 8 ] <- semd( dt[[i]], dt[[1]], dir=TRUE, a=a )
            
            # histogram emd, if possible
            if( !dec ){
                detscores[ j, 9 ] <- hemd( hists[i,], hists[1,]  )
            }
        }
    }
                    
    # maybe calculate the probabilistic scores as well
    if( n > 2 ){
        ensscores <- list()
        # observed mean spectrum
        dto <- c( dt[[1]] )
        
        #re-format predicted mean spectra
        dtf <- c( dt[[2]] )
        for( i in 3:n ) dtf <- cbind( dtf, c(dt[[i]]) )
        # get the energy score
        ensscores$en <- energyscore( dto, dtf )
        if( !dec ){
            ensscores$he <- hemd( hists[1,], colMeans( hists[-1,] ), mids=mi )
        }
    }
                    
    # save the settings
    settings <- data.frame( dec=dec, xmin=xmin, log=log, a=a, nbr=nbr, rsm=rsm, Nx=Nx, Ny=Ny, J=J, boundaries=boundaries )
                    
    # create the results
    res <- list( settings=settings, centres=cen, detscores=detscores )
    if( !dec )  res$hists <- hists
    if( n > 2 ) res$ensscores <- ensscores
    if( return_specs ) res$dt <- dt
    
    # add the time
    res$time <- unname( proc.time() - t0  )[3]
                    
    class( res ) <- "sadverif"
    
    return( res )
}


#' @rdname sadverif
#' @export
plot.sadverif <- function(x,  ...){
    # define colors
    cols <- 1:nrow( x$cen )
    
    # prepare the hexagon for our plot
    phi <- pi*60/180*0:5
    xx   <- cos( phi )
    yy   <- sin( phi )
    # get the x and y component you want to plot
    xy  <- apply(x$cen[,1:3], 1, cen2xy)
    
    # plot the centre in the x-y-plane
    graphics::plot( xx, yy, type="n" )
    graphics::abline( h=0, v=0, lty=3 )
    graphics::points( c(xx, xx[1]), c(yy, yy[1]), type="b", pch=NA, cex=1.5 )
    graphics::text( xx, yy, labels=paste0(seq(15,,30,6),"\u00B0") )
    graphics::text( xy[ 1, ], xy[2,], col=cols, labels=rownames(x$cen) )
    graphics::title( main="anisotropy and direction" )
    
    
    # dotchart for the central scales
    graphics::dotchart( xy[3,], xlim=c(1,x$settings$J), col=cols )
    graphics::abline(v=xy[3,1], lty=2)
    graphics::title( main="central scale" )
    
    # get the scores
    sc  <- x$detscores 
    
    # find the pareto-set in the dz-dxy-plane
    winners <- getpareto( abs( sc[ ,c(3,6),drop=FALSE ] ) )
    fonts   <- rep( 1, nrow(sc) )
    fonts[winners] <- 2
    lab <- rownames(sc)
    lab[winners] <- paste0( lab[winners], "*" ) 
    # plot dz against dxy
    if( any( !is.na( sc[,3] + sc[,6] ) ) ){
        graphics::plot( abs(sc[,3]), sc[,6], type="n", xlab="|dz|", ylab="dxy",
                        xlim=range(abs(sc[,3]), na.rm=TRUE)+c(-.1,.1), 
                        ylim=range(abs(sc[,6]), na.rm=TRUE)+c(-.1,.1)
                       )
        graphics::grid()
        graphics::text( abs(sc[,3]), sc[,6], col=cols[-1], labels=lab, font=fonts )
        graphics::title( main="scale- and directional error" )
    }
    
    # plot the histograms if we have them
    if( !is.null( x$hists ) ){    
        hists <- x$hists
        mi    <- colnames(hists)
        graphics::plot( mi, hists[1,], col=cols[1], type="l", ylim=range(hists),
                        xlab="scale", ylab="density" )
        for( i in 2:nrow(hists) ) graphics::points( mi, hists[i,], col=cols[i], type="l" )
        graphics::legend( x="topleft", lty=1, col=cols, legend=rownames( x$cen ), bty="n" )
        graphics::title( main="histograms of central scales" )
    }
}

#' @rdname sadverif
#' @export
summary.sadverif <- function(object, ...){
    with(object,{
        nam <- rownames(detscores)
        cat( paste0("\nVerification of ", paste(nam[-1], collapse=", "), " against ", nam[1], ".\n\n" ) )
        cat( "Deterministic verification of the individual forecasts:\n")
        print( detscores, digits=2 )
        if( nrow(detscores) > 1 ){
            winners <- getpareto( abs(detscores[,c(3,6)]) )
            cat( "\nMembers of the pareto set (based on dxy and dz):\n" )
            cat( nam[winners], "\n" )
            cat( "\nEnsemble scores:\n")
            print( unlist( ensscores), digits=2  )
        }
        cat( "\n settings:\n" )
        print(settings) 
        cat( "\n took ", round(time, 2), " seconds.\n" )
    })
}
