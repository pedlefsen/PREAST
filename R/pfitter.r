
# poisson fitter derived from code downloaded from LANL by Paul Edlefsen.


# Paul Edlefsen downloaded this on Oct 4, 2015 from
# "http://www.hiv.lanl.gov/tmp/POISSON_FITTER/pp429jC8Td/PFitter.R".
# (followed link 'Download R Script' from PoissonFitter results page
# http://www.hiv.lanl.gov/cgi-bin/POISSON_FITTER/v2/pfitter.cgi)

#########################################################################
##########       PURPOSE: fit Poisson, calculate lambdas,      ##########
##########    U-stat standard deviation, and goodness of fit   ##########
##########      written by EEG, last modified on 5/26/15       ##########
##########          send questions to egiorgi@lanl.gov         ##########

# INPUT: 
#  pairwise hamming distances file, mutation rate, length of a sequence
#  distances file: tab-delimited 3-column file. seqname1(1st col), 
#                  seqname2(2nd) and distance between seq1 and seq2(3rd).
#                  based on large-scale formatted sequnce input, which 
#                  means every unique sequence is represented only once 
#                  and a seqname should end with _nnn wher nnn is the 
#                  the multiplicity of such sequence in the alignment.

#  example:  R CMD BATCH '--vanilla --args sample.dist  2.16e-05 2517' this_script 

# OUTPUT:
#  2 files, one with lambdas, maxhd, jackknife standard deviation 
#  and estimated days, goodness of fit p-values, and the other with the 
#  star-phylogeny estimated an dobserved numbers (if they coincide, you have a star)
#########################################################################


if (!suppressMessages(require("pacman"))) install.packages("pacman")
pacman::p_load(tools, assertthat)


pfitter <- function(dlist, epsilon, nbases) {
    # sample is a character string label used to label plot and output files.

    ### FUNCTIONS ###
    iseven <- function(c) {
	c1 <- c/2-as.integer(c/2)
	if(c1==0){ ev <- TRUE } else { ev <- FALSE }
	return(ev)
    }
    
    phi <- sqrt(1+4/3)
    days <- function(l,nb,epsilon) 1.5*((phi)/(1+phi))*(l/(epsilon*nb) - (1-phi)/(phi^2))

    ### calc HD with consensus
    d1 <- dlist[-which(dlist[,1]==dlist[1,1]),]	 # PAUL moved this here, to help handle singleton clusters.
    seqnames <- unique(c( d1[,1], d1[,2] ))          # PAUL moved this here, to help handle singleton clusters.
    d0 <- dlist[which(dlist[,1]==dlist[1,1] & (dlist[,2] %in% seqnames)),] # PAUL changed this, to handle singleton clusters.
    mult0 <- as.numeric(sub('.+_(\\d+)$', '\\1', d0[,2]))
    nseq <- sum(mult0)
    yvec0 <- rep(0, (1+max(d0[,3])))
    for(i in 1:(1+max(d0[,3]))){ yvec0[i] <- sum(mult0[which(d0[,3]==(i-1))]) }

    nl0 <- length(yvec0);
    clambda <- sum((1:(nl0-1))*yvec0[-1])/sum(yvec0) #### THIS IS THE LAMBDA THAT FITS THE CONSENSUS ONLY DISTRIBUTION

    ### calc intersequence HD
    d1 <- dlist[-which(dlist[,1]==dlist[1,1]),]	
    yvec <- rep(0, (1+max(d1[,3])))
    seqnames <- unique(c( d1[,1], d1[,2] ))
    for(i in 1:length(seqnames)) {
        tmp <- d1[which(d1[,1]==seqnames[i]),,drop = FALSE]
        if( nrow( tmp ) == 0 ) {
            next;
        }
	m0 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[1,1]))
	yvec[1] <- yvec[1] + 0.5*m0*(m0-1) ## 0 bin
	for(j in 1:dim(tmp)[1]){
            m1 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[j,2]))
            val <- tmp[j,3]
            yvec[val+1] <- yvec[val+1] + m0*m1
	}
    }

    ### Fitting

    nl <- length(yvec)
    lambda <- sum((1:(nl-1))*yvec[-1])/sum(yvec)
    estdays <- days(lambda, nbases, epsilon)
    
    #### U STAT ESTIMATE OF ST DEV
    #### FORMULAE
    #### Var(HD) = (N(N-1)/2)^(-1) (2(N-2)sigma1^2 + sigma2^2)
    #### sigma1^2 = (N(N-1)(N-2)/3 -1)^(-1) sum_{i<j<l} ((Dij-mu)(Dil-mu)+(Dij-mu)(Djl-mu))
    #### sigma2^2 = (N(N-1)/2-1)^(-1) sum_{i<j} (Dij-mu)^2

    ### construct a matrix of Dij's
    ### number of unique sequences
    nuni <- dim(d0)[1]
    TX <- matrix(rep(NA,nuni^2), ncol=nuni)
    #assert_that(length(rownames(TX)) == length(seqnames))
    rownames( TX ) <- seqnames;
    colnames( TX ) <- seqnames;
    for(i in 1:(dim(d0)[1]-1)){
	useq <- d0[i,2]
	TX[d1[which(d1[,1]==useq),2],i] <- d1[which(d1[,1]==useq),3];
    }

    sigma1 <- 0
    sigma2 <- 0
    muhat <- 0
    denmu <- (sum( !is.na( TX )))^(-1)
    ## TODO: Figure out what (if any) is the right fix to the below to handle sparse distances
    den1 <- 12*(nseq*(nseq-1)*(nseq-2)*(nseq-3))^(-1)  
    den2 <- den1/4

    for(n in 1:(nuni-1)){
        for(m in (n+1):nuni){
            if( !is.na( TX[ m, n ] ) ) {
                muhat <- muhat + mult0[n]*mult0[m]*denmu*TX[m,n]
            }
        }
    }

    for(n in 1:nuni){
	dnn <- 0
	sigma1 <- sigma1 + choose(mult0[n],3)*den1*2*(dnn-muhat)^2 
	sigma2 <- sigma2 + choose(mult0[n],2)*den2*(dnn-muhat)^2
	if(n != nuni){
            for(m in (n+1):nuni){
                dnm <- TX[m,n]
                if( !is.na( dnm ) ) {
                    dmm <- 0
                    sigma2 <- sigma2 + mult0[n]*mult0[m]*(dnm - muhat)^2
                    sigma1 <- sigma1 + (2/3)*choose(mult0[n],2)*mult0[m]*(dnm-muhat)*(dnm+2*dnn-3*muhat)
                    sigma1 <- sigma1 + (2/3)*mult0[n]*choose(mult0[m],2)*(dnm-muhat)*(dnm+2*dmm-3*muhat)
                    if(m != nuni){
                        for(l in (m+1):nuni){
                            dnl <- TX[l,n]
                            dlm <- TX[l,m]
                            if( !is.na( dnl ) && !is.na( dlm ) ) {
                                sigma1 <- sigma1 + (2/3)*mult0[n]*mult0[m]*mult0[l]*((dnm-muhat)*(dnl-muhat)+(dnm-muhat)*(dlm-muhat)+(dnl-muhat)*(dlm-muhat))
                            }
                        }
                    }
                } # End if !is.na( dnm )
            } # End for( m )
 	}
    }

    ## varhd <- sqrt(denmu*(2*(nseq-2)*sigma1 + sigma2))
    A <- 8/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
    B <- 4/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
    newvarhd <- sqrt(A*sigma1 + B*sigma2)
    
    upplim <- days(lambda + 1.96*newvarhd, nbases, epsilon)
    lowlim <- days(lambda - 1.96*newvarhd, nbases, epsilon)
    uppdays <- round(upplim)
    lowdays <- round(lowlim)
    
    ### output figures
    dvec1 <- 0
    for(i in 1:length(yvec)){ dvec1 <- c(dvec1, rep((i-1),yvec[i])) }
    dvec1 <- dvec1[-1]

    meanhd <- mean(dvec1)
    maxhd <- max(dvec1)

    # construction of figures removed.

    #### FIT THE CONSENSUS ONLY HD DISTRIBUTION
    
    # construction of convolution figures removed.
    

    ### CONSTRUCT SIGMA_ij MATRIX THEN INVERT IT
    #pk <- function(x) ((nseq^2)*(2^x)*exp(-2*clambda)*(clambda^x))/factorial(x)
    pk <- function(x) (exp( ( (log(nseq)*2)+(log(2)*x)+(-2*clambda)+log(clambda)*x)-lfactorial(x) ) ) # PAUL CHANGED log( clambda^x ) to log(clambda)*x
    mui <- function(x) nseq*dpois(x, lambda=clambda)
    SIGMA.DIM.MAX <- 170; # Beyond this value, factorial stops working in R.
    sigma.dim <- min( SIGMA.DIM.MAX, (2*nl0) );
    eyvec <- 0.5*pk(0:(sigma.dim-1))
    eyvec[ !is.finite( eyvec ) ] <- 0; # Paul added this to avoid errors.  factorial(x) sometimes refuses to work (if x is too large).

    if (lambda!=0) {
	sigmaij <- matrix(nrow=sigma.dim, ncol=sigma.dim)
	coeff <- (nseq^3)*exp(-3*clambda)
	
	#### RICORDATI!!!! EYVEC[K] == E(Y_{K-1}) !!!!!!
	
	for(k in 0:(sigma.dim-1)){  
            
            for(l in 0:(sigma.dim-1)){   
                
                if(k>=l){
                    c1 <- exp( ((log( clambda )* k) - lfactorial(k)) + log( sum( exp( lchoose(k,l:0)+((log(clambda)*(0:l)) - lfactorial(0:l)) ) ) ) );
                    stopifnot( is.finite( c1 ) );
                    c2 <- exp( ((log( clambda )* l) - lfactorial(l)) + log( sum(exp( lchoose(l,l:0)+((log( clambda )*((k-l):k))-lfactorial((k-l):k))) ) ) );
                    stopifnot( is.finite( c2 ) );
                }
                if(k<l){
                    c1 <- exp( ((log( clambda )* l) - lfactorial(l)) + log( sum( exp( lchoose(l,k:0)+((log( clambda )*(0:k)) - lfactorial(0:k))) ) ) );
                    stopifnot( is.finite( c1 ) );
                    c2 <- exp( ((log( clambda )* k) - lfactorial(k)) + log( sum(exp( lchoose(k,k:0)+((log(clambda)*((l-k):l)) - lfactorial((l-k):l))) ) ) );
                    stopifnot( is.finite( c2 ) );
                }
                if( is.na( c1 ) ) {
                    c1 <- 0;
                }
                if( is.na( c2 ) ) {
                    c2 <- 0;
                }
                sigmaij[k+1,l+1] <- 0.5*coeff*(c1+c2);
                if(k==l){
                    if( iseven(k) ) {
                        sigmaij[k+1,l+1] <- sigmaij[k+1,l+1] - (0.25)*mui(k/2);
                    } else {
                        sigmaij[k+1,l+1] <- sigmaij[k+1,l+1] + (0.5)*pk(k) ;
                    }
                }
                stopifnot( is.finite( sigmaij[k+1,l+1] ) );
            }
	}
        
        #print( "about to run La.svd" );
	sdec <- La.svd(sigmaij)
        #print( "ran La.svd" );
	diag <- ifelse(sdec$d>1e-4,sdec$d,0)
	diagmat <- matrix(rep(0,sigma.dim^2), ncol=sigma.dim)
	for(ii in 1:sigma.dim){diagmat[ii,ii]<-ifelse(diag[ii]==0,0,1/diag[ii])}
	sigmainv <- sdec$u%*%diagmat%*%sdec$vt
	
	h <- hist(dvec1[ dvec1 < sigma.dim ], breaks=seq(-1,(sigma.dim-1),1), plot=FALSE)
	yvec <- h$counts
        
	if (length(yvec)<sigma.dim) { 
            ccvv <- sigma.dim - length(yvec) 
            yvec <- c(yvec, rep(0,ccvv)) 
	}
	
	chisq <- t(abs(yvec-eyvec))%*%sigmainv%*%(abs(yvec-eyvec))
        pval <- ifelse(chisq<0,2e-16,1-pchisq(chisq,df=nl0-1))
	if(pval==0){ pval <- 2e-16 }
	if(chisq<0){ chisq <- qchisq(p = 2e-16,df=nl0-1, lower.tail = FALSE ) } # Paul changed from "NA" - this is a high value, reflective of the low p-value.
    } else { 
	chisq <- NA
	nl0 <- NA # Paul fixed to "nl0" from "nl"
	pval <- 0
    }
  

    return(list(lambda=lambda,
                stdev=newvarhd,
                nseq=nseq,
                nbases=nbases,
                meanhd=meanhd,
                maxhd=maxhd,
                days=paste(round(estdays), " (", lowdays, ", ", uppdays, ")", sep=""),
                chi2=as.numeric(chisq),
                df=nl0-1,
                goodness.pval=as.numeric(pval)
                ))
}


# add a consensus sequence.
# Useful for formatting a fasta file in a form acceptable
# to the PFitter tool at http://www.hiv.lanl.gov/content/sequence/POISSON_FITTER/pfitter.html
#
# Typical usage:
# 	read.dna( fasta.file, format = "fasta" ) %>%
#     	    add.consensus( label=basename(file_path_sans_ext(fasta.file))) %>%
#     	    write.dna('output2.fa', format = "fasta", nbcol=-1, colsep = "",indent = 0,blocksep = 0)
#
add.consensus <- function ( in.fasta, label='sample' ) {
    # in.fasta <- read.dna( fasta.file, format = "fasta" );
    
    # Add the consensus.
    .consensus.mat <- matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 );
    consensus <- as.DNAbin( .consensus.mat );
    rownames( consensus ) <- paste0( label, ".CONSENSUS" );
    
    rbind( consensus, in.fasta );
}

# calculate Hamming distances and prepare a dataframe for pfitter.
# Written by Paul Edlefsen.  Adapted for use within R by cswarth@gmail.com
#
# typical usage,
#     seq <- read.dna(path, format = "fasta")
#     d <- prep.distances(seq)
#     r <- pfitter(d$distances, 2.16e-05, d$seq.length)

prep.distances <- function ( in.fasta, include.gaps.in.Hamming=FALSE ) {
    stopifnot(inherits(in.fasta, 'DNAbin'))
    
    # Add the consensus.
    .consensus.mat <- matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 );
    consensus <- as.DNAbin( .consensus.mat );
    rownames( consensus ) <- paste( "CONSENSUS" );
    
    fasta.with.consensus <- rbind( consensus, in.fasta );
    
    # Remove any columns with a consensus that is a gap, which means
    # that over half of seqs have gaps.  This needs to be removed
    # because it becomes not sensible to consider poisson rates of
    # insertions.  We do however consider rates of deletions, treating
    # them as point mutations (by including the indel counts in the
    # Hamming distance calculation).
    fasta.with.consensus <- fasta.with.consensus[ , .consensus.mat[ 1, ] != "-" ];
          
    # The pairwise.deletion = TRUE argument is necessary so that columns with any gaps are not removed.
    # The optional second call adds in the count of sites with gap differences
    # (gap in one sequence but not the other), which are not included
    # in the first call's "raw" count. Adding them back in, we are
    # effectively treating those as differences.
    fasta.with.consensus.dist <- dist.dna( fasta.with.consensus, model = "N", pairwise.deletion = TRUE );
    fasta.with.consensus.dist[ is.nan( fasta.with.consensus.dist ) ] <- 0;
    if( include.gaps.in.Hamming ) {
        fasta.with.consensus.dist <- fasta.with.consensus.dist + dist.dna( fasta.with.consensus, model = "indel", pairwise.deletion = TRUE );
    }
    
    if( any( is.null( fasta.with.consensus.dist ) ) || any( is.na( fasta.with.consensus.dist ) ) || any( !is.finite( fasta.with.consensus.dist ) ) ) {
        ## TODO: REMOVE
        warning( "UH OH got illegal distance value" );
        print( "UH OH got illegal distance value" );
        print( fasta.with.consensus.dist );
    }
    
    pairwise.distances.as.matrix <- as.matrix( fasta.with.consensus.dist );

    # prepare a dataframe in a format accesptable to the PFitter routine.
    #  seqname1(1st col), seqname2(2nd) and distance between seq1 and seq2(3rd).
    #                  based on large-scale formatted sequnce input, which 
    #                  means every unique sequence is represented only once 
    #                  and a seqname should end with _nnn wher nnn is the 
    #                  the multiplicity of such sequence in the alignment
    pairwise.distances.as.matrix.flat <- matrix( "", nrow = ( ( nrow( pairwise.distances.as.matrix ) * ( ncol( pairwise.distances.as.matrix ) - 1 ) ) / 2 ), ncol = 3 );
    line.i <- 1;
    for( row.i in 1:( nrow( pairwise.distances.as.matrix ) - 1 ) ) {
        for( col.i in ( row.i + 1 ):ncol( pairwise.distances.as.matrix ) ) {
            if( row.i == 1 ) { # consensus, no _1 (multiplicity of the observed sequence)
                pairwise.distances.as.matrix.flat[ line.i, ] <-
                    c( rownames( pairwise.distances.as.matrix )[ row.i ], paste( colnames( pairwise.distances.as.matrix )[ col.i ], "_1", sep = "" ), pairwise.distances.as.matrix[ row.i, col.i ] );
            } else {
                pairwise.distances.as.matrix.flat[ line.i, ] <-
                    c( paste( rownames( pairwise.distances.as.matrix )[ row.i ], "_1", sep = "" ), paste( colnames( pairwise.distances.as.matrix )[ col.i ], "_1", sep = "" ), pairwise.distances.as.matrix[ row.i, col.i ] );
            }
            line.i <- line.i + 1;
        }
    }

    # convert matrix to data frame
    pairwise.distances.flat <- data.frame(pairwise.distances.as.matrix.flat, stringsAsFactors=F)
    colnames(pairwise.distances.flat) <- c('seqname1', 'seqname2','distance')

    # fixup the distance to integer as it comes out of the matrix as character
    pairwise.distances.flat$distance <- as.integer(pairwise.distances.flat$distance)

    return( list(distances=pairwise.distances.flat, seq.length=ncol( fasta.with.consensus )) );
} 
