#!/home1/wim/R-3.2.3/bin/Rscript --no-save --no-restore
#
# W.M. Otte (w.m.otte@umcutrecht.nl)
###############################################################
#
# i) Determine phase lag index for 14 x 14 channel combinations
# ii) Determine network characteristics, including mst

library( 'reshape2' )
library( 'ggplot2' )
library( 'igraph' )

#################################
# FUNCTIONS
#################################

# ggplot function
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

###
# Get pli matrix
##
phaselagindex <- function( EEGdata, frequency = 128 )
{
	## return wave list
	inputw2 <- function( wave, f, channel = 1 )
	{ 
		if(is.data.frame(wave)) {f<-f ; wave <- as.matrix(wave[,channel])}
		if(is.vector(wave)) {f<-f ; wave <- as.matrix(wave)}
		if(is.matrix(wave) && !stats::is.mts(wave)) {f<-f ; wave <- wave[,channel,drop=FALSE]}  
		if(is.ts(wave)) {f<-frequency(wave) ; wave <- as.matrix(wave)} 
		if(stats::is.mts(wave)) {f<-frequency(wave) ; wave <- as.matrix(wave[, channel])} 
		if(class(wave)=="Sample") {f<-wave$rate ; wave <- as.matrix(wave$sound[channel, ])}
		if(class(wave)=="audioSample"){f<-wave$rate ; wave <- as.matrix(wave)}
		if(class(wave)=="Wave")
		{    
			f <- wave@samp.rate
			if(channel==1) {wave <- as.matrix(wave@left)}   
			if(channel==2) {wave <- as.matrix(wave@right)}     
		}    
		return( list( w = wave,f = f ) )
	} 

	## hilbert transform
	hilbert <- function( wave, f )
	{ 
		wave <- inputw2( wave = wave, f = f )$w
		n <- nrow( wave )
		ff <- fft( wave )
		h <- rep( 0, n )
		if( n > 0 & 2 * floor( n / 2 ) == n ) { h[ c(1, n / 2 + 1 ) ] <- 1; h[ 2:( n / 2 ) ] <- 2 }
		else{ if( n > 0 ) { h[ 1 ] <-1; h[ 2:( ( n + 1 ) / 2 ) ] <- 2 } }
		ht <- fft( ff * h, inverse = TRUE ) / length( ff )
		return( ht )
	} 

	## phase lag index between two channels
	pli <- function( chan1, chan2, f )
	{ 
		return( abs( mean( sign( Arg( hilbert( chan1, f = f ) ) - Arg( hilbert( chan2, f = f ) ) ) ) ) )
	}

	nchan <- ncol( EEGdata )
	plimatrix <- matrix( NA, nchan, nchan )  

	# create matrix for all pairs of channel combinations
	for (i in 1:nchan )
	{
		for (j in 1:nchan )
		{
			plimatrix[ i, j ] = pli( EEGdata[ , i ], EEGdata[ , j ], f = frequency )
		}
	}

	# convert matrix to data.frame to set column and row names
	plimatrix <- as.data.frame( plimatrix )

	# get labels
	labels <- colnames( EEGdata )

	# chop labels to channel name only (e.g. F7, T7, ...)
	for( i in 1:length( labels ) )
		labels[ i ] <- gsub( "corrected_", "", unlist( strsplit( as.character( labels[ i ]  ), '.level'))[ 1 ] )

	colnames( plimatrix ) <- labels
	rownames( plimatrix ) <- labels

	return( plimatrix )
}

###
# Plot traces.
##
get.plot.eeg <- function( df )
{
	df$time <- 1:nrow( df ) / sf
	df$subject.id <- NULL
	df$session <- NULL

	plot.eeg <- ggplot( data = melt( df, id = 'time' ), aes( x = time, color = variable, group = variable, y = value ) ) + 
				geom_line() + 
				facet_wrap( ~variable, scale = 'free_y', ncol = 2 ) + 
				theme_minimal( base_size = 12 ) + theme( legend.position = 'none' ) + xlab( "Time (s)" ) +
				scale_x_continuous( breaks = number_ticks( 25 ) ) +
				scale_y_continuous( breaks = number_ticks( 6 ) ) 

	return( plot.eeg )
}

###
# matrix to long list
##
matrix2list <- function( plimatrix, subject.id, session, Eyes, epoch, level )
{
	# add channel label
	plimatrix$channel <- rownames( plimatrix )

	# wide -> long
	m <- melt( plimatrix, id = 'channel' )
	m$subject.id <- subject.id
	m$session <- session	
	m$Eyes <- Eyes
	m$epoch <- epoch
	m$level <- level

	# remove self-connectivity
	m <- m[ m$channel != m$variable, ]

	# remove edges with NA
	m <- m[ ! is.na( m$value ), ]

	return( m )
}

###
# Convert df to weighted igraph object
##
matrix2graph <- function( df )
{
	return( igraph::graph.adjacency( as.matrix( df ), mode = "undirected", weighted = TRUE, diag = FALSE ) )
}

###
# Get vector with several graph and weighted mst properties.
##
get.network.properties <- function( plimatrix, subject.id, session, Eyes, epoch, level )
{
	# get graph
	g <- matrix2graph( plimatrix )

	# invert matrix for g inverted (required for mst (high values -> short paths))
	plimatrix.inv <- 1 / as.matrix( plimatrix )
	diag( plimatrix.inv ) <- 0
	g.inv <- matrix2graph( plimatrix.inv )

	# get minimum spanning tree using Prim's algorithm (based on weights)	
	g.mst <- igraph::minimum.spanning.tree( g.inv, algorithm = 'prim' )

	# set values back to normal values
	E( g.mst )$weight  <- 1 / E( g.mst )$weight 

	# get strength values
	( g.strength.max <- max( igraph::strength( g ) ) )
	( g.strength.mean <- mean( igraph::strength( g ) ) )
	( g.mst.strength.max <- max( igraph::strength( g.mst ) ) )
	( g.mst.strength.mean <- mean( igraph::strength( g.mst ) ) )

	# get degree values
	( g.mst.degree.max <- max( igraph::degree( g.mst ) ) )

	# get betweenness centrality values
	( g.bc.max <- max( igraph::betweenness( g ) ) )
	( g.bc.median <- median( igraph::betweenness( g ) ) )
	( g.mst.bc.max <- max( igraph::betweenness( g.mst ) ) )
	( g.mst.bc.median <- median( igraph::betweenness( g.mst ) ) )

	# get closeness centrality values
	( g.cc.max <- max( igraph::closeness( g ) ) )
	( g.cc.median <- median( igraph::closeness( g ) ) )
	( g.mst.cc.max <- max( igraph::closeness( g.mst ) ) )
	( g.mst.cc.median <- median( igraph::closeness( g.mst ) ) )

	# get leaf-number (e.g. nodes with degree 1 )
	( g.mst.leaf <- sum( degree( g.mst ) == 1 ) )

	# get m (links)
	( g.mst.m <- length( V( g.mst ) ) - 1 )

	# get diameter
	( g.mst.diameter <- g.mst.m - g.mst.leaf + 2 )

	# mean eccentricity
	( g.mst.ecc <- mean( igraph::eccentricity( g.mst ) ) )

	# smallest eccentricity
	( g.mst.radius <- igraph::radius( g.mst ) )

	# tree hierarchy
	( g.mst.Th <- g.mst.leaf / ( 2 * g.mst.m * g.mst.bc.max ) )

	# kappa
	( g.mst.kappa <- mean( igraph::degree( g.mst )^2 ) / mean( igraph::degree( g.mst ) ) )

    # store in output frame
	out <- data.frame( subject.id = subject.id,
						session = session,
						Eyes = Eyes,
						epoch = epoch,
						level = level,
						g.strength.max = g.strength.max, 
						g.strength.mean = g.strength.mean, 
						g.mst.strength.max = g.mst.strength.max, 
						g.mst.strength.mean = g.mst.strength.mean, 
						g.mst.degree.max = g.mst.degree.max, 
						g.bc.max = g.bc.max,
						g.bc.median = g.bc.median,
						g.mst.bc.max = g.mst.bc.max,
						g.mst.bc.median = g.mst.bc.median,
						g.cc.max = g.cc.max, 
						g.cc.median = g.cc.median, 
						g.mst.cc.max = g.mst.cc.max, 
						g.mst.cc.median = g.mst.cc.median,
						g.mst.leaf = g.mst.leaf,
						g.mst.m = g.mst.m,
						g.mst.diameter = g.mst.diameter,
						g.mst.ecc = g.mst.ecc,
						g.mst.radius = g.mst.radius,
						g.mst.Th = g.mst.Th,
						g.mst.kappa = g.mst.kappa )

	return( out )
}

#######################
# END FUNCTIONS
#######################

# to be sure to replicate stuff consistantly.
set.seed( 123 )

# create output directory
indir <- 'epochs_nigeria'
outdir <- 'fc_epochs_nigeria'
dir.create( outdir, showWarnings = FALSE )

# sample freq
sf <- 128

# select levels and epoch
n.levels <- 7

# get all processed eeg time-series csv.gz files
infiles <- dir( indir )[ grep( "processed-eeg-", dir( indir ) ) ] 

# loop over infiles
for( infile in infiles )
{
	# output files
	outfile.pli <- paste0( outdir, '/pli_', infile )
	outfile.network <- paste0( outdir, '/network_', infile )

	# read all wavelet eeg data
	eeg <- read.csv( paste0( indir, '/', infile ), row.names = 1 )

	# container
	all.pli <- NULL
	all.network <- NULL

	for( epoch in unique( eeg$epochs ) )
	{
		for( level in 1:n.levels )
		{
			print( paste( infile, epoch, level ) )

			head( df <- eeg[ eeg$epochs == epoch, c( "subject.id", "session", "Eyes", "epochs", colnames( eeg )[ grep( paste0( "level.", level ), colnames( eeg ) ) ] ) ] )

			# eyes open/closed
			( Eyes <- as.character( unique( df$Eyes ) ) )
		
			# subject number
			( subject.id <- unique( df$subject.id ) )

			# session
			( session <- as.character( unique( df$session ) ) )

			df$Eyes <- NULL
			df$epochs <- NULL

			#plot.eeg <- get.plot.eeg( df )

			# get 14 x 14 matrix with functional connectivity (symmetric)
			plimatrix <- phaselagindex( df[ , !colnames( df ) %in% c( 'Eyes', 'epochs', 'session', 'subject.id' ) ] )

			# get list from functional matrix
			subject.pli <- matrix2list( plimatrix, subject.id, session, Eyes, epoch, level )
			all.pli <- rbind( all.pli, subject.pli )

			# get list of functional network properties
			subject.network <- get.network.properties( plimatrix, subject.id, session, Eyes, epoch, level )
			all.network <- rbind( all.network, subject.network )
		}
	}

	# set frequency bands
	#1: 32-64 Hz	gamma
	#2: 16-32 Hz 	beta <- change between 'open' - 'close'
	#3: 8-16 Hz 	alpha <- much change between 'open' - 'close'
	#4: 4-8 Hz 		theta <- change between 'open' - 'close'
	#5: 2-4 Hz 		delta
	#6: 1-2 Hz 		delta
	#7: 0.5-1 Hz 	delta
	all.pli$bands <- NA
	all.pli[ all.pli$level == 1, 'bands' ] <- '5_gamma'
	all.pli[ all.pli$level == 2, 'bands' ] <- '4_beta'
	all.pli[ all.pli$level == 3, 'bands' ] <- '3_alpha'
	all.pli[ all.pli$level == 4, 'bands' ] <- '2_theta'
	all.pli[ all.pli$level >= 5, 'bands' ] <- '1_delta'

	all.network$bands <- NA
	all.network[ all.network$level == 1, 'bands' ] <- '5_gamma'
	all.network[ all.network$level == 2, 'bands' ] <- '4_beta'
	all.network[ all.network$level == 3, 'bands' ] <- '3_alpha'
	all.network[ all.network$level == 4, 'bands' ] <- '2_theta'
	all.network[ all.network$level >= 5, 'bands' ] <- '1_delta'

	# write subject data to output
	write.csv( file = outfile.pli, all.pli )
	write.csv( file = outfile.network, all.network )
} 


