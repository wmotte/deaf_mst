#!/home1/wim/R-3.2.3/bin/Rscript --no-save --no-restore

library( 'reshape2' )
library( 'ggplot2' )
library( 'dplyr' )

#################################
# FUNCTIONS
#################################

# ggplot function
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

###
# Plot traces.
##
get.plot.eeg <- function( df )
{
	df$time <- 1:nrow( df ) / sf

	plot.eeg <- ggplot( data = melt( df, id = 'time' ), aes( x = time, color = variable, group = variable, y = value ) ) + 
				geom_line() + 
				facet_wrap( ~variable, scale = 'free_y', ncol = 2 ) + 
				theme_minimal( base_size = 12 ) + theme( legend.position = 'none' ) + xlab( "Time (s)" ) +
				scale_x_continuous( breaks = number_ticks( 25 ) ) +
				scale_y_continuous( breaks = number_ticks( 6 ) ) 

	return( plot.eeg )
}

###
# 95% confidence interval for the mean.
# https://www.r-bloggers.com/standard-deviation-vs-standard-error/
##
ci.dev <- function( x )
{
	#computation of the standard error of the mean
	sem <- sd( x, na.rm = TRUE ) / sqrt( length( x ) )
	
	#95% confidence intervals of the mean
	return( 2 * sem )
}

###
# get network data
##
get.network.data <- function( indir, metadir )
{
	# get meta data
	head( meta <- read.csv( paste0( metadir, '/', 'demographics.csv' ), row.names = 1 ) )

	# keep relevant columns only
	head( meta <- meta[, c( "subject.id", "group", "deafness", "asl_duration", "sex", "age" ) ] ) 

	head( df <- read.csv( paste0( indir, '/', 'all.network.csv' ), row.names = 1 ) )

	# merge with meta-data
	head( mdata <- merge( df, meta ) )

	# rename bands
	levels( mdata$bands ) <- c( 'delta', 'theta', 'alpha', 'beta', 'gamma' )

	# only keep deaf and control
	mdata <- mdata[ mdata$group %in% c( 'control', 'deaf' ), ]
	mdata$group <- as.factor( as.character( mdata$group ) )

	# return data
	return( mdata )
}


#######################
# END FUNCTIONS
#######################

# to be sure to replicate stuff (like mst calculation).
set.seed( 123 )

# create output directory
metadir <- 'demographics'
indir <- 'merged_fc_epochs_nigeria'
outdir <- 'plot_subset-alpha-beta'
dir.create( outdir, showWarnings = FALSE )

# get network data
df <- get.network.data( indir, metadir )

# write to file
levels( df$deafness ) <- c( "04. unknown", "01. birth", "03. post-lingual", "02. pre-lingual" )
write.csv( df, file = paste0( outdir, '/data.csv' ) )


metrics <- c( "g.mst.strength.max", "g.mst.strength.mean", "g.mst.degree.max", "g.mst.bc.max", "g.mst.bc.median", "g.mst.cc.max", 
				"g.mst.cc.median", "g.mst.leaf", "g.mst.diameter", "g.mst.ecc", "g.mst.radius", "g.mst.Th", "g.mst.kappa" )


# only keep

# average across epochs to get single value for condition A and single value for condition B, per subject
for( metric in metrics )
{
	print( metric )
	df$tmp <- NA
	df[ , 'tmp' ] <- df[, metric ]

	s <- NULL
	s <- as.data.frame( dplyr::summarise( grouped <- dplyr::group_by( df, group, Eyes, bands ), mean = mean( tmp ), ci.dev = ci.dev( tmp ) ) )

	# only select alpha and beta bands
	s <- s[ s$bands %in% c( 'alpha', 'beta' ), ]

	metric.label <- "NOT DEFINED"
	if( metric == "g.mst.strength.max" ) metric.label <- 'Strength (max)'
	if( metric == "g.mst.strength.mean" ) metric.label <- 'Strength (mean)'
	if( metric == "g.mst.degree.max" ) metric.label <- 'Degree (max)'

	if( metric == "g.mst.bc.max" ) metric.label <- 'BC (max)'
	if( metric == "g.mst.bc.mean" ) metric.label <- 'BC (mean)'
	if( metric == "g.mst.bc.median" ) metric.label <- 'BC (median)'

	if( metric == "g.mst.cc.max" ) metric.label <- 'Closeness centrality (max)'
	if( metric == "g.mst.cc.mean" ) metric.label <- 'Closeness centrality (mean)'
	if( metric == "g.mst.cc.median" ) metric.label <- 'Closeness centrality (median)'

	if( metric == "g.mst.leaf" ) metric.label <- 'Leaf number'
	if( metric == "g.mst.diameter" ) metric.label <- 'Diameter'
	if( metric == "g.mst.ecc" ) metric.label <- 'Eccentricity'
	if( metric == "g.mst.radius" ) metric.label <- 'Radius'
	if( metric == "g.mst.kappa" ) metric.label <- 'Kappa'
	if( metric == "g.mst.Th" ) metric.label <- 'Tree-hierarchy'

	# plot
	dodge <- position_dodge( 1 )
	p <- ggplot( s, aes( x = Eyes, y = mean, ymin = mean - ci.dev, ymax = mean + ci.dev, group = group, fill = group, colour = group ) ) +
			geom_line( alpha = 0.5, position = dodge, size = 1 ) +
			geom_errorbar( colour = 'gray50', position = dodge, width = 0.5 ) +
			geom_point( alpha = 1, stat = 'identity', position = dodge, size = 3 ) +
			facet_wrap( ~bands, nrow = 1 ) +
			scale_y_continuous( breaks = number_ticks( 8 ) ) +
			scale_fill_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
			scale_color_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
			xlab( "Eyes condition" ) +
			ylab( metric.label ) +
			theme_classic( base_size = 15 ) 
            
            p2 <- p + theme( legend.position = 'none', axis.title = element_text( face = "bold" ) )
            
            if( metric == 'g.mst.diameter' | metric == 'g.mst.cc.max' | metric == 'g.mst.cc.median' )
                p2 <- p + theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )

	# save to file
	ggsave( file = paste0( outdir, '/network-metrics__', metric, '.png' ), plot = p2, dpi = 200, height = 5, width = 4 )
}

stop( "..." )

## for type of deafness
for( metric in metrics )
{
	print( metric )
	df$tmp <- NA
	df[ , 'tmp' ] <- df[, metric ]

	s <- NULL
	s <- as.data.frame( dplyr::summarise( grouped <- dplyr::group_by( df, group, deafness, Eyes, bands ), mean = mean( tmp ), ci.dev = ci.dev( tmp ) ) )
	s$deafness <- ordered( as.character( s$deafness ) )
	s <- s[ s$deafness != "04. unknown", ]
	levels( s$deafness ) <- c( "congenital", "pre-lingual", "post-lingual", "unknown" )

	metric.label <- paste0( "Network metrics [", metric, "]" )

	if( metric == 'g.mst.strength.mean' )
		metric.label <- 'Mean strength'

	# plot
	dodge <- position_dodge( 1 )
	p <- ggplot( s, aes( x = Eyes, y = mean, ymin = mean - ci.dev, ymax = mean + ci.dev, group = deafness, fill = deafness, colour = deafness ) ) +
			geom_errorbar( colour = 'gray90', position = dodge, width = 0.5 ) +
			geom_line( alpha = 0.5, position = dodge, size = 1 ) +
			geom_point( alpha = 1, stat = 'identity', position = dodge, size = 3 ) +
			facet_wrap( ~bands, nrow = 1 ) +
			scale_y_continuous( breaks = number_ticks( 8 ) ) +
			scale_fill_manual( values = c( "red", "#E69F00", "#1f78b4", "gray50" ) ) +
			scale_color_manual( values = c("red", "#E69F00", "#1f78b4", "gray50" ) ) +
			xlab( "Eyes condition" ) +
			ylab( metric.label ) +
			theme_classic( base_size = 14 ) +
			theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )

	# save to file
	ggsave( file = paste0( outdir, '/deafness___network-metrics__', metric, '.png' ), plot = p, dpi = 200, height = 6, width = 8 )
}


