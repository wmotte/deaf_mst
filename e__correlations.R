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
	head( meta <- meta[, c( "subject.id", "group", "deafness", "asl_duration", "sex", "age" ) ] ) # "time_of_deafness", "hears", "years_education"

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
outdir <- 'plot_correlations'
dir.create( outdir, showWarnings = FALSE )

# get network data
df <- get.network.data( indir, metadir )

# write to file
levels( df$deafness ) <- c( "04. unknown", "01. birth", "03. post-lingual", "02. pre-lingual" )
write.csv( df, file = paste0( outdir, '/data.csv' ) )


metrics <- c( "g.mst.strength.mean", "g.mst.strength.max",
				"g.mst.degree.max", 
				"g.mst.bc.max", "g.mst.bc.median", 
				"g.mst.cc.max", "g.mst.cc.median", 
				"g.mst.leaf", "g.mst.diameter", 
				"g.mst.ecc", "g.mst.radius", "g.mst.Th", "g.mst.kappa" )

# selection groups 
for( group in c( 'control', 'deaf' ) )
{
	# average across epochs to get single value for condition A and single value for condition B, per subject
	for( metric in metrics )
	{
		print( metric )
		df$tmp <- NA
		df[ , 'tmp' ] <- df[, metric ]

		# selection groups
		df.select <- df[ df$group == group, ]

		s <- NULL
		s <- as.data.frame( dplyr::summarise( grouped <- dplyr::group_by( df.select, subject.id, group, Eyes, bands ), 
										mean = mean( tmp ), asl = mean( asl_duration, na.rm = TRUE ) ) )

		p <- ggplot( df.select, aes( x = asl_duration, y = tmp, group = Eyes, fill = Eyes, colour = Eyes ) ) + 
				geom_point( data = s, aes( x = asl, y = mean ), alpha = 0.3 ) +
				geom_smooth( method = 'lm' ) +
				facet_wrap( ~bands, nrow = 1 ) +
				scale_y_continuous( breaks = number_ticks( 8 ) ) +
				scale_fill_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
				scale_color_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
				xlab( "Sign language experience (years)" ) +
				ylab( paste0( "Network metrics [", metric, "]" ) ) +
				theme_classic( base_size = 14 ) +
				theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )

		# save to file
		ggsave( file = paste0( outdir, '/correlation-metrics__', group, '__', metric, '.png' ), plot = p, dpi = 200, height = 6, width = 8 )
	}
}

