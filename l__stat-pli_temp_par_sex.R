#!/home1/wim/R-3.2.3/bin/Rscript --no-save --no-restore

library( 'reshape2' )
library( 'ggplot2' )
library( 'dplyr' )
library( 'BayesFactor' )


#################################
# FUNCTIONS
#################################

###
# Label bayes factors.
# https://www.r-bloggers.com/what-does-a-bayes-factor-feel-like/
##
label.bf <- function( df )
{
	df$classification <- NA
	
	df[ df$bf.with.asl > 100, 'classification' ] <- "01. Extreme evidence for H1"
	df[ df$bf.with.asl > 30 & df$bf.with.asl <= 100, 'classification' ] <- "02. Very strong evidence for H1"
	df[ df$bf.with.asl > 10 & df$bf.with.asl <= 30, 'classification' ] <- "03. Strong evidence for H1"
	df[ df$bf.with.asl > 3 & df$bf.with.asl <= 10, 'classification' ] <- "04. Moderate evidence for H1"
	df[ df$bf.with.asl > 1 & df$bf.with.asl <= 3, 'classification' ] <- "05. Anecdotal evidence for H1"

	df[ df$bf.with.asl == 1, 'classification' ] <- "No evidence"

	return( df )
}

###
# Return bayes factor for model with interaction term included.
##
get.bf <- function( data, freq )
{
	# get subset
	data <- data[ data$bands == freq, ]

	# remove missings
	data <- data[ !is.na( data$value ), ]

	# bayes
	withGroup <- lmBF( value ~ group + Eyes + sex, data = data )
	noGroup <- lmBF( value ~ Eyes + sex, data = data )

	# Bayes Factor in favor of model WITH asl_duration
	bf <- withGroup / noGroup

	# get MCMC proportional error < 2%
	# bf <- recompute( bf, iterations = 1500000 )

	# bf is stored as log10 value
	bf.with.group <- exp( bf@bayesFactor$bf )

	# plot bfs
	#plot( allBFs <- c( withAsl, noAsl ) )

	out <- data.frame( freq = freq, metric = 'pli', bf.with.asl = bf.with.group )
	
	return( out )
}

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
# get functional connectivity data
##
get.functional.connectivity.data <- function( indir, metadir )
{
	# get meta data
	head( meta <- read.csv( paste0( metadir, '/', 'demographics.csv' ), row.names = 1 ) )

	# keep relevant columns only
	head( meta <- meta[, c( "subject.id", "group", "deafness", "asl_duration", "sex", "age" ) ] ) # "time_of_deafness", "hears", "years_education"

	head( df <- read.csv( paste0( indir, '/', 'all.pli.csv' ), row.names = 1 ) )

	# merge O1 and O2 and P1 and P2
	df$meta.channel <- NA
	df$meta.variable <- NA

	df[ df$channel %in% c( "P7", "P8" ), 'meta.channel' ] <- 'Parietal'
	df[ df$channel %in% c( "T7", "T8" ), 'meta.channel' ] <- 'Temporal'

	df[ df$variable %in% c( "P7", "P8" ), 'meta.variable' ] <- 'Parietal'
	df[ df$variable %in% c( "T7", "T8" ), 'meta.variable' ] <- 'Temporal'

	df$meta.channel <- as.factor( df$meta.channel )
	df$meta.variable <- as.factor( df$meta.variable )

	# only keep Parietal and Temporal channel data
	head( df.sub <- df[ df$meta.channel == 'Parietal' | df$meta.channel == 'Temporal', ] )
	head( df.sub <- df.sub[ ! is.na( df.sub$meta.channel ), ] )
	head( df.sub <- df.sub[ ! is.na( df.sub$meta.variable ), ] )

	# remove self-connectivity
	s <- df.sub[ df.sub$meta.channel != df.sub$meta.variable, ]

	# remove duplicates by sorting channel names and combine this with fc-value
	s$name <- NA
	for( i in 1:nrow( s ) )
		s[ i, 'name' ] <- paste0( 
								s[ i, 'bands' ],
								s[ i, 'subject.id' ], 
								s[ i, 'session' ], 
								sort( c( as.character( s[ i, 'channel' ] ), as.character( s[ i, 'variable' ] ) ) ), 
								s[ i, 'value' ], collapse = '-' )

	s <- s[ !duplicated( s$name ), ]

	# merge with meta-data
	mdata <- merge( s, meta )

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
metadir <- '/home1/wim/projects/willem_buitenhuis/3_health_res/demographics'
indir <- '/home1/wim/projects/willem_buitenhuis/3_health_res/merged_fc_epochs_nigeria'
outdir <- 'stats_pli_nigeria__temp_par_sex'
dir.create( outdir, showWarnings = FALSE )

# get functional connectivity
head( fc <- get.functional.connectivity.data( indir, metadir ) )

# save FC data for sharing online
write.csv( fc, file = paste0( outdir, '/epoch_data.csv' ) )

# average across epochs to get single value for condition A and single value for condition B, per subject
grouped <- dplyr::group_by( fc, group, Eyes, bands, meta.channel, meta.variable )
( s <- as.data.frame( dplyr::summarise( grouped, mean = mean( value ), ci.dev = ci.dev( value ) ) ) )

# write raw data for stats
write.csv( fc, file = paste0( outdir, '/pli-parietal-temporal.csv' ) )

# frequency bands to compute bayes factors for
freqs <- c( 'delta', 'theta', 'alpha', 'beta', 'gamma' )

# storage container
all <- NULL

for( freq in freqs )
{
	print( paste( freq ) )
	all <- rbind( all, get.bf( fc, freq ) )
}

# classify H1 bayes factors with labels (anecdotal, moderate, strong, extreme support) 
res <- label.bf( all )

# write output
write.csv( res, file = paste0( outdir, '/group_differences_in_temp-par_pli-FC_sex__bayes-factors.csv' ) )

