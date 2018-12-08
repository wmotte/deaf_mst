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
	withGroup <- lmBF( value ~ group + Eyes, data = data )
	noGroup <- lmBF( value ~ Eyes, data = data )

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

	df[ df$channel %in% c( "O1", "O2" ), 'meta.channel' ] <- 'Occipital'
	df[ df$channel %in% c( "P7", "P8" ), 'meta.channel' ] <- 'Parietal'

	df[ df$variable %in% c( "O1", "O2" ), 'meta.variable' ] <- 'Occipital'
	df[ df$variable %in% c( "P7", "P8" ), 'meta.variable' ] <- 'Parietal'

	df$meta.channel <- as.factor( df$meta.channel )
	df$meta.variable <- as.factor( df$meta.variable )

	# only keep Parietal and Temporal channel data
	head( df.sub <- df[ df$meta.channel == 'Occipital' | df$meta.channel == 'Parietal', ] )
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
outdir <- 'violinplot_pli_nigeria__occ_par'
dir.create( outdir, showWarnings = FALSE )

# get functional connectivity
head( fc <- get.functional.connectivity.data( indir, metadir ) )

# save FC data for sharing online
write.csv( fc, file = paste0( outdir, '/epoch_data.csv' ) )

# average across epochs to get single value for condition A and single value for condition B, per subject
#grouped <- dplyr::group_by( fc, group, Eyes, bands, meta.channel, meta.variable )
#( s <- as.data.frame( dplyr::summarise( grouped, mean = mean( value ), ci.dev = ci.dev( value ) ) ) )

# average per individual
grouped <- dplyr::group_by( fc, subject.id, group, Eyes, bands, meta.channel, meta.variable )
( s <- as.data.frame( dplyr::summarise( grouped, mean = mean( value ), ci.dev = ci.dev( value ) ) ) )

# FC between parietal - temporal cortex
p <- ggplot( data = s, aes( x = Eyes, y = mean, color = group, fill = group ) ) + 
			facet_wrap( ~bands, nrow = 1, scale = 'free_y' ) + 
			geom_violin(position = position_dodge(1) ) +		
			geom_boxplot(width=0.1, position = position_dodge(1), color = "black" ) +			
			#geom_bar( stat = 'identity', position = position_dodge(), alpha = 0.8, color = 'gray30' ) + 
			#geom_errorbar( stat = 'identity', position = position_dodge( 0.9), color = 'gray70', width = 0.4 ) + 
            scale_y_continuous( breaks = number_ticks( 4 ) ) + 
            scale_fill_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) + 
            scale_color_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) + 
			ylab( "Functional connectivity (PLI)" ) +
			xlab( "Eyes condition" ) +            
            theme_classic( base_size = 14 ) + 
            theme( legend.position = 'top', axis.title = element_text( face = "bold" ) ) 

# save
ggsave( file = paste0( outdir, '/pli-parietal-occipital-violin.png' ), plot = p, dpi = 200, height = 8, width = 12 )

# write raw data for stats
write.csv( fc, file = paste0( outdir, '/pli-parietal-occipital-violin.csv' ) )


# BFs are already calculated in previous script for bar graph plots


###################
# CREATE Delta's
###################

# average across epochs to get single value for condition A and single value for condition B, per subject
grouped <- dplyr::group_by( fc, subject.id, group, Eyes, bands, meta.channel, meta.variable )
s <- as.data.frame( dplyr::summarise( grouped, mean = mean( value ) ) )

dim( s.open <- s[ s$Eyes == 'open', ] )
s.open$value.open <- s.open$mean
s.open$mean <- s.open$Eyes <- NULL

dim( s.closed <- s[ s$Eyes == 'closed', ] )
s.closed$value.closed <- s.closed$mean
s.closed$mean <- s.closed$Eyes <- NULL

# get %delta relative to open condition -> (open - closed) / open
dim( s.m <- merge( s.open, s.closed ) )
s.m$delta <- 100 * ( ( s.m$value.open - s.m$value.closed ) / s.m$value.open )

# average per subject
grouped <- dplyr::group_by( s.m, subject.id, group, bands, meta.channel, meta.variable )

s.m.mean <- as.data.frame( dplyr::summarise( grouped, mean = mean( delta ) ) )
s.m.ci <- as.data.frame( dplyr::summarise( grouped, ci.dev = ci.dev( delta ) ) )
s.m.avg <- merge( s.m.mean, s.m.ci )
s.m.avg$Eyes <- as.factor( 'Delta' )

# write to file
write.csv( s.m.avg, file = paste0( outdir, '/delta_FC_values_used_for_delta_plot_violinplot.csv' ) )

# FC between parietal - temporal cortex
p.delta <- ggplot( data = s.m.avg, aes( x = Eyes, y = mean, color = group, fill = group ) ) + 
			facet_wrap( ~bands, nrow = 1 ) + #, scale = 'free_y' ) + 
			geom_violin(position = position_dodge(1) ) +		
			geom_boxplot(width=0.1, position = position_dodge(1), color = "black" ) +	
			#geom_bar( stat = 'identity', position = position_dodge(), alpha = 0.8, color = 'gray30' ) + 
			#geom_errorbar( stat = 'identity', position = position_dodge( 0.9), color = 'gray70', width = 0.4 ) + 
            scale_y_continuous( breaks = number_ticks( 8 ) ) + 
            scale_fill_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) + 
            scale_color_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) + 
			ylab( expression( bold( Delta~"Functional connectivity (%)"  ) ) ) +
			xlab( "Eyes condition" ) +            
            theme_classic( base_size = 14 ) + 
            theme( legend.position = 'top', axis.title = element_text( face = "bold" ) ) 
# save
ggsave( file = paste0( outdir, '/pli-occipital-parietal__delta_violin.png' ), plot = p.delta, dpi = 200, height = 8, width = 5)

