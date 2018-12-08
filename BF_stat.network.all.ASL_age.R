#!/home1/wim/R-3.2.3/bin/Rscript --no-save --no-restore

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

	#df[ df$bf.with.interaction >= 1/3 & df$bf.with.interaction < 1, 'classification' ] <- "06. Anecdotal evidence for H0"
	#df[ df$bf.with.interaction >= 1/10 & df$bf.with.interaction < 1/3, 'classification' ] <- "07. Moderate evidence for H0"
	#df[ df$bf.with.interaction >= 1/30 & df$bf.with.interaction < 1/10, 'classification' ] <- "08. Strong evidence for H0"
	#df[ df$bf.with.interaction >= 1/100 & df$bf.with.interaction < 1/30, 'classification' ] <- "09. Very strong evidence for H0"
	#df[ df$bf.with.interaction < 1/100, 'classification' ] <- "10. Extreme evidence for H0"

	return( df )
}

###
# Return bayes factor for model with interaction term included.
##
get.bf <- function( data, freq, metric )
{
	# get subset
	data <- data[ data$bands == freq, ]
	data$y <- data[, metric ]

	# remove missings
	data <- data[ !is.na( data$asl_duration ), ]

	# bayes
	withAsl <- lmBF( y ~ asl_duration + Eyes + age, data = data )
	noAsl <- lmBF( y ~ Eyes + age, data = data )

	# Bayes Factor in favor of model WITH asl_duration
	bf <- withAsl / noAsl

	# get MCMC proportional error < 2%
	# bf <- recompute( bf, iterations = 1500000 )

	# bf is stored as log10 value
	bf.with.asl <- exp( bf@bayesFactor$bf )

	# plot bfs
	#plot( allBFs <- c( withAsl, noAsl ) )

	out <- data.frame( freq = freq, metric = metric, bf.with.asl = bf.with.asl )
	
	return( out )
}


#######################
# END FUNCTIONS
#######################

# to be sure to replicate stuff (like mst calculation).
set.seed( 123 )

# create output directory
metadir <- '/home1/wim/projects/willem_buitenhuis/2_edf2csv/demographics'
indir <- '/home1/wim/projects/willem_buitenhuis/3_health_res/plot_correlations'
outdir <- 'stat_ASL_age_all_correlations_nigeria'
dir.create( outdir, showWarnings = FALSE )

# metrics to compute bayes factors for
metrics <- c( "g.mst.strength.max", "g.mst.strength.mean", 
				"g.mst.degree.max", 
				"g.mst.bc.max", "g.mst.bc.median", 
				"g.mst.cc.max", "g.mst.cc.median", 
				"g.mst.leaf", "g.mst.diameter", 
				"g.mst.ecc", "g.mst.radius", "g.mst.Th", "g.mst.kappa" )

# frequency bands to compute bayes factors for
freqs <- c( 'delta', 'theta', 'alpha', 'beta', 'gamma' )

# get network data
data <- read.csv( paste0( indir, '/data.csv' ), row.names = 1 )
data$tmp <- NULL

# storage container
all <- NULL

for( freq in freqs )
{
	for( metric in metrics )
	{
		print( paste( freq, metric ) )
		all <- rbind( all, get.bf( data, freq, metric ) )
	}
}

# classify H1 bayes factors with labels (anecdotal, moderate, strong, extreme support) 
res <- label.bf( all )

# write output
write.csv( res, file = paste0( outdir, '/bayes_factors_asl_age.csv' ) )



