#!/home1/wim/R-3.2.3/bin/Rscript --no-save --no-restore

library( 'reshape2' )
library( 'ggplot2' )
library( 'dplyr' )
require( 'reshape2' )

#################################
# FUNCTIONS
#################################

# ggplot function
number_ticks <- function( n ) { function( limits ) pretty( limits, n ) }

#######################
# END FUNCTIONS
#######################

# to be sure to replicate stuff (like mst calculation).
set.seed( 123 )

# create output directory
metadir <- 'demographics'

indir <- 'plot_correlations'
outdir <- 'plot_correlations_extra'

dir.create( outdir, showWarnings = FALSE )

# get network data
head( df <- read.csv( paste0( indir, '/data.csv' ), row.names = 1 ) )

# only keep theta
df <- df[ df$bands == 'theta', ]
df$bands <- NULL
df$deafness <- NULL
df$sex <- NULL
df$age <- NULL

df$g.strength.max <- NULL
df$g.strength.mean <- NULL
df$g.bc.max <- NULL
df$g.bc.median <- NULL
df$g.cc.max <- NULL
df$g.cc.median <- NULL
df$g.mst.m <- NULL

df$g.mst.strength.max <- NULL
df$g.mst.strength.mean <- NULL
df$g.mst.bc.max <- NULL
df$g.mst.cc.max <- NULL
df$g.mst.cc.median <- NULL

# from wide to long
mdata <- melt( df, id = c( "subject.id", "session", "Eyes", "epoch", "level", "group", "asl_duration" ) )

levels( mdata$variable ) <- c( "Degree (max)", "BC (median)", "Leaf", "Diameter", 
									"Eccentricity", "Radius", "Tree-hierarchy", "Kappa" ) 

s <- as.data.frame( dplyr::summarise( grouped <- dplyr::group_by( mdata, subject.id, group, Eyes, variable ), 
										mean = mean( value ), asl = mean( asl_duration, na.rm = TRUE ) ) )

# all subjects together
p <- ggplot( mdata, aes( x = asl_duration, y = value, group = Eyes, fill = Eyes, colour = Eyes ) ) + 
				geom_point( data = s, aes( x = asl, y = mean ), alpha = 0.3 ) +
				geom_smooth( method = 'lm' ) +
				facet_wrap( ~variable, ncol = 4, scale = 'free_y' ) +
				scale_y_continuous( breaks = number_ticks( 8 ) ) +
				scale_fill_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
				scale_color_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
				xlab( "Sign language experience (years)" ) +
				ylab( "Metrics value" ) +
				theme_classic( base_size = 12 ) +
				theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )

# save to file
ggsave( file = paste0( outdir, '/correlation-metrics__all.png' ), plot = p, dpi = 200, height = 6, width = 7 )


# all subjects together
p2 <- ggplot( mdata[ mdata$group == 'deaf', ], aes( x = asl_duration, y = value, group = Eyes, fill = Eyes, colour = Eyes ) ) + 
				geom_point( data = s[ s$group == 'deaf', ], aes( x = asl, y = mean ), alpha = 0.3 ) +
				geom_smooth( method = 'lm' ) +
				facet_wrap( ~variable, ncol = 4, scale = 'free_y' ) +
				scale_y_continuous( breaks = number_ticks( 8 ) ) +
				scale_fill_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
				scale_color_manual( values = c("#1f78b4", "#E69F00", "gray50" ) ) +
				xlab( "Sign language experience (years)" ) +
				ylab( "Metrics value" ) +
				theme_classic( base_size = 12 ) +
				theme( legend.position = 'top', axis.title = element_text( face = "bold" ) )

# save to file
ggsave( file = paste0( outdir, '/correlation-metrics__deaf.png' ), plot = p2, dpi = 200, height = 6, width = 7 )

