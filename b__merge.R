#!/home1/wim/R-3.2.3/bin/Rscript --no-save --no-restore

library( 'reshape2' )


# only process EEGs with raw data file size > min.file.size (in Mb)
min.file.size <- 1.0

# to be sure to replicate stuff (like mst calculation).
set.seed( 123 )

# create output directory
indir <- 'fc_epochs_nigeria'
outdir <- 'merged_fc_epochs_nigeria'
dir.create( outdir, showWarnings = FALSE )

# get all processed eeg time-series csv.gz files
infiles <- dir( indir )[ grep( "pli_processed-eeg", dir( indir ) ) ]

#processed-eeg-34_session_first.csv.gz

# columns to select
#sel.cols <- c( 'subject.id', 'sex', 'age', 'asl', 'asl_years', 'group', 'session', 'Eyes', 'epoch', 'bands', 'channel', 'variable', 'value' )

# container
all.pli.1 <- NULL
all.pli.2 <- NULL
all.pli.3 <- NULL

i <- 1
#infile <- infiles[i]

# loop over infiles
for( infile in infiles )
{
	print( paste0( infile, " -> ", i, " from ", length( infiles ) ) )
	# read epochs
	head( df <- read.csv( paste0( indir, '/', infile ), row.names = 1 ) )
	#df$sex <- as.character( df$sex )
	#df[ df$sex == 'FALSE', 'sex' ] <- 'female' # 'F' is converted to 'FALSE'...
	#df[ df$sex == 'M', 'sex' ] <- 'male'

	# check minimal raw EEG file size
	#if( unique( df$size.mb ) > min.file.size )
	#{
		#df.sel <- df[ , sel.cols ] # subset
		df.sel <- df # subset

		if( i < 50 ) {
			all.pli.1 <- rbind( all.pli.1, df.sel ) # store in container
		} else if( i > 100 ) {
			all.pli.3 <- rbind( all.pli.3, df.sel ) # store in container
		} else{
			all.pli.2 <- rbind( all.pli.2, df.sel ) # store in container
		}
	#}
	i <- i + 1
}

# merge the containers
all.pli <- rbind( rbind( all.pli.1, all.pli.2 ), all.pli.3 )

# remove the epilepsy group
#all.pli <- all.pli[ all.pli$group != 'epilepsy', ]

# write merged data to file
write.csv( all.pli, file = paste0( outdir, '/all.pli.csv' ) )


####################### NETWORK ################################


# get all processed eeg time-series csv.gz files
infiles <- dir( indir )[ grep( "network_processed-eeg", dir( indir ) ) ]

# columns to select 
#sel.cols <- c( 'no', 'sex', 'age', 'asl', 'asl_years', 'group', 'Eyes', 'epoch', 'bands' )
#sel.cols.metrics <- c( "g.strength.max", "g.strength.mean", "g.mst.strength.max", "g.mst.strength.mean", "g.mst.degree.max", 
#						"g.mst.degree.mean", "g.bc.max", "g.bc.median", "g.mst.bc.max", "g.mst.bc.median", "g.cc.max",           
#						"g.cc.median", "g.mst.cc.max", "g.mst.cc.median", "g.mst.leaf", "g.mst.m", "g.mst.diameter", "g.mst.ecc", 
#						"g.mst.radius", "g.mst.Th", "g.mst.kappa" )

#sel.cols <- c( sel.cols, sel.cols.metrics )

# container
all.network.1 <- NULL
all.network.2 <- NULL
all.network.3 <- NULL

i <- 1
# loop over infiles
for( infile in infiles )
{
	print( paste0( infile, " -> ", i, " from ", length( infiles ) ) )
	# read epochs
	head( df <- read.csv( paste0( indir, '/', infile ), row.names = 1 ) )
	#df$sex <- as.character( df$sex )
	#df[ df$sex == 'FALSE', 'sex' ] <- 'female' # 'F' is converted to 'FALSE'...
	#df[ df$sex == 'M', 'sex' ] <- 'male'

	# check minimal raw EEG file size
	#if( unique( df$size.mb ) > min.file.size & unique( df$group ) != 'epilepsy' )
	#{
		#df.sel <- df[ , sel.cols ] # subset
		df.sel <- df
	
		if( i < 50 ) {
			all.network.1 <- rbind( all.network.1, df.sel ) # store in container
		} else if( i > 100 ) {
			all.network.3 <- rbind( all.network.3, df.sel ) # store in container
		} else{
			all.network.2 <- rbind( all.network.2, df.sel ) # store in container
		}
	#}
	i <- i + 1
}

# merge the containers
all.network <- rbind( rbind( all.network.1, all.network.2 ), all.network.3 )

# remove the epilepsy group
#all.network <- all.network[ all.network$group != 'epilepsy', ]

# write merged data to file
write.csv( all.network, file = paste0( outdir, '/all.network.csv' ) )
 
