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

# container
all.pli.1 <- NULL
all.pli.2 <- NULL
all.pli.3 <- NULL

i <- 1

# loop over infiles
for( infile in infiles )
{
	print( paste0( infile, " -> ", i, " from ", length( infiles ) ) )
	# read epochs
	head( df <- read.csv( paste0( indir, '/', infile ), row.names = 1 ) )

		df.sel <- df # subset

		if( i < 50 ) {
			all.pli.1 <- rbind( all.pli.1, df.sel ) # store in container
		} else if( i > 100 ) {
			all.pli.3 <- rbind( all.pli.3, df.sel ) # store in container
		} else{
			all.pli.2 <- rbind( all.pli.2, df.sel ) # store in container
		}
	i <- i + 1
}

# merge the containers
all.pli <- rbind( rbind( all.pli.1, all.pli.2 ), all.pli.3 )

# write merged data to file
write.csv( all.pli, file = paste0( outdir, '/all.pli.csv' ) )


####################### NETWORK ################################


# get all processed eeg time-series csv.gz files
infiles <- dir( indir )[ grep( "network_processed-eeg", dir( indir ) ) ]

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

		df.sel <- df
	
		if( i < 50 ) {
			all.network.1 <- rbind( all.network.1, df.sel ) # store in container
		} else if( i > 100 ) {
			all.network.3 <- rbind( all.network.3, df.sel ) # store in container
		} else{
			all.network.2 <- rbind( all.network.2, df.sel ) # store in container
		}
	
    i <- i + 1
}

# merge the containers
all.network <- rbind( rbind( all.network.1, all.network.2 ), all.network.3 )

# write merged data to file
write.csv( all.network, file = paste0( outdir, '/all.network.csv' ) )
