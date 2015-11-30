updateCalls <- function( calls ) {
    calls[ calls==-1 ] <- NA;
    calls[ calls==0 ] <- "AA";
    calls[ calls==1 ] <- "AB";
    calls[ calls==2 ] <- "BB";

    return( calls );
}
CalcCenters <- function( ii, alleleA, alleleB, geno ) {
    alleleA <- alleleA[ , ii ];
    alleleB <- alleleB[ , ii ];
    geno <- geno[ , ii ];

    AmAA <- median( alleleA[ geno==0 ] );
    AmAB <- median( alleleA[ geno==1 ] );
    AmBB <- median( alleleA[ geno==2 ] );

    BmAA <- median( alleleB[ geno==0 ] );
    BmAB <- median( alleleB[ geno==1 ] );
    BmBB <- median( alleleB[ geno==2 ] );

    c(  atan2( BmAA, AmAA ) * 2 /pi,
        atan2( BmAB, AmAB ) * 2 /pi,
        atan2( BmBB, AmBB ) * 2 /pi
        );
}

CalcBafCor <- function( ii, theta, centers, geno ) {
    theta <- theta[ , ii ];
    centers <- centers[ , ii ];
    geno <- geno[ , ii ];

    lessAA <- theta < centers[ "AA" ] & geno != 1;
    lessAB <- theta < centers[ "AB" ] & geno != 1;
    lessBB <- theta < centers[ "BB" ] & geno != 1;
    isAB <- geno == 1;

    down  <- lessAA 
    up    <- !lessBB
    mdown <- !lessAA & lessAB
    mup   <- !lessAB & lessBB

    M <- rep( NA, length( theta ) );
    M[ mdown ] <- 0.5 * ( theta[ mdown ] - centers[ "AA" ] ) / ( centers[ "AB" ] - centers[ "AA" ] )
    M[ mup   ] <- 0.5 * ( theta[ mup   ] - centers[ "AB" ] ) / ( centers[ "BB" ] - centers[ "AB" ] ) + 0.5
    M[ down  ] <- 0
    M[ up    ] <- 1
    M[ isAB  ] = theta[ isAB ]; 

    M[ M < 0 ] <- 0
    M[ M > 1 ] <- 1
    M;
}

cal_LRR_BAF<-function( calls.file = "AxiomGT1.calls.txt", summary.file = "AxiomGT1.summary.txt", annotation.file = "Axiom_GW_GT_Chicken.na35.annot.csv", LRR.output.name = "LRR.txt", BAF.output.name = "BAF.txt" ){
	# read "AxiomGT1.calls.txt"  <calls.file>
	calls <- read.table( calls.file, comment.char="#", header=TRUE )
	rownames( calls ) <- calls$probeset_id
    calls <- as( calls[ , -1 ], "matrix" )
    calls <- calls[ order( rownames( calls ) ), ]
	
	# Loading the summary file with the normalized intensities
	intensities <- read.table( summary.file, comment.char="#", header=TRUE )
	rownames( intensities ) <- intensities$probeset_id
    intensities <- as( intensities[ ,-1 ], "matrix" )
	
	# Loading the annotation file
	ann <- read.csv ( annotation.file, comment.char="#", header=TRUE )
    rownames( ann ) <- ann$Probe.Set.ID;
    ann <- ann[ , c( "Probe.Set.ID", "Chromosome", "Physical.Position" ) ]
    ann <- ann[ order( rownames( ann ) ), ]
    ann <- ann[ rownames( ann ) %in% rownames( calls ), ]
	
	#Filtering intensities
	names <- gsub( "\\-B", "", gsub( "\\-A", "", as.character( rownames( intensities ) ) ) )
    intensities <- intensities[ names %in% rownames( calls ), ]
	alleleA <- intensities[ seq( 1, nrow( intensities ), 2 ), ]
    alleleB <- intensities[ seq( 2, nrow( intensities ), 2 ), ]
	
	# Calculate LRR
	R <- as( alleleA + alleleB, "matrix" )
    medians <- apply( R, MARGIN=1, FUN=median)
    LRR <- log2( sweep( R, 1, medians, FUN="/" ) )/4
    rownames( LRR ) <- rownames( calls )
    colnames( LRR ) <- colnames( intensities )
	
	# Calculate BAF
	## Calculate centers
	centers <- sapply( 1:length( colnames( alleleA ) ), CalcCenters, alleleA=alleleA, alleleB=alleleB, geno=calls )
    colnames( centers ) <- colnames( calls )
    rownames( centers ) <- c( "AA", "AB", "BB" )
	
	## Calculate thita
	theta <- atan2( alleleB, alleleA ) * 2 / pi
    rownames( theta ) <- rownames( calls )
	
	## Calculate BAF
	BAF <- sapply( 1:length( colnames( theta ) ), CalcBafCor, theta=theta, centers=centers, geno=calls )
    colnames( BAF ) <- colnames( theta )
    rownames( BAF ) <- rownames( theta )
	
    total.LRR <- cbind( ann, LRR )
	total.BAF <- cbind( ann, BAF )
	
	colnames( total.LRR ) <- c( "Name", "Chr", "Position", paste0(gsub(".CEL", "", colnames( intensities ) ), "Log.R.Ratio" ) )
	colnames( total.BAF ) <- c( "Name", "Chr", "Position", paste0(gsub(".CEL", "", colnames( intensities ) ), "B.Allele.Freq" ) )
		
	write.table( total.LRR, file=LRR.output.name, quote=FALSE, row.names=FALSE, sep="\t", na="NA" )
	write.table( total.BAF, file=BAF.output.name, quote=FALSE, row.names=FALSE, sep="\t", na="NA" )
	}
	
	
	
	
	
	
	