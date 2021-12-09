addReproduciblePeakSet_win <- function(
	ArchRProj = NULL,
	groupBy = "Clusters",
	peakMethod = "Macs2",
	reproducibility = "2",
	peaksPerCell = 500,
	maxPeaks = 150000,
	minCells = 25,
	excludeChr = c("chrM","chrY"),
	pathToMacs2 = if(tolower(peakMethod)=="macs2") findMacs2() else NULL,
	genomeSize = NULL, 
	shift = -75, 
	extsize = 150, 
	method = if(tolower(peakMethod)=="macs2") "q" else "p", #P-Method for Tiles Results Better Agree w/ Macs2
	cutOff = 0.1, 
	additionalParams = "--nomodel --nolambda",
	extendSummits = 250,
	promoterRegion = c(2000, 100),
	genomeAnnotation = getGenomeAnnotation(ArchRProj),
	geneAnnotation = getGeneAnnotation(ArchRProj),
  plot = TRUE,
	threads = getArchRThreads(),
	parallelParam = NULL,
	force = FALSE,
	verbose = TRUE,
	logFile = createLogFile("addReproduciblePeakSet"),
	...
	){
	
	ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
	ArchR:::.validInput(input = groupBy, name = "groupBy", valid = c("character"))
	ArchR:::.validInput(input = reproducibility, name = "reproducibility", valid = c("character"))
	ArchR:::.validInput(input = peaksPerCell, name = "peaksPerCell", valid = c("integer"))
	ArchR:::.validInput(input = maxPeaks, name = "maxPeaks", valid = c("integer"))
	ArchR:::.validInput(input = minCells, name = "minCells", valid = c("integer"))
	ArchR:::.validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))

	if(tolower(peakMethod) == "macs2"){
		#ArchR:::.validInput(input = pathToMacs2, name = "pathToMacs2", valid = c("character"))
		ArchR:::.validInput(input = genomeSize, name = "genomeSize", valid = c("character", "numeric", "null"))
		ArchR:::.validInput(input = shift, name = "shift", valid = c("integer"))
		ArchR:::.validInput(input = extsize, name = "extsize", valid = c("integer"))
		ArchR:::.validInput(input = method, name = "method", valid = c("character"))
		ArchR:::.validInput(input = additionalParams, name = "additionalParams", valid = c("character"))
		ArchR:::.validInput(input = extendSummits, name = "extendSummits", valid = c("integer"))
		#.checkMacs2Options(pathToMacs2) #Check Macs2 Version		
	}else if(tolower(peakMethod) == "tiles"){

	}else{
		stop("peakMethod not recognized! Supported peakMethods are Macs2 or Tiles!")
	}

	ArchR:::.validInput(input = cutOff, name = "cutOff", valid = c("numeric"))
	ArchR:::.validInput(input = promoterRegion, name = "promoterRegion", valid = c("integer"))
	geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)
	genomeAnnotation <- ArchR:::.validGenomeAnnotation(genomeAnnotation)
	geneAnnotation <- ArchR:::.validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation)
  ArchR:::.validInput(input = plot, name = "plot", valid = c("boolean"))
	ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
	ArchR:::.validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
	ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
	ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
	ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

	tstart <- Sys.time()
	ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "ReproduciblePeakSet Args", logFile=logFile)

	if(tolower(peakMethod) == "macs2"){

		ArchR:::.logMessage("Calling Peaks with Macs2", logFile = logFile)

		#utility <- ArchR:::.checkPath(pathToMacs2)

		coverageMetadata <- ArchR:::.getCoverageMetadata(ArchRProj = ArchRProj, groupBy = groupBy, minCells = minCells)
		coverageParams <- ArchR:::.getCoverageParams(ArchRProj = ArchRProj, groupBy = groupBy)

		#####################################################
		# Peak Calling Summary
		#####################################################
		tableGroups <- table(getCellColData(ArchRProj, groupBy, drop = TRUE))
		groupSummary <- lapply(seq_along(coverageParams$cellGroups), function(y){
			x <- coverageParams$cellGroups[[y]]
			uniq <- unique(unlist(x))
			n <- lapply(x, length) %>% unlist %>% sum
			nmin <- lapply(x, length) %>% unlist %>% min
			nmax <- lapply(x, length) %>% unlist %>% max
			data.frame(
			  Group=names(coverageParams$cellGroups)[y], 
			  nCells=tableGroups[names(coverageParams$cellGroups)[y]], 
			  nCellsUsed=length(uniq), 
			  nReplicates=length(x), 
			  nMin=nmin, 
			  nMax=nmax, 
			  maxPeaks = min(maxPeaks, length(uniq) * peaksPerCell)
			)
		}) %>% Reduce("rbind",.)

		ArchR:::.logDiffTime("Peak Calling Parameters!", tstart, verbose = verbose, logFile = logFile)
		ArchR:::.logThis(groupSummary, "PeakCallSummary", logFile = logFile)
		if(verbose) print(groupSummary)

		#####################################################
		# Create Output Directory
		#####################################################
		outDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls")
		outSubDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls", "ReplicateCalls")
		outBedDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls", "InsertionBeds")
		dir.create(outDir, showWarnings = FALSE)
		dir.create(outSubDir, showWarnings = FALSE)
		dir.create(outBedDir, showWarnings = FALSE)

		#####################################################
		# Genome Size Presets
		#####################################################
		if(is.null(genomeSize)){
			if(grepl("hg19|hg38", getGenome(ArchRProj), ignore.case = TRUE)){
				genomeSize <- 2.7e9
			}else if(grepl("mm9|mm10", getGenome(ArchRProj), ignore.case = TRUE)){
				genomeSize <- 1.87e9
			}
		}

		#####################################################
		# Arguments for Peak Calling
		#####################################################
		coverageFiles <- coverageMetadata$File
		names(coverageFiles) <- coverageMetadata$Name
		args <- list()
		args$X <- seq_len(nrow(coverageMetadata))
		args$FUN <- callSummitsOnCoverages_win
		args$coverageFiles <- coverageFiles
		args$outFiles <- file.path(outSubDir, paste0(make.names(coverageMetadata$Name),"-summits.rds"))
		args$bedDir <- outBedDir
		args$excludeChr <- excludeChr
		args$peakParams <- list(
				pathToMacs2 = pathToMacs2,
				genomeSize = genomeSize, 
				shift = shift, 
				extsize = extsize, 
				cutOff = cutOff, 
				method = method,
				additionalParams = additionalParams
			)
		args$parallelParam <- parallelParam
		args$threads <- threads
		args$logFile <- logFile
		args$registryDir <- file.path(outDir, "batchRegistry")

		#####################################################
		# Batch Call Peaks
		#####################################################
		ArchR:::.logDiffTime("Batching Peak Calls!", tstart, verbose = verbose, logFile = logFile)
		ArchR:::.logThis(args, "PeakArgs", logFile = logFile)

		#back lapply
		outSummmits <- unlist(ArchR:::.batchlapply(args))
		
		#Summarize Output
		outSummitList <- split(outSummmits, coverageMetadata$Group)
		summitNamesList <- split(coverageMetadata$Name, coverageMetadata$Group)

		#####################################################
		# BSgenome for Add Nucleotide Frequencies!
		#####################################################
		ArchR:::.requirePackage(genomeAnnotation$genome)
		ArchR:::.requirePackage("Biostrings",source="bioc")
		BSgenome <- eval(parse(text = genomeAnnotation$genome))
		BSgenome <- validBSgenome(BSgenome)

		#####################################################
		# Identify Reproducible Peaks!
		#####################################################
		ArchR:::.logDiffTime("Identifying Reproducible Peaks!", tstart, verbose = verbose, logFile = logFile)
		groupPeaks <- ArchR:::.safelapply(seq_along(outSummitList), function(i){
			prefix <- sprintf("Rep. Peaks Group (%s of %s) :", i, length(outSummitList))
			ArchR:::.logDiffTime(sprintf("%s Creating Reproducible Peaks", prefix), tstart, verbose = FALSE, logFile = logFile)
			peaks <- suppressMessages(ArchR:::.identifyReproduciblePeaks(
				summitFiles = outSummitList[[i]],
				summitNames = summitNamesList[[i]],
				reproducibility = reproducibility,
	    	extendSummits = extendSummits,
	    	blacklist = genomeAnnotation$blacklist,
	    	prefix = prefix,
	    	logFile = logFile
			))
			ArchR:::.logDiffTime(sprintf("%s Annotating and Filtering Peaks", prefix), tstart, verbose = FALSE, logFile = logFile)
			peaks <- sort(sortSeqlevels(peaks))
			peaks <- subsetByOverlaps(peaks, genomeAnnotation$chromSizes, type = "within")
			peaks <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion, logFile = logFile)
			peaks <- peaks[which(mcols(peaks)$N < 0.001)] #Remove N Containing Peaks
			peaks <- peaks[order(peaks$groupScoreQuantile, decreasing = TRUE)]
			peaks <- head(peaks, groupSummary[names(outSummitList)[i],"maxPeaks"])
			mcols(peaks)$N <- NULL #Remove N Column
			print(file.path(outDir, paste0(make.names(names(outSummitList)[i]), "-reproduciblePeaks.gr.rds")))
			saveRDS(peaks, file.path(outDir, paste0(make.names(names(outSummitList)[i]), "-reproduciblePeaks.gr.rds")))
			return(peaks)
		}, threads = threads) %>% GRangesList()
		names(groupPeaks) <- names(outSummitList)

		#Construct Union Peak Set
		ArchR:::.logDiffTime("Creating Union Peak Set!", tstart, verbose = verbose, logFile = logFile)
		unionPeaks <- unlist(groupPeaks)
		unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)

		#Summarize Output
		peakDF <- lapply(seq_along(groupPeaks), function(x){
			data.frame(Group = names(groupPeaks)[x], table(groupPeaks[[x]]$peakType))
		}) %>% Reduce("rbind", .)
		peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
		peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
		peakDF$Freq <- peakDF$Freq / 1000
		metadata(unionPeaks)$PeakCallSummary <- peakDF

	}else if(tolower(peakMethod) == "tiles"){

		ArchR:::.logMessage("Calling Peaks with TileMatrix. We recommend using the Macs2 Version.\nThis method is still under development.", logFile = logFile)

		useMatrix <- "TileMatrix"
		
		cellGroups <- addGroupCoverages(
			ArchRProj = ArchRProj, 
			groupBy = groupBy,
			returnGroups = TRUE,
			minCells = minCells,
			logFile = logFile,
			...
		)[[1]]$cellGroups
    if(verbose) print(cellGroups)

		#####################################################
		# Peak Calling Summary
		#####################################################
		tableGroups <- table(getCellColData(ArchRProj, groupBy, drop = TRUE))
		groupSummary <- lapply(seq_along(cellGroups), function(y){
			x <- cellGroups[[y]]
			uniq <- unique(unlist(x))
			n <- lapply(x, length) %>% unlist %>% sum
			nmin <- lapply(x, length) %>% unlist %>% min
			nmax <- lapply(x, length) %>% unlist %>% max
			data.frame(
			  Group=names(cellGroups)[y], 
			  nCells=tableGroups[names(cellGroups)[y]], 
			  nCellsUsed=length(uniq), 
			  nReplicates=length(x), 
			  nMin=nmin, 
			  nMax=nmax, 
			  maxPeaks = min(maxPeaks, length(uniq) * peaksPerCell)
			)
		}) %>% Reduce("rbind",.)

		ArchR:::.logDiffTime("Peak Calling Parameters!", tstart, verbose = verbose, logFile = logFile)
		ArchR:::.logThis(groupSummary, "PeakCallSummary", logFile = logFile)
		#printSummary <- groupSummary
		#rownames(printSummary) <- NULL
		if(verbose) print(groupSummary)

		#####################################################
		# Peak Calling from TileMatrix
		#####################################################

		#MatrixFiles
		ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]
		chrToRun <- ArchR:::.availableSeqnames(ArrowFiles, subGroup = useMatrix)
    featureDF <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix)
    
    #Determine Resolution
    d1 <- featureDF$start[2] - featureDF$start[1]
    d2 <- featureDF$start[3] - featureDF$start[2]
    if(d1 != d2){
    	ArchR:::.logMessage("Something is wrong with TileMatrix, Could not determine a resolution!", logFile = logFile)
    	stop("Something is wrong with TileMatrix, Could not determine a resolution!")
    }else{
    	res <- d1
    }

		#Compute Row Sums Across All Samples
		ArchR:::.logDiffTime("Computing Total Accessibility Across All Features", tstart, addHeader = FALSE, verbose = verbose)
		totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
  	ArchR:::.logThis(totalAcc, "PeakCallTiles-totalAcc", logFile=logFile)		
		nTiles <- nrow(totalAcc)
		gc()

		#Pre-Filter 0s
		topFeatures <- totalAcc[which(totalAcc$rowSums != 0), ]
  	ArchR:::.logThis(topFeatures, "PeakCallTiles-topFeatures", logFile=logFile)		

		#Group Matrix
		#Consider reading in group-wise if this is getting too large?
		ArchR:::.logDiffTime("Computing Pseudo-Grouped Tile Matrix", tstart, addHeader = FALSE, verbose = verbose)
    groupMat <- ArchR:::.getGroupMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = topFeatures,
      useMatrix = useMatrix, 
      threads = threads,
      groupList = unlist(cellGroups),
      useIndex = FALSE,
      asSparse = TRUE,
      verbose = FALSE
    )
  	ArchR:::.logThis(groupMat, "PeakCallTiles-groupMat", logFile=logFile)		

		ArchR:::.logDiffTime(sprintf("Created Pseudo-Grouped Tile Matrix (%s GB)", round(object.size(groupMat) / 10^9, 3)), tstart, addHeader = FALSE, verbose = verbose)
    expectation <- Matrix::colSums(groupMat) / nTiles
    ArchR:::.logMessage(paste0("colSums = ", Matrix::colSums(groupMat)), logFile = logFile)
    ArchR:::.logMessage(paste0("nTiles = ", nTiles), logFile = logFile)
    ArchR:::.logMessage(paste0("Expectation = ", expectation), logFile = logFile)

		#####################################################
		# Peak Calling Per Group
		#####################################################

    groupPeaks <- ArchR:::.safelapply(seq_along(cellGroups), function(i){

    	ArchR:::.logDiffTime(sprintf("Computing Peaks from Tile Matrix Group (%s of %s)", i, length(cellGroups)), tstart, addHeader = FALSE, verbose = verbose)

    	gx <- grep(paste0(names(cellGroups)[i],"."), colnames(groupMat))
	    
	    pMat <- lapply(seq_along(gx), function(x){
		    pval <- ppois(q = groupMat[,gx[x]]-1, lambda = expectation[gx[x]], lower.tail = FALSE, log=FALSE)
		    if(tolower(method) == "q"){
		    	pval <- p.adjust(pval, method = "fdr", n = nTiles)
		    }else if(tolower(method)=="p"){
		    }else{
		    	ArchR:::.logMessage("method should be either p for p-value or q for adjusted p-value!", logFile = logFile)
		    	stop("method should be either p for p-value or q for adjusted p-value!")
		    }
		    as(as.matrix(-log10(pmax(pval, 10^-250))), "dgCMatrix")
	    }) %>% Reduce("cbind", .)

	    n <- ncol(pMat)
	    passPeaks <- Matrix::rowSums(pMat >= -log10(cutOff)) >= eval(parse(text=reproducibility))
	    mlogp <- Matrix::rowSums(Matrix::t(Matrix::t(pMat) / Matrix::colSums(pMat)) * 10^6) / ncol(pMat)

	    rm(pMat)
	    gc()

	    passPeaks <- passPeaks[order(mlogp, decreasing=TRUE)]
	    mlogp <- mlogp[order(mlogp, decreasing=TRUE)]

	    nMax <- groupSummary[names(cellGroups)[i], "maxPeaks"]
	    passPeaks <- head(passPeaks, nMax)
	    mlogp <- head(mlogp, nMax)

	    mlogp <- mlogp[which(passPeaks)]
	    passPeaks <- passPeaks[which(passPeaks)]

	    if(length(passPeaks) > 0){
		    DataFrame(
		    	Group = Rle(names(cellGroups)[i]), 
		    	peaks = names(passPeaks), 
		    	mlog10p = mlogp, 
		    	normmlogp = 10^6 * mlogp / sum(mlogp)
		    )
	    }else{
	    	NULL
	    }

    }, threads = threads) %>% Reduce("rbind", .)

  	ArchR:::.logThis(groupPeaks, "PeakCallTiles-groupPeaks", logFile=logFile)

    groupPeaks <- groupPeaks[order(groupPeaks$normmlogp, decreasing=TRUE), ]

		#####################################################
		# BSgenome for Add Nucleotide Frequencies!
		#####################################################
		ArchR:::.requirePackage(genomeAnnotation$genome)
		ArchR:::.requirePackage("Biostrings",source="bioc")
		BSgenome <- eval(parse(text = genomeAnnotation$genome))
		BSgenome <- validBSgenome(BSgenome)
		outDir <- file.path(getOutputDirectory(ArchRProj), "PeakCalls")

		#####################################################
		# Create Group Peaks
		#####################################################
		ArchR:::.logDiffTime("Creating Group Peak Sets with Annotations!", tstart, verbose = verbose, logFile = logFile)

		peakDF <- ArchR:::.safelapply(seq_along(cellGroups), function(i){
			
			groupPeaksi <- groupPeaks[groupPeaks$Group == names(cellGroups)[i], ]

			if(nrow(groupPeaksi) == 0){
				return(NULL)
			}
			
			groupPeaksi <- cbind(
				topFeatures[as.integer(gsub("f", "", groupPeaksi$peaks)),],
				groupPeaksi
			)
			groupPeaksi <- groupPeaksi[order(groupPeaksi$seqnames, groupPeaksi$idx), ]
			groupPeaksi$peaks <- NULL
			groupPeaksi$normmlogp <- NULL
			groupPeaksi$mlog10p <- round(groupPeaksi$mlog10p, 3)
			groupPeaksGRi <- GRanges(groupPeaksi$seqnames, 
				IRanges(
					start = pmax((groupPeaksi$idx - 1) * res, 1),
					end = (groupPeaksi$idx) * res - 1
				)
			)
			mcols(groupPeaksGRi) <- groupPeaksi[,c("mlog10p", "Group")]

			groupPeaksGRi <- subsetByOverlaps(groupPeaksGRi, genomeAnnotation$chromSizes, type = "within")
			groupPeaksGRi <- subsetByOverlaps(groupPeaksGRi, genomeAnnotation$blacklist, invert = TRUE)
			groupPeaksGRi <- .fastAnnoPeaks(
				groupPeaksGRi, 
				BSgenome = BSgenome, 
				geneAnnotation = geneAnnotation, 
				promoterRegion = promoterRegion
			)
			groupPeaksGRi <- groupPeaksGRi[which(mcols(groupPeaksGRi)$N < 0.001)] #Remove N Containing Peaks
			mcols(groupPeaksGRi)$N <- NULL #Remove N Column

			#Table Peak Types
			tabPT <- data.frame(Group = names(cellGroups)[i], table(groupPeaksGRi$peakType))

			#Save
			saveRDS(groupPeaksGRi, file.path(outDir, paste0(make.names(names(cellGroups)[i]), "-reproduciblePeaks.gr.rds")))

			#Remove
			rm(groupPeaksGRi)
			gc()

			tabPT

		}, threads = threads) %>% Reduce("rbind", .)

		#####################################################
		# Call Union Peaks
		#####################################################
		ArchR:::.logDiffTime("Creating Union Peak Set with Annotations!", tstart, verbose = verbose, logFile = logFile)
	    unionPeaks <- groupPeaks[!duplicated(groupPeaks$peaks), ]
	    rm(groupPeaks)
	    gc()

	    unionPeaks <- cbind(
	    	topFeatures[as.integer(gsub("f", "", unionPeaks$peaks)),],
	    	unionPeaks
	    )
	    unionPeaks <- unionPeaks[order(unionPeaks$seqnames, unionPeaks$idx), ]
	    unionPeaks$peaks <- NULL
	    unionPeaks$normmlogp <- NULL
	    unionPeaks$mlog10p <- round(unionPeaks$mlog10p, 3)

	    unionPeaksGR <- GRanges(unionPeaks$seqnames, 
	    	IRanges(
	    		start = pmax((unionPeaks$idx - 1) * res, 1),
	    		end = (unionPeaks$idx) * res - 1
	    	)
	    )
	    mcols(unionPeaksGR) <- unionPeaks[,c("mlog10p", "Group")]

		unionPeaksGR <- subsetByOverlaps(unionPeaksGR, genomeAnnotation$chromSizes, type = "within")
		unionPeaksGR <- subsetByOverlaps(unionPeaksGR, genomeAnnotation$blacklist, invert = TRUE)
		unionPeaksGR <- ArchR:::.fastAnnoPeaks(
			unionPeaksGR, 
			BSgenome = BSgenome, 
			geneAnnotation = geneAnnotation, 
			promoterRegion = promoterRegion
		)
		unionPeaksGR <- unionPeaksGR[which(mcols(unionPeaksGR)$N < 0.001)] #Remove N Containing Peaks
		mcols(unionPeaksGR)$N <- NULL #Remove N Column

		#Set To unionPeaks
		unionPeaks <- unionPeaksGR
		rm(unionPeaksGR)
		gc()

		#Summarize Output
		peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
		peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
		peakDF$Freq <- peakDF$Freq / 1000
		metadata(unionPeaks)$PeakCallSummary <- peakDF


	}else{
		ArchR:::.logMessage("method not recognized! Supported methods are Macs2 or Tiles!", logFile = logFile)
		stop("method not recognized! Supported methods are Macs2 or Tiles!")
	}	

	#Add Peak Set
	ArchRProj <- addPeakSet(ArchRProj, unionPeaks, force = TRUE)

  if(plot){
    plotPDF(ArchR:::.plotPeakCallSummary(ArchRProj), name = "Peak-Call-Summary", width = 8, height = 5, ArchRProj = ArchRProj, addDOC = FALSE)
  }

	ArchR:::.logDiffTime(sprintf("Finished Creating Union Peak Set (%s)!", length(unionPeaks)), tstart, verbose = verbose, logFile = logFile)

	return(ArchRProj)

}


callSummitsOnCoverages_win <- function(
	i = NULL,
	coverageFiles = NULL,
	outFiles = NULL,
	peakParams = NULL,
	bedDir = NULL,
	excludeChr = NULL,
	subThreads = 1,
	tstart = NULL,
	logFile = NULL
	){
	
	ArchR:::.logDiffTime(sprintf("Group %s of %s, Calling Peaks with MACS2!", i, length(coverageFiles)), tstart, verbose = TRUE, logFile = logFile)

	################
	# Create Bed File from Coverage File
	################
	bedFile <- file.path(bedDir, paste0(make.names(basename(names(coverageFiles)[i])),"-",i,".insertions.bed"))
	o <- ArchR:::.writeCoverageToBed(coverageFiles[i], bedFile, excludeChr = excludeChr, logFile = logFile)
	peakParams$bedFile <- bedFile
	
	################
	# MACS2 Peak-Calling Leave Room For Other Options?
	################
	peakParams$logFile <- logFile
	summits <- do.call(callSummitsMACS2_win, peakParams)
	rmf <- file.remove(bedFile)
	
	################
	# Save output
	################
	saveRDS(summits, outFiles[i])

	#.logDiffTime(sprintf("Group %s of %s, Finished Calling Peaks with MACS2!", i, length(coverageFiles)), tstart, verbose = verbose, logFile = logFile)

	outFiles[i]

}

callSummitsMACS2_win <- function(
	i = NULL,
	coverageFiles = NULL,
	outFiles = NULL,
	peakParams = NULL,
	bedDir = NULL,
	excludeChr = NULL,
	subThreads = 1,
	tstart = NULL,
	logFile = NULL
	){
	
	ArchR:::.logDiffTime(sprintf("Group %s of %s, Calling Peaks with MACS2!", i, length(coverageFiles)), tstart, verbose = TRUE, logFile = logFile)

	################
	# Create Bed File from Coverage File
	################
	bedFile <- file.path(bedDir, paste0(make.names(basename(names(coverageFiles)[i])),"-",i,".insertions.bed"))
	o <- ArchR:::.writeCoverageToBed(coverageFiles[i], bedFile, excludeChr = excludeChr, logFile = logFile)
	peakParams$bedFile <- bedFile
	
	################
	# MACS2 Peak-Calling Leave Room For Other Options?
	################
	peakParams$logFile <- logFile
	summits <- do.call(callSummitsMACS2_win, peakParams)
	rmf <- file.remove(bedFile)
	
	################
	# Save output
	################
	saveRDS(summits, outFiles[i])

	#.logDiffTime(sprintf("Group %s of %s, Finished Calling Peaks with MACS2!", i, length(coverageFiles)), tstart, verbose = verbose, logFile = logFile)

	outFiles[i]

}

callSummitsMACS2_win <- function(
	bedFile = NULL,
	pathToMacs2 = "macs2",
	genomeSize = 2.7e9, 
	shift = -75, 
	extsize = 150, 
	cutOff = 0.05, 
	method = "q",
	additionalParams = "--nomodel --nolambda",
	logFile = NULL
	){

	stopifnot(tolower(method) %in% c("p","q"))
	stopifnot(!is.null(genomeSize))
	#utility <- .checkPath(pathToMacs2)

	#Output Files
	bedName <- gsub("\\.insertions.bed", "", bedFile)
	summitsFile <- paste0(bedName, "_summits.bed")
	narrowPeaksFile <- paste0(bedName, "_peaks.narrowPeak")
	xlsFile <- paste0(bedName, "_peaks.xls")

	#Create MACS2 Command
	#cmd <- sprintf("callpeak -g %s --name %s --treatment %s --outdir %s --format BED --call-summits --keep-dup all %s", 
	#	genomeSize, basename(bedName), bedFile, dirname(bedName), additionalParams)

	#if(!is.null(shift) & !is.null(extsize)){
	#	cmd <- sprintf("%s --shift %s --extsize %s", cmd , shift, extsize)
	#}

	#if(tolower(method) == "p"){
	#	cmd <- sprintf("%s -p %s", cmd , cutOff)
	#}else{
	#	cmd <- sprintf("%s -q %s", cmd , cutOff)
	#}

	#ArchR:::.logMessage(paste0("Running Macs2 with Params : macs2 ", cmd), logFile = logFile)

	#run <- system2(pathToMacs2, cmd, wait=TRUE, stdout=NULL, stderr=NULL)
	bedFile2<-gsub("\\\\","/",bedFile)
	mount_drive<-substring(bedFile2, 1, 1)
	bedFile2<-gsub(sprintf("*%s:",toupper(mount_drive)),sprintf("/mnt/%s",tolower(mount_drive)),bedFile2)
	bedName2<-gsub("\\\\","/",bedName)
	bedName2<-gsub(sprintf("*%s:",toupper(mount_drive)),sprintf("/mnt/%s",tolower(mount_drive)),bedName2)
	
	pathToMacs2_temp<-system("bash -c 'find /home -name macs2'", intern = TRUE)
	run <- system(sprintf('bash -c "%s callpeak -t %s -g %s --format BED --nomodel --nolambda --extsize 150 --shift -75 --call-summits --keep-dup all -q 0.05 --name %s --outdir %s" ',pathToMacs2_temp, bedFile2, genomeSize, basename(bedName2), dirname(bedName2)))

	#Read Summits!
	out <- data.table::fread(summitsFile, select = c(1,2,3,5))
	out <- GRanges(out$V1, IRanges(out$V2 + 1, out$V3), score = out$V5)

	#Remove Files
	r2 <- suppressWarnings(file.remove(summitsFile, narrowPeaksFile, xlsFile))

	return(out)

}

