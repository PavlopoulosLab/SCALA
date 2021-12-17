addPeak2GeneLinks_shiny <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  addEmpiricalPval = FALSE,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
  ){

  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  ArchR:::.validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  ArchR:::.validInput(input = k, name = "k", valid = c("integer"))
  ArchR:::.validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  ArchR:::.validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  ArchR:::.validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "addPeak2GeneLinks Input-Parameters", logFile = logFile)

  ArchR:::.logDiffTime(main="Getting Available Matrices", t1=tstart, verbose=verbose, logFile=logFile)
  AvailableMatrices <- getAvailableMatrices(ArchRProj)

  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }

  if(useMatrix %ni% AvailableMatrices){
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }

  ArrowFiles <- getArrowFiles(ArchRProj)

  tstart <- Sys.time()

  dfAll <- ArchR:::.safelapply(seq_along(ArrowFiles), function(x){
    cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames")))
    pSx <- tryCatch({
      h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    }, error = function(e){
      if(getArchRVerbose()) message("No predictionScore found. Continuing without predictionScore!")
      rep(9999999, length(cNx))
    })
    DataFrame(
      cellNames = cNx,
      predictionScore = pSx
    )
  }, threads = threads) %>% Reduce("rbind", .)

  ArchR:::.logDiffTime(
    sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
    sum(dfAll[,2] < predictionCutoff), 
    nrow(dfAll), 
    round(sum(dfAll[,2] < predictionCutoff) / nrow(dfAll), 3)
    ), t1=tstart, verbose=verbose, logFile=logFile)

  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  ArchR:::.logThis(peakSet, "peakSet", logFile = logFile)

  #Gene Info
  geneSet <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  ArchR:::.logThis(geneStart, "geneStart", logFile = logFile)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  ArchR:::.logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  ArchR:::.logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  ArchR:::.logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)

  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)

  #Group Matrix RNA
  ArchR:::.logDiffTime(main="Getting Group RNA Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatRNA <- ArchR:::.getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  rawMatRNA <- groupMatRNA
  ArchR:::.logThis(groupMatRNA, "groupMatRNA", logFile = logFile)

  #Group Matrix ATAC
  ArchR:::.logDiffTime(main="Getting Group ATAC Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatATAC <- ArchR:::.getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  rawMatATAC <- groupMatATAC
  ArchR:::.logThis(groupMatATAC, "groupMatATAC", logFile = logFile)

  ArchR:::.logDiffTime(main="Normalizing Group Matrices", t1=tstart, verbose=verbose, logFile=logFile)

  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo

  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }

  names(geneStart) <- NULL

  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj
  ArchR:::.logThis(seRNA, "seRNA", logFile = logFile)

  names(peakSet) <- NULL

  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
    rowRanges = peakSet
  )
  metadata(seATAC)$KNNList <- knnObj
  ArchR:::.logThis(seATAC, "seATAC", logFile = logFile)

  rm(groupMatRNA, groupMatATAC)
  gc()

  #Overlaps
  ArchR:::.logDiffTime(main="Finding Peak Gene Pairings", t1=tstart, verbose=verbose, logFile=logFile)
  o <- DataFrame(
    findOverlaps(
      ArchR:::.suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), 
      ignore.strand = TRUE
    )
  )

  #Get Distance from Fixed point A B 
  o$distance <- GenomicRanges::distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] )
  colnames(o) <- c("B", "A", "distance")

  #Null Correlations
  if(addEmpiricalPval){
    ArchR:::.logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
    nullCor <- ArchR:::.getNullCorrelations(seATAC, seRNA, o, 1000)
  }

  ArchR:::.logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")  
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart

  if(addEmpiricalPval){
    out$EmpPval <- 2*pnorm(-abs(((out$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }
  
  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  ArchR:::.safeSaveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  ArchR:::.safeSaveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA

  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out

  ArchR:::.logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
  ArchR:::.endLogging(logFile = logFile)

  ArchRProj

}
determineOverlapCpp <- function(m, overlapCut) {    .Call('_ArchR_determineOverlapCpp', PACKAGE = 'ArchR', m, overlapCut)}



