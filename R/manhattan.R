#' Creates a manhattan plot using ggplot2 and returns the plot object.
#' 
#' @param x A data frame with association result data
#' @param y A data frame with annotation data
#' @param z A data frame with best hits data
#' @param snp SNP identifier column in data frame
#' @param chr Chromosome column in data frame
#' @param bp SNP position column in data frame
#' @param maf MAF column in data frame (ignored if NA)
#' @param p P-value column in data frame
#' @param typed the column in the data frame indicating whether the markers are genotyped or imputed (ignored if NA)
#' @param annotation the column in x data frame to use for annotation with text labels (ignored if NA)
#' @param category the column in y indicating the markers category to highlight in color
#' @param categoryColors list of the colors to use for the categories (ignored if NA)
#' @param categoryFlanking the flanking size in kbp (10 by default)
#' @param categoryMinP the worse log transformed p-value to consider for category annotation (5 by default)
#' @param thresholdLow the low threshold value (log10)
#' @param thresholdHigh the high threshold value (log10)
#' @param thresholdLowColor the color of the low threshold
#' @param thresholdHighColor the color of the high threshold
#' @param mafColor the color of the low maf values
#' @param build What build to use for plotting ('b37' or 'b38', default is 'b37')
#' @param title Title of plot (date by default, ignored if NA)
#' 
#' @return A manhattan plot (ggplot2 object)
#' 
#' @import ggplot2
#' @import ggrepel
#' 
#' @export

#devtools::use_package("ggplot2", "Suggests")
#devtools::use_package("ggrepel", "Suggests")

manhattan <- function(x, y = NA, z = NA, snp='SNP', chr='CHR', bp='BP', p='P', maf = NA, typed = NA, annotation = NA, 
                      category = "label", categoryColors = NA, categoryFlanking = 10, categoryMinP = 5, thresholdLow = 5, thresholdHigh = -log10(5e-8), 
                      thresholdLowColor = "blue", thresholdHighColor = "red", mafColor = "black", build = 'b37', title = Sys.time()){
  
  
  # Build specific variables
  
  chromosomeLength38 <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468) # bp length from Ensembl GRCh38.p10
  chromosomeLength37 <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566) # bp length from Ensembl GRCh37.p13
  if (build == 'b37') {
    chromosomeLength <- chromosomeLength37
  } else if (build == 'b38') {
    chromosomeLength <- chromosomeLength38
  } else {
    stop("Build not supported")
  }
  genomeLength <- sum(chromosomeLength)
  
  
  # Make a data frame for the plot
  
  manhattanData <- data.frame(snp = x[[snp]], chr = x[[chr]], bp = x[[bp]], p = x[[p]], stringsAsFactors = F)
  if (!is.na(maf)) {
    manhattanData$maf <- x[[maf]]
  }
  if (!is.na(typed)) {
    manhattanData$typed <- x[[typed]]
    if (!is.factor(manhattanData$typed)) {
      manhattanData$typed <- as.factor(manhattanData$typed)
    }
  }
  if (!is.na(annotation)) {
    manhattanData$annotation <- x[[annotation]]
  }
  manhattanData <- manhattanData[!is.na(manhattanData$p) & manhattanData$p > 0, ]
  manhattanData$logP <- -log10(manhattanData$p)
  
  
  # Create data frame for the annotation plot values
  
  if (length(y) > 1 || !is.na(y)) {
    annotationDataFrame <- data.frame(snp = y[[snp]], chr = y[[chr]], bp = y[[bp]], category = y[[category]], stringsAsFactors = F)
  }
  
  
  # Create data frame for the best hits plot values
  
  if (length(z) > 1 || !is.na(z)) {
    bestHitsDataFrame <- data.frame(chr = z[[chr]], bp = z[[bp]], stringsAsFactors = F)
  }
  
  # Sec color by category if available
  
  xValues <- c() # The position of every SNP on the x axis
  xBreak <- c() # The center of every chromosome
  xBreakLabels <- c() # The labels to use for every chromosome
  chrStart <- c() # The start of every second chromosome
  chrEnd <- c() # The end of every second chromosome
  if (length(z) > 1 || !is.na(z)) {
    bestHitsDataFrame$x <- 0 # The position of every gene
  }
  
  manhattanData <- manhattanData[order(manhattanData$chr, manhattanData$bp), ]
  
  yMax <- max(round(max(manhattanData$logP)+1), thresholdHigh)
  
  xOffset <- 0
  start <- T
  for (chromosomeNumber in 1:22) {
    
    bpTemp <- manhattanData$bp[manhattanData$chr == chromosomeNumber]
    xTemp <- bpTemp + xOffset
    breakValue <- xTemp[1] + (xTemp[length(xTemp)] - xTemp[1]) / 2
    xBreak <- c(xBreak, breakValue)
    if (chromosomeNumber < 12 || chromosomeNumber %% 2 == 0) {
      xBreakLabels <- c(xBreakLabels, chromosomeNumber)
    } else {
      xBreakLabels <- c(xBreakLabels, "")
    }
    xValues <- c(xValues, xTemp)
    
    if (length(y) > 1 || !is.na(y)) {
      
      bpTemp <- annotationDataFrame$bp[annotationDataFrame$chr == chromosomeNumber]
      xTemp <- bpTemp + xOffset
      startTemp <- xTemp - categoryFlanking * 1000
      annotationDataFrame$xStart[annotationDataFrame$chr == chromosomeNumber] <- ifelse(startTemp < 0, 0, startTemp)
      endTemp <- xTemp + categoryFlanking * 1000
      annotationDataFrame$xEnd[annotationDataFrame$chr == chromosomeNumber] <- ifelse(endTemp > genomeLength, genomeLength, endTemp)
    
    }
    
    if (length(z) > 1 || !is.na(z)) {
      bestHitsDataFrame$x[bestHitsDataFrame$chr == chromosomeNumber] <- bestHitsDataFrame$bp[bestHitsDataFrame$chr == chromosomeNumber] + xOffset
    }
    
    
    xOffset <- xOffset + chromosomeLength[chromosomeNumber]
    if (start) {
      chrStart <- c(chrStart, xOffset)
    } else {
      chrEnd <- c(chrEnd, xOffset)
    }
    start <- !start
    
  }
  manhattanData$xValues <- xValues
  
  
  if (length(y) > 1 || !is.na(y)) {
    
    manhattanData$category <- ""
    
    for (i in 1:nrow(annotationDataFrame)) {
      
      iChr <- annotationDataFrame$chr[i]
      iStart <- annotationDataFrame$xStart[i]
      iEnd <- annotationDataFrame$xEnd[i]
      iCategory <- annotationDataFrame$category[i]
      
      if (sum(is.na(manhattanData$logP)) > 0) {
        stop("logP NA")
      }
      if (sum(is.na(manhattanData$chr)) > 0) {
        stop("Chr NA")
      }
      if (sum(is.na(manhattanData$xValues)) > 0) {
        stop("x NA")
      }
      if (sum(is.na(iStart)) > 0) {
        stop("iStart NA")
      }
      if (sum(is.na(iEnd)) > 0) {
        stop("iEnd NA")
      }
      
      snpInWindow <- manhattanData$chr == iChr & manhattanData$xValues >= iStart & manhattanData$xValues <= iEnd
      
      if (sum(is.na(snpInWindow)) > 0) {
        stop("snp NA")
      }
      
      snpOverThreshold <- snpInWindow & manhattanData$logP >= categoryMinP
      
      if (sum(is.na(snpOverThreshold)) > 0) {
        stop("snp NA")
      }
      
      if (sum(snpOverThreshold) > 0) {
        
        manhattanData$category[snpInWindow] <- iCategory
        
      }
      
    }
    
    manhattanData$category <- factor(manhattanData$category)
    
  }
  
  
  # Order by category and maf to see common markers in front
  
  if (!is.na(maf)) {
    manhattanData <- manhattanData[order(manhattanData$maf), ]
  }
  if (length(y) > 1 || !is.na(y)) {
    manhattanData <- manhattanData[order(manhattanData$category, na.last = F), ]
  }
  
  
  # Make a ggplot object
  
  manhattanPlot <- ggplot()
  
  
  # Add background rectangle for every second chromosome
  
  manhattanPlot <- manhattanPlot + geom_rect(aes(xmin = chrStart, xmax = chrEnd, ymin = 0, ymax = yMax), alpha = 0.2)
  
  
  # Add line for every best hit
  
  if (length(z) > 1 || !is.na(z)) {
    manhattanPlot <- manhattanPlot + geom_vline(aes(xintercept = bestHitsDataFrame$x), col = "black", linetype = "dotted")
  }
  
  
  # Plot all markers
  
  if (!is.na(typed)) {
    if (!is.na(maf)) {
      if (length(y) > 1 || !is.na(y)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf, size = typed, col = category), shape = 21)
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf, size = typed), col = "black", shape = 21)
      }
    } else {
      if (length(y) > 1 || !is.na(y)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, size = typed, col = category))
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, size = typed), col = "black")
      }
    }
    manhattanPlot <- manhattanPlot + scale_size_manual(name = "", values = c(2, 1))
  } else {
    if (!is.na(maf)) {
      if (length(y) > 1 || !is.na(y)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf, col = category), shape = 21)
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf), col = "black", shape = 21)
      }
    } else {
      if (length(y) > 1 || !is.na(y)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, col = category))
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP), col = "black")
      }
    }
  }
  if (!is.na(maf)) {
    
    manhattanPlot <- manhattanPlot + scale_fill_gradientn(name = "MAF", colors = c("white", mafColor))
    
  }
  if (length(y) > 1 || !is.na(y)) {
    
    categoryColorsTemp <- levels(manhattanData$category)
    
    if (length(categoryColors) > 1 || !is.na(categoryColors)) {
      
      if ("" %in% categoryColorsTemp) {
        
        categoryColorsTemp <- c("black", categoryColors)
        
      } else {
        
        categoryColorsTemp <- categoryColors
        
      }
      
      
    } else {
      
      if ("" %in% categoryColorsTemp) {
        
        categoryColorsTemp <- c("black")
        
        colorOffset <- 1
        
      } else {
        
        categoryColorsTemp <- c()
        
        colorOffset <- 0
        
      }
      
      if (length(categoryColorsTemp) > colorOffset) {
        
        categoryColorsTemp <- c(categoryColorsTemp, scales::hue_pal()(length(categoryLevels)-colorOffset))
        
      } else {
        
        manhattanPlot <- manhattanPlot + guides(color = 'none')
        
      }
      
    }
    
    
    manhattanPlot <- manhattanPlot + scale_color_manual(name = "", values = categoryColorsTemp)
    
  }
  
  
  # Plot thresholds
  
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdLow), col = thresholdLowColor)
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdHigh), col = thresholdHighColor)
  
  
  # Set axes labels 
  
  manhattanPlot <- manhattanPlot + scale_y_continuous(name = "-log10(p)", breaks = 0:yMax, limits = c(0, yMax), expand = c(0, 0))
  manhattanPlot <- manhattanPlot + scale_x_continuous(name = NULL, breaks = xBreak, label = xBreakLabels, expand = c(0.01, 0), limits = c(0, genomeLength))
  
  
  # Format background and grid
  
  manhattanPlot <- manhattanPlot + theme(panel.background = element_rect(fill = NA, colour = "grey50"),
                                         panel.grid.major.y = element_line(colour = "grey50", linetype = "dotted"),
                                         panel.grid.minor.y = element_blank(),
                                         panel.grid.major.x = element_blank(),
                                         panel.grid.minor.x = element_blank(),
                                         legend.box.background = element_rect(fill = NA, colour = NA))
  
  # Add annotation if provided
  
  if (!is.na(annotation)) {
    
    manhattanData$annotation[is.na(manhattanData$annotation)] <- ""
    
    nudge_p <- yMax - manhattanData$logP - 0.5
    
    manhattanPlot <- manhattanPlot + geom_text_repel(data = manhattanData, aes(x = xValues, y = logP, label=annotation), nudge_y = nudge_p)
  
  }
  
  # Add title to plot if provided
  
  if (!is.na(title)) {
    manhattanPlot <-manhattanPlot + ggtitle(title)
  }
  
  # Return the plot
  
  return(manhattanPlot)
}

