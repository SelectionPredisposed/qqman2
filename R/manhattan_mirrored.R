#' Creates a mirrored manhattan plot using ggplot2 and returns the plot object.
#' 
#' @param x1 A data frame with association result data to be plotted on the upper side
#' @param x2 A data frame with association result data to be plotted on the lower side
#' @param x1Name the name of the first association result
#' @param x2Name the name of the second association result
#' @param y A data frame with annotation data to color markers in an area and annotate the top SNP (ignored if NA)
#' @param z A data frame with reference hits to be annotated using vertical lines (ignored if NA)
#' @param x.snp SNP identifier column in the x1 and x2 data frames
#' @param x.chr Chromosome column in the x1 and x2 data frames
#' @param x.bp SNP position column in the x1 and x2 data frames
#' @param x.maf MAF column in the x1 and x2 data frames (ignored if NA)
#' @param x.p P-value column in the x1 and x2 data frames
#' @param x.typed the column in the x1 and x2 data frames indicating whether the markers are genotyped or imputed (ignored if NA)
#' @param x.annotation the column in the x1 and x2 data frames indicating used to annotate markers with text (ignored if NA)
#' @param x.color the color used to annotate high maf values (black by default, ignored if x.maf is NA)
#' @param y.snp SNP identifier column in the y data frame
#' @param y.chr Chromosome column in the y data frame
#' @param y.bp SNP position column in the y data frame
#' @param y.category the column in y indicating the markers category to highlight in color
#' @param y.name the column in y indicating the names to use to annotate the best hit in the colored region (ignored if NA)
#' @param y.colors the colors to use to annotate the different categories
#' @param y.flanking the flanking to use for coloring in kbp (50 by default)
#' @param y.minP the the minimal p-value a category has to reach to be included in the coloring (5 by default)
#' @param z.chr Chromosome column in the z data frame
#' @param z.bp SNP position column in the z data frame
#' @param z.name the column in z indicating the names to use to annotate the markers of z with text (ignored if NA)
#' @param thresholdLow the low threshold value (log10)
#' @param thresholdHigh the high threshold value (log10)
#' @param thresholdLowColor the color of the low threshold
#' @param thresholdHighColor the color of the high threshold
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

manhattan_mirrored <- function(x1, x2, x1.name, x2.name, y = NA, z = NA, 
                               x.snp='SNP', x.chr='CHR', x.bp='BP', x.p='P', x.maf = NA, x.typed = NA, x.annotation = NA, x.color = "black", 
                               y.snp='SNP', y.chr='CHR', y.bp='BP', y.category = "category", y.name = NA, y.colors = NA, y.flanking = 50, y.minP = 5, 
                               z.chr='CHR', z.bp='BP', z.name = "name",  
                      thresholdLow = 5, thresholdHigh = -log10(5e-8), thresholdLowColor = "blue", thresholdHighColor = "red", build = 'b37', title = Sys.time()){
  
  
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
  
  manhattanData1 <- data.frame(snp = x1[[x.snp]], chr = x1[[x.chr]], bp = x1[[x.bp]], p = x1[[x.p]], stringsAsFactors = F)
  if (!is.na(x.maf)) {
    manhattanData1$maf <- x1[[x.maf]]
  }
  if (!is.na(x.typed)) {
    manhattanData1$typed <- x1[[x.typed]]
    if (!is.factor(manhattanData1$typed)) {
      manhattanData1$typed <- as.factor(manhattanData1$typed)
    }
  }
  if (!is.na(x.annotation)) {
    manhattanData1$annotation <- x1[[x.annotation]]
  }
  manhattanData1 <- manhattanData1[!is.na(manhattanData1$p) & manhattanData1$p > 0, ]
  manhattanData1$logP <- -log10(manhattanData1$p)
  
  manhattanData2 <- data.frame(snp = x2[[x.snp]], chr = x2[[x.chr]], bp = x2[[x.bp]], p = x2[[x.p]], stringsAsFactors = F)
  if (!is.na(x.maf)) {
    manhattanData2$maf <- x2[[x.maf]]
  }
  if (!is.na(x.typed)) {
    manhattanData2$typed <- x2[[x.typed]]
    if (!is.factor(manhattanData2$typed)) {
      manhattanData2$typed <- as.factor(manhattanData2$typed)
    }
  }
  if (!is.na(x.annotation)) {
    manhattanData2$annotation <- x2[[x.annotation]]
  }
  manhattanData2 <- manhattanData2[!is.na(manhattanData2$p) & manhattanData2$p > 0, ]
  manhattanData2$logP <- log10(manhattanData2$p)
  
  manhattanData <- rbind(manhattanData1, manhattanData2)
  
  
  # Create data frame for the annotation plot values
  
  if (length(y) > 1 || !is.na(y)) {
    
    annotationDataFrame <- data.frame(snp = y[[y.snp]], chr = y[[y.chr]], bp = y[[y.bp]], category = y[[y.category]], stringsAsFactors = F)
    
    if (!is.na(y.name)) {
      
      annotationDataFrame$id <- y[[y.name]]
      annotationDataFrame$annotationX <- 0
      annotationDataFrame$annotationP <- 0
      
    }
    
  }
  
  
  # Create data frame for the best hits plot values
  
  if (length(z) > 1 || !is.na(z)) {
    
    bestHitsDataFrame <- data.frame(chr = z[[z.chr]], bp = z[[z.bp]], stringsAsFactors = F)
    
    if (!is.na(z.name)) {
      
      bestHitsDataFrame$id <- z[[z.name]]
      
    }
  }
  
  
  # Map chromosomic coordinates to the x axis
  
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
  yPMax <- yMax
  
  if (length(z) > 1 || !is.na(z)) {
    
    yAnnotation <- yMax +0.5
    
    yMax <- yMax + 1
    
  }
  
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
      startTemp <- xTemp - y.flanking * 1000
      annotationDataFrame$xStart[annotationDataFrame$chr == chromosomeNumber] <- ifelse(startTemp < 0, 0, startTemp)
      endTemp <- xTemp + y.flanking * 1000
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
      
      snpInWindow <- manhattanData$chr == iChr & manhattanData$xValues >= iStart & manhattanData$xValues <= iEnd
      
      snpOverThreshold <- snpInWindow & abs(manhattanData$logP) >= y.minP
      
      if (sum(snpOverThreshold) > 0) {
        
        manhattanData$category[snpInWindow] <- iCategory
        
        if (!is.na(y.name)) {
          
          windowData <- manhattanData[snpInWindow, ]
          
          maxP <- max(abs(windowData$logP))
          if (maxP != max(windowData$logP)) {
            maxP <- -maxP
          }
          
          j <- round(median(which(windowData$logP == maxP)))
          annotationDataFrame$annotationX[i] <- windowData$xValues[j]
          annotationDataFrame$annotationP[i] <- maxP
          
        }
        
      } else if (!is.na(y.name)) {
        
        annotationDataFrame$id[i] <- ""
        
      }
      
    }
    
    manhattanData$category <- factor(manhattanData$category)
    
  }
  
  
  # Order by category and maf to see common markers in front
  
  if (!is.na(x.maf)) {
    manhattanData <- manhattanData[order(manhattanData$maf), ]
  }
  if (length(y) > 1 || !is.na(y)) {
    manhattanData <- manhattanData[order(manhattanData$category, na.last = F), ]
  }
  
  
  # Make a ggplot object
  
  manhattanPlot <- ggplot()
  
  
  # Add background rectangle for every second chromosome
  
  manhattanPlot <- manhattanPlot + geom_rect(aes(xmin = chrStart, xmax = chrEnd, ymin = -yPMax, ymax = yMax), alpha = 0.2)
  
  
  # Add line for every best hit
  
  if (length(z) > 1 || !is.na(z)) {
    
    manhattanPlot <- manhattanPlot + geom_vline(aes(xintercept = bestHitsDataFrame$x), col = "black", linetype = "dotted")
    
  }
  
  
  # Plot all markers
  
  if (!is.na(x.typed)) {
    if (!is.na(x.maf)) {
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
    if (!is.na(x.maf)) {
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
  if (!is.na(x.maf)) {
    manhattanPlot <- manhattanPlot + scale_fill_gradientn(name = "MAF", colors = c("white", x.color))
  }
  if (length(y) > 1 || !is.na(y)) {
    
    categoryColorsTemp <- levels(manhattanData$category)
    
    if (length(y.colors) > 1 || !is.na(y.colors)) {
      
      if ("" %in% categoryColorsTemp) {
        
        categoryColorsTemp <- c("black", y.colors)
        
      } else {
        
        categoryColorsTemp <- y.colors
        
      }
      
      
    } else {
      
      if ("" %in% categoryColorsTemp) {
        
        categoryColorsTemp <- c("black")
        
        colorOffset <- 1
        
      } else {
        
        categoryColorsTemp <- c()
        
        colorOffset <- 0
        
      }
      
      if (length(categoryLevels) > colorOffset) {
        
        categoryColorsTemp <- c(categoryColorsTemp, scales::hue_pal()(length(categoryLevels)-colorOffset))
        
      }
      
    }
    
    
    manhattanPlot <- manhattanPlot + scale_color_manual(name = "", values = categoryColorsTemp)
    
  }
  
  
  # Plot thresholds
  
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdLow), col = thresholdLowColor)
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdHigh), col = thresholdHighColor)
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = -thresholdLow), col = thresholdLowColor)
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = -thresholdHigh), col = thresholdHighColor)
  
  
  # Set axes labels 
  
  axisLabel <- paste(x2.name, " | ", x1.name, " [-log10(p)]")
  manhattanPlot <- manhattanPlot + scale_y_continuous(name = axisLabel, breaks = -yPMax:yPMax, limits = c(-yPMax, yMax), expand = c(0, 0))
  manhattanPlot <- manhattanPlot + scale_x_continuous(name = NULL, breaks = xBreak, label = xBreakLabels, expand = c(0.01, 0), limits = c(0, genomeLength))
  
  
  # Format background and grid
  
  manhattanPlot <- manhattanPlot + theme(panel.background = element_rect(fill = NA, colour = "grey50"),
                                         panel.grid.major.y = element_line(colour = "grey50", linetype = "dotted"),
                                         panel.grid.minor.y = element_blank(),
                                         panel.grid.major.x = element_blank(),
                                         panel.grid.minor.x = element_blank(),
                                         legend.box.background = element_rect(fill = NA, colour = NA))
  
  # Add annotation if provided
  
  if (!is.na(x.annotation)) {
    
    tempDataFrame <- manhattanData[!is.na(manhattanData$annotation) | manhattanData$logP > quantile(manhattanData$logP, 0.9), ]
    
    tempDataFrame$annotation[is.na(tempDataFrame$annotation)] <- ""
    
    nudge_p <- yMax - tempDataFrame$logP - 0.5
    
    manhattanPlot <- manhattanPlot + geom_text_repel(data = tempDataFrame, aes(x = xValues, y = logP, label=annotation), nudge_y = nudge_p)
    
  }
  
  if ((length(y) > 1 || !is.na(y)) && !is.na(y.name)) {
    
    if (sum(is.na(annotationDataFrame$annotationX)) > 0) {
      
      stop("Null in y genomic coordinates mapping.")
      
    }
    
    if (sum(is.na(annotationDataFrame$annotationP)) > 0) {
      
      stop("Null in y best P.")
      
    }
    
    manhattanPlot <- manhattanPlot + geom_text_repel(data = annotationDataFrame, aes(x = annotationX, y = annotationP, label=id), nudge_y = 1)
    
  }
  
  if ((length(z) > 1 || !is.na(z)) && !is.na(z.name)) {
    
    if (sum(is.na(bestHitsDataFrame$x)) > 0) {
      
      stop("Null in z genomic coordinates mapping.")
      
    }
    
    if (sum(is.na(yAnnotation)) > 0) {
      
      stop("Null in z y value.")
      
    }
    
    if (sum(is.na(bestHitsDataFrame$id)) > 0) {
      
      stop("Null in z id value.")
      
    }
    
    manhattanPlot <- manhattanPlot + geom_text_repel(data = bestHitsDataFrame, aes(x = x, y = yAnnotation, label = id), nudge_y = 0, point.padding = NA)
    
  }
  
  # Add title to plot if provided
  
  if (!is.na(title)) {
    
    manhattanPlot <-manhattanPlot + ggtitle(title)
    
  }
  
  # Return the plot
  
  return(manhattanPlot)
}

