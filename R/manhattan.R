#' Creates a manhattan plot using ggplot2 and returns the plot object.
#' 
#' @param x A data frame with result data
#' @param snp SNP column in data frame
#' @param chr Chromosome column in data frame
#' @param pb SNP position column in data frame
#' @param maf MAF column in data frame (ignored if NA)
#' @param category the column in the data frame indicating the markers category (ignored if NA)
#' @param p P-value column in data frame
#' @param typed the column in the data frame indicating whether the markers are genotyped or imputed (ignored if NA)
#' @param annotation a vector of annotation (ignored if NA)
#' @param thresholdLow the low threshold value (log10)
#' @param thresholdHigh the high threshold value (log10)
#' @param thresholdLowColor the color of the low threshold
#' @param thresholdHighColor the color of the high threshold
#' @param mafColor the color of the low maf values
#' @param categoryColor the color of the category provided (ggplot default if NA)
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

manhattan <- function(x, snp='SNP', chr='CHR', bp='BP', p='P', maf = 'MAF', category = NA, typed = NA, annotation = NA, 
                      thresholdLow = 5, thresholdHigh = -log10(5e-8), thresholdLowColor = "blue", thresholdHighColor = "red", 
                      mafColor = "black", categoryColors = NA, build = 'b37', title = Sys.time()){
  
  
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
  
  
  # Make a data frame for the plot
  
  manhattanData <- data.frame(snp = x[[snp]], chr = x[[chr]], bp = x[[bp]], p = x[[p]], stringsAsFactors = F)
  manhattanData$logP <- -log10(manhattanData$p)
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
  if (!is.na(category)) {
    manhattanData$category <- x[[category]]
    if (sum(is.na(manhattanData$category)) > 0) {
      manhattanData$category[is.na(manhattanData$category)] <- ""
    }
    if (!is.factor(manhattanData$category)) {
      manhattanData$category <- as.factor(manhattanData$category)
    }
    if (length(categoryColors) > 1 || !is.na(categoryColors)) {
      nColorsFound <- length(categoryColors)
      nColorsNeeded <- length(levels(manhattanData$category))
      if (nColorsFound == nColorsNeeded-1) {
        categoryColors <- c("black", categoryColors)
      }
      if (length(categoryColors) != length(levels(manhattanData$category))) {
        stop(paste(nColorsFound, " colors provided for ", nColorsNeeded, " category.", sep = ""))
      }
    }
  }
  
  
  # create vectors for horizontal annotation
  
  xValues <- c() # The position of every SNP on the x axis
  xBreak <- c() # The center of every chromosome
  xBreakLabels <- c() # The labels to use for every chromosome
  xStart <- c() # The start of every second chromosome
  xEnd <- c() # The end of every second chromosome
  
  manhattanData <- manhattanData[order(manhattanData$chr, manhattanData$bp), ]
  
  yMax <- max(round(max(manhattanData$logP)+1), thresholdHigh)
  
  xOffset <- 0
  start <- T
  for (chr in 1:22) {
    bpTemp <- manhattanData$bp[manhattanData$chr == chr]
    xTemp <- bpTemp + xOffset
    breakValue <- xTemp[1] + (xTemp[length(xTemp)] - xTemp[1]) / 2
    xBreak <- c(xBreak, breakValue)
    if (chr < 12 || chr %% 2 == 0) {
      xBreakLabels <- c(xBreakLabels, chr)
    } else {
      xBreakLabels <- c(xBreakLabels, "")
    }
    xValues <- c(xValues, xTemp)
    xOffset <- xOffset + chromosomeLength[chr]
    if (start) {
      xStart <- c(xStart, xOffset)
    } else {
      xEnd <- c(xEnd, xOffset)
    }
    start <- !start
  }
  manhattanData$xValues <- xValues
  
  
  # Order by maf to see common markers in front
  
  if (!is.na(maf)) {
    manhattanData <- manhattanData[order(manhattanData$maf), ]
  }
  
  
  # Make a ggplot object
  
  manhattanPlot <- ggplot()
  
  
  # Add background rectangle for every second chromosome
  
  manhattanPlot <- manhattanPlot + geom_rect(aes(xmin = xStart, xmax = xEnd, ymin = 0, ymax = yMax), alpha = 0.2)
  
  
  # Plot all markers
  
  if (!is.na(typed)) {
    if (!is.na(maf)) {
      if (!is.na(category)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf, size = typed, col = category), shape = 21)
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf, size = typed), col = "black", shape = 21)
      }
    } else {
      if (!is.na(category)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, size = typed, col = category))
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, size = typed), col = "black")
      }
    }
    manhattanPlot <- manhattanPlot + scale_size_manual(name = "", values = c(2, 1))
  } else {
    if (!is.na(maf)) {
      if (!is.na(category)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf, col = category), shape = 21)
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, fill = maf), col = "black", shape = 21)
      }
    } else {
      if (!is.na(category)) {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP, col = category))
      } else {
        manhattanPlot <- manhattanPlot + geom_point(data = manhattanData, aes(x = xValues, y = logP), col = "black")
      }
    }
  }
  if (!is.na(maf)) {
    manhattanPlot <- manhattanPlot + scale_fill_gradientn(name = "Maf", colors = c("white", mafColor))
  }
  if (!is.na(category) && length(categoryColors) > 1 || !is.na(categoryColors)) {
    manhattanPlot <- manhattanPlot + scale_color_manual(values = categoryColors)
  }
  
  
  # Plot thresholds
  
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdLow), col = thresholdLowColor)
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdHigh), col = thresholdHighColor)
  
  
  # Set axes labels 
  
  manhattanPlot <- manhattanPlot + scale_y_continuous(name = "-log10(p)", breaks = 0:yMax, limits = c(0, yMax), expand = c(0, 0))
  manhattanPlot <- manhattanPlot + scale_x_continuous(name = NULL, breaks = xBreak, label = xBreakLabels, expand = c(0.01, 0.01))
  
  
  # Format background and grid
  
  manhattanPlot <- manhattanPlot + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                                         panel.grid.major.y = element_line(colour = "grey50", linetype = "dotted"),
                                         panel.grid.minor.y = element_blank(),
                                         panel.grid.major.x = element_blank(),
                                         panel.grid.minor.x = element_blank())
  
  # Add annotation if provided
  
  if (!is.na(annotation)) {
    manhattanPlot <- manhattanPlot + geom_text_repel(data = manhattanData[!is.na(manhattanData$annotation), ], aes(x = x, y = logP, label=annotation))
  }
  
  # Add title to plot if provided
  
  if (!is.na(title)) {
    manhattanPlot <-manhattanPlot + ggtitle(title)
  }
  
  # Return the plot
  
  return(manhattanPlot)
}

