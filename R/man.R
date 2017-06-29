#' Creates a manhattan plot
#' @param x A data frame with result data
#' @param SNP SNP column in data frame
#' @param CHR Chromosome column in data frame
#' @param BP SNP position column in data frame
#' @param P P-value column in data frame
#' @param build What build to use for plotting (default to 'b37')
#' @param highlight Whether to highlight SNPs (default to FALSE)
#' @param highlight.col Vector containing SNPs to label
#' @param title Title of plot
#' 
#' @return A manhattan plot (ggplot2 object)
#' 
#' @import ggplot2
#' @import ggrepel
#' 
#' @export

#devtools::use_package("ggplot2", "Suggests")
#devtools::use_package("data.table", "Suggests")
#devtools::use_package("ggrepel", "Suggests")

#res <- fread('/home/oyvind/Documents/manhattan-tool/prunedman', header = T, stringsAsFactors = F, data.table = F)
#hl <- c('rs72921490','rs4478530')
manh <- function(x, SNP='SNP', CHR='CHR', BP='BP', P='P', 
                 build='b37', highlight=F, highlight.snps=NA, 
                 highlight.col='SNP', title=Sys.time()){
  
  # Load data
  manhattanData <- x
  
  # Plot settings
  thresholdLow <- 5
  thresholdHigh <- -log10(5e-8)
  chromosomeLength38 <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468) # bp length from Ensembl GRCh38.p10
  chromosomeLength37 <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566) # bp length from Ensembl GRCh37.p13
  #yMax <- max(round(max(manhattanData[[P]])+1), thresholdHigh)
  yMax <- max(round(-log10(min(res[['frequentist_add_pvalue']])) +1), thresholdHigh)
  
  if (build == 'b37') {
    chromosomeLength <- chromosomeLength37
  } else if (build == 'b38') {
    chromosomeLength <- chromosomeLength38
  } else {
    stop("Build not supported")
  }
  
  xOffset <- 0
  x <- c()
  xBreak <- c()
  xBreakLabels <- c()
  for (chr in 1:22) {
    bpTemp <- manhattanData[[BP]][manhattanData[[CHR]] == chr]
    xTemp <- bpTemp + xOffset
    breakValue <- xTemp[1] + (xTemp[length(xTemp)] - xTemp[1]) / 2
    xBreak <- c(xBreak, breakValue)
    if (chr < 12 || chr %% 2 == 0) {
      xBreakLabels <- c(xBreakLabels, chr)
    } else {
      xBreakLabels <- c(xBreakLabels, "")
    }
    x <- c(x, xTemp)
    xOffset <- xOffset + chromosomeLength[chr]
  }
  manhattanData$x <- x
  
  manhattanData$col <- as.factor(manhattanData$chr %% 2)

  # Create manhattan plot
  
  manhattanPlot <- ggplot()
  manhattanPlot <- manhattanPlot + geom_point(aes(x = manhattanData$x, y = -log10(manhattanData[[P]]), col = manhattanData$col))
  manhattanPlot <- manhattanPlot + scale_color_manual(values = c("slategray4", "midnightblue"))
  
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdLow), col = "blue")
  manhattanPlot <- manhattanPlot + geom_hline(aes(yintercept = thresholdHigh), col = "green")
   
  manhattanPlot <- manhattanPlot + scale_y_continuous(name = "-log10(p)", breaks = 0:yMax, limits = c(0, yMax), expand = c(0, 0))
  manhattanPlot <- manhattanPlot + scale_x_continuous(name = NULL, breaks = xBreak, label = xBreakLabels, expand = c(0.01, 0.01))
  manhattanPlot <- manhattanPlot + theme(legend.position="none", 
                                          panel.background = element_rect(fill = "white", colour = "gray50"),
                                          panel.grid.major.y = element_line(colour = "grey50", linetype = "dotted"),
                                          panel.grid.minor.y = element_blank(),
                                          panel.grid.major.x = element_blank(),
                                          panel.grid.minor.x = element_blank())
  
  # Add labels to plot if specified
  if(highlight==T) {
    manhattanPlot <-manhattanPlot + geom_text_repel(data= subset(manhattanData, get(SNP) %in% highlight.snps), aes(x = x, y = -log10(get(P)), label=get(highlight.col)))
  }
  
  # Add title to plot (uses system time if none in specified)
  manhattanPlot <-manhattanPlot + ggtitle(title)
  
  return(manhattanPlot)
} 

#manh(res,SNP = 'rsid', CHR = 'chromosome',BP = 'position',P = 'frequentist_add_pvalue', highlight = T, highlight.col = 'position', highlight.snps = hl)
#manh(res,SNP = 'rsid', CHR = 'chromosome',BP = 'position',P = 'frequentist_add_pvalue')