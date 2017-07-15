#' Plots the SNPs in the TopLevelPathways in Manhattan plot style
#'
#' @param pathway the pathways associated to the snp. One SNP can be associated to multiple (equivalent of chr)
#' @param rsid the SNPs (equivalent of bp)
#' @param p the p-values
#' @param maf the maf values
#' @param thresholdLow the low threshold value (log10)
#' @param thresholdHigh the high threshold value (log10)
#' @param thresholdLowColor the color of the low threshold
#' @param thresholdHighColor the color of the high threshold
#' @param mafColor the color of the low maf values
pathwayManhattanPlot <- function(pathway, rsid, p, maf, thresholdLow = 5, thresholdHigh = -log10(5e-8), 
                                 thresholdLowColor = "blue", thresholdHighColor = "green", mafColor = "black") {
  
  data <- data.frame(pathway, rsid, p = -log10(p))
  data  <- unique(data)
  data <- data[order(data$pathway, data$rsid), ]
  pathwaySet <- unique(data$pathway)
  
  plot <- ggplot(data, aes(pathway, p, colour = pathway)) + geom_point()
  plot <- plot + geom_hline(aes(yintercept = thresholdLow), col = thresholdLowColor)
  plot <- plot + geom_hline(aes(yintercept = thresholdHigh), col = thresholdHighColor)
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot <- plot + geom_jitter(width = 0.3, height = 0.1)
  plot <- plot + theme(legend.position="none")
  
  plot(plot)
}

printTime <- function(startTime, endTime, message = ""){
  diff <- endTime - startTime
  diffMin <- round(diff[3]/60)
  print(paste(Sys.time(), " ", message, " (", diffMin, " min)", sep = ""))
}


## Main script

startTimeAll <- proc.time()

# Load Libraries

# if(!require("ggplot2")){
#     install.packages("ggplot2")
    #  library("ggplot2")
# }

# Parse command line arguments

args <- commandArgs(TRUE)
if (length(args) >= 3) {
  pathwayFile <- args[1]
  associationFile <- args[2]
  outputFile <- args[3]
} else {
  pathwayFile <- paste(getwd(), "\\resources\\pathwayMap\\tlp0.csv", sep = "")
  associationFile <- paste(getwd(), "\\resources\\pathwayMap\\zBMI0-autosome-maf-above-0-005.result", sep = "")
  outputFile <- paste(getwd(), "\\resources\\pathwayMap\\pathMan0.png", sep = "")
}
print("Parsing arguments completed")

imageFile <- file.path(outputFile)

# Load pathway mapping data

print(paste(Sys.time(), " Loading pathway mapping data from ", pathwayFile, sep = ""))
startTime <- proc.time()

pathwayData <- read.csv(pathwayFile, header = TRUE, sep = ",")

printTime(startTime, proc.time(), "Loading pathway mapping data completed")

# Load association data
# Original source at: /media/local-disk2/helgeland/erc-analyses-results/snptest_unique/singlepoint/merged/

print(paste(Sys.time(), " Loading association data from ", associationFile, sep = ""))
startTime <- proc.time()

associationData <- read.table(associationFile, stringsAsFactors = F, header = T)
associationData <- associationData[!is.na(associationData$frequentist_add_pvalue) & associationData$all_maf > 0.05, c("rsid", "frequentist_add_pvalue", "all_maf")]

printTime(startTime, proc.time(), "Loading data completed")

# Merging the data

print(paste(Sys.time(), " Merging data", sep = ""))
startTime <- proc.time()

data <- merge(associationData, pathwayData, by.x = "rsid")

printTime(startTime, proc.time(), "Merging data completed")

# Make Manhattan plot

print(paste(Sys.time(), " Writing Pathway Manhattan plot to ", imageFile, sep = ""))
startTime <- proc.time()

png(filename = imageFile, width = 1920, height = 1080)
pathwayManhattanPlot(pathway = data$TopLevelPathwayName, rsid = data$rsid, p = data$frequentist_add_pvalue, maf = data$all_maf)
dummy <- dev.off()

printTime(startTime, proc.time(), "Plotting completed")


# Make random manhattan plot
if(length(args) == 4){
    newName <- paste("pathManRandom", sep = "")
    imageFile <- gsub("pathMan", newName, imageFile)
    for(run in 1:args[4]){
      newName <- gsub(".png", paste("_", run, ".png", sep = ""), imageFile)
      newName <- file.path(newName)
      print(paste(Sys.time(), " Writing Random plot ", run, " to ", newName, sep = ""))
      startTime <- proc.time()
      png(filename = imageFile, width = 1920, height = 1080)
      randomData <- data[, c("rsid", "frequentist_add_pvalue", "TopLevelPathwayName")]
      randomData$frequentist_add_pvalue <- -log10(runif(nrow(data)))
      pathwayManhattanPlot(pathway = randomData$TopLevelPathwayName, rsid = randomData$rsid, p = randomData$frequentist_add_pvalue, maf = data$all_maf)
      dummy <- dev.off()
    }  
}

# End of process

printTime(startTimeAll, proc.time(), "Process completed")
