
library(ggplot2)

# Load data

associationData <- read.table("C:\\Github\\post-association\\tmp\\zBMI11_nial2-autosome.result.gz", header = T, stringsAsFactors = F, sep = " ")

typedVector <- scan("C:\\Github\\post-association\\resources\\genotype\\typed_vector.gz")
associationData$typed <- typedVector
associationData$typed[associationData$typed == 0] <- "Imputed"
associationData$typed[associationData$typed == 1] <- "Genotyped"
associationData$typed <- factor(associationData$typed, levels = c("Genotyped", "Imputed"))

associationData <- associationData[associationData$all_maf > 0.005, ]


# Load category data

birthTable <- read.table("C:\\Projects\\ERC\\correlation\\zBMI0_nial2-egg-bw-comparison-table", header = T, stringsAsFactors = F, sep = " ")
childhoodTable <- read.table("C:\\Projects\\ERC\\correlation\\zBMI0_nial2-egg-comparison-table", header = T, stringsAsFactors = F, sep = " ")
adulthoodTable <- read.table("C:\\Projects\\ERC\\correlation\\zBMI0_nial2-GIANT-comparison-table", header = T, stringsAsFactors = F, sep = " ")

annotationTable <- associationData[, c("rsid", "chromosome", "position")]
annotationTable <- annotationTable[annotationTable$rsid %in% birthTable$rsid
                                   | annotationTable$rsid %in% childhoodTable$rsid
                                   | annotationTable$rsid %in% adulthoodTable$rsid, ]
annotationTable$birth <- ifelse(annotationTable$rsid %in% birthTable$rsid, 0, 1)
annotationTable$childhood <- ifelse(annotationTable$rsid %in% childhoodTable$rsid, 0, 1)
annotationTable$adulthood <- ifelse(annotationTable$rsid %in% adulthoodTable$rsid, 0, 1)

annotationTableBirth <- birthTable[, c("rsid", "chromosome", "position")]
annotationTableBirth$label <- "Birth"

annotationTableAdulthood <- adulthoodTable[, c("rsid", "chromosome", "position")]
annotationTableAdulthood$label <- "Adulthood"


# Best hits to annotate

bestHits <- data.frame(chromosome = c(1, 7, 16), position = c(65991203, 127860163, 53830055), name = c("LEPR", "LEP", "FTO"), stringsAsFactors = F)


# Trim to smaller data frame

# annotatedIndexes <- which(associationData$rsid %in% birthTable$rsid | associationData$rsid %in% childhoodTable$rsid | associationData$rsid %in% adulthoodTable$rsid)
# trimmedIndexes <- unique(c(annotatedIndexes, sample(1:nrow(associationData), 10000)))

# associationDataSmall <- associationData[trimmedIndexes, ]
# associationDataBackUp <- associationData
# associationData <- associationDataSmall


# Make MH

mhPlot <- manhattan(x = associationData, y = annotationTableAdulthood, snp = "rsid", chr = "chromosome", bp = "position", p = "frequentist_add_pvalue", typed = "typed", category = "label", categoryColors = c("red"), categoryFlanking = 20, categoryMinP = 5)

png("C:\\Github\\qqman2\\tmp\\mh_11_adulthood.png", width = 800, height = 600)
plot(mhPlot)
dummy <- dev.off()

# Make MH with MAF

mhPlot <- manhattan(associationData, y = annotationTableAdulthood, snp = "rsid", chr = "chromosome", bp = "position", p = "frequentist_add_pvalue", maf = "all_maf", typed = "typed", category = "label", categoryColors = c("red"))

png("C:\\Github\\qqman2\\tmp\\mh_11_maf_adulthood.png", width = 800, height = 600)
plot(mhPlot)
dummy <- dev.off()

# Make MH with best hits

mhPlot <- manhattan(associationData, y = annotationTableAdulthood, z = bestHits, snp = "rsid", chr = "chromosome", bp = "position", p = "frequentist_add_pvalue", maf = "all_maf", typed = "typed", category = "label", categoryColors = c("red"))
mhGrob <- ggplotGrob(mhPlot)

# Categories settings

categories <- c("birth", "childhood", "adulthood")

colorList <- list(
  c(NA, "darkgreen"),
  c(NA, "darkorange"),
  c(NA, "darkred")
)

# Merge with categories

mergedPlot <- add_reference(mhGrob = mhGrob, y = annotationTable, z = bestHits, snp = "rsid", chr = "chromosome", bp = "position", categories = categories, categoryColors = colorList, flanking = 3)

png("C:\\Github\\qqman2\\tmp\\mh_11_adulthood_annotated.png", width = 800, height = 600)
grid.draw(mergedPlot)
dummy <- dev.off()
