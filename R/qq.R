#' Creates a QQ plot using ggplot2 and returns the plot object.
#' 
#' @param x A data frame with result data
#' @param p P-value column in data frame
#' @param maf MAF column in data frame (ignored if NA)
#' @param typed the column in the data frame indicating whether the markers are genotyped or imputed (ignored if NA)
#' @param thresholdLow the low threshold value (log10)
#' @param thresholdHigh the high threshold value (log10)
#' @param thresholdLowColor the color of the low threshold
#' @param thresholdHighColor the color of the high threshold
#' @param title Title of plot (date by default, ignored if NA)
#' 
#' @return A qq plot (ggplot2 object)
#' 
#' @import ggplot2
#' 
#' @export

#devtools::use_package("ggplot2", "Suggests")

qq <- function(x, p='P', maf = 'MAF', typed = NA, 
                      thresholdLow = 5, thresholdHigh = -log10(5e-8), 
               thresholdLowColor = "blue", thresholdHighColor = "red", title=Sys.time()){
  
  # Make a data frame for the plot
  
  qqData <- data.frame(p = x[[p]])
  if (!is.na(maf)) {
    qqData$maf <- x[[maf]]
  }
  qqData$logP <- -log10(qqData$p)
  if (!is.na(typed)) {
    qqData$typed <- x[[typed]]
  }
  
  
  # Get the expected p-values
  
  qqData <- qqData[order(qqData$logP), ]
  qqData$expectedP <- rev(-log10(ppoints(n = length(qqData$p))))
  
  
  # Get lambda
  
  chisq <- qchisq(1-qqData$p,1)
  lambda <- median(chisq)/qchisq(0.5,1)
  lambdaLabel <- paste("lambda:", round(lambda, digits = 3), sep = " ")
  
  
  # Get the maximal p-value to plot
  
  yMax <- max(round(max(qqData$logP)+1), round(max(qqData$expectedP)+1), thresholdHigh)
  
  
  # Order by maf to have the common variants in front
  
  if (!is.na(maf)) {
    qqData <- qqData[order(qqData$maf), ]
  }
  
  
  # Create the plot
  
  qqPlot <- ggplot()
  
  
  # Add the points
  
  if (!is.na(maf)) {
    if (length(typed) > 1 || !is.na(typed)) {
      qqPlot <- qqPlot + geom_point(data = qqData, aes(x = expectedP, y = logP, col = maf, size = typed))
      qqPlot <- qqPlot + scale_size_manual(name = "", values = c(2, 1))
    } else {
      qqPlot <- qqPlot + geom_point(data = qqData, aes(x = expectedP, y = logP, col = maf))
    }
  } else {
    if (length(typed) > 1 || !is.na(typed)) {
      qqPlot <- qqPlot + geom_point(data = qqData, aes(x = expectedP, y = logP, size = typed), col = "black")
      qqPlot <- qqPlot + scale_size_manual(name = "", values = c(2, 1))
    } else {
      qqPlot <- qqPlot + geom_point(data = qqData, aes(x = expectedP, y = logP), col = "black")
    }
  }
  
  
  # Add the thresholds
  
  qqPlot <- qqPlot + geom_hline(aes(yintercept = thresholdLow), col = "blue")
  qqPlot <- qqPlot + geom_hline(aes(yintercept = thresholdHigh), col = "green")
  
  
  # Add the annotation
  
  qqPlot <- qqPlot + geom_abline(intercept = 0, slope = 1, col = "grey40", size = 1, linetype = "dotted")
  qqPlot <- qqPlot + annotate("text", x = yMax/3, y = yMax-1.5, label = lambdaLabel)
  
  
  # Format plot
  
  qqPlot <- qqPlot + scale_y_continuous(name = "-log10(p) Observed", breaks = 0:yMax, limits = c(0, yMax), expand = c(0, 0))
  qqPlot <- qqPlot + scale_x_continuous(name = "-log10(p) Expected", breaks = 0:yMax, limits = c(0, yMax), expand = c(0, 0))
  if (!is.na(maf)) {
    qqPlot <- qqPlot + scale_color_gradientn(name = "MAF", colors = c("darkred", "chocolate", "darkgreen"))
  }
  qqPlot <- qqPlot + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                           panel.grid.major = element_line(colour = "grey50", linetype = "dotted"),
                           panel.grid.minor = element_blank())
  
  # Return the plot
  
  return(qqPlot)
}

