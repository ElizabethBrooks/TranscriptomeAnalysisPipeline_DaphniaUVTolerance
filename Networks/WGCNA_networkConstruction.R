#Retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# allow multi-threading within WGCNA
# caution: skip this line if you run RStudio or other third-party R environments
# see note above
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets){
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, 
                                                     powerVector=powers,
                                                     verbose = 2)[[2]])
}
collectGarbage()

# Plot the results:
colors = c("black", "red", "blue", "green")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", 
             "Mean connectivity", 
             "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
#sizeGrWindow(8, 6)
#par(mfcol = c(4,1))
#par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets){
  if (set==1){
    plot(powerTables[[set]]$data[,1], 
         -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",
         ylab=colNames[col],
         type="n", 
         ylim = ylim[, col],
         main = colNames[col])
    addGrid()
  }
  if (col==1){
    text(powerTables[[set]]$data[,1], 
         -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,
         cex=cex1,
         col=colors[set])
  } else
    text(powerTables[[set]]$data[,1], 
         powerTables[[set]]$data[,plotCols[col]],
         labels=powers,
         cex=cex1,
         col=colors[set])
  if (col==1){
    legend("bottomright", 
           legend = setLabels, 
           col = colors, 
           pch = 20)
  } else
    legend("topright", 
           legend = setLabels, 
           col = colors, 
           pch = 20)
}
dev.off()

# set soft thresholding power
softPower = 8
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate adjacencies in each individual data set
for (set in 1:nSets){
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
}

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate TOMs in each individual data set
for (set in 1:nSets){
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
}

# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list()
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets){
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8)
  # Scale the male TOM
  if (set>1){
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
    TOM[set, ,] = TOM[set, ,]^scalePowers[set]
  }
}

# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list()
for (set in 1:nSets){
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
}
# Open a suitably sized graphics window
sizeGrWindow(6,6)
#pdf(file = "Plots/TOMScaling-QQPlot.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()