################################################################################
################################################################################
###
###  Function:  stackBarPlot.R
###  --------------------
###  analyse a vector of observations and return barPlot
###  
###  momo.sander@ebi.ac.uk
###  
################################################################################
################################################################################ 

vector <- gsub("-","", commandArgs()[8])
vector <- strsplit(vector, ',')
vector <- lapply(vector, function(z) as.numeric(z))
vector <- unlist(vector)
mm <- matrix(vector, ncol = 3)
pdf(file = 'visual/stackBarplot.pdf')
barplot(mm)
dev.off()
