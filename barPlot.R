################################################################################
################################################################################
###
###  Function:  barPlot.R
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
pdf(file = 'visual/barplot.pdf')
barplot(vector, horiz = TRUE)
dev.off()
