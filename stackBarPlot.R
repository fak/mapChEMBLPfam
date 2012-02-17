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

outfile <- gsub("-","", commandArgs()[9])
ncol <- gsub("-","", commandArgs()[10])
ncol <- as.numeric(ncol)

vector <- gsub("-","", commandArgs()[8])
print(vector)
vector <- strsplit(vector, ',')
vector <- lapply(vector, function(z) as.numeric(z))
vector <- unlist(vector)
mm <- matrix(vector, ncol = ncol, byrow = T)
pdf(file = sprintf('visual/%s',outfile))
barplot(mm)
dev.off()
