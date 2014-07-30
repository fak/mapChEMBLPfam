################################################################################
################################################################################
###
###  Function:  ecdf.R
###  --------------------
###  plots an ecdf
###  
###  momo.sander@ebi.ac.uk
###  
################################################################################
################################################################################
library(ggplot2)

key <- gsub("-","", commandArgs()[8])
filename <- gsub("-","", commandArgs()[9])
release <- gsub("-","", commandArgs()[10])


infile <- sprintf('data/%s_%s_%s.tab', key , filename, release)
intable <- read.table(infile, sep = '\t', header = T)
mainFrame <- as.data.frame(intable)

## @knitr ecdf
pp <- ggplot(mainFrame)+geom_step(aes(unique(pPfam), ecdf(pPfam)(unique(pPfam))*length(pPfam)))

## @knitr save_ecdf
ggsave(pp, file = sprintf('visual/%s_%s_%s.pdf',key,  release, filename))
