################################################################################
################################################################################
###
###  Function:  plotPfamStat.R
###  --------------------
###  analyse a vector of observations and return barPlot
###  
###  momo.sander@ebi.ac.uk
###  
################################################################################
################################################################################
library(ggplot2)


release<- gsub("-","", commandArgs()[8])
print(release)
infile <- sprintf('data/pfamTable_%s.tab', release)
intable <- read.table(infile, sep = '\t', header = T)
mainFrame <- as.data.frame(intable)
mainFrame$source <- as.factor(mainFrame$source)


pp <- ggplot(mainFrame)+
  geom_bar(aes(x=nDomains, y = ..density..), binwidth = 1, position = 'stack')+
  scale_x_continuous(limits = c(-1, 12), breaks = seq(-.5,11.5, by = 2))+
  facet_wrap( ~ source, ncol = 3)
  
ggsave(pp, file= sprintf('visual/histogram_%s.pdf', release))

pp <- ggplot(mainFrame)+
  geom_boxplot(aes(x=source, y = pPfam, fill = source))+
  scale_y_continuous(limits = c(0,1))
  
ggsave(pp, file= sprintf('visual/boxplot_%s.pdf', release))

pp <- ggplot(mainFrame)+
  geom_density(aes(y=..density.., x = pPfam, col = source))+
  scale_x_continuous(limits = c(0,1))

  t.test(mainFrame$pPfam[mainFrame$source == 1], mainFrame$pPfam[mainFrame$source == 0])
  
  
  

