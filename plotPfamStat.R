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


## @knitr prep
intable <- read.table(infile, sep = '\t', header = T)
mainFrame <- as.data.frame(intable)
mainFrame$source <- as.factor(mainFrame$source)



## @knitr key_figures
  sum_fig <- ddply(mainFrame, 'source', nrow) 
  sum_fig$V2 <- ddply(mainFrame[mainFrame$nDomains == 0,], 'source', nrow)$V1
  sum_fig$V3 <- sum_fig$V2/sum_fig$V1

  quants <- ddply(mainFrame, 'source', summarize, q_5= median(pPfam), q_25 = quantile(pPfam, 0.25, na.rm = TRUE), q_75 = quantile(pPfam, 0.75, na.rm = TRUE))
  quants$iqr <- quants$q_75 - quants$q_25  

## @knitr plot_histogram
  ggplot(mainFrame[complete.cases(mainFrame),])+
  geom_bar(aes(x=nDomains, y = ..density.., fill = source), col = 'black', binwidth = 1, position = 'stack')+
  scale_x_continuous(limits = c(-1, 12), breaks = seq(-.5,11.5, by = 2), labels = seq(0,12, by = 2))+
  facet_wrap( ~ source, ncol = 3)
  

## @knitr plot_boxplot
  ggplot(mainFrame[complete.cases(mainFrame),])+
  geom_boxplot(aes(x=source, y = pPfam, fill = source))+
  scale_y_continuous(limits = c(0,1))

## @knitr plot_dens 
  ggplot(mainFrame)+
  geom_density(aes(y=..density.., x = pPfam, col = source))+
  scale_x_continuous(limits = c(0,1))

## @knitr remains
  t.test(mainFrame$pPfam[mainFrame$source == 1], mainFrame$pPfam[mainFrame$source == 0])
  
  
## @knitr plot_ecdf
  counts <- ecdf(mainFrame$pPfam)(unique(mainFrame$pPfam))*length(mainFrame$pPfam)
  values <- unique(mainFrame$pPfam)
  plotFrame <- data.frame(c = counts, v = values)
  ggplot(plotFrame)+
  geom_hline(yintercept = 0.5, linetype = 'longdash')+
  geom_step(aes(c,v), direction='vh', col = 'blue', size = 1)+ 
  theme_bw()

## @knitr plot_ecdf_pred
  singleFrame <- mainFrame[mainFrame$mapType == 'single',]
  counts <- ecdf(singleFrame$pPfam)(unique(singleFrame$pPfam))
  values <- unique(singleFrame$pPfam)
  singleFrame <- data.frame(c = counts, v = values, t= 'single')

  multiFrame <- mainFrame[mainFrame$mapType == 'multi',]
  counts <- ecdf(multiFrame$pPfam)(unique(multiFrame$pPfam))
  values <- unique(multiFrame$pPfam)
  multiFrame <- data.frame(c = counts, v = values, t= 'multi')
  
  plotFrame = rbind(multiFrame, singleFrame)

  ggplot(plotFrame)+
  geom_hline(yintercept = 0.5, linetype = 'longdash')+
  geom_step(aes(c,v, col = t), direction='vh',  size = 1)+
  theme_bw()+
  theme(legend.position = "none")



## @knitr sum_up
  sum_up <- function(data){
	tt <- ddply(data, .(target_class_L2,mapType),summarize, tt = mean(pPfam), cc= length(pPfam))
        return(tt)
  }





