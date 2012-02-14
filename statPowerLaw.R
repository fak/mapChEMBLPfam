################################################################################
################################################################################
###
###  Function:  statPowerlaw.R
###  --------------------
###  analyse a vector of observations and return results
###  
###  momo.sander@ebi.ac.uk
###  
################################################################################
################################################################################ 


### Source the reference power law functions.
source('data/powerlaw/pareto.R')
#source('~/R/powerlaw/zeta.R') at the moment, this function is not available...
#source('~/R/powerlaw/discexp.R') ...and as a consequence these aren't loaded either.
#source('~/R/powerlaw/discweib.R')
#source('~/R/powerlaw/disclnorm.R')
#source('~/R/powerlaw/poisson.R')

### Source the functions to run comparisons against.
source('data/powerlaw/exp.R')
source('data/powerlaw/lnorm.R')
source('data/powerlaw/weibull.R')
source('data/powerlaw/yule.R')


### Load the datasets.
inFile <- gsub("-","", commandArgs()[8])
print(inFile)
inpath = (sprintf("data/%s", inFile))
intable <- read.table(inpath, header = T, sep = '\t')
genFrame <- as.data.frame(intable)
genFreq <- genFrame$freq

### 1. Estimate the values for xmin and scaling coefficient alpha.
source('data/powerlaw/NewmanPowerLawFunctions.R')
newmanModel <- plfit(genFreq)
minx <- newmanModel$xmin
al <- newmanModel$alpha

### 2. Calculate goodness of fit using Rick Wash's script. 
source('data/powerlaw/WashPowerLawFunctions.R')
washModel <- plm(genFreq, xmin = minx)
pp <- summary(washModel)
pVal <- pp$p

x <- c(inFile, 'alpha:',al, 'xmin:',minx,'p-value:', pVal)

write(x, file = sprintf("data/powerLawLog%s" ,inFile), append = TRUE, ncolumns = 7,  sep = '\t')
### 3. Compare the likelihood of power-law distribution vs alternative models.
source('data/powerlaw/power-law-test.R')

mm <- pareto.fit(genFreq, minx)
di <- c("lognormal", "exp", "weibull")


c1 <- lnorm.fit(genFreq)
vo <- vuong(pareto.lnorm.llr(genFreq,mm,c1)) 
p1 <- vo$p.one.sided
p2 <- vo$p.two.sided
vu <- vo$Vuong
x <- c(di[1],vu,p1, p2)
write(x, file = sprintf("data/powerLawLog%s" ,inFile) , append = TRUE,ncolumns = 4, sep = '\t')
 

c2 <- exp.fit(genFreq)
vo <- vuong(pareto.exp.llr(genFreq,mm,c2))
p1 <- vo$p.one.sided
p2 <- vo$p.two.sided
vu <- vo$Vuong
x <- c(di[2],vu,p1, p2)
write(x, file = sprintf("data/powerLawLog%s" ,inFile), append = TRUE,ncolumns = 4, sep = '\t')
 
c3 <- weibull.fit(genFreq)
vo <- vuong(pareto.weibull.llr(genFreq,mm,c3))
p1 <- vo$p.one.sided
p2 <- vo$p.two.sided
vu <- vo$Vuong
x <- c(di[3],vu,p1, p2)
write(x, file = sprintf("data/powerLawLog%s" ,inFile), append = TRUE,ncolumns = 4, sep = '\t')
 
