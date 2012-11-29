library(ggplot2)
library(MASS)

infile <- gsub("-","", commandArgs()[8])
intable <- read.table(infile,  header = T, sep = '\t')
mainFrame <- as.data.frame(intable)
mainFrame$rtb <- as.numeric(mainFrame$rtb)
mainFrame$HBA <- as.numeric(mainFrame$HBA)
mainFrame$HBD <- as.numeric(mainFrame$HBD)
mainFrame$acd_most_bpka <- NULL
mainFrame$acd_most_apka <- NULL

completeFrame <- mainFrame[complete.cases(mainFrame),]

bmolweight <- quantile(completeFrame$molweight, 0.01, na.rm = TRUE)
tmolweight <- quantile(completeFrame$molweight, 0.99, na.rm = TRUE)

blogP <- quantile(completeFrame$logP, 0.01, na.rm = TRUE)
tlogP <- quantile(completeFrame$logP, 0.99, na.rm = TRUE)

bPSA <- quantile(completeFrame$PSA, 0.01, na.rm = TRUE)
tPSA <- quantile(completeFrame$PSA, 0.99, na.rm = TRUE)
                                    
bHBA <- quantile(completeFrame$HBA, 0.01, na.rm = TRUE)
tHBA <- quantile(completeFrame$HBA, 0.99, na.rm = TRUE) 
                         
bHBD <- quantile(completeFrame$HBD, 0.01, na.rm = TRUE)
tHBD <- quantile(completeFrame$HBD, 0.99, na.rm = TRUE)   

brtb <- quantile(completeFrame$rtb, 0.01, na.rm = TRUE)   
trtb <- quantile(completeFrame$rtb, 0.99, na.rm = TRUE)   
                   
cleanFrame <- completeFrame[completeFrame$molweight >= bmolweight & completeFrame$molweight <= tmolweight & completeFrame$logP >= blogP & completeFrame$logP <= tlogP & completeFrame$PSA >= bPSA & completeFrame$PSA <= tPSA & completeFrame$HBA >= bHBA & completeFrame$HBA <= tHBA & completeFrame$HBD >= bHBD & completeFrame$HBD <= tHBD & completeFrame$rtb >= brtb & completeFrame$rtb <= trtb,] 

pca <- prcomp(~ logP + molweight+ HBD + HBA + rtb + PSA, scale. =T, retx = T, data = cleanFrame )

comp1 <- as.data.frame(pca$x[,1])
comp2 <- as.data.frame(pca$x[,2])               

cleanFrame$component1 <- comp1[,1]
cleanFrame$component2 <- comp2[,1]


ggplot()+geom_boxplot(aes(cleanFrame$domain, cleanFrame$component1, fill = cleanFrame$domain))
ggsave('visual/pcaBoxTop6.pdf', useDingbats = F)

dens <- kde2d(cleanFrame$component1, cleanFrame$component2, n=50)

ggplot()+
facet_wrap(~factor(cleanFrame$domain),ncol = 2)+
geom_density2d(aes(cleanFrame$component1, cleanFrame$component2),size= .75, breaks = c(max(dens$z),max(dens$z)/2,max(dens$z)/4,max(dens$z)/10))

ggsave('visual/pcaTop6.pdf', useDingbats = F)

dens <- kde2d(cleanFrame$molweight, cleanFrame$logP, n=10)

ggplot()+
facet_wrap(~factor(cleanFrame$domain),ncol = 2)+
geom_point(aes(cleanFrame$molweight, cleanFrame$logP), col = 'darkgrey')+
geom_density2d(aes(cleanFrame$molweight, cleanFrame$logP, col = ..level..),size= .75, breaks = c(max(dens$z),max(dens$z)/2,max(dens$z)/4,max(dens$z)/10))+
coord_cartesian(ylim=c(-2,8), xlim = c(150,1020))

ggsave('visual/molweightLogpTop6.pdf', width=4.4, height = 5,  useDingbats = F)

# Plot rtb distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$rtb))+
geom_hline(yintercept=brtb, col = 'red')+
geom_hline(yintercept=trtb, col = 'red')
ggsave('visual/rtbTop6.pdf',width = 2, height = 6,  useDingbats = F)


# Plot hbd distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$HBD))+
geom_hline(yintercept=bHBD, col = 'red')+
geom_hline(yintercept=tHBD, col = 'red') 
ggsave('visual/hbdTop6.pdf',width = 2, height = 6,  useDingbats = F)


# Plot PSA distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$PSA))+
geom_hline(yintercept=bPSA, col = 'red')+
geom_hline(yintercept=tPSA, col = 'red')  
ggsave('visual/psaTop6.pdf', width= 2, height = 6, useDingbats = F)


# Plot molweight distributions
ggplot()+
geom_boxplot(aes('all',mainFrame$molweight))+
geom_hline(yintercept=bmolweight, col = 'red')+
geom_hline(yintercept=tmolweight, col = 'red') 

ggsave('visual/molweightTop6.pdf',width = 2, height = 6, useDingbats = F)



# Plot logP distributions
ggplot()+
geom_boxplot(aes('all',mainFrame$logP))+
geom_hline(yintercept=blogP, col = 'red')+
geom_hline(yintercept=tlogP, col = 'red')  

ggsave('visual/logPTop6.pdf',width = 2, height = 6, useDingbats = F)


# Plot hba distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$HBA))+
geom_hline(yintercept=bHBA, col = 'red')+
geom_hline(yintercept=tHBA, col = 'red') 
ggsave('visual/hbaTop6.pdf',width = 2, height = 6,  useDingbats = F)  

         
