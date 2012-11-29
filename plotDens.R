library(ggplot2)
library(MASS)

infile <- gsub("-","", commandArgs()[8])
intable <- read.table(infile, header = T, sep = '\t')
mainFrame <- as.data.frame(intable)

# Plot molweight distributions
ggplot()+
geom_boxplot(aes('all',mainFrame$molweight))
ggsave('visual/molweightTop6.pdf',width = 2, height = 6, useDingbats = F)

# Plot logP distributions
ggplot()+
geom_boxplot(aes('all',mainFrame$logP))
ggsave('visual/logPTop.pdf',width = 2, height = 6, useDingbats = F)


# Plot PSA distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$PSA))
ggsave('visual/psaTop6.pdf', width= 2, height = 6, useDingbats = F)

# Plot hba distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$HBA))
ggsave('visual/hbaTop6.pdf',width = 2, height = 6,  useDingbats = F)


# Plot hbd distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$HBD))
ggsave('visual/hbdTop6.pdf',width = 2, height = 6,  useDingbats = F)


# Plot rtb distributions
ggplot()+
geom_boxplot(aes('all', mainFrame$rtb))
ggsave('visual/rtbTop6.pdf',width = 2, height = 6,  useDingbats = F)


stripFrame <- mainFrame[complete.cases(mainFrame$molweight),]
stripFrame <- stripFrame[complete.cases(stripFrame$logP),]


dens <- kde2d(stripFrame$molweight, stripFrame$logP, n=50)

ggplot()+
facet_wrap(~factor(stripFrame$domain),ncol = 2)+
geom_density2d(aes(stripFrame$molweight, stripFrame$logP),size= .75, breaks = c(max(dens$z),max(dens$z)/2,max(dens$z)/4,max(dens$z)/10))+
scale_x_continuous(limits = c(-10,1000))+
scale_y_continuous(limits = c(-5, 10))

ggsave('visual/molweightLogpTop6.pdf', useDingbats = F)


