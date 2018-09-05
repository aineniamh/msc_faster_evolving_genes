cc <- read.csv('duplicability/all_codeml_results.txt', header = T, sep=',')
sp <- cc[which(cc$type=='singleton'), 1:9]
dp <- cc[which(cc$type=='duplicate'), 1:9]
sing_high_dnds <- sp[which(sp$dnds_ML>2), 1:9]
length(which(cc$type== 'singleton'))
length(which(cc$type== 'duplicate'))
unique(cc$type)
nrow(d)
dim(s)
head(s)
names(d)
str(s)
levels(cc$type)

length(which(s$type == 'duplicate'))
s_dn <- s[which(s$)]

isPositive <- function(x) x>=0 & x<10
isInRange <- function(x) x>0.01 & x<2
islow <- function(x) x<2

s_dnds_NG <- Filter(isPositive, sp[ ,'dnds_NG'])
d_dnds_NG <- Filter(isPositive, dp[ ,'dnds_NG'])
s_dnds_ML <- Filter(isPositive, sp[ ,'dnds_ML'])
d_dnds_ML <- Filter(isPositive, dp[ ,'dnds_ML'])

s_dn_ML <- Filter(islow, sp[ ,'X2ML.dN'])
d_dn_ML <- Filter(islow, dp[ ,'X2ML.dN'])
s_dn_NG <- Filter(islow, sp[ ,'X2NG.dN'])
d_dn_NG <- Filter(islow, dp[ ,'X2NG.dN'])

s_ds_ML <- Filter(isInRange, sp[ ,'X2ML.dS'])
d_ds_ML <- Filter(isInRange, dp[ ,'X2ML.dS'])
s_ds_NG <- Filter(isInRange, sp[ ,'X2NG.dS'])
d_ds_NG <- Filter(isInRange, dp[ ,'X2NG.dS'])

stest <- shapiro.test(s_dnds_NG)
stest$p.value
qqnorm(s_dnds_ML)
#definitely not a normal distribution for any of them
wilcox.test(d_dn_NG,s_dn_NG)
wilcox.test(d_dn_ML,s_dn_ML)
wilcox.test(d_ds_ML,s_ds_ML)
wilcox.test(d_dnds_ML,s_dnds_ML)
wilcox.test(d_dnds_NG,s_dnds_NG)
boxplot(d_dn_ML,s_dn_ML,
        names = c('Duplicates','Singletons'),
        col = c('pink','purple'),ylim = c(0,0.13))


par(mfrow=c(1, 2))
boxplot(d_dn_NG,s_dn_NG,d_ds_NG,s_ds_NG,
        names = c('Duplicates dn','Singletons dn','Duplicates ds','Singletons ds'),
        col = c('orange','purple','orange','purple'),outline = FALSE)
title(main = list("Nei-Gojobori dn and ds for duplicate and singleton gene sets"), sub = "*Outliers not shown")

boxplot(d_dnds_NG,s_dnds_NG,
        names = c('Duplicates','Singletons'),
        col = c('orange','purple'), outline = FALSE)
title(main = list("Nei-Gojobori dn/ds for duplicate and singleton gene sets"), ylab = 'dn/ds', sub = "*Outliers not shown")


par(mfrow=c(1, 2))
boxplot(d_dn_ML,s_dn_ML,d_ds_ML,s_ds_ML,
        names = c('Duplicates dn','Singletons dn','Duplicates ds','Singletons ds'),
        col = c('orange','purple','orange','purple'))
title(main = list("PAML dn and ds for duplicate and singleton gene sets"))
boxplot(d_dnds_ML,s_dnds_ML,
        names = c('Duplicates','Singletons'),
        col = c('orange','purple'))
title(main = list("PAML dn/ds for duplicate and singleton gene sets"), ylab = 'dn/ds')



boxplot(d_dn_NG,d_ds_NG,s_dn_NG,s_ds_NG,
        names = c('Duplicates dn','Duplicates ds','Singletons dn','Singletons ds'),
        col = c('orange','purple'), outline = FALSE)
median(d_dnds_ML)
Median_dnds <- c(median(d_dnds_ML), median(s_dnds_ML), median(d_dnds_NG), median(s_dnds_NG))
Median_dn <- c(median(d_dn_ML), median(s_dn_ML), median(d_dn_NG), median(s_dn_NG))
median_df <- data.frame(Median_dnds, Median_dn)
rownames(median_df) = c('D_ML','S_ML','D_NG','S_NG')