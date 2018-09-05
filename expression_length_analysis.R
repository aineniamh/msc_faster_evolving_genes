res <- read.csv('expression_length_data_for_sing_dup_muscle_pairwise.txt', header = T, sep=',')
res <- read.csv('genomic_length.txt', header = T, sep=',')


res$log_expr <- log10(res$median_expression + 1)
res$log_len <- log10(res$length)
res$log_genomic_len <- log10(res$genomic_len)
hist(res$log_genomic_len)
sing <- res[which(res$type=='Singleton'), ]
dup <- res[which(res$type=='Duplicate'), ]
hist(res$log_expr)
res$logdn <-log10(res$dN)
res$logds <-log10(res$dS)
res$logo <-log10(res$omega)
hist(res$logdn)
hist(res$logds)
hist(res$logo)

rcorr(res$log_len, res$dN, type= 'pearson')
rcorr(res$log_len, res$omega, type= 'pearson')
rcorr(res$log_genomic_len, res$omega, type= 'pearson')

pe <- ggplot(res, aes(x=median_expression)) + 
  geom_histogram(colour="black", fill="purple") 
pe
pl <- ggplot(res, aes(x=length)) + 
  geom_histogram(colour="black", fill="purple") 
pl

pl <- ggplot(res, aes(x=log_genomic_len)) + 
  geom_histogram(colour="black", fill="purple") 
pl

p1a <- ggplot(res, aes(type, median_expression))
p1 <- p1a + geom_boxplot(aes(fill = type)) +   
  scale_y_continuous(limits = c(0,25)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "Expression")
p1

p2a <- ggplot(res, aes(type, log_genomic_len))
p2 <- p2a + geom_boxplot(aes(fill = type)) +   
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "CDS Length")
p2

grid.arrange(pe,pl,ncol=2)

p3 <- ggplot(res, aes(x=median_expression, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20)) +
  scale_x_continuous(limits = c(0,75)) +
  labs(x = "Expression Level",y="Density")
p3

p4 <- ggplot(res, aes(x=length, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20)) +
  labs(x = "CDS Length",y="Density")
p4


l1 <- loess(formula = res$omega ~ res$median_expression)
summary(l1)
h <- predict(l1)
plot(res$omega ~ res$median_expression)
lines(res$median_expression[order(res$median_expression)], h[order(h)], col="red")
(r_sq_loess <- cor(res$omega, h)^2)
plot(l1$fitted, l1$residuals, main="LOESS")

l2 <- loess(formula = res$omega ~ res$length)
summary(l2)
h <- predict(l2)
plot(res$omega ~ res$length)
lines(res$length[order(res$length)], h[order(h)], col="red")
(r_sq_loess <- cor(res$omega, h)^2)
plot(l2$fitted, l2$residuals, main="LOESS")

pr <- ggplot(aes(x=l1$residuals, y=l2$residuals))
pr
fit1 <- loess(median_expression ~ omega, data = res)
fit2 <- loess(length ~ omega, data = res)

# Look at the TWO residual plots to assess linearity and constant spread
qplot(fit2$residuals,fit1$residuals)  #Residuals Vs. fitted values

pr <- ggplot(data=res, aes(x = loess(median_expression ~ omega)$residuals, y= loess(length ~ omega)$residuals))
pr
p5 <- ggplot(res, aes(x=median_expression, y=omega)) +
  geom_point(aes(colour=type), shape=1) +
  theme(text = element_text(size=20), legend.title=element_blank()) +
  scale_x_continuous(limits = c(0,75)) +
  labs(x = "Expression Level",y="dnds") +
  scale_color_manual(values=c("purple", "maroon")) +
  geom_smooth(method= "loess",
              se=TRUE, # method=lm,  #Add linear regression line
              fullrange=TRUE,
              colour = 'black')
p5
smooth_vals_e = predict(loess(median_expression~omega,res), res$omega)
smooth_vals_l = predict(loess(length~omega,res), res$omega)

qplot(smooth_vals_l,smooth_vals_e)

p6 <- ggplot(res, aes(x=log_len, y=dN, colour=type)) +
  geom_point(aes(colour=type), shape=1) + 
  theme(text = element_text(size=20), legend.title=element_blank()) +
  #scale_x_continuous(limits = c(0,8000)) +
  labs(x = "CDS Length",y="dnds") +
  scale_color_manual(values=c("purple", "maroon")) +
  geom_smooth(method= "loess",
              se=TRUE, # method=lm,  #Add linear regression line
              fullrange=TRUE,
              colour = res$type)

p6

library(hmisc)


grid.arrange(p5,p6,ncol=2)
grid.arrange(p1, p2, p3, p4, p5, p6,ncol=2)

p7 <- ggplot(res, aes(x=length, y=dS, color=type)) +
  scale_color_manual(values=c("purple", "maroon")) +
  geom_point(shape=1) + 
  theme(text = element_text(size=20), legend.title=element_blank()) +
  scale_x_continuous(limits = c(0,8000)) +
  scale_y_continuous(limits = c(0,0.5)) +
  labs(x = "CDS Length",y="dS") +
  geom_smooth(method=lm,
              se=TRUE, # method=lm,  #Add linear regression line
              fullrange=TRUE) 
p7
p8 <- ggplot(res, aes(x=length, y=dN, color=type)) +
  scale_color_manual(values=c("purple", "maroon")) +
  geom_point(shape=1) + 
  theme(text = element_text(size=20), legend.title=element_blank()) +
  scale_x_continuous(limits = c(0,8000)) +
  scale_y_continuous(limits = c(0,0.5)) +
  labs(x = "CDS Length",y="dN") +
  geom_smooth(method=lm,
              se=TRUE, # method=lm,  #Add linear regression line
              fullrange=TRUE) 
p8
grid.arrange(p7,p8,ncol=2)

mean(res$length)
median(res$length)
dup_long <- dup[which(dup$length > 1623), ]
sing_long <-sing[which(sing$length > 1623), ]
res_long <-res[which(res$length > 1623), ]

dup_sh <- dup[which(dup$length < 1623), ]
sing_sh <-sing[which(sing$length < 1623), ]
res_sh <-res[which(res$length < 1623), ]

p7 <- ggplot(res_sh, aes(x=length, y=dS, color=type)) +
  scale_color_manual(values=c("purple", "maroon")) +
  geom_point(shape=1) + 
  theme(text = element_text(size=20), legend.title=element_blank()) +
  scale_x_continuous(limits = c(0,1623)) +
  scale_y_continuous(limits = c(0,0.5)) +
  labs(x = "CDS Length",y="dS") +
  geom_smooth(method=lm,
              se=TRUE, # method=lm,  #Add linear regression line
              fullrange=TRUE) 
p7
p8 <- ggplot(res_sh, aes(x=length, y=dN, color=type)) +
  scale_color_manual(values=c("purple", "maroon")) +
  geom_point(shape=1) + 
  theme(text = element_text(size=20), legend.title=element_blank()) +
  scale_x_continuous(limits = c(0,1623)) +
  scale_y_continuous(limits = c(0,0.5)) +
  labs(x = "CDS Length",y="dN") +
  geom_smooth(method=lm,
              se=TRUE, # method=lm,  #Add linear regression line
              fullrange=TRUE) 
p8
grid.arrange(p7,p8,ncol=2)

p2a <- ggplot(res_long, aes(type, length))
p2 <- p2a + geom_boxplot(aes(fill = type)) +   
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "CDS Length")
p2

wilcox.test(dup_long$length,sing_long$length)
wilcox.test(dup_sh$omega,sing_sh$omega)


fit_s_e <- lm(sing$dS~sing$length)
fit_d_e <- lm(dup$dS~dup$length)
fit_s_e
fit_d_e

fit_s_e <- lm(sing$omega~sing$median_expression)
fit_d_e <- lm(dup$omega~dup$median_expression)

fit_s_l <- lm(sing$omega~sing$length)
fit_s_l
fit_d_l <- lm(dup$omega~dup$length)
fit_d_l

rcorr(sing$omega,sing$median_expression,type="spearman")
rcorr(dup$omega,dup$median_expression,type="spearman")
rcorr(res$omega,res$median_expression,type="spearman")

modele <- lm(res$omega~res$median_expression*type,data=res)
modell <- lm(res$omega~res$length*type,data=res)


model2 <- lm (dup$dN ~ dup$dS)
rcorr(sing$omega,sing$length,type="spearman")
rcorr(dup$omega,dup$length,type="spearman")
wilcox.test(sing$median_expression, dup$median_expression)

ks.test(dup$median_expression, sing$median_expression, alternative = "greater")

input <-res[c("type","omega","median_expression","length")]
result1<- aov(median_expression~omega*type,data=input)
print(summary(result1))
result2<- aov(median_expression~omega+type,data=input)
print(summary(result2))

print(anova(result1,result2))

result3<- aov(length~omega*type,data=input)
print(summary(result3))
result4<- aov(length~omega+type,data=input)
print(summary(result4))

print(anova(result3,result4))

twocol <- input[c("median_expression","length")]
mat <- as.matrix(twocol)
ab <- with(res, ancova.np <- sm.ancova(twocol, res$omega, res$type, model="equal"))


p1a <- ggplot(resdel, aes(type, delta_dnds))
p1 <- p1a + geom_boxplot(aes(fill = type)) + 
  scale_y_continuous(limits = c(-1,1.2)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "∆dn/ds")
p1

lessthan2 <- res[which(res$omega < 2), ]
par(mfrow=c(0,2))
a1 <- with(lessthan2, ancova.np <- sm.ancova(lessthan2$median_expression, lessthan2$omega, lessthan2$type, model="equal"))
a2 <- with(lessthan2, ancova.np <- sm.ancova(lessthan2$length, lessthan2$omega, lessthan2$type, model="equal"))
preda1 <- predict(a1)
par(mfrow=c(1,2))


par(mfrow=c(4, 2))
plot(log_expr~ type, data=res, ylab='log(median expression + 1)', xlab='Gene set')
plot(log(length) ~ type, data=res, xlab='Gene set')
#plot(median_expression ~ omega, data=res)
#plot(length ~ omega, data=res)
#qqnorm(res$length)
#qqnorm(res$median_expression)

plot(res$omega[res$type == 'singleton'], res$log_expr[res$type == 'singleton'], xlab='omega',
     ylab='log(median expression + 1)', pch=19, col='maroon')
points(res$omega[res$type == 'duplicate'], res$log_expr[res$type == 'duplicate'], pch=19, col='dark blue')
legend("topright", c("singleton", "duplicate"), pch=(c(19,19)), col=c('maroon','dark blue'))

plot(res$omega[res$type == 'singleton'], res$length[res$type == 'singleton'], xlab='omega',
     ylab='length', pch=19, col='maroon')
points(res$omega[res$type == 'duplicate'], res$length[res$type == 'duplicate'], pch=19, col='dark blue')
legend("topright", c("singleton", "duplicate"), pch=(c(19,19)), col=c('maroon','dark blue'))

plot(res$dN[res$type == 'singleton'], res$log_expr[res$type == 'singleton'], xlab='dN',
     ylab='log(median expression + 1)', pch=19, col='maroon')
points(res$dN[res$type == 'duplicate'], res$log_expr[res$type == 'duplicate'], pch=19, col='dark blue')
legend("topright", c("singleton", "duplicate"), pch=(c(19,19)), col=c('maroon','dark blue'))

plot(res$dN[res$type == 'singleton'], res$length[res$type == 'singleton'], xlab='dN',
     ylab='length', pch=19, col='maroon')
points(res$dN[res$type == 'duplicate'], res$length[res$type == 'duplicate'], pch=19, col='dark blue')
abline
legend("topright", c("singleton", "duplicate"), pch=(c(19,19)), col=c('maroon','dark blue'))
#interaction.plot(res$median_expression, res$type, res$omega)

plot(res$dS[res$type == 'singleton'], res$log_expr[res$type == 'singleton'], xlab='dS',
     ylab='log(median expression + 1)', pch=19, col='maroon')
points(res$dS[res$type == 'duplicate'], res$log_expr[res$type == 'duplicate'], pch=19, col='dark blue')
legend("topright", c("singleton", "duplicate"), pch=(c(19,19)), col=c('maroon','dark blue'))

plot(res$dS[res$type == 'singleton'], res$length[res$type == 'singleton'], xlab='dS',
     ylab='length', pch=19, col='maroon')
points(res$dS[res$type == 'duplicate'], res$length[res$type == 'duplicate'], pch=19, col='dark blue')
legend("topright", c("singleton", "duplicate"), pch=(c(19,19)), col=c('maroon','dark blue'))
#interaction.plot(res$median_expression, res$type, res$omega)

rcorr(res$dN,res$median_expression,type='spearman');  #monotonic   
rcorr(res$dN,res$median_expression,type='pearson');  #linear  
rcorr(res$dS,res$median_expression,type='spearman');  #monotonic   
rcorr(res$dS,res$median_expression,type='pearson');  #linear   
rcorr(res$omega,res$median_expression,type='spearman');  #monotonic   
rcorr(res$omega,res$median_expression,type='pearson');  #linear 

rcorr(res$dN,res$length,type='spearman');  #monotonic   
rcorr(res$dN,res$length,type='pearson');  #linear 
rcorr(res$dS,res$length,type='spearman');  #monotonic   
rcorr(res$dS,res$length,type='pearson');  #linear   
rcorr(res$omega,res$length,type='spearman');  #monotonic   
rcorr(res$omega,res$length,type='pearson');  #linear 

rcorr(sing$omega,sing$median_expression,type='spearman');  #monotonic 
rcorr(dup$omega,dup$median_expression,type='spearman');  #monotonic 


rcorr(res$median_expression,res$length,type='spearman');  #monotonic 

#hist(res$median_expression)
#hist(res$log_exp)
#hist(res$length)
loglen <- log(res$length)
#hist(loglen)
shapiro.test(res$log_exp)
#qqnorm(res$log_exp)
#qqnorm(loglen)
d_expr <- res$median_expression[which(res$type=='Duplicate')]
s_expr <- res$median_expression[which(res$type=='Singleton')]
wilcox.test(res$length[which(res$type=='Singleton')],res$length[which(res$type=='Duplicate')])
wilcox.test(res$median_expression[which(res$type=='Singleton')],res$median_expression[which(res$type=='Duplicate')])

p1a <- ggplot(res, aes(type, length))
p1 <- p1a + geom_boxplot(outlier.shape=NA, aes(fill = type)) + 
      scale_y_continuous(limits = c(0,4900)) +
      scale_fill_manual(values=c("purple", "maroon")) +
      theme(text = element_text(size=20), axis.title.x = element_blank()) +
      ylab("CDS Length (nt)") +
      guides(fill=guide_legend(title=NULL))

p2a <- ggplot(res, aes(type, median_expression))
p2 <- p2a + geom_boxplot(outlier.shape=NA, aes(fill = type)) + 
      scale_y_continuous(limits = c(0,15)) +
      scale_fill_manual(values=c("purple", "maroon")) +
      theme(text = element_text(size=20), axis.title.x = element_blank()) +
      ylab("Median Expression") +
      guides(fill=guide_legend(title=NULL))

grid.arrange(p1,p2,ncol=2)

covariate1 <- matrix(res$median_expression)
covariate2 <- matrix(res$length)
covariates <-cbind(covariate1, covariate2)
mfrow=c(0, 1)
an <- sm.ancova(res$median_expression, res$omega, res$type)
print(an)
res$logdn <- log(res$dN)
res$logds <- log(res$dS)
xyplot(res$dN ~ res$dS, data=res, groups=type, 
       key = list(text=list(levels(res$type)), space="right",
                  points=list(pch=c(1,1), col=c("blue","maroon"))
                  ),
       pch=c(1,1),
       panel = function(x, y, groups, yscale.components, aspect, type, key, col,...) {
         panel.xyplot(x, y, groups, aspect="iso", type=c("p"), key, col=c("blue","maroon"),...)
         panel.abline(0, 1, col="black")
         yscale.components = yscale.components.log10ticks(lim=c(0.0001,1), logsc = TRUE, ...)
        },
       xlab="Substitutions/Silent Site",
       ylab="Substitutions/Replacement Site",
       main="Substitution values for Macaque and Gibbon Singleton and Duplicate Genes")

a1 <- with(res, ancova.np <- sm.ancova(res$median_expression, res$omega, res$type, model="equal"))
a2 <- with(res, ancova.np <- sm.ancova(res$log_expr, res$omega, res$type, model="equal"))
mfrow=c(0, 1)

p1 <- ggplot(res, aes(x=median_expression, y=omega, color=type)) +
  geom_point(shape=1) + 
  scale_y_continuous(limits = c(0,2)) +
  guides(colour=FALSE) +
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE,
              fullrange=TRUE) 
p1
p2 <- ggplot(res, aes(x=length, y=omega, color=type)) +
  geom_point(shape=1) + 
  scale_y_continuous(limits = c(0,2)) +
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE,
              fullrange=TRUE) 
p2
p3 <- ggplot(res, aes(x=log_expr, y=omega, color=type)) +
  geom_point(shape=1) + 
  scale_y_continuous(limits = c(0,2)) +
  guides(colour=FALSE) +
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE,
              fullrange=TRUE) 
p3
p4 <- ggplot(res, aes(x=log_len, y=omega, color=type)) +
  geom_point(shape=1) + 
  scale_y_continuous(limits = c(0,2)) +
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE,
              fullrange=TRUE) 
p4
grid.arrange(p1,p2,p3,p4,ncol=2)
plot1 <- xyplot(res$omega ~ res$log_expr, data=res, groups=type, 
          key = list(text=list(levels(res$type)), space="right",
                  points=list(pch=c(19,19), col=c("blue","maroon"))
          ),
          pch=c(19,19),
          panel = function(x, y, groups, aspect, type, key, col,...) {
          panel.xyplot(x, y, groups, aspect="iso", type=c("p","r"), key, col=c("blue","maroon"),...)
          }
          )
plot2 <- xyplot(res$omega ~ res$log_len, data=res, groups=type, 
          key = list(text=list(levels(res$type)), space="right",
                  points=list(pch=c(19,19), col=c("blue","maroon"))
          ),
          pch=c(19,19),
          panel = function(x, y, groups, aspect, type, key, col,...) {
          panel.xyplot(x, y, groups, aspect="iso", type=c("p","r"), key, col=c("blue","maroon"),...)
          }
          )
grid.arrange(plot1,plot2, ncol=2)

plot3 <- xyplot(res$dN ~ res$log_expr, data=res, groups=type, 
                key = list(text=list(levels(res$type)), space="right",
                           points=list(pch=c(19,19), col=c("blue","maroon"))
                ),
                pch=c(19,19),
                panel = function(x, y, groups, aspect, type, key, col,...) {
                  panel.xyplot(x, y, groups, aspect="iso", type=c("p","r"), key, col=c("blue","maroon"),...)
                }
)
plot4 <- xyplot(res$dN ~ res$log_len, data=res, groups=type, 
                key = list(text=list(levels(res$type)), space="right",
                           points=list(pch=c(19,19), col=c("blue","maroon"))
                ),
                pch=c(19,19),
                panel = function(x, y, groups, aspect, type, key, col,...) {
                  panel.xyplot(x, y, groups, aspect="iso", type=c("p","r"), key, col=c("blue","maroon"),...)
                }
)
grid.arrange(plot1,plot2,plot3,plot4, ncol=2)

model <- lm(res$omega~res$median_expression*res$type)



