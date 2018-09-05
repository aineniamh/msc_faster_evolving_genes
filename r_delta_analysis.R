resdel <- read.csv('all_delta_codeml_parsed_results.txt', header = T, sep=',')
resdel$Post.duplication<- resdel$dN_2/resdel$dS_2
colnames(resdel)[which(names(resdel) == "Post.duplication")] <- "Post-duplication"
colnames(resdel)[which(names(resdel) == "Pre.duplication")] <- "Pre-duplication"
resdel$delta_dnds <- resdel$"Post-duplication" - resdel$"Pre-duplication"

resdel <-resdel[resdel$"Pre-duplication" < 10, ]
resdel <-resdel[resdel$"Post-duplication" < 10, ] 
resdel <-resdel[resdel$"dS_2" >= 0.01, ]

omegas <-subset(resdel, select=c("Pre-duplication", "Post-duplication"))

d<- subset(resdel, type=='Duplicate')
d<- na.omit(d)

d$Post.duplication<- is.finite(d$Post.duplication)
d<- subset(d, dS_2 != 0)

write.csv(d, file = "duplicate_delta_analysis.csv",quote=FALSE,row.names = FALSE)

(sd(sing$delta_dnds)*1.96)+mean(sing$delta_dnds)      

h<-subset(dup, delta_dnds>0.575)
       
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Mode(sing$delta_dnds)
Mode(dup$delta_dnds)
less_than_zero <- subset(resdel, delta_dnds <= 0)
less_than_zero <- subset(dup, delta_dnds <= 0)
less_than_zero <- subset(sing, delta_dnds <= 0)

wilcox.test(resdel$delta_dnds[which(resdel$type=='Singleton')], resdel$delta_dnds[which(resdel$type=='Duplicate')])
median(resdel$delta_dnds[which(resdel$type=='Singleton')])
median(resdel$delta_dnds[which(resdel$type=='Duplicate')])

dup <- subset(resdel, type=='Duplicate')
sing <- subset(resdel, type=='Singleton')

ggplot(dup, aes(delta_dnds))+geom_histogram()
ggplot(dup, aes(dS_1))+geom_histogram()
ggplot(dup, aes(dS_2))+geom_histogram()
ggplot(dup, aes(dN_1))+geom_histogram()
ggplot(dup, aes(dN_2))+geom_histogram()

ggplot(dup, aes(dN_1, dN_2))+geom_point()

ggplot(resdel, aes(dS_2, delta_dnds, colour=type))+geom_point()

melted<-melt(omegas, id=c('Pre', 'Post'))

sorted_df <- resdel[order(resdel$"dS_2"), ]
dup_sorted<- sorted_df[which(sorted_df$type == "Duplicate"), ]

1496/2

rcorr(resdel$delta_dnds, resdel$dS_2)
rcorr(lowhalf$residuals, lowhalf$dS_2)


lfit_expr <- lowess(resdel$delta_dnds, resdel$dS_2, f=0.3)
lfun_expr <- approxfun(lfit_expr)
fitted_expr <- lfun_expr(resdel$dS_2)
resdel$resid_expr <- resdel$delta_dnds-fitted_expr

r1b <- ggplot(data = resdel, aes(dS_2, resid_expr)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none") +
  labs(x = "dS",y="LOWESS Residuals") +
  geom_hline(yintercept = 0, col="black", linetype="dotted") 
r1b

median(lowhalf$delta_dnds[which(lowhalf$type=='Singleton')])
median(lowhalf$delta_dnds[which(lowhalf$type=='Duplicate')])

median(resdel$delta_dnds[which(resdel$type=='Singleton')])
median(resdel$delta_dnds[which(resdel$type=='Duplicate')])


wilcox.test(resdel$resid_expr[which((resdel$type=='Singleton'))], resdel$resid_expr[which((resdel$type=='Duplicate'))])
########
#LOW GROUP
rcorr(lowhalf$delta_dnds, lowhalf$dS_2)


lfit_expr <- lowess(lowhalf$delta_dnds, lowhalf$dS_2, f=0.3)
lfun_expr <- approxfun(lfit_expr)
fitted_expr <- lfun_expr(lowhalf$dS_2)
lowhalf$resid_expr <- lowhalf$delta_dnds-fitted_expr

r1b <- ggplot(data = lowhalf, aes(dS_2, residuals)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none") +
  labs(x = "dS",y="LOWESS Residuals") +
  geom_hline(yintercept = 0, col="black", linetype="dotted") 
r1b

wilcox.test(lowhalf$resid_expr[which((lowhalf$type=='Singleton'))], lowhalf$resid_expr[which((lowhalf$type=='Duplicate'))])

sorted_df <- resdel[order(resdel$"dS_2"), ]
lowhalf<- sorted_df[1:742, ]
write.csv(lowhalf, file = "lowhalf.csv",quote=FALSE,row.names = FALSE)

highhalf<- sorted_df[743:1484, ]

dup <- subset(lowhalf, type=='Duplicate')
sing <- subset(lowhalf, type=='Singleton')
wilcox.test(dup$delta_dnds, sing$delta_dnds)

model1 <- lm(resdel$delta_dnds ~ resdel$dS_2)
hist(model1$residuals)
resdel$residuals<- model1$residuals
dup <- subset(resdel, type=='Duplicate')
sing <- subset(resdel, type=='Singleton')

wilcox.test(dup$residuals,sing$residuals)
rcorr(resdel$residuals, resdel$dS_2)

model2 <- lm(lowhalf$delta_dnds ~ lowhalf$dS_2)
hist(model2$residuals)
lowhalf$residuals<- model2$residuals
dup <- subset(lowhalf, type=='Duplicate')
sing <- subset(lowhalf, type=='Singleton')

wilcox.test(dup$residuals,sing$residuals)
rcorr(lowhalf$delta_dnds, lowhalf$dS_2 )

wilcox.test(resdel$delta_dnds, resdel$dS_2)

wilcox.test(lowhalf$delta_dnds, lowhalf$dS_2)



model3 <- lm(formula = highhalf$delta_dnds ~ highhalf$dS_2 )
hist(model3$residuals)
highhalf$residuals <- model3$residuals

dup <- subset(highhalf, type=='Duplicate')
sing <- subset(highhalf, type=='Singleton')

wilcox.test(dup$residuals,sing$residuals)
rcorr(highhalf$residuals, highhalf$dS_2)


model1 <- lm(formula = lowhalf$delta_dnds ~ lowhalf$dS_2)
hist(model1$residuals)
lowhalf$residuals<- model1$residuals
dup <- subset(lowhalf, type=='Duplicate')
sing <- subset(lowhalf, type=='Singleton')

wilcox.test(dup$residuals,sing$residuals)



d<-subset(lowhalf, type=='Duplicate')
s<-subset(lowhalf, type=='Singleton')
wilcox.test(d$delta_dnds, s$delta_dnds)

d<-subset(highhalf, type=='Duplicate')
s<-subset(highhalf, type=='Singleton')
wilcox.test(d$delta_dnds, s$delta_dnds)


ggplot(resdel, aes(delta_dnds, dS_2,  colour=type))+
  geom_point()+
  scale_color_manual(values=c("#E69F00","#56B4E9")) 

ggplot(resdel, aes(type, dS_2,  fill=type))+
  geom_boxplot()+
  annotate("text", x = -Inf, y = Inf, 
           label = "{p}-value==0.04",parse = TRUE, hjust=-0.1, vjust=1, size=10)+
  theme(legend.position = "none", text = element_text(size=35))+
  labs(x = expression(Delta*"d"["N"]*"d"["S"]), y= 'Density')+
  scale_fill_manual(values=c("#E69F00","#56B4E9")) 


lfit_expr <- lowess(resdel$dS_2,resdel$delta_dnds,f=0.3)
lfun_expr <- approxfun(lfit_expr)
fitted_expr <- lfun_expr(resdel$dS_2)
resdel$resid_delta <- resdel$omega-fitted_expr




d_low <- dup_sorted[1:16,]
d_med <- dup_sorted[17:32, ]
d_high <- dup_sorted[33:49, ]

as_low <- sorted_df[1:498, ]
as_med <- sorted_df[499:997, ]
as_high <- sorted_df[998:1496, ]
dl <- subset(as_high, type=='Singleton')

#as_low <- sorted_df[which(sorted_df$"dS_2">=0 & sorted_df$"dS_2"<=0.0699), ]
#as_med <- sorted_df[which(sorted_df$"dS_2">=0.0417 & sorted_df$"dS_2"<=0.0699), ]
#as_high <- sorted_df[which(sorted_df$"dS_2">=0.0701 & sorted_df$"dS_2"<=0.7483), ]

t_high <- sorted_df[which(sorted_df$"dS_2"<=10 & sorted_df$"dS_2">=0.15), ]
t_med <- sorted_df[which(sorted_df$"dS_2">0.05 & sorted_df$"dS_2"<0.15), ]
t_low <- sorted_df[which(sorted_df$"dS_2">=0.00 & sorted_df$"dS_2"<=0.05), ]


s_high <- sorted_df[which(sorted_df$"dS_2"<=0.0527 & sorted_df$"dS_2">=0.0324), ]
s_med <- sorted_df[which(sorted_df$"dS_2">=0.0127 & sorted_df$"dS_2"<=0.06324), ]
s_low <- sorted_df[which(sorted_df$"dS_2">=0.0026 & sorted_df$"dS_2"<=0.0125), ]

a <- ggplot(resdel, aes(x=delta_dnds, fill = type)) + geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(-1.85, 1.85)) + 
  #ggtitle("Low ds
          #ds < 0.05") +
  scale_fill_manual(values=c("#E69F00","#56B4E9")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=25), legend.position="none", plot.title = element_text(size = 10)) +
  #annotate("text", x = Inf, y = Inf, label = "n = 698
           #p = 0.03725", hjust=1, vjust=1, size=3)+
  labs(x = expression(paste(Delta, "d"["N"]*"d"["S"])), y="Density")
a



l <- ggplot(lowhalf, aes(x=type, y = delta_dnds)) + geom_boxplot() +
  #scale_x_continuous(limits = c(-1.85, 1.85)) + 
  ggtitle("Low ds
  ds < 0.05") +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=15), legend.position="none", plot.title = element_text(size = 10)) +
  annotate("text", x = Inf, y = Inf, label = "n = 698
  p = 0.03725", hjust=1, vjust=1, size=3)+
  labs(x = "∆dn/ds",y="Density")
l
m <- ggplot(as_med, aes(x=delta_dnds, fill=type)) + geom_density(alpha=.5) +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  ggtitle("Medium ds
  0.05 < ds < 0.15") +  
  scale_fill_manual(values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) +
  annotate("text", x = Inf, y = Inf, label = "n = 679
  p = 0.6079", hjust=1, vjust=1, size=3)+
  theme(text = element_text(size=15), legend.position='none', plot.title = element_text(size = 10)) +
  labs(x = "∆dn/ds", y=" ")
m
h <- ggplot(as_high, aes(x=delta_dnds, fill=type)) + geom_density(alpha=.5) +
  scale_x_continuous(limits = c(-1, 1)) + 
  ggtitle("High ds
  0.15 < ds") +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  annotate("text", x = Inf, y = Inf, label = "n = 119
  p = 0.7282", hjust=1, vjust=1, size=3)+
  theme(text = element_text(size=15), legend.position="none", plot.title = element_text(size = 10)) +
  labs(x = "∆dn/ds",y=" ")
h

grid.arrange(l,m,h,ncol=3)


wilcox.test(t_low$delta_dnds[which(t_low$type == 'Singleton')],t_low$delta_dnds[which(t_low$type == 'Duplicate')])
wilcox.test(t_med$delta_dnds[which(t_med$type == 'Singleton')],t_med$delta_dnds[which(t_med$type == 'Duplicate')])
wilcox.test(t_high$delta_dnds[which(t_high$type == 'Singleton')],t_high$delta_dnds[which(t_high$type == 'Duplicate')])

sing_delta <- resdel[which(resdel$type=='Singleton'), ]
dup_delta <- resdel[which(resdel$type=='Duplicate'), ]

ps <- ggplot(sing_delta, aes(x=delta_dnds)) + 
  geom_histogram(binwidth = 0.05,colour="black", fill="maroon")
ps
pd <- ggplot(dup_delta, aes(x=delta_dnds)) + 
  geom_histogram(binwidth = 0.2,colour="black", fill="purple") + 
  scale_x_continuous(limits = c(-1.85, 1.85))
pd
grid.arrange(ps,pd,ncol=2)
dat <- data.frame(dens = c(resdel$"Pre-duplication", resdel$"Post-duplication")
                  , lines = rep(c("Pre", "b"), each = 100))

median(resdel$delta_dnds[which(resdel$type=="Duplicate")])


#Plot.
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)

ps <- ggplot(resdel, aes(x=delta_dnds, fill=type)) + geom_density(alpha=.5) +
  scale_x_continuous(limits = c(-1.85, 1.85)) + 
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20)) +
  labs(x = "∆dn/ds",y="Density")
ps

p1a <- ggplot(resdel, aes(type, delta_dnds))
p1 <- p1a + geom_boxplot(aes(fill = type)) + 
  scale_y_continuous(limits = c(-1,1.2)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "∆dn/ds")
p1
r<- resdel[which(resdel$"type"=='Duplicate'), ]
p1 <- qplot(r, aes(x=r$"dN_2", y=r$"Post-duplication")) +
  theme(text = element_text(size=20), legend.title=element_blank()) 


p1



pdn <- ggplot(resdel, aes(x=resdel$"Post-duplication", fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20)) +
  labs(x = "Post-Duplication dNdS",y="Density")
pdn

pds <- ggplot(resdel, aes(x=resdel$dS_2, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20)) +
  labs(x = "Post Dup ds",y="Density")
pds


grid.arrange(p1,ps,ncol=2)

ks.test(d_del,s_del,alternative = "less")

d_del <- dup_delta[ ,'delta_dnds']
s_del <- sing_delta[ ,'delta_dnds']

p1 <- ggplot(resdel, aes(x=omega1, y=omega2)) +
  geom_point(shape=1) 
p1

wilcox.test(d_del, s_del)
boxplot(d_del, s_del)


hist(resdel$deltadnds)

wilcox.test(dup_delta$"Pre-duplication", dup_delta$"Post-duplication")

median(dup_delta$"Pre-duplication")

p1a <- ggplot(resdel, aes(dn., omega2))
p1 <- p1a + geom_boxplot(outlier.shape=NA)+ 
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) 
p1

ggplot(data = melt(omegas), 
       aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) 

####
  
p1a <- ggplot(resdel, aes(type, delta_dnds))
p1 <- p1a + geom_boxplot(aes(fill = type)) + 
  scale_y_continuous(limits = c(-1,1.2)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=20), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "∆dn/ds")
p1

p2a <- ggplot(data = melt(omegas), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_y_continuous(limits = c(0,1.5)) +
  theme(text = element_text(size=20), axis.title.y = element_blank(), axis.title.x = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  scale_x_discrete(labels=c("Pre-duplication","Post-duplication")) + 
  guides(fill=guide_legend(title='dN/dS', values= c("Pre-duplication","Post-duplication")),label=FALSE)
p2a
grid.arrange(p2a,p1,ncol=2)

ks.test(d_del,s_del,alternative = "two.sided")

sig$delta<- sig$omega-sig$background_omega
