#I'm trying to plot the residuals of one loess regression (omega and expression) against another (omega and cds length). 
# The problem is, I want them to also be a different colour for singleton vs duplicate. 
resall <- read.csv('expression_length_data_for_sing_dup_muscle_pairwise.txt', header = T, sep=',')
res <- read.csv('expression_length_data_for_sing_dup_muscle_pairwise.txt', header = T, sep=',')
res$log_expr <- log10(res$median_expression + 1)
res$log_len <- log10(res$length)


qplot(res$log_expr)
qplot(res$log_len)

res <- res[which(is.finite(res$log_expr)), ]

res$crt_omega <- (res$omega)^(1/3)

qplot(res$crt_omega)
qplot(res$median_expression, res$omega)
qplot(res$log_expr,res$crt_omega)

sing <- res[which(res$type=='Singleton'), ]
dup <- res[which(res$type=='Duplicate'), ]

r1 <- ggplot(data = res, aes(res$log_expr, res$omega)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("blue", "maroon")) +
  theme(text = element_text(size=15), legend.title=element_blank(), legend.position="none") +
  labs(x = "Expression level",y="dN/dS") +
  geom_line(aes(x = lowess(res$omega ~ res$log_expr, f=0.3)$x, y = lowess(res$omega ~ res$log_expr, f=0.3)$y))
r1
lfit <- lowess(res$log_expr,res$omega,f=0.3)
lfun <- approxfun(lfit)
fitted <- lfun(res$log_expr)
resid <- res$omega-fitted
r1b <- ggplot(data = res, aes(res$log_expr, resid)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("blue", "maroon")) +
  theme(text = element_text(size=15), legend.title=element_blank(), legend.position="none") +
  labs(x = "Expression level",y="LOWESS Residuals") +
  geom_hline(yintercept = 0, col="black", linetype="dotted") 
r1b
r2 <- ggplot(data = res, aes(res$log_len, res$omega)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("blue", "maroon")) +
  theme(text = element_text(size=15), legend.title=element_blank(), legend.position="none") +
  labs(x = "CDS Length",y="dN/dS") +
  geom_line(aes(x = lowess(res$omega ~ res$log_len, f=0.3)$x, y = lowess(res$omega ~ res$log_len, f=0.3)$y))
r2
rcorr(res$omega, res$log_expr, type= "spearman")
rcorr(res$omega, res$log_len, type= "spearman")

grid.arrange(r2, r3, r4, ncol= 3)
lfit <- lowess(res$log_len,res$omega,f=0.3)
lfun <- approxfun(lfit)
fitted <- lfun(res$log_len)
resid <- res$omega-fitted
r2b <- ggplot(data = res, aes(res$log_len, resid)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("blue", "maroon")) +
  theme(text = element_text(size=15), legend.title=element_blank(), legend.position="none") +
  labs(x = "Expression level",y="LOWESS Residuals") +
  geom_hline(yintercept = 0, col="black", linetype="dotted") 
r2b
grid.arrange(r1,r1b, r2, r2b)

lfit_expr <- lowess(res$log_expr,res$omega,f=0.3)
lfun_expr <- approxfun(lfit_expr)
fitted_expr <- lfun_expr(res$log_expr)
res$resid_expr <- res$omega-fitted_expr
lfit_len <- lowess(res$log_len,res$omega,f=0.3)
lfun_len <- approxfun(lfit_len)
fitted_len <- lfun(res$log_len)
res$resid_len <- res$omega-fitted_len

wilcox.test(res$resid_expr[which((res$type=='Singleton'))], res$resid_expr[which((res$type=='Duplicate'))])
wilcox.test(res$resid_len[which((res$type=='Singleton'))], res$resid_len[which((res$type=='Duplicate'))])

median_len <- median(res$resid_len)
median_expr <- median(res$resid_expr)
median_len_S <- median(res$resid_len[which(res$type=='Singleton')])
median_len_D <- median(res$resid_len[which(res$type=='Duplicate')])
median_expr_S <- median(res$resid_expr[which(res$type=='Singleton')])
median_expr_D <- median(res$resid_expr[which(res$type=='Duplicate')])
q_len <- quantile(res$resid_len, c(.05, .95))
q_expr <- quantile(res$resid_expr, c(.05, .95))

rr <- ggplot(data = res, aes(resid_len, resid_expr)) +
      geom_point(aes(colour=type), shape=1) +
      scale_color_manual(values=c("blue", "maroon")) +
      theme(text = element_text(size=15), legend.title=element_blank(), legend.position="none") +
      labs(x = "Residuals on CDS Length",y="Residuals on Expression Level") +
      annotate(geom = "point", x = median_len_S, y= median_expr_S, size = 3, colour= "black") +
      annotate(geom = "point", x = median_len_D, y= median_expr_D, size = 3, colour= "yellow") +
      geom_vline(xintercept = median_len, col="gray30", linetype="dotted") +
      geom_vline(xintercept = q_len[1], col="gray30", linetype="dotted") +
      geom_vline(xintercept = q_len[2], col="gray30", linetype="dotted") +
      geom_hline(yintercept = q_expr[1], col="gray30", linetype="dotted") +
      geom_hline(yintercept = q_expr[2], col="gray30", linetype="dotted") +
      geom_hline(yintercept = median_expr, col="gray30", linetype="dotted") 
rr

p_omegas <- wilcox.test(resall$omega[which(res$type=="Singleton")],resall$omega[which(res$type=="Duplicate")])
p_expr <- wilcox.test(res$log_expr[which(res$type=="Singleton")],res$log_expr[which(res$type=="Duplicate")])
p_len <- wilcox.test(res$log_len[which(res$type=="Singleton")],res$log_len[which(res$type=="Duplicate")])


pbp1 <- ggplot(res, aes(type, omega)) + 
  geom_boxplot(aes(fill = type)) + 
  scale_y_continuous(limits = c(0,2.25)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=15), legend.position="none") +
  guides(fill=FALSE) +
  labs(y = "dN/dS", x = " ")
pbp1
pbp2 <- ggplot(res, aes(type, log_expr)) + 
  geom_boxplot(aes(fill = type)) +   
  scale_y_continuous(limits = c(0,2.5)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=15), legend.position="none") +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "Expression Level", x = " ")
pbp2
pbp3 <- ggplot(res, aes(type, log_len)) + 
  geom_boxplot(aes(fill = type)) +  
  scale_y_continuous(limits = c(2.25,4.5)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=15), legend.position="none") +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "CDS Length", x = " ")
pbp3

pd1 <- ggplot(res, aes(x=log_expr, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=15), legend.position="none") +
  labs(x = "Expression level",y="Density")
pd1
pd2 <- ggplot(res, aes(x=log_len, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=15), legend.position="none") +
  labs(x = "CDS Length",y="Density")
pd2


grid.arrange(pbp2, pbp3, r1, r2, rr, ncol = 4, nrow = 2,
             layout_matrix = cbind(c(1,2),c(3,4),c(5,5),c(5,5)))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

dfsing <- data.frame(res$resid_len[which(res$type=="Singleton")],res$resid_expr[which(res$type=="Singleton")])
colnames(dfsing) <- c("resid_len", "resid_expr")
dfdup<- data.frame(res$resid_len[which(res$type=="Duplicate")],res$resid_expr[which(res$type=="Duplicate")])
colnames(dfdup) <- c("resid_len", "resid_expr")

pbp3 <- ggplot(res, aes(type, resid_expr)) + 
  geom_boxplot(aes(fill = type)) +  
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values=c("purple", "maroon")) +
  theme(text = element_text(size=15), legend.position="none") +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "Residuals Expression Level", x = " ")
pbp3


mp1 <- ggplot(res, aes(x=resid_expr, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20), legend.position="none") +
  labs(x = "Residuals Expression Level",y="Density")
mp1
mp2 <- ggplot(res, aes(x=resid_len, fill=type)) + geom_density(alpha=.5) +
  scale_fill_manual( values = c("purple","red")) +
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20), legend.position="none") +
  labs(x = "Residuals CDS Length",y="Density")
mp2

rr

kde.test(x1=as.matrix(dfsing), x2=as.matrix(dfdup))$pvalue

legend = gtable_filter(ggplotGrob(p3), "guide-box") 

grid.arrange(bp, vp, legend, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.2))
print(s) 
print(d)


wilcox.test(l5$residuals, l6$residuals)

s <- c(median(l3$residuals), median(l5$residuals))
d <- c(median(l2$residuals), median(l6$residuals))
sm <- c(mean(l3$residuals), mean(l5$residuals))
dm <- c(mean(l2$residuals), mean(l6$residuals))
