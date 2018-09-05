#lowess regression code starts on line 94

library(Hmisc)
library(ggplot2)
library(grid)
library(gridExtra)
##########################


#read in file
res <- read.csv('expression_length_data_for_sing_dup.txt', header = T, sep=',')


#log transform the expression and length data
res$log_expr <- log10(res$median_expression + 1)
res$log_len <- log10(res$length)
qplot(res$log_expr)
qplot(res$log_len)
##########################

#correlation between omega and expression/length
rcorr(res$omega, res$log_expr, type= "spearman")
rcorr(res$omega, res$log_len, type= "spearman")

#create singleton and duplicate subsets of the data
sing<- subset(res, type=='Singleton')
dup<- subset(res, type=='Duplicate')

#mwu tests between dups and singletons for length and expression
wilcox.test(sing$log_len, dup$log_len)
wilcox.test(sing$log_expr, dup$log_expr)
##########################


#BOXPLOTS

#to create the expression boxplot
pbp2 <- ggplot(res, aes(type, log_expr)) + 
  geom_boxplot(aes(fill = type)) +  
  scale_x_discrete(labels=c("Duplicable","Singleton")) +
  annotate("text", x = -Inf, y = Inf, 
           label = "{p}-value==0.237",parse = TRUE, hjust=-0.1, vjust=2, size=6)+
  scale_y_continuous(limits = c(0,2.5)) +
  ggtitle("a") +
  scale_fill_manual(values=c("#E69F00","#56B4E9")) +
  theme(text = element_text(size=20), legend.position="none") +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "Expression Level", x = " ")
pbp2

#to create the length boxplot
pbp3 <- ggplot(res, aes(type, log_len)) + 
  geom_boxplot(aes(fill = type)) + 
  scale_x_discrete(labels=c("Duplicable","Singleton")) +
  
  annotate("text", x = -Inf, y = Inf, 
           label = "{p}-value==0.001",parse = TRUE, hjust=-0.1, vjust=2, size=6)+
  ggtitle("c") +
  scale_y_continuous(limits = c(2.25,4.5)) +
  scale_fill_manual(values=c("#E69F00","#56B4E9")) +
  theme(text = element_text(size=20), legend.position="none") +
  guides(fill=guide_legend(title=NULL)) +
  labs(y = "CDS Length", x = " ")
pbp3
##########################

#SCATTERPLOTS

#omega vs expression level scatter plot with lowess regression line
r1 <- ggplot(data = res, aes(res$log_expr, res$omega,colour=type)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c('black',"#E69F00","#56B4E9")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none",  axis.title.y=element_text(face="italic")) +
  ggtitle("b") +
  geom_point(data = dup, aes(x = log_expr, y = omega), size=2) +
  labs(x = "Expression level",y=expression("d"["N"]*"d"["S"])) +
  geom_line(data = sing, aes(colour = 'black', x = lowess(omega ~ log_expr, f=0.2)$x, y = lowess(omega ~ log_expr, f=0.2)$y), size=1.5)+
  geom_line(data = dup, aes(x = lowess(omega ~ log_expr, f=0.4)$x, y = lowess(omega ~ log_expr, f=0.4)$y), size=1.5)
r1

#omega vs cds length scatter plot with lowess regression line
r2 <- ggplot(data = res, aes(res$log_len, res$omega,colour=type)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c('black',"#E69F00","#56B4E9")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none") +
  ggtitle("d") +
  geom_point(data = dup, aes(x = log_len, y = omega), size=2) +
  labs(x = "CDS Length",y=expression("d"["N"]*"d"["S"])) +
  geom_line(data = sing, aes(colour = 'black', x = lowess(omega ~ log_len, f=0.2)$x, y = lowess(omega ~ log_len, f=0.2)$y), size=1.5)+
  geom_line(data = dup, aes(x = lowess(omega ~ log_len, f=0.4)$x, y = lowess(omega ~ log_len, f=0.4)$y), size=1.5)
r2
##########################

#the regression calculations
lfit_expr <- lowess(res$log_expr,res$omega,f=0.3)
lfun_expr <- approxfun(lfit_expr)
fitted_expr <- lfun_expr(res$log_expr)
res$resid_expr <- res$omega-fitted_expr
lfit_len <- lowess(res$log_len,res$omega,f=0.3)
lfun_len <- approxfun(lfit_len)
fitted_len <- lfun_len(res$log_len)
res$resid_len <- res$omega-fitted_len
##########################

#examining the spread of the residuals
r1b <- ggplot(data = res, aes(res$log_expr, resid_expr)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none") +
  labs(x = "Expression level",y="LOWESS Residuals") +
  geom_hline(yintercept = 0, col="black", linetype="dotted") 
r1b
r2b <- ggplot(data = res, aes(res$log_len, resid_len)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none") +
  labs(x = "CDS Length",y="LOWESS Residuals") +
  geom_hline(yintercept = 0, col="black", linetype="dotted") 
r2b
##########################

#mwu testing between the residuals of singleton and duplicate sets
wilcox.test(res$resid_expr[which((res$type=='Singleton'))], res$resid_expr[which((res$type=='Duplicate'))])
wilcox.test(res$resid_len[which((res$type=='Singleton'))], res$resid_len[which((res$type=='Duplicate'))])
##########################

#median and 95 confidence interval calculation
median_len <- median(res$resid_len)
median_expr <- median(res$resid_expr)
median_len_S <- median(res$resid_len[which(res$type=='Singleton')])
median_len_D <- median(res$resid_len[which(res$type=='Duplicate')])
median_expr_S <- median(res$resid_expr[which(res$type=='Singleton')])
median_expr_D <- median(res$resid_expr[which(res$type=='Duplicate')])
q_len <- quantile(res$resid_len, c(.05, .95))
q_expr <- quantile(res$resid_expr, c(.05, .95))
#add these values to the sing and dup subset dataframes
sing<- subset(res, type=='Singleton')
dup<- subset(res, type=='Duplicate')
##########################

#residual vs residual plot
rr <- ggplot(data = res, aes(resid_len, resid_expr,colour=type)) +
  geom_point(aes(colour=type), shape=1) +
  scale_color_manual(values=c("#E69F00","#56B4E9")) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position="none") +
  labs(x = "Residuals on CDS Length",y="Residuals on Expression Level") +
  ggtitle("e") +
  geom_vline(xintercept = median_len, col="gray30", linetype="dotted") +
  geom_vline(xintercept = q_len[1], col="gray30", linetype="dotted") +
  geom_vline(xintercept = q_len[2], col="gray30", linetype="dotted") +
  geom_hline(yintercept = q_expr[1], col="gray30", linetype="dotted") +
  geom_hline(yintercept = q_expr[2], col="gray30", linetype="dotted") +
  geom_hline(yintercept = median_expr, col="gray30", linetype="dotted") +
  geom_point(data = dup, aes(x =resid_len, y =resid_expr), size=2) +
  annotate(geom = "point", x = median_len_S, y= median_expr_S, size = 4, colour= "black") +
  annotate(geom = "point", x = median_len_D, y= median_expr_D, size = 4, colour= "white") 
rr
##########################

#creating the final figure
grid.arrange(pbp2, pbp3, r1, r2, rr, ncol = 4, nrow = 2,
             layout_matrix = cbind(c(1,2),c(3,4),c(5,5),c(5,5)))
##########################

