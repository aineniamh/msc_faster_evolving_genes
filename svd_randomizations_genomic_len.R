library(Hmisc)


res <- read.csv('genomic_length.txt', header = T)

#log transform the expression and length data
res$log_expr <- log10(res$median_expression + 1)
res$log_len <- log10(res$genomic_len)


#correlation between omega and expression/length
rcorr(res$omega, res$log_expr, type= "spearman")
rcorr(res$omega, res$log_len, type= "spearman")

#create singleton and duplicate subsets of the data
sing<- subset(res, type=='Singleton')
dup<- subset(res, type=='Duplicate')

#mwu tests between dups and singletons for length and expression
wt.svd.len<- wilcox.test(sing$log_len, dup$log_len)
wt.svd.expr <- wilcox.test(sing$log_expr, dup$log_expr)

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

#mwu testing between the residuals of singleton and duplicate sets
wt.svd.expr.res<- wilcox.test(res$resid_expr[which((res$type=='Singleton'))], res$resid_expr[which((res$type=='Duplicate'))])

P.real.svd.expr <- wt.svd.expr.res$p.value
wt.svd.len.res<- wilcox.test(res$resid_len[which((res$type=='Singleton'))], res$resid_len[which((res$type=='Duplicate'))])
P.real.svd.len <- wt.svd.len.res$p.value
##########################

randP.expr <- c()
randP.len <- c()


for (i in c(1:10000)) {



r.expr <- sample(res$log_expr)
r.len <- sample(res$log_len)


lfit_expr <- lowess(r.expr,res$omega,f=0.3)
lfun_expr <- approxfun(lfit_expr)
fitted_expr <- lfun_expr(r.expr)
res$resid_expr.r <- res$omega-fitted_expr


lfit_len <- lowess(r.len,res$omega,f=0.3)
lfun_len <- approxfun(lfit_len)
fitted_len <- lfun_len(r.len)
res$resid_len.r <- res$omega-fitted_len

##########################

#mwu testing between the residuals of singleton and duplicate sets - randomized
wt.svd.expr.res.rand <- wilcox.test(res$resid_expr.r[which((res$type=='Singleton'))], res$resid_expr.r[which((res$type=='Duplicate'))])
P.rand.svd.expr <- wt.svd.expr.res.rand$p.value
randP.expr <- c(randP.expr, P.rand.svd.expr)

wt.svd.len.res.rand <- wilcox.test(res$resid_len.r[which((res$type=='Singleton'))], res$resid_len.r[which((res$type=='Duplicate'))])
P.rand.svd.len <- wt.svd.len.res.rand$p.value
randP.len <- c(randP.len, P.rand.svd.len)


res$resid_expr.r <- c()
res$resid_len.r <- c()
}

hist(randP.len)
abline(v=0.031, col='orange')
median(randP.len)
abline(v=0.025, col='blue')

lesssig.len <- randP.len[randP.len >= P.real.svd.len]
overallP.len <- length(lesssig.len)/length(randP.len)



moresig.expr <- randP.expr[randP.expr <= P.real.svd.expr]
overallP.expr <- length(moresig.expr)/length(randP.expr)
