
rm(list=ls())
library(MASS)
setwd('/Users/jiachengli/Desktop/Human Trafifcking/Academia paper/SAE/clean code')

################################################################################
## input data
## The data consists of 
## m: total number of chiefdoms
## yi: CT rate for each chiefdom
## X: design matrix of covariates
## n: sample size for each chiefdom
################################################################################
dat_all = read.csv('dat.csv')
y = dat_all[,2]
X = cbind(1,dat_all[,c(3,4)])
n = dat_all[,5]
m = dim(X)[1]
y_bar = sum(round(dat_all[,2]*n))/sum(n) 
DD.raw = (y_bar * (1-y_bar))/n


################################################################################
## Bayesian method using RStan
## Input parameters
## m: number of small areas
## y: target variables
## X: covariate matrix, include coefficient 1
## DD.raw: sd for each small areas
## p: number of coefficients
bayes_dat = list(m=m,X=X,y=y,D=DD.raw,I=rep(1,m),p=dim(X)[2])
################################################################################
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = 6)

bayes_result = stan(file = "bayes.stan",dat = bayes_dat,
                    pars = c("sig2v","beta","theta"),
                    iter=25000,warmup=5000,
                    control = list(adapt_delta=0.99,max_treedepth = 15),
                    set.seed(30606))
bayes_result
bayes_extract = extract(bayes_result, permuted = TRUE, inc_warmup = FALSE, include = TRUE)
thin = seq(1,80000,by=10)
sig2v_bayes = bayes_extract$sig2v[thin]
beta_bayes = bayes_extract$beta[thin,]
theta_bayes = bayes_extract$theta[thin,]
theta.ci.lower.bayes = apply(theta_bayes, 2, quantile, c(0.025, 0.975))[1,]
theta.ci.upper.bayes = apply(theta_bayes, 2, quantile, c(0.025, 0.975))[2,]
credible.interval.bayes = theta.ci.upper.bayes - theta.ci.lower.bayes
# posterior SD of theta
theta.sd.bayes = apply(theta_bayes, 2, sd)
theta.mean.bayes = colMeans(theta_bayes)
mean(sig2v_bayes)
bayes_result  = matrix(NA,nrow=m,ncol=3)
bayes_result[,1] = round(theta.mean.bayes,3)
bayes_result[,2] = round(theta.sd.bayes,4)
bayes_result[,3] = round(credible.interval.bayes,3)
colnames(bayes_result) = c('theta_mean_bayes','theta_sd_bayes','CI_bayes')
bayes_result
## save the bayesian results
write.csv(bayes_result, 'bayes_result.csv')
bayes_result = read.csv('bayes_result.csv')[,2:4]

################################################################################
## EBLUPs
################################################################################
library(smallarea)
library(sae)
X = as.matrix(X)
# ML----------------------------------------------------------------------------
EBLUP.fun = function(y,X,DD.raw,method){
  sae = mseFH(formula=y~X[,-1], vardir=DD.raw, method = method)
  theta.est = sae$est$eblup
  a.est = sae$est$fit$refvar
  beta.est = sae$est$fit$estcoef[,1]
  mspe = sae$mse

  rmspe = sqrt(mspe)
  z.025 = qnorm(0.975)

  quant = cbind(theta.est - z.025 * rmspe, theta.est + z.025 * rmspe)
  theta.ci.lower = quant[,1]
  theta.ci.upper = quant[,2]
  al = 2 * z.025 * rmspe
  a = a.est
  output = list(
    theta.est = theta.est,
    rmspe = rmspe,
    al = al,
    a = a
  )

}

ml_res = EBLUP.fun(y,X,DD.raw,method='ML')
reml_res = EBLUP.fun(y,X,DD.raw,method='REML')
fh_res = EBLUP.fun(y,X,DD.raw,method='FH')

# PR----------------------------------------------------------------------------
sae.pr = smallareafit(y~DD.raw+X[,-1],method='PR')
theta.est.pr = sae.pr$smallmean.est
a.est.pr = sae.pr$var.comp
beta.est.pr = sae.pr$est.coef
mspe.pr = sae.pr$smallmean.mse

rmspe.pr = sqrt(mspe.pr)
z.025 = qnorm(0.975)

quant.pr = cbind(theta.est.pr - z.025 * rmspe.pr, theta.est.pr + z.025 * rmspe.pr)
theta.ci.lower = quant.pr[,1]
theta.ci.upper = quant.pr[,2]
al = 2 * z.025 * rmspe.pr
a = a.est.pr

table1 = matrix(NA,nrow=m,ncol=9)
colnames(table1) = c('y','Bayes','EBLUP_ML',
                     'sqrt_D','Post SDs','RMSE',
                     'length CrI','length CI',
                     'n')
table1[,1] = y
table1[,2] = round(bayes_result[,1],3)
table1[,3] = round(ml_res$theta.est,3)

table1[,4] = round(sqrt(DD.raw),3)
table1[,5] = round(bayes_result[,2],3)
table1[,6] = round(ml_res$rmspe,3)

table1[,7] = round(bayes_result[,3],3)
table1[,8] = round(ml_res$al,3)

table1[,9] = dat_all$n
table1
library(xtable)
xtable(table1)
write.csv(table1, 'table1.csv')
# 


table2 = matrix(NA,nrow=m,ncol=8)
colnames(table2) = c('x1_age','x2_goods','y',
                     'ML','REML','FH','PR',
                     'n')
table2[,1] = X[,2]
table2[,2] = X[,3]
table2[,3] = y
table2[,4] = round(ml_res$theta.est,3)
table2[,5] = round(reml_res$theta.est,3)
table2[,6] = round(fh_res$theta.est,3)
table2[,7] = round(theta.est.pr,3)
table2[,8] = dat_all$n
table2
library(xtable)
xtable(table2)
write.csv(table2, 'table2_eblup.csv')

