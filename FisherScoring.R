#---------------------------------------------------------------#
# Fisher scoring iteration with Wald, Likelihood ratio test,   #
# & deviance statistics & person Chi-square test                #
#---------------------------------------------------------------#



###################################
#  Fisher scoring iterative method  #
#    with binomial distribution     #
###################################

# start: p_dim from uniform
# make x as study matrix

# create datasets    
x <- matrix(c(rep(1,10), log(c(1/128, 1/64, 1/32, 1/16, 1/8,1/4,1/2,1,2,4))), nrow = 10, ncol = 2) # x is on log scale. 
y <- matrix(c(0,0,0.4,0.4,0.6,0.8,1,1,1,1), ncol = 1, byrow = FALSE)   # y is proportion
m <- 5 # m is equal for each x_i

# fisher scoring function to get MLE of regression coefficients. 
# num of iteration: 20
fisher.scoring <- function(y, x, start = runif(ncol(x))) {
  m <- 5
  p <- ncol(x)
  beta <- matrix(rep(0,20*p),ncol = 20)
  score <- rep(0, p)
  theta <- matrix(rep(0,20*10), ncol =20)
  mu <- matrix(rep(0,20*10), ncol = 20)
  beta[,1] <- start
  
  for (i in 1:19) {
    # cannonical link function
    theta[,i] <- x %*% as.matrix(beta[,i])
    mu[,i] <- exp(as.matrix(theta[,i]))/(1+exp(as.matrix(theta[,i])))
    # w 
    w <- diag(as.vector(mu[,i] * (1-mu[,i])))
    #beta
    beta[,i+1] <- as.matrix(beta[,i]) + solve(t(x)%*% w %*% x) %*% t(x) %*% (y-as.matrix(mu[,i]))
  } 
return (beta)
}
# Call fisher scoring function to get iteration values. Check for convergence.
beta <- fisher.scoring(y,x)

# Get MLE of beta
beta_MLE <- as.matrix(beta[,20])

#################################
#  Wald Test; H0: beta_i = 0;  #
#################################
# to get standard error of MLE
theta.mle <- x %*% as.matrix(beta[,20])
mu.mle <- exp(as.matrix(theta.mle))/(1+exp(as.matrix(theta.mle)))
w.mle <- diag(as.vector(mu.mle * (1-mu.mle)))
se.mle <- sqrt(solve(m*t(x)%*%w.mle%*%x))  # square root of covariance matrix
cov.mle <- se.mle^2      # variance-covariance matrix of estimated coefficent

# Z test statistics for beta0 and beta1
z.0 <- (beta_MLE[1,1] - 0) / se.mle[1,1] # square root of Chisq test
z.1 <- (beta_MLE[2,1] - 0) / se.mle[2,2] # square root of Chisq test

# Chisq test statistics for beta0 and beta1
Chisq.0 <- z.0^2
Chisq.1 <- z.1^2
# CI of beta0 and beta1 for Wald test
CI.0 <- beta_MLE[1,1] +c(1,-1)*1.96* se.mle[1,1]
CI.1 <- beta_MLE[2,1] +c(1,-1)*1.96* se.mle[2,2]

Chisq.0
CI.0
Chisq.1
CI.1


#################################
#  Likelihood Ratio Test;       #
#    H0: beta_i = 0;            #
#################################
# under H1: beta_MLE
pi <- exp(theta.mle)/(1+exp(theta.mle))
loglik_1 <- sum(m*(y*theta.mle - log(1+exp(theta.mle))))

# under H0:beta_i = 0, testing for beta1 
beta1.h0 <- 0 # slope is fixed as 0. 
beta0 <- fisher.scoring(y,as.matrix(x[,1]))  # if beta1 = 0, only intercept in model, hence we can use fisher scoring function well here. 
alpha <- beta0[,20] # assign MLE of intercept as alpha. 
theta0.h0 <- as.matrix(x[,1]) * alpha + 0
mu0.h0 <- exp(as.matrix(theta0.h0))/(1+exp(as.matrix(theta0.h0)))
w <- diag(as.vector(mu0.h0 * (1-mu0.h0)))
loglik_0 <- sum(m*(y*theta0.h0 - log(1+exp(theta0.h0))))
# LR for beta1=0
2*(loglik_1-loglik_0)

# Compute CI for likelihood ratio test
# function to find MLE for intercept, when slope is fixed as j. 
GetAlphaFun <- function(y,x,j, start = 0) {
  p <- 1
  alpha <- matrix(rep(0,20*p),ncol = 20)
  score <- rep(0, p)
  theta <- matrix(rep(0,20*10), ncol =20)
  mu <- matrix(rep(0,20*10), ncol = 20)
  alpha[,1] <- start
for (i in 1:19) {
    # cannonical link function
    theta[,i] <- as.matrix(x[,1] * alpha[,i]) + as.matrix(x[,2] * j) 
    mu[,i] <- exp(as.matrix(theta[,i]))/(1+exp(as.matrix(theta[,i])))
    # w 
    w <- diag(as.vector(mu[,i] * (1-mu[,i])))
    #beta
    alpha[,i+1] <- as.matrix(alpha[,i]) + sum(y-mu[,i])/sum(mu[,i]*(1-mu[,i]))
    }
return(alpha)
  }

ran <- seq(0,3, by = 0.01)
nran <- length(ran)
TS <- rep(0,nran)
for (i in 1:nran) {
  j <- ran[i]
  alpha <- GetAlphaFun(y,x,j)[,20]
  beta.full <- as.matrix(c(alpha,j))
  theta0.h0 <- x %*% beta.full
  mu0.h0 <- exp(as.matrix(theta0.h0))/(1+exp(as.matrix(theta0.h0)))
  w <- diag(as.vector(mu0.h0 * (1-mu0.h0)))
  loglik_0 <- sum(m*(y*theta0.h0-log(1+exp(theta0.h0))))
  
  TS[i] <- 2*(loglik_1-loglik_0)
}
as.matrix(sign(TS - 3.84))

c(ran[92], ran[276])



###########################################
#  Goodness of test for deviance statistic #
###########################################

s <- m*y
s.hat <- m*mu.mle
# Be sure to remove 0 multipliers in equation. 
GOF <- 2*(sum(s[c(3:10),]*log(s[c(3:10),]/s.hat[c(3:10),])) +sum((m-s)[c(1:6)]*log((m-s)[c(1:6)]/(m-s.hat)[c(1:6)])))
GOF


#################################
#  Person Chi-Square Test        #
#################################
var.est <- mu.mle*(1-mu.mle)/5
chi2 <- sum((y-mu.mle)^2/var.est)
chi2





