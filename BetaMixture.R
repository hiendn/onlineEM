##############################################################
#
# Online EM for mixtures of beta distributions
#
##############################################################

library(pracma) # for the psi (digamma) function


# Mixture of Beta distributions density
dmixbeta <- function(x, pis, alphas, betas, log = F) {
  K = length(pis)
  dens = 0
  for (k in 1:K) {
    dens = dens + pis[k]*dbeta(x, alphas[k], betas[k])
  }
  if (log) {return(log(dens))}
  else {return(dens)}
}

##############################################################
# Simulated data: Mixture of Beta distributions
# parametrization is in alpha beta 

# N number of samples and iterations
N<-3*10^5

# Example for a 3-component mixture (K=3)
# True values for K=3 
pitrue<-c(.5, .25,.25)
atrue<- c(3,2,9)
btrue<-c(1,2,1)

K<-length(pitrue)

# Initial parameter values for K=3, to be changed accordingly
piinit<-rep(1/K, K)
ainit<-c(2,2,10)
binit<-c(1,1,2)


# Simulated data: Y is a vector of length N
# Home made function to simulate a beta mixture
rmixbeta<-function(Nsize,pis, alphas, betas){
  datasim<-rep(0, Nsize) 
  Kdim = length(pis)
  classes<-rmultinom(Nsize,1,pis)  #K x N
  betacomp<-matrix(0,Kdim, Nsize)
  for (k in 1:Kdim) {
    betacomp[k,]<-rbeta(Nsize, alphas[k], betas[k])
  }
  for(ii in 1:Nsize){
    datasim[ii]<-betacomp[classes[,ii]==1,ii]
  }
  return(datasim)
}

Y <- rmixbeta(N,pitrue,atrue, btrue)
# check
hist(Y)

# Set learning rate as indicated in the paper (gamma_i sequence)
learning <- (1-10^-10)*(1:N)^(-6/10)


# Set up sequences of  parameters: pivect, alphavect, betavect
pivect <- matrix(rep(piinit, N), ncol=N)
alphavect<-matrix(rep(ainit, N), ncol=N)
betavect<-matrix(rep(binit, N), ncol=N)

# Define functions riK3
# create a 3K-dim vector where the first 3 elements are ri1, the following 3 are
# ri2 etc... multiplied resp by 1, log y, log(1-y) for the beta distribution
betariK3_fun <- function(y,piv,av,bv) {
  Kdim<-length(piv)
  riK3<-rep(0,3*Kdim)
  cste<-dmixbeta(y, piv,av,bv)
  for(k in 1:Kdim){
    riK3[((k-1)*3+1):(3*k)]<-(piv[k]* dbeta(y,av[k], bv[k])/cste)
    } 
  return(riK3*rep(c(1,log(y), log(1-y) ),Kdim))
}

# sufficent statistics
svect<-matrix(0, 3*K, N)
# sbar just for testing 
# sbarvect<-matrix(0, 3*K, N)
# Initialize the sufficient statistics
svect[,1]<-betariK3_fun(Y[1],pivect[,1],alphavect[,1], betavect[,1])


###############################
# Main loop: Online EM
###############################
for(ii in 2:N){
  
  # Update sufficient statistics
  svect[,ii]<-svect[,ii-1] + learning[ii-1]*(betariK3_fun(Y[ii],pivect[,ii-1],alphavect[,ii-1], betavect[,ii-1]) - svect[,ii-1])
  ## test for check
  ##sbarvect[,ii]<-betariK3_fun(Y[ii],pivect[,ii-1],alphavect[,ii-1], betavect[,ii-1])
  
  # Update pi, alpha, theta
  if(ii>500){
    for(k in 1:K){
      
      # unormalized pik
      pivect[k,ii]<-svect[3*k-2,ii] 
      
      # Update alpha and beta
      # requires optimization (minimization ie -Q): Attention optim is not to find zeros!!
      # For beta the optim is jointly on the 2 parameters
      ob <- function(x) {
        a<-x[1]
        b<-x[2]
        -(a*svect[3*k-1,ii]/svect[3*k-2,ii] + b*svect[3*k,ii]/svect[3*k-2,ii]  -lgamma(a) - lgamma(b) + lgamma(a+b))
      }
      
      res<-optim(c(alphavect[k,ii-1],betavect[k,ii-1]),ob)$par
      alphavect[k,ii] <- res[1]
      betavect[k,ii]<- res[2]
  
    } 
    pivect[,ii]<-pivect[,ii]/sum(pivect[,ii])
  } # end if
} # end for main loop

##########################
# Plot the sequences
# First component
plot(501:N,alphavect[1,-c(1:500)],type='l',xlab='i',ylab='alpha1', ylim=c(1.1, 3.1))
#plot(501:N,alphavect[1,-c(1:500)],type='l',xlab='i',ylab='alpha1')
abline(h=3,col='red',lty=2,lwd=2)
grid()
plot(501:N,betavect[1,-c(1:500)],type='l',xlab='i',ylab='beta1')
abline(h=1,col='red',lty=2,lwd=2)
grid()
#plot(501:N,pivect[1,-c(1:500)],type='l',xlab='i',ylab='pi1')
plot(501:N,pivect[1,-c(1:500)],type='l',xlab='i',ylab='pi1', ylim=c(0.355, 0.505))
abline(h=.5,col='red',lty=2,lwd=2)
grid()

# Second component
plot(501:N,alphavect[2,-c(1:500)],type='l',xlab='i',ylab='alpha2')
abline(h=2,col='red',lty=2,lwd=2)
grid()
#plot(501:N,betavect[2,-c(1:500)],type='l',xlab='i',ylab='beta2')
plot(501:N,betavect[2,-c(1:500)],type='l',xlab='i',ylab='beta2',ylim=c(0.69,2.02))
abline(h=2,col='red',lty=2,lwd=2)
grid()
#plot(501:N,pivect[2,-c(1:500)],type='l',xlab='i',ylab='pi2')
plot(501:N,pivect[2,-c(1:500)],type='l',xlab='i',ylab='pi2', ylim=c(0.24, 0.415))
abline(h=.25,col='red',lty=2,lwd=2)
grid()

# Third component
plot(501:N,alphavect[3,-c(1:500)],type='l',xlab='i',ylab='alpha3')
abline(h=9,col='red',lty=2,lwd=2)
grid()
plot(501:N,betavect[3,-c(1:500)],type='l',xlab='i',ylab='beta3')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[3,-c(1:500)],type='l',xlab='i',ylab='pi3')
abline(h=.25,col='red',lty=2,lwd=2)
grid()


## CHECK for non identifiability issue
# Estimations:  means over the last 50 iterations over 300000 here
# eg for pi_k
rowMeans(pivect[,300000-50:300000])
#[1] 0.4011144 0.4011144 0.1977712

# Log-likelihoods
# True loglike
sum(dmixbeta(Y,pitrue, atrue,btrue, log=T))
#[1] 100526.1
# log-like for last 50 means
sum(dmixbeta(Y,rowMeans(pivect[,300000-50:300000]), rowMeans(alphavect[,300000-50:300000]),rowMeans(betavect[,300000-50:300000]), log=T))
#[1] 100521.4
# same after truncation to 2 digits
sum(dmixbeta(Y,c(0.4,0.4,0.2), c(1.99, 1.99,10.40),c(0.93,0.93,1.12), log=T))
#[1] 100519.1
# For the pdf comparaison:  plot
plot(xval,dmixbeta(xval,pitrue, atrue,btrue), type="l", col="red", ylab="Beta mixture pdf", xlab="y")
points(xval,dmixbeta(xval,c(0.4,0.4,0.2), rowMeans(alphavect[,300000-50:300000]),rowMeans(betavect[,300000-50:300000])), type="l")
