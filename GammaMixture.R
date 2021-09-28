
##############################################################
#
# Online EM for mixtures of gamma distributions
#
##############################################################

# Libraries
library(pracma) # for the computation of the psi (digamma) function
library(bmixture) # to simulate mixtures of Gamma distributions 

##############################################################
# Simulated data: Mixture of Gamma distributions
# parametrization is in (alpha, beta) = (k, 1/theta)

# N number of samples and iterations
N<-10^5

# Example for a 3-component mixture (K=3)
# True parameter values for K=3
pitrue<-c(0.3, 0.2,0.5) # mixture weights
atrue<- c(9,20,1) # alpha or k
ttrue<-c(0.5,0.5,1) # theta or 1/beta

K<-length(pitrue)

# Initial parameter values 
piinit<-rep(1/K, K)
# Initial values for K=3
ainit<-c(5,9,2)
tinit<-c(1,1,1)

#simulated data: Y is a vector of length N
Y <- rmixgamma(N,weight=pitrue,alpha=atrue, beta=1/ttrue)

# To check the separability of the mixture visually
hist(Y)

# Set learning rate as indicated in the paper (gamma_i sequence)
learning <- (1-10^-10)*(1:N)^(-6/10)

# Set up sequences of parameters: pivect, alphavect, thetavect
# Remark: it is convenient to put the initial values in each column for the
# "burnin" period which does not update the parameters but only the
# sufficient statistics
pivect <- matrix(rep(piinit, N), ncol=N)  # size 3 x N
alphavect<-matrix(rep(ainit, N), ncol=N)  # size 3 x N
thetavect<-matrix(rep(tinit, N), ncol=N)  # size 3 x N

# Compute \bar{s}: in the mixture case it involves the responsibilities (rik) and depending
# on the mixture, the "\bar{s}" of the component distribution, here a gamma
# 1) compute the rik
# 2) compute the gamma specific statistics and add 1, ie (1, stat1, stat2, etc...)

# Define function riK3 (3 because 1+2 Gamma specific statistics)
# create a 3K-dim vector where the first 3 elements are ri1, the following 3 are
# ri2 etc... multiplied respectively by 1, log y, y (gamma specific statistics)
riK3_fun <- function(y,piv,av,thetav) {
 Kdim<-length(piv)
 riK3<-rep(0,3*Kdim)
 cste<-dmixgamma(y, piv,av,beta=(1/thetav))
 for(k in 1:Kdim){riK3[((k-1)*3+1):(3*k)]<-(piv[k]* dgamma(y,av[k], scale=thetav[k])/cste)} 
 return(riK3*rep(c(1,log(y), y),Kdim))
}

# sequence of "s" sufficent statistics
svect<-matrix(0, 3*K, N)
# Initialize the sufficient statistics with \bar{s} 
svect[,1]<-riK3_fun(Y[1],pivect[,1],alphavect[,1], thetavect[,1])

#############################
# Main loop: Online EM
#############################
for(ii in 2:N){

# Update sufficient statistics
svect[,ii]<-svect[,ii-1] + learning[ii-1]*(riK3_fun(Y[ii],pivect[,ii-1],alphavect[,ii-1], thetavect[,ii-1]) - svect[,ii-1])

# After 500 iterations, start updating parameters: pi, alpha, theta
  if(ii>500){
for(k in 1:K){

# unormalized pik
pivect[k,ii]<-svect[3*k-2,ii] 

# Update alpha and theta
# alpha (= k) requires optimization: Attention function optim is not to find zeros!!!

# ob is actually just -Q as given in the paper for a single component but where each
# statistic has to be divided by "s0", the statistic associated to "1"
 ob <- function(x) {
   -((x-1)*svect[3*k-1,ii]/svect[3*k-2,ii] - x + x*log(svect[3*k-2,ii]*x/svect[3*k,ii]) - lgamma(x))
 }

     alphavect[k,ii] <- optim(alphavect[k,ii-1],ob)$par

# theta can be deduced from alpha (see paper)
thetavect[k,ii] <- svect[3*k,ii]/(svect[3*k-2,ii]*alphavect[k,ii])
} 

pivect[,ii]<-pivect[,ii]/sum(pivect[,ii])
  } # end if
} # end for main loop

#############################
# Plot the sequences
# To save a  pdf
#pdf(file = "name.pdf", width = 5.5, height = 4)
#[...]
#dev.off()

# First component
plot(501:N,alphavect[1,-c(1:500)],type='l',xlab='i',ylab='k1')
abline(h=9,col='red',lty=2,lwd=2)
grid()
plot(501:N,thetavect[1,-c(1:500)],type='l',xlab='i',ylab='theta1')
abline(h=0.5,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[1,-c(1:500)],type='l',xlab='i',ylab='pi1')
abline(h=0.3,col='red',lty=2,lwd=2)
grid()

# Second component
plot(501:N,alphavect[2,-c(1:500)],type='l',xlab='i',ylab='k2')
abline(h=20,col='red',lty=2,lwd=2)
grid()
plot(501:N,thetavect[2,-c(1:500)],type='l',xlab='i',ylab='theta2')
abline(h=0.5,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[2,-c(1:500)],type='l',xlab='i',ylab='pi2')
abline(h=0.2,col='red',lty=2,lwd=2)
grid()

# Third component
plot(501:N,alphavect[3,-c(1:500)],type='l',xlab='i',ylab='k3')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,thetavect[3,-c(1:500)],type='l',xlab='i',ylab='theta3')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[3,-c(1:500)],type='l',xlab='i',ylab='pi3')
abline(h=.5,col='red',lty=2,lwd=2)
grid()
