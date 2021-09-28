##############################################################
#
# Online EM for mixtures of univariate Student distributions
#
##############################################################


# Libraries
library(pracma) # for the psi (digamma) function
library(bmixture) # to simulate mixtures of univariate Student distributions
library(extraDistr) # for the Student distribution 

##############################################################
# Simulated data: Mixture of Student distributions
# parametrization is in mu, sigma2, nu

# N number of samples and iterations
N<-5*10^5

# Example for a 3-component mixture (K=3)
# True paramter values for K=3
pitrue<-c(.3, .5, .2) # mixture weights
mutrue<- c(0,3,6)     # mixture means
sigma2true<-c(1,1,1)  # mixture variances
nutrue<-c(3,2,1)      # degree-of-freedom parameters

# Remark: it seems that lower dof are easier to estimate, maybe the likelihood is less
# flat than for higher dof ?
# Then a lower dof (eg 1) maybe ok with a small weight while higher dof may need
# a bigger weight...

K<-length(pitrue)

# Initial parameter values 
piinit<-rep(1/K, K)
# initial values for K=3
muinit<-c(1,4,7)
sigma2init<-c(2,2,2)
nuinit<-c(4,4,4)

# Simulated data: Y is a vector of length N
Y <- rmixt(N,weight=pitrue, df= nutrue, mean=mutrue, sd=sqrt(sigma2true))

# To check the separability of the mixture visually
hist(Y) # hist(Y, nclass=500, xlim=c(-10,10))

# Set learning rate as indicated in the paper (gamma_i sequence)
learning <- (1-10^-10)*(1:N)^(-6/10)

# Set up sequences of parameters: pivect, muvect, sigma2vect, nuvect
# Remark: it is convenient to put the initial values in each column for the
# "burnin" period which does not update the parameters but only the
# sufficient statistics
pivect <- matrix(rep(piinit, N), ncol=N)
muvect<-matrix(rep(muinit, N), ncol=N)
sigma2vect<-matrix(rep(sigma2init, N), ncol=N)
nuvect<-matrix(rep(nuinit, N), ncol=N)

# Define functions  to compute U and U_tilde
U_fun <- function(y,nu,mu,sigma2) {
  (nu + 1)/(nu + (y-mu)^2/sigma2)
}

Util_fun <- function(y,nu,mu,sigma2) {
  digamma(nu/2 + 1/2) - log(nu/2 + 0.5*(y-mu)^2/sigma2)
}


# Compute \bar{s}: in the mixture case it involves the responsibilities (rik) and depending
# on the mixture, the "\bar{s}" of the component distribution, here a Student
# 1) compute the rik
# 2) compute the Student specific statistics and add 1, ie (1, stat1, stat2, etc...)

# Define function riK5 (5 because 1+4 Student specific statistics)
# create a 5K-dim vector where the first 5 elements are ri1, the following 5 are
# ri2 etc... multiplied resp by 1, u y, u y^2, u, tildeu (Student specific)
StudentriK5_fun <- function(y,piv,muv,sigma2v, nuv) {
  Kdim<-length(piv)
  riK5<-rep(0,5*Kdim)
  cste<-dmixt(y, piv,nuv, muv, sd=sqrt(sigma2v) )
  tstatv<-rep(0,5*Kdim)
  for(k in 1:Kdim){
    # student pdf with dlst (extrDistr)
    riK5[((k-1)*5+1):(5*k)]<-(piv[k]*dlst(y, df=nuv[k],mu=muv[k], sigma=sqrt(sigma2v[k])  )/cste)
    # compute u and tildeu for each k
    u<-U_fun(y,nuv[k],muv[k],sigma2v[k])
    tildeu<-Util_fun(y,nuv[k],muv[k],sigma2v[k])
    tstatv[((k-1)*5+1):(5*k)]<-c(1,u*y,u*y*y,u,tildeu)
    } 
  return(riK5*tstatv)
}

# sequence of "s" sufficent statistics
svect<-matrix(0, 5*K, N)
# Initialize the sufficient statistics with \bar{s} 
svect[,1]<-StudentriK5_fun(Y[1],pivect[,1],muvect[,1], sigma2vect[,1],nuvect[,1])

###############################
# Main loop: Online EM
###############################
for(ii in 2:N){
  
  # Update sufficient statistics
  svect[,ii]<-svect[,ii-1] + learning[ii-1]*(StudentriK5_fun(Y[ii],pivect[,ii-1],muvect[,ii-1], sigma2vect[,ii-1],nuvect[,ii-1]) - svect[,ii-1])
  
  # After 500 iterations, start updating parameters: pi, mu, sigma2, nu
  if(ii>500){
    for(k in 1:K){
      
      # unormalized pik
      pivect[k,ii]<-svect[5*k-4,ii] 
      
      # Update mu and sigma2: is closed form
      muvect[k,ii] <- svect[5*k-3,ii]/svect[5*k-1,ii]
      sigma2vect[k,ii] <- svect[5*k-2,ii]/svect[5*k-4,ii] - svect[5*k-3,ii]^2/(svect[5*k-1,ii]*svect[5*k-4,ii])
      
      #update nu
      ob <- function(x) {
        - (
          svect[5*k,ii]*(x-1)/svect[5*k-4,ii] - svect[5*k-1,ii]*x/svect[5*k-4,ii] - lgamma(x) + x*log(x)
        )
      }
      opti <- optim(nuvect[k,ii-1]/2,ob)
      nuvect[k,ii] <- 2*opti$par[1]
      
    } 
    
    pivect[,ii]<-pivect[,ii]/sum(pivect[,ii])
  } # end if
} # end for main loop

##########################
# Plot the sequences

# First component
plot(501:N,muvect[1,-c(1:500)],type='l',xlab='i',ylab='mu1')
abline(h=0,col='red',lty=2,lwd=2)
grid()
plot(501:N,sigma2vect[1,-c(1:500)],type='l',xlab='i',ylab='sigma21')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,nuvect[1,-c(1:500)],type='l',xlab='i',ylab='nu1')
abline(h=3,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[1,-c(1:500)],type='l',xlab='i',ylab='pi1')
abline(h=.3,col='red',lty=2,lwd=2)
grid()

# Second component
plot(501:N,muvect[2,-c(1:500)],type='l',xlab='i',ylab='mu2')
abline(h=3,col='red',lty=2,lwd=2)
grid()
plot(501:N,sigma2vect[2,-c(1:500)],type='l',xlab='i',ylab='sigma22')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,nuvect[2,-c(1:500)],type='l',xlab='i',ylab='nu2')
abline(h=2,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[2,-c(1:500)],type='l',xlab='i',ylab='pi2')
abline(h=.5,col='red',lty=2,lwd=2)
grid()

# Third component
plot(501:N,muvect[3,-c(1:500)],type='l',xlab='i',ylab='mu3')
abline(h=6,col='red',lty=2,lwd=2)
grid()
plot(501:N,sigma2vect[3,-c(1:500)],type='l',xlab='i',ylab='sigma23')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,nuvect[3,-c(1:500)],type='l',xlab='i',ylab='nu3')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[3,-c(1:500)],type='l',xlab='i',ylab='pi3')
abline(h=.2,col='red',lty=2,lwd=2)
grid()


