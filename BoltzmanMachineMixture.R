##############################################################
#
# Online EM for mixtures of Boltzman Machine distributions
#
##############################################################

library(BoltzMM)

# d=2 otherwise may be problematic to compute normalizing constants

# home made function for a mixture of Boltzmann Machine  density
dmixboltz <- function(x, pis, as, bs, log = F) {
  # as is a 2 x K matrix, 2 singleton parameters
  # bs is a vector of length K interactions parameters
  # x is a 2dim vector
  K = length(pis)
  dens = 0
  for (k in 1:K) {
    dens = dens + pis[k]*pfvbm(x, as[,k], -(diag(rep(bs[k],2))-bs[k]) )
  }
  if (log) {return(log(dens))}
  else {return(dens)}
}

##############################################################
# Simulated data: Mixture of Boltzmann machine distributions
# parametrization is in a1,a2,b12 

# N number of samples and iterations
N<-10^5

# Example for a 3-component mixture (K=3)
# True parameter values for K=3
pitrue<-c(1/3, 1/3,1/3)
atrue<- matrix(c(2,1,1,2,1,1), 2,3) # a1 is 2 1 1 and a2 is 1 2 1 for z=1,2,3
btrue<-c(-1,0,1)

K<-length(pitrue)

# initial parameter values for K=3
piinit<-rep(1/K, K)
ainit<-matrix(c(1,1,1,1,1,1), 2,3)
binit<-c(-2,1,2)

# Simulated data: Y is a matrix  of size 2 x N
# Home made function to simulate a bm mixture
rmixboltz<-function(Nsize,pis, as, bs){
  datasim<-matrix(0, 2, Nsize) 
  Kdim = length(pis)
  classes<-rmultinom(Nsize,1,pis)  #K x N
  boltzcomp<-array(0, c(2, Nsize, Kdim) )
  for (k in 1:Kdim) {
    boltzcomp[,,k]<-t(rfvbm(Nsize, as[,k], -(diag(rep(bs[k],2))-bs[k]) ) )
  }
  for(ii in 1:Nsize){
    datasim[,ii]<-boltzcomp[,ii,classes[,ii]==1]
  }
  return(datasim)
}

Y <- rmixboltz(N,pitrue,atrue, btrue)


# Set learning rate as indicated in the paper (gamma_i sequence)
learning <- (1-10^-10)*(1:N)^(-6/10)

# Set up sequences of parameters: pivect, amat (the 2 singleton parameters), bvect (1 interaction parameter)
pivect <- matrix(rep(piinit, N), ncol=N)
amat<-array(0, c(2,K,N) ) 
for(ii in 1:N){amat[,,ii]<-ainit}
bvect<-matrix(rep(binit, N), ncol=N)

# Define functions riK4
# create a 4K-dim vector where the first 4 elements are ri1, the following 4 are
# ri2 etc... multiplied resp by 1,  y1, y2, y1 y2  the statistics in the Boltzmann case
boltzriK4_fun <- function(y,piv,av,bv) {
  Kdim<-length(piv)
  riK4<-rep(0,4*Kdim)
  cste<-dmixboltz(y, piv,av,bv)
  for(k in 1:Kdim){
    riK4[((k-1)*4+1):(4*k)]<-(piv[k]* pfvbm(y,av[,k],-(diag(rep(bv[k],2)) - bv[k]) )/cste)
  } 
  return(riK4*rep(c(1,y[1], y[2], y[1]*y[2]),Kdim))
}

# sufficent statistics
svect<-matrix(0, 4*K, N)
# Initialize the sufficient statistics
svect[,1]<-boltzriK4_fun(Y[,1],pivect[,1],amat[,,1], bvect[,1])

###############################
# Main loop: Online EM
###############################
for(ii in 2:N){
  
  # Update sufficient statistics
  svect[,ii]<-svect[,ii-1] + learning[ii-1]*(boltzriK4_fun(Y[,ii],pivect[,ii-1],amat[,,ii-1], bvect[,ii-1]) - svect[,ii-1])
  
  # Update pi, alpha, theta
  if(ii>500){
    for(k in 1:K){
      
      # unormalized pik
      pivect[k,ii]<-svect[4*k-3,ii] 
      
      # Update a and b
      #  requires optimization (minimization ie -Q): Attention optim is not to find zeros!!
      ob <- function(x) {
        a1<-x[1]
        a2<-x[2]
        b12<-x[3]
        -(a1*svect[4*k-2,ii]/svect[4*k-3,ii] + a2*svect[4*k-1,ii]/svect[4*k-3,ii] + b12*svect[4*k,ii]/svect[4*k-3,ii]
        # normalizing cste is computed explicitly in this small case
            - log(exp(a1+a2+b12)+exp(-a1-a2+b12)+exp(a1-a2-b12)+exp(-a1+a2-b12)) ) 
      }
      
      res<-optim(c(amat[,k,ii-1],bvect[k,ii-1]),ob,method='BFGS')$par
      amat[,k,ii] <- c(res[1], res[2])
      bvect[k,ii]<- res[3]
      
    } 
    pivect[,ii]<-pivect[,ii]/sum(pivect[,ii])
  } # end if
} # end for main loop

#############################################################
# Plot the sequences
# First component
plot(501:N,amat[1,1,-c(1:500)],type='l',xlab='i',ylab='a11')
abline(h=2,col='red',lty=2,lwd=2)
grid()
plot(501:N,amat[2,1,-c(1:500)],type='l',xlab='i',ylab='a21')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,bvect[1,-c(1:500)],type='l',xlab='i',ylab='b121')
abline(h=-1,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[1,-c(1:500)],type='l',xlab='i',ylab='pi1')
abline(h=1/3,col='red',lty=2,lwd=2)
grid()

# Second comp
plot(501:N,amat[1,2,-c(1:500)],type='l',xlab='i',ylab='a12')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,amat[2,2,-c(1:500)],type='l',xlab='i',ylab='a22')
abline(h=2,col='red',lty=2,lwd=2)
grid()
plot(501:N,bvect[2,-c(1:500)],type='l',xlab='i',ylab='b122')
abline(h=0,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[2,-c(1:500)],type='l',xlab='i',ylab='pi2')
abline(h=1/3,col='red',lty=2,lwd=2)
grid()

# Third comp
plot(501:N,amat[1,3,-c(1:500)],type='l',xlab='i',ylab='a13')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,amat[2,3,-c(1:500)],type='l',xlab='i',ylab='a23')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,bvect[3,-c(1:500)],type='l',xlab='i',ylab='b123')
abline(h=1,col='red',lty=2,lwd=2)
grid()
plot(501:N,pivect[3,-c(1:500)],type='l',xlab='i',ylab='pi3')
abline(h=1/3,col='red',lty=2,lwd=2)
grid()


## CHECK for non identifiability issue
# Estimations:  means over the last 50 iterations/ over 100000 here
# eg for pi_k
piestim<-rowMeans(pivect[,100000-100:100000])
#0.2528333 0.3720191 0.3751476
# for a
toto<-matrix(0,2,3)
for(i in (100000-100):100000){toto<-toto+amat[,,i]}
aestim<-toto/100
#[,1]     [,2]     [,3]
#[1,] 0.6777772 1.490487 1.598868
#[2,] 0.2501806 1.123918 1.479846
bestim<-rowMeans(bvect[,100000-100:100000])
# [1] -1.2152052  0.9759693  1.6695007

## true pdf
dmixboltz(c(1,1), pitrue,atrue,btrue)
#[1] 0.7602125 -- 0.76 vs 0.76
## estimated pdf
dmixboltz(c(1,1), piestim,aestim,bestim)
#[1] 0.7571592

dmixboltz(c(-1,1), pitrue,atrue,btrue)
#[1] 0.06590868 -- 0.07 vs 0.07
## estimated pdf
dmixboltz(c(-1,1), piestim,aestim,bestim)
#0.07071415

dmixboltz(c(-1,-1), pitrue,atrue,btrue)
#[1] 0.006888386  -- 0.01 vs 0.01
## estimated pdf
dmixboltz(c(-1,-1), piestim,aestim,bestim)
#0.006634842

dmixboltz(c(1,-1), pitrue,atrue,btrue)
#[1]  0.1669905 -- 0.17 vs 0.17
## estimated pdf
dmixboltz(c(1,-1), piestim,aestim,bestim)
# 0.1654918
