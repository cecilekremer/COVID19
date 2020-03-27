
##########
## DATA ##
##########

rm(list = ls())

### DATA UPDATED 27/02/2020

Cluster <- c(rep(1,8), rep(5,23), rep(2,8), rep(3,3), rep(4,5), rep(6,4), rep(7,3))
Case.ID <- c(8,9,31,33,38,83,90,91,66,68,70,71,80,84,88,48,49,51,54,57,58,60,53,61,62,63,67,73,74,78,81,19,20,21,24,25,27,34,40,30,36,39,42,47,52,56,69,2,13,26,44,50,55,77)
Case <- c(1:length(Case.ID))
Time <- c(4,4,3,10,14,8,20,3,9,10,14,12,15,15,27,12,14,15,21,22,21,19,21,17,20,21,20,20,23,20,27,9,5,13,10,4,12,7,10,1,4,9,12,17,18,23,22,1,8,8,11,18,10,21)
vi <- c(0,0,rep(NA,6), rep(NA,18),24,18,rep(NA,3),0,0,32,0,0,32,rep(NA,8),44,NA,0,0,0,NA,NA,0,52)

data <- data.frame(Cluster,Case,Time,vi)

NCases <- length(Case)

# Which cases are possible?
vi.list <- list()

vi.list[[3]] <- c(1,2,4:8)
vi.list[[4]] <- c(1,2,3,5:8)
vi.list[[5]] <- c(1:4,6,7,8)
vi.list[[6]] <- c(1,2)
vi.list[[7]] <- c(1:6,8)
vi.list[[8]] <- c(1,2)

vi.list[[9]] = vi.list[[10]] = vi.list[[11]] = vi.list[[12]] = vi.list[[13]] = vi.list[[14]] = vi.list[[15]] <- c(6,8)
vi.list[[16]] <- c(9, 17:22)
vi.list[[17]] <- c(9, 16,18:22)
vi.list[[18]] <- c(9, 16,17,19:22)
vi.list[[19]] <- c(9, 16:18,20:22)
vi.list[[20]] <- c(9, 16:19,21,22)
vi.list[[21]] <- c(9, 16:20,22)
vi.list[[22]] <- c(9, 16:21)
vi.list[[23]] <- c(9, 16:22, 24:31)
vi.list[[24]] <- c(9, 27)
vi.list[[25]] <- c(9, 16:22, 23,24,26:31)
vi.list[[26]] <- c(9, 16:22, 23:25,27:31)
vi.list[[29]] <- c(9, 16:22, 23:28,30,31)
vi.list[[30]] <- c(9, 16:22,23:29,31)
vi.list[[31]] <- c(9, 16:22, 23:30)

vi.list[[38]] <- c(32,33)
vi.list[[39]] <- c(32,33)

vi.list[[40]] <- c(41,42)
vi.list[[41]] <- c(40,42)
vi.list[[42]] <- c(40,41)

vi.list[[43]] <- c(44,45)
vi.list[[44]] <- c(43,45,46)
vi.list[[45]] <- c(43,44)
vi.list[[47]] <- c(43:46)

vi.list[[51]] <- c(49,50)

vi.list[[52]] <- c(53,54)
vi.list[[53]] <- c(0)
vi.list[[54]] <- c(52)

for(i in 1:NCases){
  if(is.null(vi.list[[i]])) vi.list[[i]] <- vi[i] 
}
 

PossibleInfector <- matrix(nrow=NCases,ncol=max(lengths(vi.list)))

for(i in 1:NCases){
  PossibleInfector[i,1:lengths(vi.list)[i]] <- vi.list[[i]]
}
NPossibleInfector    <- rowSums(!is.na(PossibleInfector))
IsNotContributorToLikel <- c(which(vi==0)) # index cases
IsContributorToLikel <- Case[!Case%in%IsNotContributorToLikel]

###############
# likelihood
likelihood <- function(){
alpha<-c(); Beta<-c() 
# alpha[1] & Beta[1] are shape & rate of the serial interval distriution  
alpha[1] <- theta[1]^2/theta[2] # shape = mean^2/var
Beta[1]  <- theta[1]/theta[2]; 	# rate = mean/var
# alpha[2] & Beta[2] are shape & rate of the distribution of sum of 2 independent gamma variables
alpha[2] <- theta[3]^2/theta[4] # shape = mean^2/var
Beta[2]  <- theta[3]/theta[4]   # rate = mean/var
# evaluate density
f_Z <- function(Z){
monteCarloN <- 300
delta_i     <- rgamma(monteCarloN, shape=alpha[2], rate=Beta[2])
delta_j     <- rgamma(monteCarloN, shape=alpha[2], rate=Beta[2])
Y           <- delta_i-delta_j  
Z_Y         <- Z - Y
return(mean(dgamma(Z_Y, shape=alpha[1], rate=Beta[1], log=F)))
}
SerialInterval   <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]] # serial interval
return(sum(log(1e-50+sapply(SerialInterval, function(x) f_Z(x)))))
}

# prior
prior <- function(){
alpha<-c(); Beta<-c()   
alpha[1] <- theta[1]^2/theta[2] # shape = mean^2/var
Beta[1]  <- theta[1]/theta[2]; 	# rate = mean/var  
alpha.prior <-  dunif(alpha[1], 0, 30)
Beta.prior  <-  dunif(Beta[1], 0, 20)
return(log(1e-50+alpha.prior)+log(1e-50+Beta.prior))
}

# posterior
posterior <- function(){
return (likelihood()+prior())
}

#--------------------------------------------------
#----MCMC Algorithm--------------------------------
#--------------------------------------------------
Network <- numeric(NCases)+0
Update <- IsContributorToLikel
Draw <- round(runif(length(Update),min=0.5,max=NPossibleInfector[Update]+0.5))
for(i in 1:length(IsContributorToLikel)){
  Network[Update[i]] <- PossibleInfector[Update[i],Draw[i]]
}
AcceptedNetwork <- Network

AcceptedTheta=theta <- c(1, 1, 5.2, 2.8^2) 

SerialInterval   <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]] # serial interval

P <- posterior()
AcceptedP <- P

NRuns <- 3000000
NUpdate <- length(IsContributorToLikel)
Burnin <- 500000
Thinning<-200
SaveP <- numeric()
SaveNetwork <- matrix(nrow=NCases,ncol=(NRuns-Burnin)/Thinning)
Savetheta <- matrix(nrow=(NRuns-Burnin)/Thinning,ncol=(1+length(theta))) 

anetwork=asd<-0
tuning <- c(0.5,0.5)
a <- 0
sd = 0.5

progressbar <- txtProgressBar(min = 0, max = NRuns, style = 3)
ptm <- proc.time()
for(b in 1:NRuns){
  
  if(b%%2 != 0){
    theta <- AcceptedTheta
    Update <- IsContributorToLikel
    Draw <- round(runif(length(Update),min=0.5,max=NPossibleInfector[Update]+0.5))
    for(i in 1:NUpdate){ 
      Network[Update[i]] <- PossibleInfector[Update[i],Draw[i]]
    }
  }
  
  
  if(b%%2 == 0){
    Network <- AcceptedNetwork
    theta[1] <- runif(1, (AcceptedTheta[1]-tuning[1]), (AcceptedTheta[1]+tuning[1]))
    if(theta[1]<0){theta[1] <- AcceptedTheta[1]}
    
    theta[2] <- runif(1, (AcceptedTheta[2]-tuning[2]), (AcceptedTheta[2]+tuning[2]))
    if(theta[2]<0){ theta[2] <- AcceptedTheta[2]}
    
  }
  
  P <- posterior()
  
  AcceptYN <- runif(1,min=0,max=1) <= exp(P-AcceptedP)
  if(AcceptYN){
    if(b%%2 != 0){
      anetwork <- anetwork + 1
      AcceptedNetwork <- Network
    }
    if(b%%2 == 0){
      asd <- asd+1
      AcceptedTheta <- theta
    }
    AcceptedP <- P
  }
  if(b%%Thinning == 0 & b>Burnin){
    a <- a + 1
    Savetheta[a,] <- c(AcceptedTheta, (AcceptedTheta[2]+2*AcceptedTheta[4]))
    SaveNetwork[,a] <- AcceptedNetwork
    SaveP[a] <- AcceptedP
  }
  setTxtProgressBar(progressbar, b)
}
close(progressbar)
proc.time() - ptm


#####################
## SERIAL INTERVAL ##

par(mfrow=c(1,1))
n = 1000000
inc1=rgamma(n,shape=5.2^2/(2.8^2),rate = 5.2/(2.8^2))
inc2=rgamma(n,shape=5.2^2/(2.8^2),rate = 5.2/(2.8^2))
hist(inc1); hist(inc2)
xg=rgamma(n,shape=5.20^2/1.72^2,rate=5.20/1.72^2)
hist(xg)
mean(xg)
xs=xg+inc1-inc2
hist(xs)
mean(xs); quantile(xs, c(0.025,0.5,0.975))
sd(xs)

# Proportion pre-symptomatic transmission

n=1000

inc2=rgamma(n,shape=5.2^2/(2.8^2),rate = 5.2/(2.8^2))

ci.fxn <- function(mean.gi.sample, var.gi.sample){
  xg.sample <- rgamma(n,shape=mean.gi.sample^2/var.gi.sample,rate=mean.gi.sample/var.gi.sample)
  mean(xg.sample<inc2)  
}

distr.p <- c()
Nresamples <- nrow(Savetheta)
progressbar <- txtProgressBar(min = 0, max = Nresamples, style = 3)
for (k in 1:Nresamples){
  set.seed(k*245+98)
  m <- Savetheta[k,1] 
  v <- Savetheta[k,2] 
  distr.p <- c(distr.p, ci.fxn(m, v))
  setTxtProgressBar(progressbar, k)
}
close(progressbar)

quantile(distr.p, c(0.025, 0.5, 0.975))

## Calculating R0

Time <- c(1,1,3,4,5,4,4,4,1,7,6,8,8,10,3,10,10,9,5,13,10,4,8,12,16,8,1,3,13,10,7,10,4,10,14,9,10,12,12,10,11,16,17,12,14,18,15,18,21,21,10,23,22,21,18,19,17,20,21,14,9,20,10,22,14,12,21,20,23,21,20,23,15,27,20,8,15,25,25,27,14,20,3)

EpiCurve <- data.frame(table(Time))
I        <- EpiCurve$Freq[1:which(EpiCurve$Freq==max(EpiCurve$Freq))] # incidence curve (exponential phase)
plot(I, type="l")

t = 1:length(I)

growth.rate <- lm(log(I) ~ t)
r <- coef(growth.rate)[2]

distr.r <- c()
Nresamples <- nrow(Savetheta)
progressbar <- txtProgressBar(min = 0, max = Nresamples, style = 3)
for (j in 1:Nresamples){
  set.seed(j*245+98)
  m <- Savetheta[j,1] 
  v <- Savetheta[j,2]    
  distr.r <- c(distr.r, exp(r*m-0.5*r^2*v))
  setTxtProgressBar(progressbar, j)
}
close(progressbar)


quantile(distr.r, c(0.025, 0.5, 0.975))
