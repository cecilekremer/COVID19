
##########
## DATA ##
##########

rm(list=ls())
data <- read.csv("G:/My Drive/UHasselt/PhD/HAV outbreak data/2019nCoV/Tianjin/FINAL VERSION/Tianjin cluster cases (updated Feb 29).csv", header=T)
data$symptomOnset <- as.Date(data$symptom_onset, "%d/%m/%Y")

Case <- data$CaseID
Time <- data$symptomOnset - min(data$symptomOnset, na.rm=T) 

# Which cases are possible?
NCases  <- nrow(data)
vi.list <- list()
for (i in 1:NCases){
  vi.list[[i]] <- as.numeric(unlist(strsplit(as.character(data$Source[i]),",")))
}

PossibleInfector <- matrix(nrow=NCases,ncol=max(lengths(strsplit(as.character(data$Source),","))))

for(i in 1:NCases){
  PossibleInfector[i,1:lengths(vi.list)[i]] <- vi.list[[i]]
}

NPossibleInfector <- rowSums(!is.na(PossibleInfector))
IsNotContributorToLikel <-  c(which(PossibleInfector[,1]==0))
IsContributorToLikel <- c(which(PossibleInfector[,1]!=0))

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

## Proportion asymptomatic transmission
ci.fxn <- function(mean.gi.sample, var.gi.sample){
  xg.sample <- rgamma(n,shape=mean.gi.sample^2/var.gi.sample,rate=mean.gi.sample/var.gi.sample)
  mean(xg.sample<inc2)  
}

n=1000000

distr.p <- c()
Nresamples <- 500
progressbar <- txtProgressBar(min = 0, max = Nresamples, style = 3)
for (k in 1:Nresamples){
  set.seed(k*245+98)
  m <- runif(1, min(Savetheta[,1]), max(Savetheta[,1]))    
  v <- runif(1, min(Savetheta[,2]), max(Savetheta[,2]))  
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
Nresamples <- 5000
progressbar <- txtProgressBar(min = 0, max = Nresamples, style = 3)
for (j in 1:Nresamples){
  set.seed(j*245+98)
  m <- runif(1, min(Savetheta[,1]), max(Savetheta[,1]))  
  v <- runif(1, min(Savetheta[,2]), max(Savetheta[,2]))     
  distr.r <- c(distr.r, exp(r*m-0.5*r^2*v))
  setTxtProgressBar(progressbar, j)
}
close(progressbar)

quantile(distr.r, c(0.025, 0.5, 0.975))


