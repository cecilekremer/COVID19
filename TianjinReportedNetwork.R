
setwd("/Users/tapiwaganyani/Google Drive/RESEARCH/Corona/Codes/R/Tianjin/Final")

rm(list=ls())
data <- read.csv("Tianjin cluster cases (updated Feb 29).csv", header=T)
data$symptomOnset <- as.Date(data$symptom_onset, "%d/%m/%Y")

Case <- data$CaseID
Time <- data$symptomOnset - min(data$symptomOnset, na.rm=T) 

# Which cases are possible?
reportedNegSI <- rep(F, length(data$CaseID))
for (i in data$CaseID){
all.sources <- as.numeric(unlist(strsplit(as.character(data$Source[i]),","))) 
if(!all.sources[1]==0){
all.sources.neg <- c()
for (k in all.sources){
if((data$symptomOnset[i]-data$symptomOnset[k])<0){
all.sources.neg <- c(all.sources.neg,k)  
}  
}  
if(!is.null(all.sources.neg)){
reportedNegSI[i] <- T 
}
}
}

vi.list <- list()
for(i in which(!reportedNegSI==T)){
all.sources <- as.numeric(unlist(strsplit(as.character(data$Source2[i]),",")))   
if(all.sources[1]==0){
vi.list[[i]] <- all.sources  
}else{
all.sources.pos <- c()
for (k in all.sources){
if((data$symptomOnset[i]-data$symptomOnset[k])>=0){
all.sources.pos <- c(all.sources.pos,k)  
}  
}  
if(!is.null(all.sources.pos)){
vi.list[[i]] <- all.sources.pos 
}else{
vi.list[[i]] <- 0 
}
}
}

for(i in which(reportedNegSI==T)){
vi.list[[i]] <- as.numeric(unlist(strsplit(as.character(data$Source[i]),","))) 
}



NCases  <- nrow(data)
PossibleInfector <- matrix(nrow=NCases,ncol=max(lengths(vi.list)))

for(i in 1:NCases){
PossibleInfector[i,1:lengths(vi.list)[i]] <- vi.list[[i]]
}

NPossibleInfector <- rowSums(!is.na(PossibleInfector))
IsNotContributorToLikel <-  c(which(PossibleInfector[,1]==0))
IsContributorToLikel <- c(which(PossibleInfector[,1]!=0))


##################
SI <- c()
for(b in 1:5000){
Network <- numeric(NCases)+0
Update <- IsContributorToLikel
Draw <- round(runif(length(Update),min=0.5,max=NPossibleInfector[Update]+0.5))
for(i in 1:length(IsContributorToLikel)){
Network[Update[i]] <- PossibleInfector[Update[i],Draw[i]]
}
SI<- c(SI, Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]])
}
hist(SI, ylab = "Frequency", main = "", yaxt='n')


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
delta_i     <- rgamma(monteCarloN, shape = alpha[2], rate = Beta[2])
delta_j     <- rgamma(monteCarloN, shape = alpha[2], rate = Beta[2])
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
#----MCMC---------------------
#--------------------------------------------------
Network <- numeric(NCases)+0
Update <- IsContributorToLikel
Draw <- round(runif(length(Update),min=0.5,max=NPossibleInfector[Update]+0.5))
for(i in 1:length(IsContributorToLikel)){
Network[Update[i]] <- PossibleInfector[Update[i],Draw[i]]
}
AcceptedNetwork <- Network

m.ll <- 5.2; m.ul <- 5.2 # lower and upper limit for mean 
v.ll <- 2.8^2; v.ul <- 2.8^2 # lower and upper limit for variance

AcceptedTheta=theta <- c(1, 1, m.ll, v.ll) 

P <- posterior()
AcceptedP <- P

NRuns <- 3000000
NUpdate <- length(IsContributorToLikel)
Burnin <- 100000
Thinning<-200
SaveP <- numeric()
SaveNetwork <- matrix(nrow=NCases,ncol=(NRuns-Burnin)/Thinning)
Savetheta <- matrix(nrow=(NRuns-Burnin)/Thinning,ncol=(2+length(theta)))

anetwork=asd<-0
tuning <- c(0.5, 0.5)
a <- 0

progressbar <- txtProgressBar(min = 0, max = NRuns, style = 3)
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

theta[3] <- runif(1, m.ll, m.ul)
theta[4] <- runif(1, v.ll, v.ul)
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
monteCarloN <- 300
delta_i     <- rgamma(monteCarloN, shape = AcceptedTheta[1]^2/AcceptedTheta[2], rate = AcceptedTheta[1]/AcceptedTheta[2])
delta_j     <- rgamma(monteCarloN, shape = AcceptedTheta[1]^2/AcceptedTheta[2], rate = AcceptedTheta[1]/AcceptedTheta[2])
Y           <- delta_i-delta_j  
Savetheta[a,] <- c(AcceptedTheta, AcceptedTheta[1]+mean(Y), (AcceptedTheta[2]+2*AcceptedTheta[4]))
SaveNetwork[,a] <- AcceptedNetwork
SaveP[a] <- AcceptedP
}
setTxtProgressBar(progressbar, b)
}
close(progressbar)

save.image("29FebDataFixedPriorReportedNetwork.RData")

# some trace plots
par(mfrow=c(2,1))
plot(Savetheta[,1], type="l", col=3, main="mean of generation interval distr.")
plot(Savetheta[,2], type="l", col=3, main="var of generation interval distr.")
plot(Savetheta[,3], type="l", col=3, main="mean of incubation period distr.")
plot(Savetheta[,4], type="l", col=3, main="var of incubation period distr.")

plot(SaveP, type="l")
hist(Savetheta[,1])
hist(Savetheta[,2])
plot(Savetheta[,5], type="l")
hist(Savetheta[,5])
plot(SaveNetwork[5,])
plot(SaveNetwork[6,])
plot(SaveNetwork[26,])
plot(SaveNetwork[34,])
plot(SaveNetwork[37,])
plot(SaveNetwork[39,])


# serial interval distribution parameter estimates
mean(Savetheta[,1], na.rm = T); quantile(Savetheta[,1], c(0.025, 0.5, 0.975))#; sd(Savetheta[,1], na.rm = T)
mean(Savetheta[,2]^0.5, na.rm = T); quantile(Savetheta[,2]^0.5, c(0.025, 0.5, 0.975))#; sd(Savetheta[,2], na.rm = T)
mean(Savetheta[,5], na.rm = T); quantile(Savetheta[,5], c(0.025, 0.5, 0.975))#; sd(Savetheta[,5], na.rm = T)
mean(Savetheta[,6]^0.5, na.rm = T); quantile(Savetheta[,6]^0.5, c(0.025, 0.5, 0.975))#; sd(Savetheta[,5], na.rm = T)

