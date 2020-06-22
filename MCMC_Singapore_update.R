
rm(list = ls())

Cluster <- c(rep(1,8), rep(5,23), rep(2,8), rep(3,3), rep(4,5), rep(6,4), rep(7,3))
Case.ID <- c(8,9,31,33,38,83,90,91,66,68,70,71,80,84,88,48,49,51,54,57,58,60,53,61,62,63,67,73,74,78,81,19,20,21,24,25,27,34,40,30,36,39,42,47,52,56,69,2,13,26,44,50,55,77)
Case <- c(1:length(Case.ID))
Time <- c(4,4,3,10,14,8,20,3,9,10,14,12,15,15,27,12,14,15,21,22,21,19,21,17,20,21,20,20,23,20,27,9,5,13,10,4,12,7,10,1,4,9,12,17,18,23,22,1,8,8,11,18,10,21)
vi <- c(0,0,rep(NA,6), rep(NA,18),24,18,rep(NA,3),0,0,32,0,0,32,rep(NA,8),44,NA,0,0,0,NA,NA,0,52)

data <- data.frame(Cluster,Case,Case.ID,Time,vi)

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

#vi.list[[40]] <- c(41,42,0)
vi.list[[40]] <- c(0) # for SI > -3
vi.list[[41]] <- c(40,42,0)
vi.list[[42]] <- c(40,41,0)

vi.list[[43]] <- c(44,45,0)
vi.list[[44]] <- c(43,45,46,0)
vi.list[[45]] <- c(43,44,0)
vi.list[[47]] <- c(43:46,0)

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
IsNotContributorToLikelorg <- c(which(vi==0)) # known index cases
IsContributorToLikelorg <- Case[!Case%in%IsNotContributorToLikelorg]

## Constrain neg SI to be max 3 days
x = -3

inf.mat = matrix(nrow=NCases, ncol=NCases)
for(i in 1:NCases){
  for(j in 1:NCases){
    inf.mat[i,j] = ifelse(j%in%PossibleInfector[i,],1,0)
  }
}
# all possible SI
time.mat = matrix(nrow=NCases, ncol=NCases)
for(i in 1:NCases){
  for(j in 1:NCases){
    time.mat[i,j] = Time[i] - Time[j]
  }
}
SI.mat = inf.mat*time.mat # min(SI.mat) = -17 (lowest SI)
# weigh by number of infectors
w.SI <- NULL
for(i in IsContributorToLikelorg){
  w.SI[i] = sum(SI.mat[i,]) / NPossibleInfector[i]
}

# Select SI >= x
for(i in 1:NCases){
  for(j in 1:NCases){
    inf.mat[i,j]=ifelse(time.mat[i,j] < x, 0, inf.mat[i,j])
  }
}
SI.mat2 = inf.mat*time.mat
# weigh by number of infectors
w.SI2 <- NULL
for(i in IsContributorToLikelorg){
  w.SI2[i] = sum(SI.mat2[i,]) / sum(inf.mat[i,])
}

# change vi.list
vi.list2 <- rep(list(NULL),NCases)
#k=1
for(i in Case){
  for(j in Case){
    if(inf.mat[i,j]==1) vi.list2[[i]] <- c(vi.list2[[i]],Case[j])
    # k=k+1
  }
}

for(i in 1:NCases){
  if(is.null(vi.list2[[i]])) vi.list2[[i]]=0
}

PossibleInfector <- matrix(nrow=NCases,ncol=max(lengths(vi.list2)))
for(i in 1:NCases){
  PossibleInfector[i,1:lengths(vi.list2)[i]] <- vi.list2[[i]]
}
NPossibleInfector    <- rowSums(!is.na(PossibleInfector))

IsNotContributorToLikelorg <- c(which(vi==0 | Case==40 | Case==43)) # when allowing SI >= -3
IsContributorToLikelorg <- Case[!Case%in%IsNotContributorToLikelorg]


library(igraph)
# Function to find cycles in network
FindCycles = function(g) {
  Cycles = NULL
  for(v1 in V(g)) {
    if(degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}

##############-----MCMC-----------------------
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
#----MCMC---------------------
#--------------------------------------------------
Network <- numeric(NCases)+0
Update <- IsContributorToLikelorg

# Resample network until there are no cylces
n.cycles=1
while(n.cycles>0){
  Draw <- round(runif(length(Update),min=0.5,max=NPossibleInfector[Update]+0.5))
  for(i in 1:length(IsContributorToLikelorg)){
    Network[Update[i]] <- PossibleInfector[Update[i],Draw[i]]
  }
  AcceptedNetwork <- Network
  Infector=Network[IsContributorToLikelorg]
  Infectee=Case[IsContributorToLikelorg]
  
  tree=cbind(Infector,Infectee); tree=tree[which(Infector!=0),]
  g = graph_from_edgelist(tree)
  n.cycles = length(FindCycles(g))
}

IsContributorToLikel <- Case[Network!=0]
IsNotContributorToLikel <- Case[!Case%in%IsContributorToLikel]

AcceptedTheta=theta <- c(1, 1, 5.2, 2.8^2) # incubation: mean=5.2 days & sd=2.8 days

SerialInterval   <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]] # serial interval

P <- posterior()
AcceptedP <- P

NRuns <- 3000000
NUpdate <- length(IsContributorToLikelorg)
Burnin <- 500000
Thinning<-200
SaveP <- numeric()
SaveNetwork <- matrix(nrow=NCases,ncol=(NRuns-Burnin)/Thinning)
Savetheta <- matrix(nrow=(NRuns-Burnin)/Thinning,ncol=(2+length(theta))) # theta + var SI + mean SI

anetwork=asd<-0
tuning <- c(0.5, 0.5)
a <- 0
sd = 0.5


for(b in 1:NRuns){
  
  if(b%%100 == 0){
    theta <- AcceptedTheta
    Update <- IsContributorToLikelorg
    
    n.cycles=1
    while(n.cycles>0){
      Draw <- round(runif(length(Update),min=0.5,max=NPossibleInfector[Update]+0.5))
      for(i in 1:length(IsContributorToLikelorg)){
        Network[Update[i]] <- PossibleInfector[Update[i],Draw[i]]
      }
      AcceptedNetwork <- Network
      Infector=Network[IsContributorToLikelorg]
      Infectee=Case[IsContributorToLikelorg]
      
      tree=cbind(Infector,Infectee); tree=tree[which(Infector!=0),]
      g = graph_from_edgelist(tree)
      n.cycles = length(FindCycles(g))
    }
    
    IsContributorToLikel <- Case[Network!=0]
    IsNotContributorToLikel <- Case[!Case%in%IsContributorToLikel]
    
    meanSI = mean(Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]])
  }
  
  
  if(b%%100 != 0){
    Network <- AcceptedNetwork
    theta[1] <- runif(1, (AcceptedTheta[1]-tuning[1]), (AcceptedTheta[1]+tuning[1]))
    if(theta[1]<0){theta[1] <- AcceptedTheta[1]}

    theta[2] <- runif(1, (AcceptedTheta[2]-tuning[2]), (AcceptedTheta[2]+tuning[2]))
    if(theta[2]<0){theta[2] <- AcceptedTheta[2]}
    
  }
  
  P <- posterior()
  
  AcceptYN <- runif(1,min=0,max=1) <= exp(P-AcceptedP)
  if(AcceptYN){
    if(b%%100 == 0){
      anetwork <- anetwork + 1
      AcceptedNetwork <- Network
    }
    if(b%%100 != 0){
      asd <- asd+1
      AcceptedTheta <- theta
    }
    AcceptedP <- P
  }
  if(b%%Thinning == 0 & b>Burnin){
    a <- a + 1
    Savetheta[a,] <- c(AcceptedTheta, (AcceptedTheta[2]+2*AcceptedTheta[4]), meanSI)
    SaveNetwork[,a] <- AcceptedNetwork
    SaveP[a] <- AcceptedP
  }
}




