library(HMM)
library(Hapi)


imputationFun1 <- function(pollens, nSPT=2) {
snpName <- rownames(pollens)
polName <- colnames(pollens)
polNum <- ncol(pollens)
snpNumTotal <- nrow(pollens)
impute.final <- matrix(rep(NA,snpNumTotal*polNum), snpNumTotal, polNum)
minNA <- 10000000
newNA <- sum(apply(pollens, 1, function(y) sum(is.na(y)))>=1)

pols <- 1:polNum
while (newNA < minNA) {
    minNA <- newNA
    for (k in pols) {
        pol.k <- pollens[,k]
        nonNAMatrix <- naAdjacentFun(pol.k)
        impute.final[,k] <- pol.k
  
        if (nrow(nonNAMatrix) == 0) {
                next
            } else {
        for (i in 1:nrow(nonNAMatrix)) {
            l <- nonNAMatrix[i,1]
            r <- nonNAMatrix[i,2]
            polPart <- pollens[l:r,]
            ll <- 1
            rr <- nrow(polPart)
            snpNum <- rr-ll+1
            polPart.k <- polPart[,k]
            impute.k <- matrix(rep(NA, snpNum*polNum), snpNum, polNum)
            impute.k[,k] <- polPart[,k]
 
            for (j in pols[pols!=k]) {
                polPart.j <- polPart[,j]
                if (is.na(polPart.j[ll]) | is.na(polPart.j[rr])) {
                                   next
                    } else if (polPart.k[ll] == polPart.j[ll] &
                               polPart.k[rr] == polPart.j[rr]) {
                               impute.k[ll:rr,j] <- polPart.j[ll:rr]
 
                    } else if (polPart.k[ll] == flipFun(polPart.j[ll]) &
                               polPart.k[rr] == flipFun(polPart.j[rr])) {
                               impute.k[ll:rr,j] <- flipFun(polPart.j[ll:rr])
  
                    } else {
                            next
  
                    }
            }
        
        impute.k <- impute.k[-c(1,nrow(impute.k)),]
        impute.k <- matrix(impute.k, ncol=polNum)
           if (is.matrix(impute.k)) {
                impute.final[(l+1):(r-1),k] <-apply(impute.k, 1, function(x) imptDetermineFun(x,nSPT=nSPT))
                impute.final[(l+1):(r-1),k][impute.final[(l+1):(r-1),k]==7] <- NA
                }else {
                impute.final[(l+1):(r-1),k] <-imptDetermineFun(impute.k)
                impute.final[(l+1):(r-1),k][impute.final[(l+1):(r-1),k]==7] <- NA
                }
 
            }
        }
  }
  pollens <- impute.final
  newNA <- sum(apply(pollens, 1, function(y) sum(is.na(y)))>=1)
 }
  
 rownames(pollens) <- snpName
 colnames(pollens) <- polName
 return (pollens)
 }





imptDetermineFun <- function(v,nSPT=2) {
   if (sum(!is.na(v))==0) {
        return (7)
   }
 
   v <- v[which(!is.na(v))]
   if (length(v) >= nSPT & length(unique(v)) == 1) {
        return (unique(v))
   } else {
        return (7)
 
   }
 
 }







naAdjacentFun <- function(gmt) {
   nonNAPos <- which(!is.na(gmt))
   nonNAPosDiff <- diff(nonNAPos)
   nonNAPosLeft <- nonNAPos[which(nonNAPosDiff != 1)]
   nonNAPosRight <- nonNAPos[which(nonNAPosDiff != 1)+1]
   nonNAMatrix <- cbind(nonNAPosLeft, nonNAPosRight)
   return (nonNAMatrix)
}






filterErrorFun <- function(gmt1, gmt2, Phasing_sample_use, hmm=NULL) {
  dComp <- gmt1 == gmt2
  idKnownPos <- which(!is.na(idComp))
  if(length(idKnownPos)<10)
  {
       return (c(0))
  }else {
  
  idComp <- as.numeric(idComp[!is.na(idComp)])
  genoSymbol <- ifelse(idComp==0,hmm$Symbols[1],hmm$Symbols[2])
  Phasing_sample_use=Phasing_sample_use
  position=row.names(Phasing_sample_use)[idKnownPos]
  position_data=c()
  j=1
  for (i in position)
  {
     position_data[j]=mean( as.numeric(strsplit(i,"_")[[1]][2]),as.numeric(strsplit(i,"_")[[1]][3])  )
     j=j+1
 
  }
 
     rec=1-exp(-diff(position_data)*1e-8)
     transProbs=list()
     for (k in 1:length(rec))
        {
            transProbs[[k]]=diag(2)*(1-rec[k])+transbasis*rec[k]
            colnames(transProbs[[k]])=c("S","D")
            row.names(transProbs[[k]])=c("S","D")

        }
 
 hmm = initHMM_our(States=c("S","D"), Symbols=c("s","d"),
 emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2),
 transProbs= transProbs,
 startProbs = c(0.5,0.5))
 
 correctGeno <- viterbi_our(hmm, genoSymbol)
 correctGeno <- ifelse(correctGeno==hmm$States[1], 0, 1)
 genoError <- correctGeno-idComp
  
 genoError <- which(genoError != 0)
 return (idKnownPos[genoError])
 }
}





viterbi_our = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  v          = array(NA,c(nStates,nObservations))
  dimnames(v)= list(states=hmm$States,index=1:nObservations)
  # Init
 for(state in hmm$States)
   {
      v[state,1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,observation[1]])
 
   }
  # Iteration
  for(k in 2:nObservations)
   {
        for(state in hmm$States)
        {
            maxi = NULL
            for(previousState in hmm$States)
                {
                    temp = v[previousState,k-1] + log(hmm$transProbs[[k-1]][previousState,state])
                    maxi = max(maxi, temp)
  
                }
            v[state,k] = log(hmm$emissionProbs[state,observation[k]]) + maxi
  
        }

   }
   
   
   
   
   
   
   
   
     ¦   ¦ # Traceback
 27     viterbiPath = rep(NA,nObservations)
 26     ¦   for(state in hmm$States)
 25     ¦   {
 24     ¦   ¦   if(max(v[,nObservations])==v[state,nObservations])
 23     ¦   ¦   {
 22     ¦   ¦   ¦   viterbiPath[nObservations] = state
 21     ¦   ¦   ¦   break
 20     ¦   ¦   ¦   ¦   ¦   ¦
 19     ¦   ¦   }
 18     ¦   ¦
 17     ¦   }
 16     ¦   for(k in (nObservations-1):1)
 15     ¦   {
 14     ¦   ¦   for(state in hmm$States)
 13     ¦   ¦   {
 12     ¦   ¦   ¦   if(max(v[,k]+log(hmm$transProbs[[k]][,viterbiPath[k+1]]))==v[state,k]+log(hmm$transProbs[[k]][state,vit    erbiPath[k+1]]))
 11     ¦   ¦   ¦   {
 10     ¦   ¦   ¦   ¦   viterbiPath[k] = state
  9     ¦   ¦   ¦   ¦   break
  8     ¦   ¦   ¦   ¦   ¦   ¦   ¦   ¦   ¦   ¦
  7     ¦   ¦   ¦   }
  6     ¦   ¦   ¦   ¦   ¦
  5     ¦   ¦   }
  4     ¦   ¦   ¦
  3     ¦   }
  2     ¦   return(viterbiPath)
  1     ¦   ¦
247 }







