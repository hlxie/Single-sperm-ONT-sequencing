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
   
   
   
   
   
    # Traceback
viterbiPath = rep(NA,nObservations)
for(state in hmm$States)
  {
    if(max(v[,nObservations])==v[state,nObservations])
        {
            iterbiPath[nObservations] = state
            break
 
        }
 
  }
for(k in (nObservations-1):1)
 {
    for(state in hmm$States)
        {
            if(max(v[,k]+log(hmm$transProbs[[k]][,viterbiPath[k+1]]))==v[state,k]+log(hmm$transProbs[[k]][state,viterbiPath[k+1]]))
                {
                    viterbiPath[k] = state
                    break
                }
  
        }
  
  }
  return(viterbiPath)

}

initHMM_our = function(States, Symbols, startProbs=NULL, transProbs=NULL, emissionProbs=NULL)
{
  nStates    = length(States)
  nSymbols   = length(Symbols)
  S          = rep(1/nStates,nStates)
  T          = 0.5*diag(nStates) + array(0.5/(nStates),c(nStates,nStates))
  E          = array(1/(nSymbols),c(nStates,nSymbols))
  names(S)   = States
  dimnames(T)= list(from=States,to=States)
  dimnames(E)= list(states=States,symbols=Symbols)
  if(!is.null(startProbs)){S[]  = startProbs[]}
  if(!is.null(emissionProbs)){E[,] = emissionProbs[,]}
  return(list(States=States,Symbols=Symbols,startProbs=S,transProbs=transProbs,emissionProbs=E))
  }
                     
                 
flipFun <- function(v){
  v2 <- ifelse(v==7,7,ifelse(v==0, 1, 0))
  return (v2)
}

                     
transbasis=rbind(c(0,1),c(1,0))
gettrans_our=function(i)
{
  return(diag(2)*(1-rec[i])+transbasis*rec[i])
 
}
 
 
cell_10_SUPP=read.table("Phasing_SV.dat")
cell_10_COV=read.table("Phasing_SV_coverage.dat")
  
cell_10_COV=cell_10_COV[row.names(cell_10_SUPP),colnames(cell_10_SUPP)]
Finall_SUPP=cell_10_SUPP/cell_10_COV
Finall_SUPP[Finall_SUPP>1]=1


############################################################################################################
Chr=c()
Star=c()
End=c()
j=1
for( i in row.names(Finall_SUPP)  )
{
  Chr[j]=strsplit(i,'_')[[1]][1]
  Star[j]=as.numeric(strsplit(i,'_')[[1]][2])
  End[j]=as.numeric(strsplit(i,'_')[[1]][3])
  j=j+1
 
}
Finall_SUPP$Chr=Chr
Finall_SUPP$Star=Star
Finall_SUPP$End=End
 
Finall_SUPP=Finall_SUPP[order(Finall_SUPP$Chr,Finall_SUPP$Star,Finall_SUPP$End),]
 
Finall_SUPP$Chr=NULL
Finall_SUPP$Star=NULL
Finall_SUPP$End=NULL
 
##############################################################################################################
  
Phasing_sample_use=Finall_SUPP
 
j=1
Chr=c()
for( i in row.names(Phasing_sample_use) )
{
 
  Chr[j]=strsplit(i,'_')[[1]][1]
  j=j+1
 
}
Phasing_sample_use$Chr=Chr
 
for(Chr_num in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"))
{
  Phasing_sample_use_chr1=Phasing_sample_use[Phasing_sample_use$Chr==Chr_num,]
  gmtDa=Phasing_sample_use_chr1[,1:1573]
  gmtDa=hapiFrameSelection(gmtDa, 5)
  
  hmm = initHMM(States=c("S","D"), Symbols=c("s","d"),
        transProbs=matrix(c(0.99,0.0001,0.01,0.99),2),
        emissionProbs=matrix(c(0.99,0.0001,0.01,0.99),2),
        startProbs = c(0.5,0.5))

  data_tmp=as.data.frame(apply(gmtDa, 2, function(y) sum(!is.na(y))))
  data_tmp$F=data_tmp[,1]
  well_cell=row.names( data_tmp[order(data_tmp$F,decreasing = T),]  )[1:100]
  gmtDa=gmtDa[,well_cell]
  gmtDa=hapiFrameSelection(gmtDa, 3)
  Phasing_sample_use_chr1=gmtDa
 
for (nn in 1:5){
  nSNP <- nrow(gmtDa)
  genoError=list()
  Cell_Error=list()
  tmp2=c()
  for (i in 1:ncol(gmtDa)) {
      genoError[[i]] <- lapply(gmtDa[,-i], function(x)
      filterErrorFun(gmtDa[,i],x,Phasing_sample_use=Phasing_sample_use_chr1,hmm=hmm))
      tmp1=as.data.frame(table(unlist(genoError[[i]])))
      tmp1=tmp1[tmp1$Var1!=0,]
  
      tmp1=tmp1[tmp1$Freq>4,]
      tmp1$cell_sup=as.data.frame(apply(Phasing_sample_use_chr1, 1, function(y) sum(!is.na(y))))[as.numeric(as.charac    ter(tmp1$Var1)),]
      tmp1$ratio=tmp1$Freq/tmp1$cell_sup
      tmp1=tmp1[tmp1$ratio>0.5,]
    
      Cell_Error[[i]]=as.character(tmp1$Var1)
      tmp2=c(tmp2,as.character(tmp1$Var1))
 
   }
 
   need_ver=as.data.frame(table(tmp2))
   need_delet_location=as.numeric(as.character(need_ver[need_ver$Freq>5,]$tmp2))
   need_ver_location=need_ver[need_ver$Freq<3,]$tmp2
   need_ver_location=as.numeric(as.character(need_ver_location))
 
   for (i in 1:ncol(gmtDa)) {
        ver_num=intersect(as.numeric(Cell_Error[[i]]),need_ver_location)
        for(j in ver_num)
            {
               gmtDa[j,i]=flipFun(gmtDa[j,i])
            }
 
   }
  
   if(length(need_delet_location)==0)
        { gmtDa <- gmtDa  }else
        { gmtDa <- gmtDa[-need_delet_location,] }
                                        
 }

                                        
                                   
gmtFrame <- hapiFrameSelection(gmtDa, 3)
 
gmtFrame <- imputationFun1(gmtFrame, nSPT=2)
gmtFrame <- imputationFun1(gmtFrame, nSPT=2)
gmtFrame <- imputationFun1(gmtFrame, nSPT=2)
gmtFrame <- imputationFun1(gmtFrame, nSPT=2)
gmtFrame <- imputationFun1(gmtFrame, nSPT=2)
 
filter <- which(apply(gmtFrame, 1, function(y) sum(is.na(y)))>2)
 
gmtFrame <- gmtFrame[-filter,]
 
filter_cell=which(apply(gmtFrame, 2, function(y) sum(is.na(y)))>0)
 
gmtFrame <- gmtFrame[,-filter_cell]
 
draftHap <- hapiPhase(gmtFrame)
 
finalDraft <- hapiBlockMPR(draftHap, gmtFrame, cvlink = 2,smallBlock=2)
  
consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)
  
file1=paste("SV phasing/",Chr_num,"/Hap1.bed",sep = "")
file2=paste("SV phasing/",Chr_num,"/Hap2.bed",sep = "")
  
  
write.table(consensusHap[consensusHap$hap1==0 ,],file1,quote=F,sep = "\t",row.names = T)
write.table(consensusHap[consensusHap$hap1==1 ,],file2,quote=F,sep = "\t",row.names = T)

                                        
                                        
                                        
