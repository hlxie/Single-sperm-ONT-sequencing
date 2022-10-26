
argv <- commandArgs(TRUE)
Snp_file <- argv[1]
sample <- argv[2]


 get_emission<-function(quam,noise=per_call_noise)
{
  ems<-do.call(rbind,lapply(1:nrow(quam),function(i){e=10^(quam[i,])+noise;return(e/sum(e))}))
  return(ems)
}

gettrans=function(i)
{
  return(diag(2)*(1-rec[i])+transbasis*rec[i])

}

viterbi=function()
{
  prob_state=matrix(nrow=nrow(emsmat),ncol=2)
  prev_state=matrix(nrow=nrow(emsmat),ncol=2)
  current_pos=1
  current_prob=c(0.5,0.5)
  current_vec=emsmat[1,]*current_prob
  prob_state[current_pos,]=current_vec
  prev_state[current_pos,]=NA
 
  for(current_pos in 1:(nrow(emsmat)-1))
     {
        for(s in 1:2)
        {
            new_vec=prob_state[current_pos,] * gettrans(current_pos)[s,]
            prob_state[current_pos+1,s] = max(new_vec) * emsmat[current_pos+1,s]
            prev_state[current_pos+1,s] = which.max(new_vec)
        }
        prob_state[current_pos+1,] = prob_state[current_pos+1,] / max(prob_state[current_pos+1,])
 
     }
   result = curr_best = which.max(prob_state[nrow(prob_state),])
   for(current_pos in nrow(prev_state):2)
        {
            curr_best=prev_state[current_pos,curr_best]
            result=c(curr_best,result)
        }
   result=result-1
   return(result)
  
    }
  
  data_raw=rbind(Snp_file)
  sample_crossover_site=list()
  w=1
  p=1
  sample_site=list()
  Chr_length=c("chr1"=195471971,"chr2"=182113224,"chr3"=160039680,"chr4"=156508116,"chr5"=151834684,"chr6"=149736546,"chr7"=145441459,"chr8"=129401213,"chr9"=124595110,"chr10"=130694993,"chr11"=122082543,"chr12"=120129022,"chr13"=120421639,    "chr14"=124902244,"chr15"=104043685,"chr16"=98207768,"chr17"=94987271,"chr18"=90702639,"chr19"=61431566)
 
  for (chr_ in chr_num)
   {
     flag=0
     data_raw_C = data_raw[data_raw$V1==chr_,]
     data_raw_C = data_raw_C[order(data_raw_C$V2),]
     if(dim(data_raw_C)[1]<100)
      { next }
     else{
       err=1
       baseline_noise=err/5;
       per_call_noise=err/5;
       num_gen=1
       corr_reads_b6=corr_reads_cast=0
       all_reads_b6=all_reads_cast=0
       corr_reads_target_test=0
       all_reads_target_test=0
       num_b6 = num_cast=0
 
       position=as.double(data_raw_C$V2)
 
       logll=strsplit(as.character(data_raw_C$V4),",",fixed=TRUE)
       logll=matrix(as.double(unlist(logll)), ncol = 3, byrow = TRUE)
       logll=logll[,c(1,3)]
 
       emsmat=get_emission(logll)
       rec=1-exp(-diff(position)*1e-9*num_gen)
 
       transbasis=rbind(c(0,1),c(1,0))
  
 ###############################################################
  
       v=viterbi()
  
 ###############################################################
     data_raw_C$HMM=v
     data_raw_C$diff=c(diff(v),0)
     data_raw_C$ID=1:dim(data_raw_C)[1]
 
 
     site=data_raw_C[data_raw_C$diff=="-1"|data_raw_C$diff=="1" ,]
 
     if (dim(site)[1]>0)
       {
        site_length=c()
        site_length[1]=data_raw_C[site[1,"ID"],"V2"]-data_raw_C[1,"V2"]
        site=rbind(site,tail(data_raw_C, n=1L))
 
        for (i in 2:dim(site)[1])
            {
       
              site_length[i]=data_raw_C[ site[i,"ID"],"V2" ] - data_raw_C[ site[i-1,"ID"]+1 ,"V2" ]
  
            }
  
        site$site_length=site_length
        site$ID_diff=c(   site$ID[1]  ,diff(site$ID) )
       }
  
     else if (dim(site)[1]==0 & Chr_length[chr_]/dim(data_raw_C)[1] <40000)
       {
        flag=1
 
       }
     else
       {
        next
       }
 
 
     finally_site = site[site$ID_diff>100 & site$site_length>50000,]
     sample_site[[p]]=site
     p=p+1
       
 ######################################################
     if( dim(finally_site)[1]>1 & flag==0   )
        {
            finally_site$Cross_site_diff=c(0,diff(finally_site$HMM))
            finally_site$Cross_site_end=finally_site$V2-finally_site$site_length
            finally_site$Cross_site_str=c(1,finally_site$V2[   1:(dim(finally_site)[1]-1)   ])
            finally_site$resolution=finally_site$Cross_site_end - finally_site$Cross_site_str
  
            tmp_site=finally_site[abs(finally_site$Cross_site_diff)==1,c("V1","Cross_site_str","Cross_site_end","resolution","HMM")]
                if ( dim(tmp_site)[1]>0   )
                    {
                        sample_crossover_site[[w]]=tmp_site
                        w=w+1

                    }
                else
                    {

                        TT=finally_site[1,"HMM"]
                        kong=as.data.frame(t(c(chr_,0,0,0,TT)))
                        colnames(kong)=c("V1","Cross_site_str","Cross_site_end","resolution","HMM")
 
                        kong$HMM=as.double(TT)
                        kong$Cross_site_str=as.double(kong$Cross_site_str)
                        kong$Cross_site_end=as.double(kong$Cross_site_end)
                        kong$resolution=as.double(kong$resolution)
                        sample_crossover_site[[w]]=kong
                        w=w+1
 
                    }
 
        }
 
        else if( dim(finally_site)[1]==1 )
        {
            if(as.numeric( as.numeric(finally_site$site_length)/Chr_length[chr_] )>0.7 & finally_site$ID_diff>500)
                {
                    TT=data_raw_C[1,"HMM"]
                    kong=as.data.frame(t(c(chr_,0,0,0,TT)))
                    colnames(kong)=c("V1","Cross_site_str","Cross_site_end","resolution","HMM")
                    kong$HMM=as.double(TT)
                    kong$Cross_site_str=as.double(kong$Cross_site_str)
                    kong$Cross_site_end=as.double(kong$Cross_site_end)
                    kong$resolution=as.double(kong$resolution)
                    sample_crossover_site[[w]]=kong
                    w=w+1
                 }
 
        }
 
 
        else if (flag==1)
        {
 
          TT=data_raw_C[1,"HMM"]
          kong=as.data.frame(t(c(chr_,0,0,0,TT)))
          colnames(kong)=c("V1","Cross_site_str","Cross_site_end","resolution","HMM")
          kong$HMM=as.double(TT)
          kong$Cross_site_str=as.double(kong$Cross_site_str)
          kong$Cross_site_end=as.double(kong$Cross_site_end)
          kong$resolution=as.double(kong$resolution)
          sample_crossover_site[[w]]=kong
          w=w+1
        }
  
     }
}

Sample_crossover_site_info=Reduce(function(x,y) merge(x,y,all=T),sample_crossover_site)
Sample_site_info=Reduce(function(x,y) merge(x,y,all=T),sample_site)

  
write.table(Sample_crossover_site_info, paste("./07.crossover_site/",sample,"_crossover_site.bed",sep=""),quote=F,row.names=F,sep="\t")

