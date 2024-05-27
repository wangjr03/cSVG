
rm(list = ls())
if(TRUE){
  args=(commandArgs(TRUE))
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}else{
  jobID=1
}

###################################################
fn_iteration<-function(iter){
library(SPARK)
  library(nnSVG)
  library(scater)
  library(scran)
  library(scry)
  library(scuttle)
  
  sample_i=6 #change it to any other sample
  load(paste0("dataset/LIBD_sample",sample_i,".RData"))
  xy_coords = as.matrix(xy_coords)
  rownames(xy_coords) = colnames(count_sub)
  
  count_df=as.matrix(count_sub) 
  loc_df=xy_coords
  
  #filtering based on spark criteria
  keep_cols=which(apply(count_df,2,sum)>10)
  percent=0.1 
  keep_rows=which( rowSums(count_df > 0) >= floor(percent*ncol(count_df)) )
  count_df=count_df[keep_rows,keep_cols]
  loc_df=loc_df[keep_cols,]
  dim(count_df) 
  dim(loc_df) 
  
  count_df=as.matrix(count_df)
  loc_df=as.matrix(loc_df)
  count_df[1:5,1:5]
  norm_df=scuttle::normalizeCounts(count_df,transform="log")
  norm_df[1:5,1:5]
  
  rm(count_df)
  
  exp_f=norm_df
  coord_df=loc_df
  genes=rownames(exp_f)
  print(paste0("dim of exp_f is:",dim(exp_f)[1]))
  print(dim(coord_df))
  

  st=(iter-1)*10+1
  en=iter*10
  if(en>(dim(exp_f)[1])){
    en=dim(exp_f)[1]
  }
genes=rownames(norm_df)
  
source("fn_main_chunk.R")
final=fn_cSVG(data_mat=exp_f,loc_mat=coord_df,method_step1="MargcorTest",thres_step1="standard",control=FALSE,st=st,en=en)
  
write.csv(final,paste0("dataset/out_step1_sample",sample_i,"_",iter,".csv"))
}
ans = fn_iteration(iter=jobID)
quit(save="no")

