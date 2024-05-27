For larger datasets, the first two steps of cSVG are typically performed using a high-performance computing cluster by breaking the data into smaller gene chunks. This makes the code runs faster and smoother. Here is the description of the process: 

Step 1: run code code_Step1.R using the .sb file step1.sb.
Based on how small the gene chunks(10 genes in the example code) are and the total no of genes in the dataset, the number N_ARRAY is decided. N_ARRAY=ceiling(total no of genes in dataset/chunk size). In this example of sample 6, N_ARRAY=ceiling(4865/10)=487. 

combining results from the chunk output:

```{r}
sample_i=6
data=matrix(0,nrow=1,ncol=12)
colnames(data)=c("X","GSP1","COS1","GSP2","COS2","GSP3","COS3","GSP4","COS4","GSP5","COS5","combined")
for(iter in 1:N_ARRAY){
    r=read.csv(paste0("dataset/out_step1_sample",sample_i,"_",iter,".csv"))
    data=rbind(data,r)
}

data=data[-1,]
write.csv(data,paste0("dataset/Results/Sample_",sample_i,"_result_step1.csv")) 
```

Step 2:


