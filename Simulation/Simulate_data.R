fn_simulation_input1<-function(n=2000,m=10,r=0.9,domain_effect=0,tau=0.2,tau1=0.5,tau2=0.1,tau3=0.05,tau4=0.01,sig_sq=0.25,kern_para=5,mean_exp=0,cov_st="CS"){
  
  if(!require(MASS)){
    stop("MASS package not installed")
  }else if(!require(Matrix)){
    stop("Matrix package not installed")
  }else if(!require(kernlab)){
    stop("kernlab package not installed")
  }else {
    if(cov_st=="CS"){
      R1=matrix(r,nrow=m,ncol=m)
      diag(R1)=1
    }else if(cov_st=="AR1") {
      R1<-toeplitz(sapply(1:m,function(i) r^(i-1)))
    }else{
      stop("Select AR1 or CS covariance structure")
    }
    d11=mvrnorm(n,mu=rep(mean_exp,m),Sigma = diag(rep(1,m)), tol = 0)
    
    sigma=R1+tau1*diag(rep(1,m))
    d1=mvrnorm(n,mu=rep(mean_exp,m),Sigma = sigma, tol = 0) #rows spots, column genes ##only coexpression(no noise)
    
    loc_df=read.csv("/Users/siktadasadhikari/Desktop/cSVG/Simulation/location_2000.csv")
    coord_df=data.frame(loc_df[1:n,2:3])
    coord_df=coord_df[order(coord_df$y),]
    rbf <- rbfdot(sigma = kern_para)
    K=kernelMatrix(kernel=rbf, x=as.matrix(coord_df), y = NULL)
    tau=tau
    K=Matrix::forceSymmetric(K)
    sig_sq=sig_sq
    V=tau*K+sig_sq*diag(1,n)
    sp=mvrnorm(1,mu=rep(mean_exp,n),Sigma = V, tol = 0) 
    co=seq(0.1,1,by=0.1)
    #co=seq(1.1,2,by=0.1)
    sigma1=tau2*diag(rep(1,m))+R1
    d2=mvrnorm(n,mu=rep(mean_exp,m),Sigma =  sigma1, tol = 0) #rows spots, column genes ##only coexpression(no noise)
    for(i in 1:m){
      d2[,i]=d2[,i]+(co[i]*sp)
    }       
    
    rbf <- rbfdot(sigma = 20)
    K=kernelMatrix(kernel=rbf, x=as.matrix(coord_df), y = NULL)
    tau=tau
    K=Matrix::forceSymmetric(K)
    sig_sq=sig_sq
    V=tau*K+sig_sq*diag(1,n)
    sp=mvrnorm(1,mu=rep(mean_exp,n),Sigma = V, tol = 0) 
    co=seq(0.1,1,by=0.1)
    sigma1=tau3*diag(rep(1,m))+R1
    d3=mvrnorm(n,mu=rep(mean_exp,m),Sigma =  sigma1, tol = 0) #rows spots, column genes ##only coexpression(no noise)
    for(i in 1:m){
      d3[,i]=d3[,i]+(co[i]*sp)
    }
    
    K=kernelMatrix(kernel=polydot(degree=2), x=as.matrix(coord_df), y = NULL)
    tau=tau
    K=Matrix::forceSymmetric(K)
    sig_sq=sig_sq
    V=tau*K+sig_sq*diag(1,n)
    sp=mvrnorm(1,mu=rep(mean_exp,n),Sigma = V, tol = 0) 
    co=seq(0.1,1,by=0.1)
    sigma1=tau4*diag(rep(1,m))+R1
    d4=mvrnorm(n,mu=rep(mean_exp,m),Sigma =  sigma1, tol = 0) #rows spots, column genes ##only coexpression(no noise)
    for(i in 1:m){
      d4[,i]=d4[,i]+(co[i]*sp)
    }
    d5=mvrnorm(n,mu=rep(0.05,3),Sigma = diag(rep(0.1,3)), tol = 0)
    d5[,1]=d5[,1]+c(rep(0,(10*n/20)),rep(2,(n/20)),rep(0,(9*n/20)))
    d5[,2]=d5[,2]+c(rep(0,(8*n/20)),rep(2,(n/20)),rep(0,(11*n/20)))
    d5[,3]=d5[,3]+c(rep(0,(12*n/20)),rep(2,(n/20)),rep(0,(7*n/20)))
    
    d=cbind(d11,d1,d2,d3,d4,d5)
    
    genes=paste0("indept",1:(5*m+3))
    genes[(m+1):(2*m)]=paste0("correlated",1:m)
    genes[(2*m+1):(3*m)]=paste0("spatial1_",1:m)
    genes[(3*m+1):(4*m)]=paste0("spatial2_",1:m)
    genes[(4*m+1):(5*m)]=paste0("spatial3_",1:m)
    genes[(5*m+1):(5*m+3)]=paste0("unique_",1:3)
    
    colnames(d)<-genes
    return(list(data_mat=t(d),loc_mat=as.matrix(coord_df)))
  }
}
