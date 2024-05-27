
fn_cSVG<-function(data_mat,loc_mat,method_step1="MargcorTest",thres_step1="standard",control=FALSE,st=1,en=10){

if (!require(SKAT)) {
    
    stop("SKAT package not installed")
    
  }else if(!require(MASS)){
    
    stop("MASS package not installed")
  }else if(!require(Matrix)){
    
    stop("Matrix package not installed")
  }else {
  }

#SKAT code #Some functions were taken from SKAT Github
Get_Davies_PVal<-function(Q, W, Q.resampling = NULL){
  
  K<-W/2
  
  Q.all<-c(Q,Q.resampling)
  
  re<-Get_PValue(K,Q.all)
  param<-list()
  param$liu_pval<-re$p.val.liu[1]
  param$Is_Converged<-re$is_converge[1]
  
  
  p.value.resampling = NULL
  if(length(Q.resampling) > 0){
    p.value.resampling<-re$p.value[-1]
    param$liu_pval.resampling<-re$p.val.liu[-1]
    param$Is_Converged.resampling<-re$is_converge[-1]
    
  }
  
  
  re<-list(p.value = re$p.value[1], param=param,p.value.resampling = p.value.resampling
           , pval.zero.msg=re$pval.zero.msg )  
  return(re)
}

###############################################################
Get_PValue<-function(K,Q){
  
  lambda<-Get_Lambda(K)
  re<-Get_PValue.Lambda(lambda,Q)
  return(re)
}

##############################################################
#only focus on positive eigenvalue
Get_Lambda<-function(K){
  
  out.s<-eigen(K,symmetric=TRUE, only.values = TRUE)
  #print(out.s$values)
  
  #out.s1<-eigen(K,symmetric=TRUE)
  #print(out.s1$values)
  
  lambda1<-out.s$values
  IDX1<-which(lambda1 >= 0)
  
  # eigenvalue bigger than sum(eigenvalues)/1000
  IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
  #cat("Lambda:", lambda1, "\n")
  #K1<<-K
  
  if(length(IDX2) == 0){
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda<-lambda1[IDX2]
  return(lambda)
  
}
################################################################
Get_PValue.Lambda<-function(lambda,Q){
  
  #print(lambda)
  n1<-length(Q)
  
  p.val<-rep(0,n1)
  p.val.liu<-rep(0,n1)
  is_converge<-rep(0,n1)
  p.val.liu<-Get_Liu_PVal.MOD.Lambda(Q, lambda)
  
  for(i in 1:n1){
    out<-SKAT_davies(Q[i],lambda,acc=10^(-6))
    
    p.val[i]<-out$Qq
    #p.val.liu[i]<-SKAT_liu(Q[i],lambda)
    
    is_converge[i]<-1
    
    # check convergence
    if(length(lambda) == 1){
      p.val[i]<-p.val.liu[i]
    } else if(out$ifault != 0){
      is_converge[i]<-0
    }
    
    # check p-value
    if(p.val[i] > 1 || p.val[i] <= 0 ){
      is_converge[i]<-0
      p.val[i]<-p.val.liu[i]
    }
  }
  
  p.val.msg = NULL
  p.val.log=NULL
  #cat(p.val[1])
  if(p.val[1] == 0){
    
    param<-Get_Liu_Params_Mod_Lambda(lambda)
    p.val.msg<-Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ, param$muX, param$sigmaQ, param$sigmaX, param$l, param$d)
    p.val.log<-Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p=TRUE)[1]
    
  }
  
  return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge, p.val.log=p.val.log, pval.zero.msg=p.val.msg))
  
}
########################################################
Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  
  return(p.value)
  
}
##############################################################
Get_Liu_PVal.MOD.Lambda.Zero<-function(Q, muQ, muX, sigmaQ, sigmaX, l, d){
  
  
  Q.Norm<-(Q - muQ)/sigmaQ
  Q.Norm1<-Q.Norm * sigmaX + muX
  
  temp<-c(0.05,10^-10, 10^-20,10^-30,10^-40,10^-50, 10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
  #qchisq(temp, df=1000000000,lower.tail=FALSE)	
  out<-qchisq(temp,df = l,ncp=d, lower.tail=FALSE)
  #cat(c(Q.Norm1,l,d, out))
  #cat("\n")
  IDX<-max(which(out < Q.Norm1))
  
  pval.msg<-sprintf("Pvalue < %e", temp[IDX])
  return(pval.msg)
  
}

##################################################################
SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  
  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}

#######################################
Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation
  
  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }
  
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2
  
  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0
  
  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a
  
  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}


SKAT.linear.Other = function(res,Z,X1, kernel, weights = NULL, s2, method,res.out,n.Resampling){

  
  n<-nrow(Z)
  m = ncol(Z) 
  if (is.matrix(kernel)) {

    K = kernel

  } else {

    #K = lskmTest.GetKernel(Z, kernel, weights,n,m)
    print("the kernel is not matrix")
  }


  Q = t(res)%*%K%*%res/(2*s2)

  Q.res = NULL
  if(n.Resampling > 0){
	Q.res<-rep(0,n.Resampling)
	for(i in 1:n.Resampling){

  		Q.res[i] = t(res.out[,i])%*%K%*%res.out[,i]/(2*s2)
  	}
  }


  W = K - X1%*%solve( t(X1)%*%X1)%*%( t(X1) %*% K)	# W = P0 K

  if(method == "davies"){
  	# P0_half = P0
	W1 = W - (W %*% X1) %*%solve( t(X1)%*%X1)%*% t(X1)
  } 


  if( method == "liu" ){

	out<-Get_Liu_PVal(Q, W,Q.res)    
	pval.zero.msg=NULL

  } else if( method == "liu.mod" ){

	out<-Get_Liu_PVal.MOD(Q, W,Q.res)   
	pval.zero.msg = NULL 

  } else if( method == "davies" ){

	out<-Get_Davies_PVal(Q, W1,Q.res)  
	pval.zero.msg = out$pval.zero.msg  

  } else {
	stop("Invalid Method!")
  }

  

  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling
  , Test.Type = method, Q = Q, param=out$param, pval.zero.msg=pval.zero.msg )  

  return(re)
}

#screening code # Some function were taken from Github: https://github.com/wwrechard/screening

.rankScreening <- function(X, Y) {
  n = dim(X)[1]
  p = dim(X)[2]
  w = rep(0,p)
  Ynew = sort(Y, index.return = T)
  for(j in 1 : n)
      for(k in j : n)
          w = w + (X[Ynew$ix[k], ] > X[Ynew$ix[j], ])
  w = w / n / (n-1)
  w = w - 1 / 4
  return (abs(w))
}

#############################################################################

.forwardRegression <- function(X, Y, num.select, family, ebic.gamma = 1) {

    # get the dimension
    n = dim(X)[1]
    p = dim(X)[2]

    # store the used variables, including good variable and bad variable
    usedVariables = NULL

    # store the selected variables
    selectedVariables = NULL

    # store the residual sum of squares
    rss = rep(0, p)

    # store the bic values
    bic = rep(Inf, num.select)

    iteration = 0
    while (iteration < num.select) {

        # to compute for each variable the deviance of the model if they were added
        for (i in setdiff(1 : p, usedVariables)) {
            activeVariables = c(selectedVariables, i)
            model = try(glm(Y ~ X[, activeVariables] - 1, family = family))
            if (inherits(model, "try-error")) {
                rss[i] = Inf
                usedVariables = c(usedVariables, i)
            } else {
                rss[i] = model$deviance
            }
        }
        if (min(rss) == Inf) {
            break
        }

        # select the variabel that gives the smallest deviance to the model
        toAdd = which.min(rss)
        selectedVariables = c(selectedVariables, toAdd)
        usedVariables = c(usedVariables, toAdd)

        # record the corresponding bic value
        iteration = iteration + 1
        bic[iteration] = .ebic(rss[toAdd], p, n, iteration, ebic.gamma)
        rss[toAdd] = Inf
    }

    return (list(select = selectedVariables, bic = bic))
}

##############################################################################

.ebic <- function(deviance, model.size, sample.size, num.select, ebic.gamma) {
    return (deviance + num.select * (log(sample.size) + 2 * ebic.gamma * log(model.size)))
}

##############################################################################

.ebicRanking <- function(X, Y, sortedVariables, family, ebic.gamma) {
    # get the dimension
    n = dim(X)[1]
    p = dim(X)[2]

    # store the currently selected variables
    selectedVariables = NULL

    # store the bic values
    bic = rep(Inf, n - 1)

    iteration = 0
    while (iteration < n - 1) {

        iteration = iteration + 1

        # to compute for each variable the deviance of the model if they were added
        i = sortedVariables[iteration]
        selectedVariables = c(selectedVariables, i)
        model = try(glm(Y ~ X[, selectedVariables] - 1, family = family))
        if (inherits(model, "try-error")) {
            rss = Inf
        } else {
            rss = model$deviance
        }

        if (rss == Inf) {
            break
        }

        # record the corresponding bic value
        bic[iteration] = .ebic(rss, p, n, iteration, ebic.gamma)
    }

    bestModel = which.min(bic)
    return (list(select = sortedVariables[1 : bestModel], bic = bic))
}



screening <- function(x, y, method = 'holp', num.select, family = 'gaussian', ebic = FALSE, ebic.gamma = 1){
#floor(dim(x)[1]/2)
    #if(num.select=="standard"){
    #  num.select=floor(dim(x)[1]/2)
    #}else{
    #  num.select=num.select
    #}  #this if-else part added by SDA #not needed.
    # standardize
    x = as.matrix(x)
    X = scale(x)
    if (family == 'gaussian'){
        Y = y - mean(y)
    }
    else {
        Y = y
    }
    if (is.null(dim(X))) {
        p = 1
        n = length(X)
    } else {
        n = dim(X)[1]
        p = dim(X)[2]
    }

    # if p is smaller than the required number of variables, return all
    if (p == 1 || (p < num.select && !ebic)) {
        selectedVariable = 1 : p
    }
    else {
        # for the linear case, it is easy to compute everything.
        if (family == 'gaussian') {
            if (method == 'holp') {
                OLS = t(X) %*% solve(X %*% t(X) + diag(n) * 1, Y)
                ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'sis') {
                OLS = t(X) %*% Y
                ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'rrcs') {
                OLS = .rankScreening(X, Y)
                ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'forward') {
                if (ebic) {
                    ranking = .forwardRegression(X, Y, n - 1, family, ebic.gamma)
                    bestModel = which.min(ranking$bic)
                    result = ranking$select[1 : bestModel]
                }
                else {
                    result = .forwardRegression(X, Y, num.select, family, ebic.gamma)$select
                }
            }
            else {
                stop('The method is unknown and not supported/')
            }
        }
        else {
            if (method == 'holp') {
                require(glmnet)
                model = glmnet(x = X, y = Y, family = family, alpha = 0, lambda = 1, intercept = FALSE)
                coefs = coef(model)[-1]
                ranking = sort(abs(coefs), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'sis') {
                residuals = rep(0, p)
                for (i in 1 : p) {
                    model = glm(Y ~ X[, i] - 1, family = family)
                    residuals[i] = model$deviance
                }
                ranking = sort(residuals, index.return = TRUE)
                if (ebic) {
                    results = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'forward') {
                if (ebic) {
                    ranking = .forwardRegression(X, Y, n - 1, family, ebic.gamma)
                    bestModel = which.min(ranking$bic)
                    result = ranking$select[1 : bestModel]
                }
                else {
                    result = .forwardRegression(X, Y, num.select, family, ebic.gamma)$select
                }
            }
            else if (method == 'rrcs') {
                stop('rrcs is not supported for GLM')
            }
            else {
                stop('The method is unknown and not supported.')
            }
        }
    }
    return (list(screen = result, method = method))
}




 ##Some code were taken from SPARK Github :https://github.com/xzhoulab/SPARK

#' Calculate Gaussian Parameter L based on Eculidean Distance
#' @param X Cell corrdinates matrix n x 2 or kernel matrix computed already
#' @param compute_distance Compute the distance matrix using generic function dist, default=TRUE
#' @export
ComputeGaussianPL <- function(X, compute_distance=TRUE){
  if(compute_distance){
    if(ncol(X)<2){stop("X has to be a coordinate matrix with number of column greater than 1")}
    D <- dist(X)
  }else{
    D <- X
  }# end fi
  
  #Dval <- unique(as.vector(D))
  Dval <- D
  Dval.nz <- Dval[Dval>1e-8]
  lmin <- min(Dval.nz)/2
  lmax <- max(Dval.nz)*2
  lrang <- 10^(seq(log10(lmin),log10(lmax),length.out=10))
  return(lrang)
}# end func



#' Summarized the multiple p-values via Cauchy combination rule
#' @param pvalues Pvalues from multiple kernels, a px10 matrix
#' @param weights The weights for combining the pvalues from multiple kernels, a px10 matrix or NULL
#' @export
CombinePValues <- function(pvalues, weights=NULL){
	if(!is.matrix(pvalues)){pvalues <- as.matrix(pvalues)}
	## to avoid extremely values
	pvalues[which(pvalues==0)] <- 5.55e-17
	pvalues[which((1-pvalues)<1e-3)] <- 0.99
	
	num_pval <- ncol(pvalues)
	num_gene <- nrow(pvalues)
	if(is.null(weights)){
		weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
	}# end fi
	if( (nrow(weights) != num_gene) || (ncol(weights) != num_pval)){
		stop("the dimensions of weights does not match that of combined pvalues")
	}# end fi
	
	Cstat <- tan((0.5 - pvalues)*pi)

	wCstat <- weights*Cstat
	Cbar <- apply(wCstat, 1, sum)
	#combined_pval <- 1.0/2.0 - atan(Cbar)/pi
	combined_pval <- 1.0 - pcauchy(Cbar)	
	combined_pval[which(combined_pval <= 0)] <- 5.55e-17
	return(combined_pval)
}# end func

####################################################################################
#Preparation

exp_f=data_mat
coord_df=loc_mat
genes=rownames(exp_f)
print("exp_f is:")
print(dim(exp_f))
print(dim(coord_df))


################

res=matrix(0,nrow=1,ncol=10)
#colnames(res)=c("combinedPval", "adjustedPval")
Z=as.matrix(coord_df)
ED <- as.matrix(dist(coord_df))
lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
l_g=list()
for(g in genes[st:en]){
    ind=which(rownames(exp_f)==g)
    y=exp_f[ind,]
    if(control!=TRUE){
    obj1<-SKAT_Null_Model(y~1, out_type = 'C')
    }else{
    dep=exp_f[-ind,]
    skip_pc=0
    if(thres_step1=="standard"){
      thres_step1=floor((dim(dep)[1])/2)
    }else{
      thres_step1=as.numeric(thres_step1)
    } 
    print("thres-step1 is:")
    print(thres_step1)
    if(method_step1=="SIS"){
    #library(devtools)
    #install_github('wwrechard/screening')
    #library(screening)
    out = screening(x=t(dep), y=as.numeric(exp_f[ind,]), method = "sis", num.select = thres_step1, family = "gaussian",  ebic=FALSE,ebic.gamma = 1)
    data_mat1= dep[out$screen,]
    }else if(method_step1=="rrcs"){
    out = screening(x=t(dep), y=as.numeric(exp_f[ind,]), method = 'rrcs', family = "gaussian", num.select = thres_step1, ebic=FALSE,ebic.gamma = 1)
    data_mat1= dep[out$screen,]
    }else if(method_step1=="MargcorTest"){
      fn_corTest<-function(x){
        model <- lm(y~x)
        #return(summary(model)$coefficients[2,4])
        if(summary(model)$coefficients[2,1]>0){   #added extra step to restrict negatively correlated genes
        return(summary(model)$coefficients[2,4])
        }else{
          return(1)
        }
      }
      res1=apply(dep,1,fn_corTest)
      res1=p.adjust(res1, method = "BY")
      ind1=which(res1<0.01)
      l_g[[ind]]= which(genes %in% names(ind1))
      if(length(ind1)>3){
        data_mat1=dep[ind1,] 
      }else{
        skip_pc=1
      }
      }else if(method_step1=="Enet"){
        require(glmnet)
        cvfit <- cv.glmnet(x=t(dep), y=y, type.measure="mse", alpha=0.2, standardize.response = TRUE)
        coeff1=coef(cvfit)[-1]
        ind1=which(coeff1>0)
        if(length(ind1)>3){
        data_mat1=dep[ind1,] 
        }else{
        skip_pc=1
      }
      }else if(method_step1=="SIS+Enet"){
        require(glmnet)
        out = screening(x=t(dep), y=as.numeric(exp_f[ind,]), method = "sis", num.select = thres_step1, family = "gaussian",  ebic=FALSE,ebic.gamma = 1)
        data_mat= dep[out$screen,]
        cvfit <- cv.glmnet(x=t(data_mat), y=y, type.measure="mse", alpha=0.2, standardize.response = TRUE)
        coeff1=coef(cvfit)[-1]
        ind1=which(coeff1>0)
        if(length(ind1)>3){
        data_mat1=data_mat[ind1,] 
        }else{
        skip_pc=1
      }
      }else{
      stop("Please specify the method used in step 1!")
      }
    
    if(skip_pc==0){
    data_mat1=t(scale(t(data_mat1))) #scales accross genes
    pc1 <- prcomp(data_mat1,
             center = FALSE,
             scale. = FALSE)
    #print(summary(pc1))
    
    prop_var=pc1$sdev^2 / sum(pc1$sdev^2)
    tot_var=cumsum(prop_var)
    n=min(2,which(tot_var>0.8)[1]) 
    pc_data1=pc1$rotation[,1:n]
    X=as.matrix(pc_data1)
  
    obj1<-SKAT_Null_Model(y~X, out_type = 'C')
    }else if(length(ind1)>1 & length(ind1)<=3){
      data_mat1=t(dep[ind1,]) 
      X=as.matrix(data_mat1)
      obj1<-SKAT_Null_Model(y~X, out_type = 'C')
    }else if(length(ind1)==1){
      X=as.matrix(dep[ind1,],nrow=length(y)) 
      obj1<-SKAT_Null_Model(y~X, out_type = 'C')
    }else{
      obj1<-SKAT_Null_Model(y~1, out_type = 'C')
    }
  }  
  p_val_all=c()
  for(ikernel in c(1:5) ){
			# Gaussian kernel
			cat(paste0("## testing Gaussian kernel: ",ikernel,"...\n"))
			kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))
      SKAT1<-SKAT.linear.Other(obj1$res,Z,obj1$X1,kernel = kernel_mat,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
      p_val_all=c(p_val_all,SKAT1$p.value)
			rm(kernel_mat)

			# Periodic kernel
			cat(paste0("## testing Periodic kernel: ",ikernel,"...\n"))
			kernel_mat <- cos(2*pi*ED/lrang[ikernel])
      SKAT1<-SKAT.linear.Other(obj1$res,Z,obj1$X1,kernel = kernel_mat,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
      p_val_all=c(p_val_all,SKAT1$p.value)
			rm(kernel_mat)
		}

      res=rbind(res,p_val_all)
      print(ind)
      dim(res)
}
l_g=l_g[st:en]
p_val_mat=res[-1,]
rownames(p_val_mat)=genes[st:en]

comb_pval=CombinePValues(p_val_mat)
final=cbind(p_val_mat,comb_pval)
colnames(final) <- c(paste0(c("GSP","COS"), rep(1:5,each=2)),"combined")

if(control==TRUE & method_step1=="MargcorTest"){
  return(list(final=final,list_g=l_g))
}else{
return(final)
}
} #end fn_cSVG

