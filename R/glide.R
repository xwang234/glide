
glide <- function(data,genotype_columns=NULL,adjusting_covariate_columns=NULL,
                  outcome_column=NULL,coeff=NULL,np=100000,qcutoff=0.2,usec=T)
{
  #check dataframe
  if(is.null(genotype_columns)) stop("column numbers of genotypes must be provided!")
  if(is.null(adjusting_covariate_columns)) 
    stop("column numbers of adjusting covariates must be provided!")
  if(is.null(outcome_column)) stop("column number of outcome must be provided!")
  if(is.null(coeff)) stop("external coefficients  must be provided!")
  if (length(coeff) != length(genotype_columns)) 
    stop("the number of genetypes are not consistent with the number of coefficients!")
  if (sum(names(coeff) %in% colnames(data)[genotype_columns]) != length(coeff)) 
    stop("names of genotypes in the dataframe are not consistant with names of coefficients!")
  print(Sys.time())
  #remvoe rows containing NAs in the related columns
  data=remove_missingdata(data[,c(genotype_columns,adjusting_covariate_columns,outcome_column)])
  nsnp=length(coeff)
  crc_coeff <- matrix(0,nsnp,3)
  for (i in 1:nsnp){
    fm=paste0(colnames(data)[outcome_column],"~",colnames(data)[which(colnames(data)==names(coeff)[i])])
    for (j in 1:length(adjusting_covariate_columns))
    {
      fm=paste0(fm,"+",colnames(data)[which(colnames(data)==colnames(data)[adjusting_covariate_columns[j]])])
    }
    fit <- glm(as.formula(fm),family=poisson,data=data)
    crc_coeff[i,1] <- fit$coef[2]
    fit <- glm(as.formula(fm),family=binomial,data=data)
    crc_coeff[i,2] <- fit$coef[2]
    #variance
    crc_coeff[i,3] <- summary(fit)$coef[2,2]^2
  }

  yy <- crc_coeff[,2]
  xx <- cbind(1,coeff)
  ww <- diag(1/crc_coeff[,3])


  bb <- solve(t(xx) %*% ww %*% xx) %*% t(xx) %*% ww %*% yy
  bbcov<- solve(t(xx) %*% ww %*% xx)


  data$grs <- 0
  for (i in 1:nsnp){
    k<- which(colnames(data)==names(coeff)[i])
    data$grs <- data$grs + coeff[i]*data[,k]
  }

  outp <- rep(0,nsnp)
  for (i in 1:nsnp) {
    k <- which(colnames(data)==names(coeff[i]))
    fm=paste0(colnames(data)[outcome_column],"~grs+",colnames(data)[which(colnames(data)==names(coeff)[i])])
    for (j in 1:length(adjusting_covariate_columns))
    {
      fm=paste0(fm,"+",colnames(data)[which(colnames(data)==colnames(data)[adjusting_covariate_columns[j]])])
    }
    fit <- glm(as.formula(fm),family=binomial,data=data)
    outp[i] <- summary(fit)$coef[3,4]
  }

  ### now work out the correlation matrix under the null hypothesis ####

  fm=paste0(colnames(data)[outcome_column],"~grs+")
  for (j in 1:length(adjusting_covariate_columns))
  {
    fm=paste0(fm,"+",colnames(data)[which(colnames(data)==colnames(data)[adjusting_covariate_columns[j]])])
  }

  fit <- glm(as.formula(fm),family=binomial,data=data,y=T)
  yfit <- fit$fitted.values
  y <-fit$y #outcome
  xmat=model.matrix(as.formula(fm),data=data)
  grs=data$grs
  #generate genotype matrix, which columns are consistent with coeff
  data_genotype=data.frame(matrix(NA,nrow=nrow(data),ncol=length(coeff)))
  colnames(data_genotype)=names(coeff)
  for (i in 1:length(coeff))
  {
    idx=which(colnames(data)==names(coeff)[i])
    data_genotype[,i]=data[,idx]
  }
  cormat <- matrix(1,nsnp,nsnp)
  # #storage.mode(x) <- storage.mode(y) <- "double"
  print(Sys.time())
  print("start to compute the correlation matrix...")
  
  if (usec==T)
  {
    result <- .C("compute_cormat",
                 nsnp=as.integer(nsnp),
                 n_subject=as.integer(nrow(data)),
                 ncol_xmat=as.integer(ncol(xmat)),
                 yfit=as.double(yfit),
                 y=as.double(y),
                 xmat=as.double(xmat),
                 data_genotype=data.matrix(data_genotype),
                 cormat=as.double(cormat),
                 PACKAGE="GLIDE"
    )
    for (j in 1:nsnp)
    {
      for (i in 1:nsnp)
      {
        cormat[i,j]=result$cormat[(j-1)*nsnp+i]
      }
    }
  }else
  {
    for (i in 1:(nsnp-1)){
      cat(i,"..")
      for (j in (i+1):nsnp) {
        ii <- which(colnames(data)==names(coeff)[i])
        jj <- which(colnames(data)==names(coeff)[j])
        fm=paste0(colnames(data)[outcome_column],"~grs+",colnames(data)[ii])
        for (k in 1:length(adjusting_covariate_columns))
        {
          fm=paste0(fm,"+",colnames(data)[which(colnames(data)==colnames(data)[adjusting_covariate_columns[k]])])
        }
        xmat1 <- model.matrix(as.formula(fm),data=data)
        fm=paste0(colnames(data)[outcome_column],"~grs+",colnames(data)[jj])
        for (k in 1:length(adjusting_covariate_columns))
        {
          fm=paste0(fm,"+",colnames(data)[which(colnames(data)==colnames(data)[adjusting_covariate_columns[k]])])
        }
        xmat2 <- model.matrix(as.formula(fm),data=data)
        
        xmat  <- cbind(xmat1,xmat2)
        bread <- matrix(0,ncol(xmat),ncol(xmat))
        bread1 <- t(xmat1)%*%(xmat1*drop(yfit*(1-yfit)))
        bread2 <- t(xmat2)%*%(xmat2*drop(yfit*(1-yfit)))
        
        bread[1:ncol(xmat1),1:ncol(xmat1)] <- bread1
        bread[(ncol(xmat2)+1):ncol(xmat),(ncol(xmat2)+1):ncol(xmat)] <- bread2
        
        score <- cbind(xmat1*drop(fit$y-yfit),xmat2*drop(fit$y-yfit))
        beef <- t(score) %*% score
        
        covmat <- solve(bread) %*% beef %*% solve(bread)
  
        cormat[i,j] <- covmat[3,ncol(xmat1)+3]/sqrt(covmat[3,3]*covmat[ncol(xmat1)+3,ncol(xmat1)+3])
        cormat[j,i] <- cormat[i,j]
      }
    }
  }
  
  
  print(Sys.time())
  zsim <- matrix(0,np,nsnp)
  print(paste0("start to null p-value..."))
  if (usec==T)
  {
    result1=.C("compute_z",
               nsnp=as.integer(nsnp),
               cormat=data.matrix(cormat),
               np=as.integer(np),
               zsim=as.double(zsim),
               PACKAGE="GLIDE")
    for (i in 1:np)
    {
      for (j in 1:nsnp)
      {
        zsim[i,j]=result1$zsim[(i-1)*nsnp+j]
      }
    }
  }else
  {
    zsim[,1] <- rnorm(np,0,1)
    for (k in 2:nsnp) {
      mu <- drop(zsim[,1:(k-1)] %*% ginv(cormat[1:(k-1),1:(k-1)]) %*% cormat[1:(k-1),k])
      if (k !=nsnp) sig<- sqrt(1- drop(cormat[1:(k-1),k] %*% ginv(cormat[1:(k-1),1:(k-1)]) %*% cormat[1:(k-1),k]))
      if (k !=nsnp) zsim[,k] <- rnorm(np,mu,sig) else zsim[,k] <- mu
    }
    for (k in 1:nsnp) zsim[,k] <- 2*(1-pnorm(abs(zsim[,k])))
    for (k in 1:np) {
      zsim[k,] <- zsim[k,order(zsim[k,])]
    }
  }
  
  print(Sys.time())  
 
  orderp <- apply(zsim,2,mean)
  
  ### compute the FWER and FDR values for each observed p-value
  
  fwer <- rep(0,length(outp))
  qval <- rep(0,length(outp))
  for (i in 1:length(outp)) {
    temp <- apply(1*(zsim <= rep(outp[i],np)),1,sum)  
    fwer[i] <- mean(temp>=1)
    qval[i] <- (sum(temp)/np)/sum(outp<=outp[i])
  }
  out <- matrix(0,nsnp,4)
  out[,1] <- outp
  out[,2] <- orderp
  out[,3] <- fwer
  out[,4] <- qval
  qval <- rep(0,length(outp))
  for (i in 1:length(outp)) {
    qval[i] <- min(out[out[,1]>= out[i,1],4])
  }
  out[,4] <- qval
  colnames(out)=c("observed_pvalue","null_pvalue","fwer","q_vlaue")
  rownames(out)=names(coeff)
  selected_idx=which(qval<=qcutoff)
  
  #draw plots
  par(mfrow=c(1,2))
  plot(coeff,crc_coeff[,2],xlab="Genetic effect for the exposure",ylab="Genetic effect for the outcome")
  abline(bb[1],bb[2],lty=2)
  text(mean(par("usr")[1:2]),mean(par("usr")[3:4]),paste0("y=",format(bb[1],digits=2),"+",format(bb[2],digits=2),"*x"),cex=1.1)
  egger_pvalue=format(2*(1-pnorm(abs(bb[1,1])/sqrt(bbcov[1,1]))),digits = 2)
  text(mean(par("usr")[1:2]),mean(c(mean(par("usr")[3:4]),par("usr")[3])),paste0("p-value=",egger_pvalue),cex=1.1)
  #text(0.07,0.01,"y=0.001+0.25*x",cex=1.1)
  if (length(selected_idx)>0)
  {
    points(coeff[selected_idx],crc_coeff[selected_idx,2],col=2,pch=16,cex=1.4)
    legend("bottomright",legend="SNP with direct effect",col=2,pch=16)
  }
  
  plot(-log(orderp,base=10),-log(outp[order(outp)],base=10),xlab="Expected null p-values (log base 10)",ylab="Observed p-values (log base 10)",type="n")
  points(-log(orderp,base=10),-log(outp[order(outp)],base=10),cex=1.4)
  
  
  if (length(selected_idx)>0)
  {
    points(-log(orderp[1:length(selected_idx)],base=10),-log(outp[order(outp)][1:length(selected_idx)],base=10),col=2,pch=16,cex=1.4)
    legend("bottomright",legend="SNP with direct effect",col=2,pch=16)
  }
  abline(0,1)
  par(mfrow=c(1,1))
  return(out)
}