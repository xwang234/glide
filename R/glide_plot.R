glide_plot <- function(out,qcutoff=0.2,xlable="Genetic effect for the exposure",ylable="Genetic effect for the outcome")
{
  out=as.data.frame(out)
  
  yy <- out$genetic_effect_outcome
  xx <- cbind(1,out$genetic_effect_exposure)
  ww <- diag(1/out$genetic_effect_outcome_variance)
  
  
  bb <- solve(t(xx) %*% ww %*% xx) %*% t(xx) %*% ww %*% yy
  bbcov<- solve(t(xx) %*% ww %*% xx)
  
  selected_idx=which(out$q_value<=qcutoff)
  #draw plots
  par(mfrow=c(1,2))
  plot(out$genetic_effect_exposure,out$genetic_effect_outcome,xlab=xlable,ylab=ylable)
  abline(bb[1],bb[2],lty=2)
  text(mean(par("usr")[1:2]),mean(par("usr")[3:4]),paste0("y=",format(bb[1],digits=2),"+",format(bb[2],digits=2),"*x"),cex=1.1)
  egger_pvalue=format(2*(1-pnorm(abs(bb[1,1])/sqrt(bbcov[1,1]))),digits = 2)
  text(mean(par("usr")[1:2]),mean(c(mean(par("usr")[3:4]),par("usr")[3])),paste0("p-value=",egger_pvalue),cex=1.1)
  #text(0.07,0.01,"y=0.001+0.25*x",cex=1.1)
  if (length(selected_idx)>0)
  {
    points(out$genetic_effect_exposure[selected_idx],out$genetic_effect_outcome[selected_idx],col=2,pch=16,cex=1.4)
    legend("bottomright",legend="SNP with direct effect",col=2,pch=16)
  }
  
  plot(-log(out$null_pvalue,base=10),-log(out$observed_pvalue[order(out$observed_pvalue)],base=10),xlab="Expected null p-values (log base 10)",ylab="Observed p-values (log base 10)",type="n")
  points(-log(out$null_pvalue,base=10),-log(out$observed_pvalue[order(out$observed_pvalue)],base=10),cex=1.4)
  
  
  if (length(selected_idx)>0)
  {
    points(-log(out$null_pvalue[1:length(selected_idx)],base=10),-log(out$observed_pvalue[order(out$observed_pvalue)][1:length(selected_idx)],base=10),col=2,pch=16,cex=1.4)
    legend("bottomright",legend="SNP with direct effect",col=2,pch=16)
  }
  abline(0,1)
  par(mfrow=c(1,1))
}
