plot.glide <- function(out,qcutoff=0.2)
{
  out=as.data.frame(out)
  
  yy <- out$geffect_outcome
  xx <- cbind(1,out$geffect_exposure)
  ww <- diag(1/out$geffect_outcome_variance)
  bb <- solve(t(xx) %*% ww %*% xx) %*% t(xx) %*% ww %*% yy
  bbcov<- solve(t(xx) %*% ww %*% xx)
  #SNPs with evidence of direct effect by GLIDE based on the selected q-value cutoff
  selected_idx=which(out$q_value<=qcutoff)
  
  plot(-log(out$expected_pvalue,base=10),-log(out$observed_pvalue[order(out$observed_pvalue)],base=10),xlab="Expected null p-values (log base 10)",ylab="Observed p-values (log base 10)",type="n")
  points(-log(out$expected_pvalue,base=10),-log(out$observed_pvalue[order(out$observed_pvalue)],base=10),cex=1.4)
  
  if (length(selected_idx)>0)
  {
    points(-log(out$expected_pvalue[1:length(selected_idx)],base=10),-log(out$observed_pvalue[order(out$observed_pvalue)][1:length(selected_idx)],base=10),col=2,pch=16,cex=1.4)
    legend("bottomright",legend="SNPs with evidence of direct effect by GLIDE",col=2,pch=16, cex=0.8)
  }
  abline(0,1)
  
}
