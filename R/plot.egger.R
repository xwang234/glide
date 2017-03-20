plot.egger <- function(out,qcutoff=0.2,xlable="Genetic association with the exposure",ylable="Genetic association with the outcome")
{
  out=as.data.frame(out)
  
  yy <- out$geffect_outcome
  xx <- cbind(1,out$geffect_exposure)
  ww <- diag(1/out$geffect_outcome_variance)
  bb <- solve(t(xx) %*% ww %*% xx) %*% t(xx) %*% ww %*% yy
  bbcov<- solve(t(xx) %*% ww %*% xx)
  #SNPs with evidence of direct effect by GLIDE based on the selected q-value cutoff
  selected_idx=which(out$q_value<=qcutoff)
  
  plot(out$geffect_exposure,out$geffect_outcome,xlab=xlable,ylab=ylable)
  abline(bb[1],bb[2],lty=2)
  text(mean(par("usr")[1:2]),mean(par("usr")[3:4]),paste0("y=",format(bb[1],digits=2),"+",format(bb[2],digits=2),"*x"),cex=1.1)
  egger_pvalue=format(2*(1-pnorm(abs(bb[1,1])/sqrt(bbcov[1,1]))),digits = 2)
  text(mean(par("usr")[1:2]),mean(c(mean(par("usr")[3:4]),par("usr")[3])),paste0("Egger intercept p-value=",egger_pvalue),cex=1.1)
  #text(0.07,0.01,"y=0.001+0.25*x",cex=1.1)
  if (length(selected_idx)>0)
  {
    points(out$geffect_exposure[selected_idx],out$geffect_outcome[selected_idx],col=2,pch=16,cex=1.4)
    legend("bottomright",legend="SNPs with evidence of direct effect by GLIDE",col=2,pch=16, cex=0.8)
  }
}
