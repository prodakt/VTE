myvolcano  = function(myx,myy,vt1,vt2,ht,myxlab="log2(fold changes)",myylab="-log10(q-values)", cex.lab=1.2, kolory, maksY, limitX,...){
  plot(myx,myy, pch=20, col=kolory, xlab=myxlab, ylab=myylab, cex.lab=cex.lab, ylim=c(0,maksY),xlim=c(-limitX,limitX),...)
  abline(v=vt1, lty=2)
  abline(v=vt2, lty=2)
  abline(h=ht, lty=2)
}


cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}
