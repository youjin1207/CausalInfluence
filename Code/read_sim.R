library(xtable)
library(gplots)
library(RColorBrewer)
###
topkl.concordance = function(k, l, centrality, tau){
  # k : top k of tau
  # l : top l of centrality
  # centrality : centrality measure
  # tau : causal influence measure
  
  n = length(centrality) # the number of nodes
  a = order(-centrality)
  b = order(-tau)
  if(k <= l){
    check = sum(b[c(1:as.integer(k))] %in% a[c(1:as.integer(l))])
    concordance = (check == as.integer(k))
  }else{
    check = sum(a[c(1:as.integer(l))] %in% b[c(1:as.integer(k))])
    concordance = (check == as.integer(l)) 
  }
  return(concordance)
}
##################  
nr = 500; n = 277
outdegree = btwn.centrality = closeout.centrality = 
  diff.centrality = random = proportion = matrix(0, nr, n)

# outdegree: nr x n matrix of degree centrality measures for n subjects from nr replicates
# btwn.centrality: nr x n matrix of betweenness centrality measures for n subjects from nr replicates 
# diff.centrality: nr x n matrix of diffusion centrality measures for n subjects from nr replicates 
# random: nr x n matrix of random selection for n subjects from nr replicates 
# proportion : nr x n matrix of causal influence measures for n subjects from nr replicates 

result.spearman = matrix(0, nr, 4)
for(r in 1:nr){
  result.spearman[r,] = c(cor.test(outdegree[r,], proportion[r,], exact = FALSE, method = 'spearman')$estimate
                          ,cor.test(btwn.centrality[r,], proportion[r,], exact = FALSE, method = 'spearman')$estimate
                          ,cor.test(diff.centrality[r,], proportion[r,], exact = FALSE, method = 'spearman')$estimate
                          ,cor.test(random[r,], proportion[r,], exact = FALSE, method = 'spearman')$estimate)
}

colMeans(result.spearman)

print.spearman = c()
for(i in 1:4){
  print.spearman[i] = paste(formatC(mean(result.spearman[,i]),2,format='f'), " (",
                           formatC(sd(result.spearman[,i]),2,format='f'), ")", sep="")
}


percentage.table1 = percentage.table2 = percentage.table3 = 
  percentage.table4 =  matrix(0, 277, 277)
for(k in 1:277){ print(k)
  for(l in 1:277){
    for(r in 1:nr){
      percentage.table1[k,l] = percentage.table1[k,l] + topkl.concordance(k, l, outdegree[r,], proportion[r,]) / nr
      percentage.table2[k,l] = percentage.table2[k,l] + topkl.concordance(k, l, btwn.centrality[r,], proportion[r,]) / nr
      percentage.table3[k,l] = percentage.table3[k,l] + topkl.concordance(k, l, diff.centrality[r,], proportion[r,]) / nr
      percentage.table4[k,l] = percentage.table4[k,l] + topkl.concordance(k, l, random[r,], proportion[r,]) / nr
    }
  }
}

percentage.table1 = percentage.table1[seq(nrow(percentage.table1), 1, -1),]
percentage.table2 = percentage.table2[seq(nrow(percentage.table2), 1, -1),]
percentage.table3 = percentage.table3[seq(nrow(percentage.table3), 1, -1),]
percentage.table4 = percentage.table4[seq(nrow(percentage.table4), 1, -1),]


#Palette = colorRampPalette(brewer.pal(9,"YlGnBu"))(20)
#Palette = colorRampPalette(brewer.pal(9,"YlGn"))(20)
Palette = colorRampPalette(brewer.pal(9,"YlOrRd"))(20)
########### figures of matrices ######
percentage.table4[277-20,50]
pdf("Figure/random_diffusion_model.pdf", width = 10, height = 8)
par(cex.main=3, cex.lab = 5, cex.axis = 5, 
    mar=c(3,4,3,10), tcl = 0.5, oma = c(0, 0, 0, 0), xpd = TRUE)
heatmap.2(percentage.table4, dendrogram='none', 
          Rowv=FALSE, Colv=FALSE,trace='none', 
          xlab = "", 
          ylab = "", 
          main = "Random",
          labCol = c(1, rep(NA, 48), "l=50", rep(NA, 49), "l=100", rep(NA,49), "l=150", rep(NA, 66), "l=277"),
          labRow = c("k=277", rep(NA, 66), "k=150", rep(NA, 49), "k=100", rep(NA,49), "k=50 ", rep(NA, 49)),
          cexRow = 2.5, cexCol = 2.5,
          key.par=list(mar=c(1,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          adjRow = c(2,1),
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)), lhei=c(0.3, 2.5, 0.2), lwid=c(0.5, 4, 1),
          key.title = "",
          col = Palette,
          srtCol=0,
          key = FALSE,
          breaks=seq(0, 1, 1/20),
          add.expr = list(rect(50, 20, 51, 21,border="blue",lwd=5), 
                          lines(c(1,50), c(20, 20), col = "blue", lty = 2),
                          lines(c(50, 50), c(1, 20), col = "blue", lty = 2)),
          offsetRow = -45, density.info = "none")
mtext("Top l of centrality", side=1, line=0, cex=2)
mtext( expression(paste("Top k of ", tau, sep="")), side=4, line=0, cex=2)
legend(x = 1, y = 1, legend=c("0%", "25%", "50%",  "75%", "100%"), 
       fill= c(Palette[1],  Palette[5], Palette[10],  Palette[15], Palette[20]), cex = 1.5, xpd = TRUE)
legend(x = 0.94, y = 0.3, "0% of top k=20 \n influential nodes \n include nodes with \n top l=50 centrality",
       col = "blue", pch = 19, cex = 1.3, bty = "n", pt.lwd = 2)
dev.off()