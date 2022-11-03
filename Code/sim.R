library(igraph)
library(parallel)
library(Matrix)
library(keyplayer)
###
source("Code/aux_functions.R")

####
n.rep = 500
homo.direct.result = traffic.model.result = diffusion.model.result = matrix(NA, n.rep, 277)
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
centrality= list()
for(r in 1:n.rep){
  set.seed(r)
  G = sample_sbm(277, pref.matrix = cbind(c(0.1, 0.01), c(0.01, 0.1)), block.sizes = c(100, 177), 
                 directed = FALSE, loops = FALSE)

  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result[r,] = homo.direct(10, 20, 50, Adj)
  traffic.model.result[r,] = traffic.model(10, 20, 0.1, traffic)
  diffusion.model.result[r,] = diff.randomwalk(p = 0.3, alpha0 = 10, alpha1 = 20, alpha2 = 0.01, time = 5, Adj)
  
  ## results
  diffusion.measure = diffusion(Adj*0.3,  T = 5)
  random.measure = sample(1:N, N, replace = FALSE)
  
  centrality[[r]] = list(rowSums(Adj), betweenness(G), diffusion.measure, random.measure)
  
  homo.direct.result.mat[r,] = c(cor(homo.direct.result[r,], rowSums(Adj), method = "spearman"),
                                 cor(homo.direct.result[r,], betweenness(G), method = "spearman"),
                                 cor(homo.direct.result[r,],  diffusion.measure, method = "spearman"),
                                 cor(homo.direct.result[r,],  random.measure, method = "spearman"))
  traffic.model.result.mat[r,] = c(cor(traffic.model.result[r,], rowSums(Adj), method = "spearman"),
                                   cor(traffic.model.result[r,], betweenness(G), method = "spearman"),
                                   cor(traffic.model.result[r,], diffusion.measure, method = "spearman"),
                                   cor(traffic.model.result[r,],  random.measure, method = "spearman"))
  diffusion.model.result.mat[r,] = c(cor(diffusion.model.result[r,], rowSums(Adj), method = "spearman"),
                                     cor(diffusion.model.result[r,], betweenness(G), method = "spearman"),
                                     cor(diffusion.model.result[r,], diffusion.measure, method = "spearman"),
                                     cor(diffusion.model.result[r,], random.measure, method = "spearman"))
}


save(homo.direct.result.mat, file = "data/homo.direct.result.mat.RData")
save(traffic.model.result.mat, file = "data/traffic.model.result.mat.RData")
save(diffusion.model.result.mat, file = "data/diffusion.model.result.mat.RData")
