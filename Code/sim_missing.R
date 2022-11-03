library(igraph)
library(parallel)
library(Matrix)
library(keyplayer)
library(gtools)
####
source("Code/aux_functions.R")

####
n.rep = 500
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
missing.p = 0.1 # 0.3, 0.5

for(r in 1:n.rep){
  set.seed(r)
  G = sample_sbm(277, pref.matrix = cbind(c(0.1, 0.01), c(0.01, 0.1)), block.sizes = c(100, 177), 
                 directed = FALSE, loops = FALSE)
  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  ## missing edges
  obs.Adj = Adj
  for(i in 1:N){
    for(j in i:N){
      obs.Adj[i,j] = ifelse(Adj[i,j] == 1, rbinom(1, 1, 1-missing.p), Adj[i,j])
    }
    for(j in 1:i){
      obs.Adj[i,j] = obs.Adj[j,i]
    }
  }
  obs.G = graph.adjacency(obs.Adj)
  
  homo.direct.result = homo.direct(10, 20, 50, Adj)
  traffic.model.result = traffic.model(10, 20, 0.1, traffic)
  diffusion.model.result = diff.randomwalk(p = 0.3, alpha0 = 10, alpha1 = 20, alpha2 = 0.01, time = 5, Adj)
  
  ## results
  diffusion.measure = diffusion(obs.Adj*0.3,  T = 5)
  random.measure = sample(1:N, N, replace = FALSE)
  homo.direct.result.mat[r,] = c(cor(homo.direct.result, rowSums(obs.Adj), method = "spearman"),
                                 cor(homo.direct.result, betweenness(obs.G), method = "spearman"),
                                 cor(homo.direct.result,  diffusion.measure, method = "spearman"),
                                 cor(homo.direct.result,  random.measure, method = "spearman"))
  traffic.model.result.mat[r,] = c(cor(traffic.model.result, rowSums(obs.Adj), method = "spearman"),
                                   cor(traffic.model.result, betweenness(obs.G), method = "spearman"),
                                   cor(traffic.model.result, diffusion.measure, method = "spearman"),
                                   cor(traffic.model.result,  random.measure, method = "spearman"))
  diffusion.model.result.mat[r,] = c(cor(diffusion.model.result, rowSums(obs.Adj), method = "spearman"),
                                     cor(diffusion.model.result, betweenness(obs.G), method = "spearman"),
                                     cor(diffusion.model.result, diffusion.measure, method = "spearman"),
                                     cor(diffusion.model.result, random.measure, method = "spearman"))
}
mat = rbind(colMeans(homo.direct.result.mat),
            colMeans(diffusion.model.result.mat),
            colMeans(traffic.model.result.mat))
colnames(mat) = c("Degree centrality", "Betweenness", "Diffusion", "Random")
print(xtable(mat, digits = 2))