library(igraph)
library(parallel)
library(Matrix)
library(keyplayer)
###
source("Code/aux_functions.R")

#### adding random errors in the outcome model ####
n.rep = 500
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
sigmas = c(0.5, 1.0, 1.5, 2.0)
for(r in 1:n.rep){
  set.seed(r)

  G = sample_sbm(277, pref.matrix = cbind(c(0.1, 0.01), c(0.01, 0.1)), block.sizes = c(100, 177), 
                 directed = FALSE, loops = FALSE)

  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result = homo.direct.sigma(10, 20, 50, Adj, sigmas[1]) 
  traffic.model.result = traffic.model.sigma(10, 20, 0.1, traffic, sigmas[1])
  diffusion.model.result = diff.randomwalk.sigma(p = 0.3, alpha0 = 10, alpha1 = 20, alpha2 = 0.01, time = 5, Adj, sigmas[1])
  
  ## results
  diffusion.measure = diffusion(Adj*0.3,  T = 5)
  random.measure = sample(1:N, N, replace = FALSE)
  homo.direct.result.mat[r,] = c(cor(homo.direct.result, rowSums(Adj), method = "spearman"),
                                 cor(homo.direct.result, betweenness(G), method = "spearman"),
                                 cor(homo.direct.result,  diffusion.measure, method = "spearman"),
                                 cor(homo.direct.result,  random.measure, method = "spearman"))
  traffic.model.result.mat[r,] = c(cor(traffic.model.result, rowSums(Adj), method = "spearman"),
                                   cor(traffic.model.result, betweenness(G), method = "spearman"),
                                   cor(traffic.model.result, diffusion.measure, method = "spearman"),
                                   cor(traffic.model.result,  random.measure, method = "spearman"))
  diffusion.model.result.mat[r,] = c(cor(diffusion.model.result, rowSums(Adj), method = "spearman"),
                                     cor(diffusion.model.result, betweenness(G), method = "spearman"),
                                     cor(diffusion.model.result, diffusion.measure, method = "spearman"),
                                     cor(diffusion.model.result, random.measure, method = "spearman"))
}

mat = rbind(colMeans(homo.direct.result.mat),
            colMeans(diffusion.model.result.mat),
            colMeans(traffic.model.result.mat))
colnames(mat) = c("Degree centrality", "Betweenness", "Diffusion", "Random")
print(xtable(mat, digits = 2))


#### Heterogeneous alpha's ####
n.rep = 500
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
alphas = c(10, 20, 30)
for(r in 1:n.rep){
  set.seed(r)
  G = sample_sbm(277, pref.matrix = cbind(c(0.1, 0.01), c(0.01, 0.1)), block.sizes = c(100, 177), 
                 directed = FALSE, loops = FALSE)
  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result = homo.direct.alphas(10, 20, 50, Adj, alphas[1])
  traffic.model.result = traffic.model.alphas(10, 20, 0.1, traffic, alphas[1])
  diffusion.model.result = diff.randomwalk.alphas(p = 0.3, alpha0 = 10, alpha1 = 20, alpha2 = 0.01, time = 5, Adj, alphas[1])
  
  ## results
  diffusion.measure = diffusion(Adj*0.3,  T = 5)
  random.measure = sample(1:N, N, replace = FALSE)
  homo.direct.result.mat[r,] = c(cor(homo.direct.result, rowSums(Adj), method = "spearman"),
                                 cor(homo.direct.result, betweenness(G), method = "spearman"),
                                 cor(homo.direct.result,  diffusion.measure, method = "spearman"),
                                 cor(homo.direct.result,  random.measure, method = "spearman"))
  traffic.model.result.mat[r,] = c(cor(traffic.model.result, rowSums(Adj), method = "spearman"),
                                   cor(traffic.model.result, betweenness(G), method = "spearman"),
                                   cor(traffic.model.result, diffusion.measure, method = "spearman"),
                                   cor(traffic.model.result,  random.measure, method = "spearman"))
  diffusion.model.result.mat[r,] = c(cor(diffusion.model.result, rowSums(Adj), method = "spearman"),
                                     cor(diffusion.model.result, betweenness(G), method = "spearman"),
                                     cor(diffusion.model.result, diffusion.measure, method = "spearman"),
                                     cor(diffusion.model.result, random.measure, method = "spearman"))
}
mat = rbind(colMeans(homo.direct.result.mat),
            colMeans(diffusion.model.result.mat),
            colMeans(traffic.model.result.mat))
colnames(mat) = c("Degree centrality", "Betweenness", "Diffusion", "Random")
print(xtable(mat, digits = 2))


#### heterogeneous beta's ####
n.rep = 500
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
betas = c(1,2,3)

for(r in 1:n.rep){
  print(r)
  set.seed(r)
  
  G = sample_sbm(277, pref.matrix = cbind(c(0.1, 0.01), c(0.01, 0.1)), block.sizes = c(100, 177), 
                 directed = FALSE, loops = FALSE)
  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result = homo.direct.betas(10, 20, 50, Adj, 50*betas[1])
  traffic.model.result = traffic.model.betas(10, 20, 0.1, traffic, 0.1*betas[1])
  diffusion.model.result = diff.randomwalk.betas(p = 0.3, alpha0 = 10, alpha1 = 20, alpha2 = 0.01, time = 5, Adj, 0.01*betas[1]*10)
  
  ## results
  diffusion.measure = diffusion(Adj*0.3,  T = 5)
  random.measure = sample(1:N, N, replace = FALSE)
  homo.direct.result.mat[r,] = c(cor(homo.direct.result, rowSums(Adj), method = "spearman"),
                                 cor(homo.direct.result, betweenness(G), method = "spearman"),
                                 cor(homo.direct.result,  diffusion.measure, method = "spearman"),
                                 cor(homo.direct.result,  random.measure, method = "spearman"))
  traffic.model.result.mat[r,] = c(cor(traffic.model.result, rowSums(Adj), method = "spearman"),
                                   cor(traffic.model.result, betweenness(G), method = "spearman"),
                                   cor(traffic.model.result, diffusion.measure, method = "spearman"),
                                   cor(traffic.model.result,  random.measure, method = "spearman"))
  diffusion.model.result.mat[r,] = c(cor(diffusion.model.result, rowSums(Adj), method = "spearman"),
                                     cor(diffusion.model.result, betweenness(G), method = "spearman"),
                                     cor(diffusion.model.result, diffusion.measure, method = "spearman"),
                                     cor(diffusion.model.result, random.measure, method = "spearman"))
}
mat = rbind(colMeans(homo.direct.result.mat),
            colMeans(diffusion.model.result.mat),
            colMeans(traffic.model.result.mat))
colnames(mat) = c("Degree centrality", "Betweenness", "Diffusion", "Random")
print(xtable(mat, digits = 2))
