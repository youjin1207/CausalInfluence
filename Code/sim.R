library(igraph)
library(parallel)
library(Matrix)
library(keyplayer)
###

# generate a
homo.direct = function(alpha0, alpha1, alpha2, Adj){
  N = nrow(Adj)
  influence.Y = rep(NA, N)  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    prop = alpha0 + alpha1*treat + alpha2*as.numeric(Adj %*% treat)/ max(rowSums(Adj))
    Y = prop
    influence.Y[i] = mean(Y) 
  }
  return(influence.Y)
}


traffic.array = function(Adj){
  traffic = array(0, dim = c(nrow(Adj), nrow(Adj), nrow(Adj)))
  G = graph.adjacency(Adj)
  for(i in 1:nrow(Adj)){
    for(k in 1:nrow(Adj)){
      tmp = all_shortest_paths(G, from = k, to = i, mode = "all")
      tmp2 = c(); if(length(tmp$res) > 0){
        for(jj in 1:length(tmp$res)){
          tmp2 = c(tmp2, as.numeric(tmp$res[[jj]]))
        }
        for(j in 1:nrow(Adj)){
          traffic[i,k,j] =  sum(tmp2 %in% j) / length(tmp$res)
          traffic[i,k,i] =  0; traffic[i,k,k] =  0
        }
      }
    }
  }
  return(traffic)
}


traffic.model = function(alpha0, alpha1, alpha2, traffic){
  N = nrow(traffic)
  influence.Y = rep(NA, N) 
  error = rnorm(N, 0, 0.0)
  for(j in 1:N){
    treat = rep(0, N); treat[j] = 1; prop = c()
    for(i in 1:N){
      prop[i] = plogis(alpha0 + alpha1*treat[i] + alpha2*sum(as.numeric(traffic[i,,] %*% treat)[-i]) + error[i])
    }
    Y = prop
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y =  rowMeans(prop))
}


diff.randomwalk = function(p, alpha0, alpha1, alpha2, time, Adj){
  # Adj : adjacency matrix
  # p: probability of adoption
  # alpha0 : probability of being active under treatment but without any influence
  # alpha1 : probability of being active under control but without any influence
  # alpha2 : probability of being active under control and having one active friends.
  # time : time spent under random walk
  N = nrow(Adj); out = matrix(0, N, time+1); influence.Y = rep(NA, N)
  
  ### define a diffusion process
  #for(i in 1:N){
  #  treat = rep(0, N); treat[i] = 1
  #  out[,1] = rep(alpha0, N)
  #  propmat = p*Adj; 
  #  out[,2] = alpha2*propmat[i, ] + alpha1*out[, 1] 
  #  for(t in 3:(time+1)){
  #    propmat = propmat + propmat %*% (p * Adj)
  #    out[,t] = alpha2*propmat[i,] + alpha1*out[, t-1] 
  #  }
  #  influence.Y[i] = mean(out[,t])
  #}
  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = alpha2*propmat[i, ] #+ alpha1*out[, 1] 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*propmat[i,] #+ alpha1*out[, t-1] 
    }
    influence.Y[i] = mean(out[,t])
  }
  
  return(influence.Y)
}


n.rep = 500
#homo.direct.result = traffic.model.result = diffusion.model.result = matrix(NA, n.rep, 217)
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)


for(r in 1:n.rep){
  print(r)
  set.seed(r)
  #G = erdos.renyi.game(217, 0.05, type = "gnp")
  G = sample_sbm(217, pref.matrix = cbind(c(0.1, 0.01), c(0.01, 0.1)), block.sizes = c(100, 117), 
                 directed = FALSE, loops = FALSE)
  #plot(G, vertex.label= NA, edge.arrow.size=0.02, vertex.size = 0.5, xlab = "Random Network: G(217,363) model")
  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result = homo.direct(0.1, 0.2, 0.3, Adj)
  traffic.model.result = traffic.model(-0.5, 1, 0.5, traffic)
  diffusion.model.result = diff.randomwalk(p = 0.3, alpha0 = 0.1, alpha1 = 0.4, alpha2 = 0.3, time = 5, Adj)
  
  ## results
  diffusion.measure = diffusion(Adj*0.3,  T = 5)
  #diffusion.measure = sample(1:N, N, replace = FALSE)
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


