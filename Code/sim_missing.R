library(igraph)
library(parallel)
library(Matrix)
library(keyplayer)
library(gtools)
###
load("Data/sample_net.RData")

# generate (1) homogeneous direct interference
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

# calculate filtering functions of the given network 
traffic.array = function(Adj){
  traffic = array(0, dim = c(nrow(Adj), nrow(Adj), nrow(Adj)))
  G = graph.adjacency(Adj)
  for(i in 1:nrow(Adj)){
    for(k in 1:nrow(Adj)){
      tmp = get.all.shortest.paths(G, from = k, to = i)
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

# generate (2) traffic-dependent process
traffic.model = function(alpha0, alpha1, alpha2, traffic){
  N = nrow(traffic)
  influence.Y = rep(NA, N) 
  for(j in 1:N){
    treat = rep(0, N); treat[j] = 1; prop = c()
    for(i in 1:N){
      prop[i] = plogis(alpha0 + alpha1*treat[i] + alpha2*sum(as.numeric(traffic[i,,] %*% treat)[-i]))
    }
    Y = prop
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y)
}

# generate (3) homogeneous diffusion process
diff.randomwalk = function(p, alpha0, alpha1, alpha2, time, Adj){
  # Adj : adjacency matrix
  # p: probability of adoption
  # alpha0 : probability of being active under treatment but without any influence
  # alpha1 : probability of being active under control but without any influence
  # alpha2 : probability of being active under control and having one active friends.
  # time : time spent under random walk
  N = nrow(Adj); out = matrix(0, N, time+1); influence.Y = rep(NA, N)
  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = alpha2*propmat[i, ] 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*propmat[i,]  
    }
    influence.Y[i] = mean(out[,t])
  }
  return(influence.Y)
}


n.rep = 500
homo.direct.result = traffic.model.result = diffusion.model.result = matrix(NA, n.rep, 277)
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
centrality = list()
missing.p = 0.1 # 0.1, 0.3, 0.5
for(r in 1:n.rep){
  set.seed(r)
  print(r)
  G = net0
  edges = as.matrix(as_edgelist(G))
  all.edges = permutations(n=217, 2)
  ind = sample(nrow(edges), round(nrow(edges)*0.2), replace = FALSE)
  ind2 = sample(nrow(all.edges), round(nrow(edges)*0.2), replace = FALSE)
  
  G =  delete_edges(G, edges = edges[ind,])
  G =  add_edges(G, edges = all.edges[ind2,] )
  Adj = as.matrix(get.adjacency(G))
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  ## add missingness 
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
  
  homo.direct.result[r,] = homo.direct(0.1, 0.2, 0.3, Adj)
  traffic.model.result[r,] = traffic.model(-0.5, 1, 0.5, traffic)
  diffusion.model.result[r,] = diff.randomwalk(p = 0.3, alpha0 = 0.1, alpha1 = 0.4, alpha2 = 0.3, time = 5, Adj)
  
  ## results
  diffusion.measure = diffusion(obs.Adj*0.3,  T = 5)
  #diffusion.measure = sample(1:N, N, replace = FALSE)
  random.measure = sample(1:N, N, replace = FALSE)
  
  centrality[[r]] = list(rowSums(Adj), betweenness(G), diffusion.measure, random.measure)
  homo.direct.result.mat[r,] = c(cor(homo.direct.result[r,], rowSums(obs.Adj), method = "spearman"),
                                 cor(homo.direct.result[r,], betweenness(obs.G), method = "spearman"),
                                 cor(homo.direct.result[r,],  diffusion.measure, method = "spearman"),
                                 cor(homo.direct.result[r,],  random.measure, method = "spearman"))
  traffic.model.result.mat[r,] = c(cor(traffic.model.result[r,], rowSums(obs.Adj), method = "spearman"),
                                   cor(traffic.model.result[r,], betweenness(obs.G), method = "spearman"),
                                   cor(traffic.model.result[r,], diffusion.measure, method = "spearman"),
                                   cor(traffic.model.result[r,],  random.measure, method = "spearman"))
  diffusion.model.result.mat[r,] = c(cor(diffusion.model.result[r,], rowSums(obs.Adj), method = "spearman"),
                                     cor(diffusion.model.result[r,], betweenness(obs.G), method = "spearman"),
                                     cor(diffusion.model.result[r,], diffusion.measure, method = "spearman"),
                                     cor(diffusion.model.result[r,], random.measure, method = "spearman"))
}



save(homo.direct.result.mat, file = "data/homo.direct.result.mat_missing.RData")
save(traffic.model.result.mat, file = "data/traffic.model.result.mat_missing.RData")
save(diffusion.model.result.mat, file = "data/diffusion.model.result.mat_missing.RData")