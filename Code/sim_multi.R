library(igraph)
library(keyplayer)
library(xtable)
library(Matrix)
library(combinat)
####
graph = make_graph(~ 1-2-5-9, 1-2-3-4-9, 1-2-6, 1-7-3-4, 1-7-6, 
                    10-1-7, 8)
plot(graph)
ind = as.character(c(1:10))
tmp.Adj = get.adjacency(graph)
Adj = as.matrix(get.adjacency(graph))[match(ind, row.names(tmp.Adj)), match(ind, row.names(tmp.Adj))]

## centrality measures for each of teh nodes
mat = cbind(degree(graph)[match(ind, row.names(tmp.Adj))],
            betweenness(graph)[match(ind, row.names(tmp.Adj))],
            diffusion(Adj*0.3,  T = 2),
            diffusion(Adj*0.3,  T = 3),
            diffusion(Adj*0.3,  T = 5))
rownames(mat) =  c(1:10)
colnames(mat) = c("Degree", "Betweenness", "Diffusion (T=2)", "Diffusion (T=3)", "Diffusion (T=5)")

print(xtable(mat, digits = 3))

## the ranking based on the five different centralities
order.mat = cbind(order(mat[,1], decreasing = TRUE),
                  order(mat[,2], decreasing = TRUE),
                  order(mat[,3], decreasing = TRUE),
                  order(mat[,4], decreasing = TRUE),
                  order(mat[,5], decreasing = TRUE))
rownames(order.mat) =  c(1:10)
colnames(order.mat) = c("Degree", "Betweenness", "Diffusion (T=2)", "Diffusion (T=3)", "Diffusion (T=5)")

print(xtable(order.mat, digits = 0))


###### selecting two most influential subjects ######
homo.direct.two = function(alpha0, alpha1, alpha2, Adj){
  N = nrow(Adj); M = N*(N-1)/2
  influence.Y = rep(NA, M)  
  ind = combn(N, 2)
  for(i in 1:M){
    treat = rep(0, N); treat[ind[,i]] = 1
    prop = alpha0 + alpha1*treat + alpha2*as.numeric(Adj %*% treat)/ max(rowSums(Adj))
    Y = prop
    influence.Y[i] = mean(Y) 
  }
  return(influence.Y)
}

traffic.model.two = function(alpha0, alpha1, alpha2, traffic){
  N = nrow(traffic); M = N*(N-1)/2
  influence.Y = rep(NA, M)
  ind = combn(N,2)  
  for(j in 1:M){
    treat = rep(0, N); treat[ind[,j]] = 1; prop = c()
    for(i in 1:N){
      prop[i] = plogis(alpha0 + alpha1*treat[i] + alpha2*sum(as.numeric(traffic[i,,] %*% treat)[-i]))
    }
    Y = prop
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y)
}


diff.randomwalk.two = function(p, alpha0, alpha1, alpha2, time, Adj){
  # Adj : adjacency matrix
  # p: probability of adoption
  # alpha0 : probability of being active under treatment but without any influence
  # alpha1 : probability of being active under control but without any influence
  # alpha2 : probability of being active under control and having one active friends.
  # time : time spent under random walk
  N = nrow(Adj); M = N*(N-1)/2; ind = combn(N,2)  
  out = matrix(0, N, time+1); influence.Y = rep(NA, M)
  
  for(i in 1:M){
    treat = rep(0, N); treat[ind[,i]] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = alpha2*(propmat[ind[1,i], ] + propmat[ind[2,i], ]) 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*(propmat[ind[1,i], ] + propmat[ind[2,i], ]) 
    }
    influence.Y[i] = mean(out[,t])
  }
  
  return(influence.Y)
}


n.rep = 500
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
homo.direct.order = traffic.model.order = diffusion.model.order = matrix(NA, n.rep, 45)

for(r in 1:n.rep){
  #print(r)
  set.seed(r)
  graph = make_graph(~ 1-2-5-9, 1-2-3-4-9, 1-2-6, 1-7-3-4, 1-7-6, 
                      10-1-7, 8)
  plot(graph)
  ind = as.character(c(1:10))
  tmp.Adj = get.adjacency(graph)
  Adj = as.matrix(get.adjacency(graph))[match(ind, row.names(tmp.Adj)), match(ind, row.names(tmp.Adj))]
  G = graph.adjacency(Adj)
  
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result = homo.direct.two(0.1, 0.2, 0.3, Adj)
  traffic.model.result = traffic.model.two(-0.5, 1, 0.5, traffic)
  diffusion.model.result = diff.randomwalk.two(p = 0.3, alpha0 = 0.1, alpha1 = 0.4, alpha2 = 0.3, time = 5, Adj)
  
  ## results 
  M = N*(N-1)/2; ind = combn(N, 2)
  result.mat = matrix(NA,M,3); row.names(result.mat) = rep("NA", M)
  for(i in 1:M){
    row.names(result.mat)[i] = paste0("(", ind[1,i], ",", ind[2,i], ")")
    result.mat[i,1] = homo.direct.result[i]
    result.mat[i,2] = traffic.model.result[i]
    result.mat[i,3] = diffusion.model.result[i]
  }
  
  homo.direct.order[r,] = order(result.mat[,1])
  traffic.model.order[r,] = order(result.mat[,2])
  diffusion.model.order[r,] = order(result.mat[,3])
  
}

summary.mat = cbind(names(result.mat[order(result.mat[,1], decreasing = TRUE) ,1])[1:7], 
                    names(result.mat[order(result.mat[,2], decreasing = TRUE) ,2])[1:7],
                    names(result.mat[order(result.mat[,3], decreasing = TRUE) ,3])[1:7])

###### selecting three most influential subjects ######
homo.direct.three = function(alpha0, alpha1, alpha2, Adj){
  N = nrow(Adj); M = N*(N-1)*(N-2)/6
  influence.Y = rep(NA, M)  
  ind = combn(N, 3)
  for(i in 1:M){
    treat = rep(0, N); treat[ind[,i]] = 1
    prop = alpha0 + alpha1*treat + alpha2*as.numeric(Adj %*% treat)/ max(rowSums(Adj))
    Y = prop
    influence.Y[i] = mean(Y) 
  }
  return(influence.Y)
}


traffic.model.three = function(alpha0, alpha1, alpha2, traffic){
  N = nrow(traffic); M = N*(N-1)*(N-2)/6
  influence.Y = rep(NA, M)
  ind = combn(N,3)  
  for(j in 1:M){
    treat = rep(0, N); treat[ind[,j]] = 1; prop = c()
    for(i in 1:N){
      prop[i] = plogis(alpha0 + alpha1*treat[i] + alpha2*sum(as.numeric(traffic[i,,] %*% treat)[-i]))
    }
    Y = prop
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y)
}


diff.randomwalk.three = function(p, alpha0, alpha1, alpha2, time, Adj){
  # Adj : adjacency matrix
  # p: probability of adoption
  # alpha0 : probability of being active under treatment but without any influence
  # alpha1 : probability of being active under control but without any influence
  # alpha2 : probability of being active under control and having one active friends.
  # time : time spent under random walk
  N = nrow(Adj); M = N*(N-1)*(N-2)/6; ind = combn(N,3)  
  out = matrix(0, N, time+1); influence.Y = rep(NA, M)
  
  for(i in 1:M){
    treat = rep(0, N); treat[ind[,i]] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = alpha2*(propmat[ind[1,i], ] + propmat[ind[2,i], ]) 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*(propmat[ind[1,i], ] + propmat[ind[2,i], ]) 
    }
    influence.Y[i] = mean(out[,t])
  }
  
  return(influence.Y)
}


n.rep = 500
homo.direct.result.mat = traffic.model.result.mat = diffusion.model.result.mat = matrix(NA, n.rep, 4)
homo.direct.order = traffic.model.order = diffusion.model.order = matrix(NA, n.rep, 120)

for(r in 1:n.rep){
  set.seed(r)
  graph <- make_graph(~ 1-2-5-9, 1-2-3-4-9, 1-2-6, 1-7-3-4, 1-7-6, 
                      10-1-7, 8)
  plot(graph)
  ind = as.character(c(1:10))
  tmp.Adj = get.adjacency(graph)
  Adj = as.matrix(get.adjacency(graph))[match(ind, row.names(tmp.Adj)), match(ind, row.names(tmp.Adj))]
  G = graph.adjacency(Adj)
  
  
  traffic = traffic.array(Adj)
  N = nrow(Adj)
  
  homo.direct.result = homo.direct.three(0.1, 0.2, 0.3, Adj)
  traffic.model.result = traffic.model.three(-0.5, 1, 0.5, traffic)
  diffusion.model.result = diff.randomwalk.three(p = 0.3, alpha0 = 0.1, alpha1 = 0.4, alpha2 = 0.3, time = 5, Adj)
  
  ## results 
  M = N*(N-1)*(N-2)/6; ind = combn(N, 3)
  result.mat = matrix(NA,M,3); row.names(result.mat) = rep("NA", M)
  for(i in 1:M){
    row.names(result.mat)[i] = paste0("(", ind[1,i], ",", ind[2,i], ",", ind[3,i], ")")
    result.mat[i,1] = homo.direct.result[i]
    result.mat[i,2] = traffic.model.result[i]
    result.mat[i,3] = diffusion.model.result[i]
  }
  
  homo.direct.order[r,] = order(result.mat[,1])
  traffic.model.order[r,] = order(result.mat[,2])
  diffusion.model.order[r,] = order(result.mat[,3])
  
}


summary.mat = cbind(names(result.mat[order(result.mat[,1], decreasing = TRUE) ,1])[1:7], 
                    names(result.mat[order(result.mat[,2], decreasing = TRUE) ,2])[1:7],
                    names(result.mat[order(result.mat[,3], decreasing = TRUE) ,3])[1:7])