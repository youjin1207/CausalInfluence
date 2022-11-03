#######################################################
############ Default setting ##########################
#######################################################

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
    out[,2] = alpha2*propmat[i, ] + alpha1*out[, 1] 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*propmat[i,] + alpha1*out[, t-1] 
    }
    influence.Y[i] = mean(out[,t])
  }
  
  return(influence.Y)
}


#######################################################
############ Sensitivity analysis######################
#######################################################

#### heterogeneous alphas ####
homo.direct.alphas = function(alpha0, alpha1, alpha2, Adj, alphas){
  N = nrow(Adj)
  influence.Y = rep(NA, N)  
  alpha.unif = runif(N, alpha1 - alphas, alpha1 + alphas)
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    prop = alpha0 + alpha.unif*treat + alpha2*as.numeric(Adj %*% treat)/ max(rowSums(Adj))
    Y = prop
    influence.Y[i] = mean(Y) 
  }
  return(influence.Y)
}

traffic.model.alphas = function(alpha0, alpha1, alpha2, traffic, alphas){
  N = nrow(traffic)
  influence.Y = rep(NA, N)  
  alpha.unif = runif(N, alpha1 - alphas, alpha1 + alphas)
  for(j in 1:N){
    treat = rep(0, N); treat[j] = 1; prop = c()
    for(i in 1:N){
      prop[i] = alpha0 + alpha.unif[i]*treat[i] + alpha2*sum(as.numeric(traffic[i,,] %*% treat)[-i])
    }
    Y = prop
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y)
}

diff.randomwalk.alphas = function(p, alpha0, alpha1, alpha2, time, Adj, alphas){
  N = nrow(Adj); out = matrix(0, N, time+1); influence.Y = rep(NA, N)
  alpha.unif = runif(N, alpha1 - alphas, alpha1 + alphas)
  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = alpha2*propmat[i, ] + (1-p)*out[, 1] 
    out[i,2] = alpha.unif[i] + (1-p)*out[i, 1] 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*propmat[i,] + (1-p)*out[, t-1] 
    }
    influence.Y[i] = mean(out[,t])
  }
  return(influence.Y)
}

#### heterogeneous betas ####
homo.direct.betas = function(alpha0, alpha1, alpha2, Adj, betas){
  N = nrow(Adj)
  influence.Y = rep(NA, N)  
  beta.unif = runif(N, alpha2 - betas, alpha2 + betas)
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    prop = alpha0 + alpha1*treat + beta.unif*as.numeric(Adj %*% treat)/ max(rowSums(Adj))
    Y = prop
    influence.Y[i] = mean(Y) 
  }
  return(influence.Y)
}

traffic.model.betas = function(alpha0, alpha1, alpha2, traffic, betas){
  N = nrow(traffic)
  influence.Y = rep(NA, N)  
  beta.unif = runif(N, alpha2 - betas, alpha2 + betas)
  for(j in 1:N){
    treat = rep(0, N); treat[j] = 1; prop = c()
    for(i in 1:N){
      prop[i] = alpha0 + alpha1*treat[i] + beta.unif[i]*sum(as.numeric(traffic[i,,] %*% treat)[-i])
    }
    Y = prop
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y)
}

diff.randomwalk.betas = function(p, alpha0, alpha1, alpha2, time, Adj, betas){
  N = nrow(Adj); out = matrix(0, N, time+1); influence.Y = rep(NA, N)
  beta.unif = runif(N, alpha2 - betas, alpha2 + betas)
  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = beta.unif*propmat[i, ] + (1-p)*out[, 1] 
    out[i,2] = alpha1 + (1-p)*out[i, 1] 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = beta.unif*propmat[i,] + (1-p)*out[, t-1] 
    }
    influence.Y[i] = mean(out[,t])
  }
  
  return(influence.Y)
}

#### adding random errors in the outcome model ####
homo.direct.sigma = function(alpha0, alpha1, alpha2, Adj, sigma){
  N = nrow(Adj)
  influence.Y = rep(NA, N)  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    prop = alpha0 + alpha1*treat + alpha2*as.numeric(Adj %*% treat)/ max(rowSums(Adj))
    Y = prop + rnorm(N, 0, sigma)
    influence.Y[i] = mean(Y) 
  }
  return(influence.Y)
}

traffic.model.sigma = function(alpha0, alpha1, alpha2, traffic, sigma){
  N = nrow(traffic)
  influence.Y = rep(NA, N)  
  for(j in 1:N){
    treat = rep(0, N); treat[j] = 1; prop = c()
    for(i in 1:N){
      prop[i] = alpha0 + alpha1*treat[i] + alpha2*sum(as.numeric(traffic[i,,] %*% treat)[-i])
    }
    Y = prop + rnorm(N, 0, sigma)
    influence.Y[j] = mean(Y) 
  }
  return(influence.Y)
}

diff.randomwalk.sigma = function(p, alpha0, alpha1, alpha2, time, Adj, sigma){
  N = nrow(Adj); out = matrix(0, N, time+1); influence.Y = rep(NA, N)
  
  for(i in 1:N){
    treat = rep(0, N); treat[i] = 1
    out[,1] = rep(alpha0, N)
    propmat = p*Adj; 
    out[,2] = alpha2*propmat[i, ] + (1-p)*out[, 1] + rnorm(N, 0, sigma)
    out[i,2] = alpha1 + (1-p)*out[i, 1] + rnorm(1, 0, sigma) 
    for(t in 3:(time+1)){
      propmat = propmat + propmat %*% (p * Adj)
      out[,t] = alpha2*propmat[i,] + (1-p)*out[, t-1] + rnorm(N, 0, sigma)
    }
    influence.Y[i] = mean(out[,t])
  }
  
  return(influence.Y)
}


