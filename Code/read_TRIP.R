pkgs <- c("igraph", "gplots", "RColorBrewer", "keyplayer")
install.packages(pkgs)
library(igraph); library(gplots); library(RColorBrewer); library(keyplayer)

#########
# load("Data/REALDATA.Rdata")
# We use the Transmission Reduction Intervention Project (TRIP) data, 
# and this data is currently not publicly available due to privacy or ethical restrictions. 
# The data are available on request from Dr.Georgios Nikolopoulos.  
head(private_edges)
head(private_nodes)
c(nrow(nodes), length(unique(nodes$EGO.ID))) # a total of 356 (unique) participants

# the number of nodes: 356
# the number of edges: 829
net = graph_from_data_frame(d = private_edges[,1:2], vertices = nodes$EGO.ID, directed = F)

sum(degree(net) == 0) # 79 isolated nodes

# the number of nodes: 277 (356-79)
# the number of edges: 829
net0 = delete.vertices(net, degree(net) == 0)

# analysis network
# the number of nodes: 277 (356-79)
# the number of edges: 542
net0 = simplify(net0, remove.multiple = T, remove.loops = T) 
n = length(V(net0))
## read covariates ##
V(net0)$egoid = names(V(net0)) # nodes id 

V(net0)$age = V(net0)$gender = V(net0)$education = V(net0)$employ = V(net0)$alert = 
  V(net0)$share.x = V(net0)$posneg = V(net0)$share.y =  V(net0)$hiv = rep(NA, n)

for(i in 1:n){
  V(net0)$age[i] = nodes$age.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$gender[i] = nodes$gender.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$education[i] = nodes$b03_1.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$employ[i] = nodes$b04.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$alert[i] = nodes$i01.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$share.x[i] = nodes$share.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$share.y[i] = nodes$share.y[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]]
  V(net0)$hiv[i] = as.integer(nodes$posneg.x[as.character(nodes$EGO.ID) %in% V(net0)$egoid[i]] == "POS")
}
V(net0)$hiv = as.numeric(V(net0)$hiv)


data277 = data.frame(id = names(V(net0)), age = V(net0)$age, gender = V(net0)$gender, 
                     education = V(net0)$education, employ = V(net0)$employ,
                     alert = V(net0)$alert, share.x = V(net0)$share.x,
                     share.y = V(net0)$share.y, hiv = V(net0)$hiv)

table(data277$alert)

data277$alert = ifelse(data277$alert == 1, 1, 0)

summary(data277$age) 

table(data277$gender)
table(data277$gender) / nrow(data277)

table(data277$hiv)
table(data277$hiv) / nrow(data277)

table(data277$education)
data277$education2 = ifelse(data277$education == 1, 2, 
                            ifelse(data277$education == 6 | data277$education == 7, 5, data277$education))


table(data277$education2)
data277$education2=factor(data277$education2, levels = 2:5, ordered = T)

table(data277$employ)
data277$employ2 = ifelse(data277$employ == 1 | data277$employ == 2 | 
                           data277$employ == 3 | data277$employ == 4 | data277$employ == 5, 7, data277$employ)
data277$employ2 =as.factor(data277$employ2)

table(data277$share.x)
table(data277$share.y)



### descriptive table ####
summary(data277$age)
mat = matrix(NA, nrow = 14, ncol = 2)
mat[1,] = c(sum(data277$gender == 0), round(mean(data277$gender == 0)*100,1))
mat[2,] = c(sum(data277$hiv == 1), round(mean(data277$hiv == 1)*100,1))
## education
mat[3,] = c(sum(data277$education2 == 2), round(mean(data277$education2 == 2)*100, 1))
mat[4,] = c(sum(data277$education2 == 3), round(mean(data277$education2 == 3)*100, 1))
mat[5,] = c(sum(data277$education2 == 4), round(mean(data277$education2 == 4)*100, 1))
mat[6,] = c(sum(data277$education2 == 5), round(mean(data277$education2 == 5)*100, 1))
## employment
mat[7,] = c(sum(data277$employ2 == 7), round(mean(data277$employ2 == 7)*100, 1))
mat[8,] = c(sum(data277$employ2 == 8), round(mean(data277$employ2 == 8)*100, 1))
mat[9,] = c(sum(data277$employ2 == 9), round(mean(data277$employ2 == 9)*100, 1))
mat[10,] = c(sum(data277$employ2 == 10), round(mean(data277$employ2 == 10)*100, 1))
## sharing
mat[11,] = c(sum(data277$share.x == 1), round(mean(data277$share.x == 1)*100, 1))
mat[12,] = c(sum(data277$share.x > 90), round(mean(data277$share.x > 90)*100, 1))
mat[13,] = c(sum(data277$share.y[!is.na(data277$share.y)] == 1), 
             round(mean(data277$share.y[!is.na(data277$share.y)] == 1)*100, 1))
mat[14,] = c(sum(is.na(data277$share.y) == 1 |  data277$share.y > 90), 
             round(mean(is.na(data277$share.y) == 1 |  data277$share.y > 90)*100, 1))
print(mat)




### Figure 1 ###
set.seed(1234)
layout = layout.fruchterman.reingold(net0)
pdf("Figure/trip_network.pdf", width = 7, height = 7)
plot(net0, vertex.size = 5, vertex.color = "azure3", layout = layout, vertex.label = NA)
dev.off()

V(net0)$color = ifelse(V(net0)$alert == 1 , "firebrick", "azure3")
pdf("Figure/trip_network_alert.pdf", width = 7, height = 7)
plot(net0, vertex.size = 4.5, vertex.color = V(net0)$color, layout = layout, vertex.label = NA)
dev.off()


### centrality ###
n = length(V(net0))
A = as.matrix(get.adjacency(net0))

## calculate centrality measures
V(net0)$degree = rowSums(A)
V(net0)$between = betweenness(net0)
V(net0)$diffusion_03_T3 = diffusion(A*0.3, T= 3)
V(net0)$diffusion_03_T5 = diffusion(A*0.3, T= 5)
V(net0)$diffusion_03_T10 = diffusion(A*0.3, T= 10)

cor(V(net0)$degree, V(net0)$between)
cor(V(net0)$degree, V(net0)$diffusion_03_T5)
cor(V(net0)$between, V(net0)$diffusion_03_T5)

## degree centrality
Palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
V(net0)$degree.color = Palette[V(net0)$degree]

pdf("Figure/trip_network_degree.pdf", width = 7, height = 7)
par(cex.main = 2)
plot(net0, vertex.size = 5, vertex.color = V(net0)$degree.color, layout = layout, 
     main = "Degree centrality", vertex.label = NA)
dev.off()

## betweenness
Palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
cuts = c(0, quantile(V(net0)$between[V(net0)$between > 0], seq(0, 1, 1/18)))
V(net0)$between.color = Palette[1]
for(k in 1:19){
  V(net0)$between.color = ifelse(V(net0)$between >= cuts[k] & V(net0)$between <= cuts[(k+1)], 
                                 Palette[k], V(net0)$between.color)
}
pdf("Figure/trip_network_betweenness.pdf", width = 7, height = 7)
par(cex.main = 2)
plot(net0, vertex.size = 4.5, vertex.color = V(net0)$between.color, layout = layout, 
     main = "Betweenness centrality", vertex.label = NA)
dev.off()

## Diffusion centrality
Palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
cuts = c(0, quantile(V(net0)$diffusion_03_T3[V(net0)$diffusion_03_T3 > 0], seq(0, 1, 1/18)))
V(net0)$diffusion_03_T3.color = Palette[1]
for(k in 1:19){
  V(net0)$diffusion_03_T3.color = ifelse(V(net0)$diffusion_03_T3 >= cuts[k] & V(net0)$diffusion_03_T3 <= cuts[(k+1)], 
                                         Palette[k], V(net0)$diffusion_03_T3.color)
}

pdf("Figure/trip_network_diffusion_03_T3.pdf", width = 7, height = 7)
par(cex.main = 2)
plot(net0, vertex.size = 4.5, vertex.color = V(net0)$diffusion_03_T3.color, layout = layout, 
     main = "Diffusion centrality (p=0.3, T=3)", vertex.label = NA)
dev.off()
##
Palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
cuts = c(0, quantile(V(net0)$diffusion_03_T5[V(net0)$diffusion_03_T5 > 0], seq(0, 1, 1/18)))
V(net0)$diffusion_03_T5.color = Palette[1]
for(k in 1:19){
  V(net0)$diffusion_03_T5.color = ifelse(V(net0)$diffusion_03_T5 >= cuts[k] & V(net0)$diffusion_03_T5 <= cuts[(k+1)], 
                                         Palette[k], V(net0)$diffusion_03_T5.color)
}
pdf("Figure/trip_network_diffusion_03_T5.pdf", width = 7, height = 7)
par(cex.main = 2)
plot(net0, vertex.size = 4.5, vertex.color = V(net0)$diffusion_03_T5.color, layout = layout, 
     main = "Diffusion centrality (p=0.3, T=5)", vertex.label = NA)
dev.off()
##
Palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
cuts = c(-1, quantile(V(net0)$diffusion_03_T10[V(net0)$diffusion_03_T10 > 0], seq(0, 1, 1/18)))
V(net0)$diffusion_03_T10.color = Palette[1]
for(k in 1:19){
  V(net0)$diffusion_03_T10.color = ifelse(V(net0)$diffusion_03_T10 > cuts[k] & V(net0)$diffusion_03_T10 <= cuts[(k+1)], 
                                          Palette[k], V(net0)$diffusion_03_T10.color)
}
pdf("Figure/trip_network_diffusion_03_T10.pdf", width = 7, height = 7)
par(cex.main = 2)
plot(net0, vertex.size = 4.5, vertex.color = V(net0)$diffusion_03_T10.color, layout = layout, 
     main = "Diffusion centrality (p=0.3, T=10)", vertex.label = NA)
dev.off()

