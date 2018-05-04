library(igraph)
library(ggnetwork)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(statnet)
library(ergm)

# --------------------------------------------------------------------------------------------

# Importing data
drug_data = read.csv("DRUGNET.csv",header=TRUE,row.names=1,check.names=FALSE)
drug_attr = read.csv("DRUGATTR.csv",header=TRUE,row.names=1,check.names=FALSE)

# --------------------------------------------------------------------------------------------

m=as.matrix(drug_data) # coerces the data set as a matrix

# Creating network
net = as.network(x = m, directed = TRUE, 
                 loops = FALSE, matrix.type = "adjacency")

Gender = drug_attr$Gender
for(i in 1:293){
  if(Gender[i] == 0) {Gender[i] = "Unknown"}
  else if (Gender[i] == 1) {Gender[i] = "Male"}
  else{Gender[i] = "Female"}
}
net = set.vertex.attribute(net, "Gender", Gender)

Ethnicity = drug_attr$Ethnicity
for(i in 1:293){
  if(Ethnicity[i] == 1) {Ethnicity[i] = "White"}
  else if (Ethnicity[i] == 2) {Ethnicity[i] = "African American"}
  else{Ethnicity[i] = "Latino"}
}
net = set.vertex.attribute(net, "Ethnicity", Ethnicity)

summary.network(net, print.adj = FALSE)

num_nodes = 293
one = get.vertex.attribute(net,"Gender")=="Male"
table(one) ["TRUE"]
# 200 (Male)
sec = get.vertex.attribute(net,"Gender")=="Female"
table(sec) ["TRUE"]
# 86 (Female)
th = get.vertex.attribute(net,"Gender")=="Unknown"
table(th) ["TRUE"]
# 7 (Unknown)

wh = get.vertex.attribute(net,"Ethnicity")=="White"
table(wh) ["TRUE"]
# 39 (White or others)
af = get.vertex.attribute(net,"Ethnicity")=="African American"
table(af) ["TRUE"]
# 99
lat = get.vertex.attribute(net,"Ethnicity")=="Latino"
table(lat) ["TRUE"]
# 155

node_names <- rep("",num_nodes)
node_title <- rep("",num_nodes)
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Gender")[i] == "Male"){
    node_names[i] <- "Male"
    node_title[i] = i
  }else if(get.vertex.attribute(net,"Gender")[i] == "Female"){
    node_names[i] <- "Female"
    node_title[i] = i
  }
  else{
    node_names[i] <- "Unknown"
    node_title[i] = i
  }
}
print(node_names)
print(node_title)

node_eth <- rep("",num_nodes)
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Ethnicity")[i] == "White"){
    node_eth[i] <- "White or other"
  }else if(get.vertex.attribute(net,"Ethnicity")[i] == "African American"){
    node_eth[i] <- "African American"
  }
  else{
    node_eth[i] <- "Latino"
  }
}
print(node_eth)


# --------------------------------------------------------------------------------------------


# Default Network

jpeg("Network_1.jpeg", width = 12, height = 8, units = 'in', res = 1200)

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       size = 6, edge.size = 0.5, edge.color = "black", node.color = "orange",
       label = node_title, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk1",0.2) ))

dev.off()

# --------------------------------------------------------------------------------------------

# Network components
## Weakly connected component
largest_weak = component.largest(net, connected = "weak")
node_lw <- rep("",num_nodes)
for(i in 1:num_nodes){
  if(largest_weak[i]){
    node_lw[i] <- "1"
    }
  else{
    node_lw[i] <- "0"
  }
}
print(node_lw)

k = 0
l = 0
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Gender")[i] == "Female" & largest_weak[i]){
    k = k + 1
  }
  else if(get.vertex.attribute(net,"Gender")[i] == "Male" & largest_weak[i]){
    l = l + 1
  }
}
# Female = 43
# Male = 148 

w = 0
a = 0
l = 0
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Ethnicity")[i] == "White" & largest_weak[i]){
    w = w + 1
  }
  else if(get.vertex.attribute(net,"Ethnicity")[i] == "African American" & largest_weak[i]){
    a = a + 1
  }
  else if(get.vertex.attribute(net,"Ethnicity")[i] == "Latino" & largest_weak[i]){
    l = l + 1
  }
}
# White = 15
# Aferican American = 69
# Latino = 109

jpeg("Network_2.jpeg", width = 12, height = 8, units = 'in', res = 1200)

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       size = 6, edge.size = 0.5, edge.color = "black", color = node_lw, palette = c("0" = "orange",
                                                                                     "1" = "maroon"),
       label = node_title, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (Weakly Connected Component)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk1",0.2) ), legend.position ="none")
dev.off()

## Strongly connected component
largest_strong = component.largest(net, connected = "strong")
node_ls <- rep("",num_nodes)
for(i in 1:num_nodes){
  if(largest_strong[i]){
    node_ls[i] <- "1"
  }
  else{
    node_ls[i] <- "0"
  }
}
print(node_ls)

k = 0
l = 0
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Gender")[i] == "Female" & largest_strong[i]){
    k = k + 1
  }
  else if(get.vertex.attribute(net,"Gender")[i] == "Male" & largest_strong[i]){
    l = l + 1
  }
}
# Female = 3
# Male = 24 

w = 0
a = 0
l = 0
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Ethnicity")[i] == "White" & largest_strong[i]){
    w = w + 1
  }
  else if(get.vertex.attribute(net,"Ethnicity")[i] == "African American" & largest_strong[i]){
    a = a + 1
  }
  else if(get.vertex.attribute(net,"Ethnicity")[i] == "Latino" & largest_strong[i]){
    l = l + 1
  }
}
# White = 0
# Aferican American = 0
# Latino = 27

jpeg("Network_3.jpeg", width = 12, height = 8, units = 'in', res = 1200)

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       size = 6, edge.size = 0.5, edge.color = "black", color = node_ls, palette = c("0" = "orange",
                                                                                     "1" = "maroon"),
       label = node_title, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (Strongly Connected Component)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk1",0.2) ), legend.position ="none")
dev.off()

# --------------------------------------------------------------------------------------------
# Network with Gender

jpeg("Network_4.jpeg", width = 12, height = 8, units = 'in', res = 1200)

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       size = 6, color = node_names, palette = c("Male" = "skyblue",
                                        "Female" = "deeppink",
                                        "Unknown" = "gold"), 
       edge.size = 0.5, edge.color = "black",
       label = node_title, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (Gender)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk1",0.2) ))
dev.off()

# --------------------------------------------------------------------------------------------

# Network with Ethnicity

cut_p = cutpoints(net)
node_cp <- rep("",num_nodes)
for(i in 1:num_nodes){
  if(i %in% cut_p){
    node_cp[i] <- "1"
  }
  else{
    node_cp[i] <- "0"
  }
}
print(node_cp)

jpeg("Network_5.jpeg", width = 12, height = 8, units = 'in', res = 1200)

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       size = 6, color = node_eth, palette = c("White or other" = "grey",
                                                 "African American" = "forestgreen",
                                                 "Latino" = "coral"), 
       edge.size = 0.5, edge.color = "black", #shape = node_cp,
       label = node_title, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (Ethnicity)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk",0.2) ))

dev.off()
  
# --------------------------------------------------------------------------------------------

# Network with Gender and Ethnicity
jpeg("Network_6.jpeg", width = 12, height = 8, units = 'in', res = 1200)

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       size = 6, color = node_eth, palette = c("White or other" = "grey",
                                               "African American" = "forestgreen",
                                               "Latino" = "coral"), 
       edge.size = 0.5, edge.color = "black", shape = node_names,
       label = node_title, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (Gender and Ethnicity)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk",0.2) ))

dev.off()

# --------------------------------------------------------------------------------------------

# Network with Ethiniicty and Gender (Indegree)
node_nt <- rep("",num_nodes)
for(i in 1:num_nodes){
  if(get.vertex.attribute(net,"Gender")[i] == "Male"){
    node_nt[i] <- "M"
  }else if(get.vertex.attribute(net,"Gender")[i] == "Female"){
    node_nt[i] <- "F"
  }
  else{
    node_nt[i] <- "U"
  }
}
print(node_nt)

jpeg("Network_7.jpeg", width = 12, height = 8, units = 'in', res = 1200) 

set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
     color = node_eth, palette = c("White or other" = "grey",
                                               "African American" = "forestgreen",
                                               "Latino" = "coral"), 
       edge.size = 0.5, edge.color = "black",size = "indegree",
       label = node_nt, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (In-degree)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk",0.2) ))

dev.off()

# Interactive network plot

library(igraph)
library(magrittr)
#install.packages("visNetwork")
library(visNetwork)
library(data.table)

graph =  graph_from_adjacency_matrix(m, mode = "directed")
graph = simplify(graph)

V(graph)$indegree = centr_degree(graph, mode = "in")$res

nodes = get.data.frame(graph, what = "vertices")
nodes = data.frame(id = nodes$name, title = nodes$name,
                   group = nodes$indegree, indegree = nodes$indegree)
setnames(nodes, "indegree", "in-degree centrality")
nodes = nodes[order(nodes$id, decreasing = F),]

edges = get.data.frame(graph, what="edges")[1:2]
plot = visNetwork(nodes, edges, height = "500px", width = "100%",
           main = "In-degree centrality") %>%
  visOptions(selectedBy = "in-degree centrality", highlightNearest = TRUE, nodesIdSelection = TRUE)%>%
  visPhysics(stabilization = FALSE) %>%
  visEdges(arrows = "to") 

visSave(plot, file = "indegree.html", selfcontained = TRUE, background = "white")

# --------------------------------------------------------------------------------------------
indeg = degree(net,cmode="indegree")
qplot(indeg,
      geom="histogram",
      binwidth = 0.5,
      main = "Indegree distribution",
      xlab = "Indegree", ylab = "Frequency")
# Network with Ethiniicty and Gender (Outdegree)
set.seed(1)
ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), 
       color = node_eth, palette = c("White or other" = "yellow",
                                     "African American" = "forestgreen",
                                     "Latino" = "coral"), 
       edge.size = 0.5, edge.color = "black",size = "outdegree",
       label = node_nt, label.color = "black", label.size = 3,
       arrow.size = 8, arrow.gap = 0.004) + 
  ggtitle("Drugnet (Outdegree)") + 
  theme(panel.background = element_rect(fill = alpha("cornsilk",0.2) ))

outdeg <- degree(net,cmode="outdegree")
qplot(outdeg,
      geom="histogram",
      binwidth = 0.5,
      main = "Outdegree distribution",
      xlab = "Outdegree", ylab = "Frequency")

library(igraph) # This loads the igraph package
g=graph.adjacency(m,mode="directed",weighted=NULL) # this will create an 'igraph object'
g 

# --------------------------------------------------------------------------------------------

### ERGM 

# Bernoulli model
set.seed(1)
m1 = ergm(net ~ edges)
summary(m1)

# m1 = ergm(net ~ edges + mutual + sender(base = 1:2) + receiver(base = 1:2))

#                                                                          
summary(m1)


# Simulation
sim = simulate(m1, burnin = 1e+6, verbose = TRUE, seed = 9)

par(mfrow=c(1,2)) 
plot(net,vertex.col = "Gender", main = "Observed Network")
plot(sim, vertex.col = "Gender", main = "Simulated Network")
dev.off()


legend("topleft", legend=levels(group), pch=16, col=unique(group))

mixingmatrix(sim, "Gender")
mixingmatrix(net, "Gender")

mixingmatrix(sim, "Ethnicity")
mixingmatrix(net, "Ethnicity")

# In-degree plot
plot(summary(net~idegree(0:10)), type = "l", lty = 1, lwd = 2,xlab = "In-Degree", ylab = "Count", 
     main= "In-degree distribution", xaxt = "n")
axis(1, at=1:11, labels=0:10)
lines(summary(sim~idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,lty = 1:2)

# Out-degree plot
plot(summary(net~odegree(0:10)), type = "l", lty = 1, lwd = 2,xlab = "Out-Degree", ylab = "Count", 
     main= "Out-degree distribution", ylim = c(0, 120), xaxt = "n")
axis(1, at=1:11, labels=0:10)
lines(summary(sim~odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,lty = 1:2)

# Selected model
set.seed(1)
m2 = ergm(net ~ edges + mutual 
          + asymmetric("Gender", diff = FALSE) + asymmetric("Ethnicity", diff = FALSE) 
          + nodematch("Ethnicity", diff = F) + nodefactor("Ethnicity") 
          + simmelianties + isolates,
          MCMCsamplesize = 1e+5, interval = 1000)
summary(net ~ edges + mutual 
 + asymmetric("Gender", diff = FALSE) + asymmetric("Ethnicity", diff = FALSE) 
 + nodematch("Ethnicity") +nodefactor("Ethnicity") 
 + simmelianties + isolates)                                                                        
summary(m2)

# Simulation
sim_nets = simulate(m2, burnin = 1e+6, verbose = TRUE, seed = 9)

par(mfrow=c(1,2)) 
plot(net,vertex.col = "Gender", main = "Observed Network")
plot(sim_nets, vertex.col = "Gender", main = "Simulated Network")
dev.off()

mixingmatrix(sim_nets, "Gender")
mixingmatrix(net, "Gender")

mixingmatrix(sim_nets, "Ethnicity")
mixingmatrix(net, "Ethnicity")

# In-degree plot
plot(summary(net~idegree(0:10)), type = "l", lty = 1, lwd = 2,xlab = "In-Degree", ylab = "Count", 
    main= "In-degree distribution",  xaxt = "n")
axis(1, at=1:11, labels=0:10)
lines(summary(sim_nets~idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,lty = 1:2)

# Out-degree plot
plot(summary(net~odegree(0:10)), type = "l", lty = 1, lwd = 2,xlab = "Out-Degree", ylab = "Count", 
     main= "Out-degree distribution", xaxt = "n")
axis(1, at=1:11, labels=0:10)
lines(summary(sim_nets~odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,lty = 1:2)

gf = gof(m2 ~ model)
gf
plot(gf)

gf_1 <- gof(m2 ~ idegree + odegree + esp + distance)
gf_1
plot(gf_1)



# ---------------------------------------------------------------------------


# Model - 1 with edges
m1 = ergm(net~edges)
summary(m1)

m1 = ergm(net~edges + mutual("Gender")+mutual("Ethnicity"))
summary(m1)


m1 = ergm(net~edges + mutual + asymmetric("Gender", diff = FALSE) + asymmetric("Ethnicity", diff = FALSE)
          + nodefactor("Ethnicity", base = 2:4)
          + nodematch("Gender") + + nodematch("Ethnicity")
          + isolates)
#nodefactor("Gender", base = 2:3)    
# + nodefactor("Gender", base = 2:3)
summary(m1)

# In-degree plot
plot(summary(net~idegree(0:10)), type = "l", lty = 1, lwd = 2,xlab = "In-Degree", ylab = "Count", 
     main= "In-degree distribution")
lines(summary(sim_nets~idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,lty = 1:2)

# Out-degree plot
plot(summary(net~odegree(0:10)), type = "l", lty = 1, lwd = 2,xlab = "Out-Degree", ylab = "Count", 
     main= "Out-degree distribution")
lines(summary(sim_nets~odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,lty = 1:2)
# Model - 2 with edges and mutuality
m2 = ergm(net ~ edges + mutual+sender +receiver)
summary(m2)

# Model - 3 with edges + mutuality + nodematch 
m3 = ergm(net ~ edges +  mutual + nodematch('Gender', diff = F)
          + nodematch('Ethnicity', diff = F))
summary(m3)

# Model - 4 with edges + mutuality + nodemix
m5 = ergm(net ~ edges +  mutual + nodemix('Gender'))
summary(m5)

# Model - 4 with edges + mutuality + nodemix
m4 = ergm(net ~ edges +  mutual + nodefactor('Gender')
          + nodefactor('Ethnicity'))
summary(m4)
