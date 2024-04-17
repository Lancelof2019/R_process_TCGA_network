##### The melanoma network and network community detection
options(stringsAsFactors = F)
library(igraph)
library(visNetwork)

library(OmnipathR)#Get interactions data
library(igraph)#network analysis
library(ggraph)#network visualization
library(RColorBrewer)#Network visualization, heat map drawing

library(factoextra)#principal component analysis
library(FactoMineR)#principal component analysis

library(DESeq2)#differential expression analysis
library(SANTA)#Node Score Calculation

library(ComplexHeatmap)#heatmap drawing
library(stringr)#regular expression

library(randomForest)#Random Forest
library(caret)#Cross-validation
library(pROC)#AUC

library(clusterProfiler)#enrichment analysis
library(org.Hs.eg.db)
library(enrichplot)#enrichment analysis
library(msigdbr)#enrichment analysis
library(dplyr)

library(survival)#survival analysis
library(survminer)#survival analysis
library(glmnet)# cox regression
library(timeROC)#AUC

library(reshape2)#ggplot2 drawing

library(coin)#permutation test

library(networkD3)#sankey diagram

library(janitor)# clean data
library(here)# home directory setting
library(ggplot2) # draw figures
library(tidyr) # data processing
if (!require("visNetwork")) install.packages("visNetwork")
if (!require("igraph")) install.packages("igraph")
library(visNetwork)
library(igraph)
### The melanoma network
 cancer_type<-"SKCM"
 setwd("C:/Users/gklizh/Documents/Workspace/fan_project/code_data")
 edgetest<-read.csv(paste0("Raw data/network_merge/edge/", cancer_type, ".csv"))
 vertextest<-read.csv(paste0("Raw data/network_merge/vertex/", cancer_type, ".csv"))


 #edge_gene<-edgetest[[1]]

 #edge_gene2<-edgetest[[2]]

 #vertex_gene<-(vertextest[[1]])




namePoint<-c("RARRES3" ,  "ACPP"  ,    "HIST1H1E" , "SEPT5"  ,   "H2AFZ"  ,   "WISP2"  ,   "WARS"  ,    "HIST1H2AI")
namechange<-c("PLAAT4",   "ACP3",      "H1-4",  "SEPTIN5" ,    "H2AZ1"  ,   "CCN5"  ,   "WARS1"   ,   "H2AC13")

name_map <- setNames(namechange, namePoint)
#gg<-name_map["RARRES3"]
#print(gg[[1]])

# index_edge_gene<-match(namePoint,edgetest[[1]])
# print(index_edge_gene)
# #print(edge_gene[index_edge_gene])
# 
# index_edge_gene2<-match(namePoint,edgetest[[2]])
# #print(edge_gene2[index_edge_gene2])
# #print(edge_gene)
# print(index_edge_gene2)
# index_vertex_gene<-match(namePoint,vertextest[[1]])
# #print(vertex_gene[index_vertex_gene])
# print(index_vertex_gene)

#######################################################

matches <- edgetest[[1]] %in% namePoint
index_edge_gene <- which(matches)
#print(index_edge_gene)


matches <- edgetest[[2]] %in% namePoint
index_edge_gene2 <- which(matches)
#print(index_edge_gene2)
# 打印出所有匹配的索引

matches <- vertextest[[1]] %in% namePoint
index_vertex_gene<-which(matches)
#print(index_vertex_gene)

#print(edgetest[[1]][index_edge_gene])

#test1<-edgetest[[1]][index_edge_gene]

#print(test1)

#print(index_edge_gene[seq_along(index_edge_gene)[2]])
#match(test1,edgetest[[1]])
#print(index_edge_gene[1])



# for(i in seq_along(index_edge_gene)){
#   
#   original_name <- edgetest[[1]][index_edge_gene[i]]
#   
#   if (original_name %in% names(name_map)) {
#     edgetest[[1]][index_edge_gene[i]] <- name_map[[original_name]]
#   }
#   
# }
# 
# for(i in seq_along(index_edge_gene)){
#   
#   original_name <- edgetest[[1]][index_edge_gene[i]]
#   
#   if (original_name %in% names(name_map)) {
#     edgetest[[1]][index_edge_gene[i]] <- name_map[[original_name]]
#   }
#   
# }

edgetest[[1]][index_edge_gene]<-unname(name_map[edgetest[[1]][index_edge_gene]])
edgetest[[2]][index_edge_gene2]<-unname(name_map[edgetest[[2]][index_edge_gene2]])
vertextest[[1]][index_vertex_gene]<-unname(name_map[vertextest[[1]][index_vertex_gene]])




###############################################################################
# valid_indices_index_edge_gene <- !is.na(index_edge_gene)
# #edgetest[[1]][index_edge_gene[valid_indices_index_edge_gene]]
# 
# valid_indices_index_edge_gene2 <- !is.na(index_edge_gene2)
# #edgetest[[2]][index_edge_gene2[valid_indices_index_edge_gene2]]
# 
# valid_indices_index_vertex_gene <- !is.na(index_vertex_gene)
# #vertextest[[1]][index_vertex_gene[valid_indices_index_vertex_gene]]
# 
# edgetest[[1]][index_edge_gene[valid_indices_index_edge_gene]]<-namechange[valid_indices_index_edge_gene]
# edgetest[[2]][index_edge_gene2[valid_indices_index_edge_gene2]]<-namechange[valid_indices_index_edge_gene2]
# vertextest[[1]][index_vertex_gene[valid_indices_index_vertex_gene]]<-namechange[valid_indices_index_vertex_gene]
# 
# 
# match(namechange,edgetest[[1]])
# 
# match(namechange,edgetest[[2]])
# 
# 
# match(namechange,vertex_gene[[1]])
# 
# 
# match(namePoint,edgetest[[1]])
# 
# match(namePoint,edgetest[[2]])
# 
# 
# match(namePoint,vertex_gene[[1]])

#edgetest<-edgetest
#vertextest<-vertextest
#write.csv(edgetest,paste0("Raw data/network_merge/edge/", cancer_type, "_edge_filter03.csv"),row.names = FALSE)
#write.csv(vertextest,paste0("Raw data/network_merge/vertex/", cancer_type, "_vertex_filter03.csv"),row.names = FALSE)

interactions_E <- edgetest
interactions_V <- vertextest

graph_comp <- graph_from_data_frame(interactions_E, directed = TRUE, vertices = interactions_V)
graph_comp <- decompose(graph_comp, min.vertices = 15)[[1]]
graph_comp <- set_vertex_attr(graph_comp, "degree", index = V(graph_comp)$name, degree(graph_comp))

print(class(degree(graph_comp)))
print(length(degree(graph_comp)))
print(class(V(graph_comp)$log2FoldChange))
print(length(V(graph_comp)$log2FoldChange))

# 检查所有源顶点是否都在顶点数据框中
# all(edgetest$source_genesymbol %in% vertextest$HGNC_ID) # 'name'应该是vertextest中包含顶点名称的列
# setdiff(edgetest$source_genesymbol, vertextest$HGNC_ID)
# 
# 
# # 检查所有目标顶点是否都在顶点数据框中
# all(edgetest$target_genesymbol %in% vertextest$HGNC_ID)
# 
# setdiff(edgetest$target_genesymbol, vertextest$HGNC_ID)
#
# g <- decompose(g, min.vertices = 15)[[1]]
# g <- set_vertex_attr(g, "degree", index = V(g)$name, degree(g))
#print(match(namePoint,edge_gene))

#print(match(namePoint,edge_gene2))


#print(match(namePoint,vertex_gene))

# network_graph <- function(data_name)
# {
#   setwd("C:/Users/gklizh/Documents/Workspace/fan_project/code_data/")
#   
#   interactions_E <- read.csv(paste0("Raw data/network_merge/edge/", data_name, ".csv"))
#   interactions_V <- read.csv(paste0("Raw data/network_merge/vertex/", data_name, ".csv"))
#   ######################################################################################
#   
#  #edgetest<-read.csv(paste0("Raw data/network_merge/edge/", data_name, ".csv"))
#  #vertextest<-read.csv(paste0("Raw data/network_merge/vertex/", data_name, ".csv"))
#   
#   
#   #edge_gene<-edgetest[[1]]
#   
#  # edge_gene2<-edgetest[[2]]
#   
#   #vertex_gene<-(vertextest[[1]])
#   
#   
#   
#   
#  # namePoint<-c("RARRES3" ,  "ACPP"  ,    "HIST1H1E" , "SEPT5"  ,   "H2AFZ"  ,   "WISP2"  ,   "WARS"  ,    "HIST1H2AI")
#  # namechange<-c("PLAAT4",   "ACP3",      "H1-4",  "SEPTIN5" ,    "H2AZ1"  ,   "CCN5"  ,   "WARS1"   ,   "H2AC13")
#   
#   #index_edge_gene<-match(namePoint,edge_gene)
#   #print(edge_gene[index_edge_gene])
#   
#  # index_edge_gene2<-match(namePoint,edge_gene2)
#   #print(edge_gene2[index_edge_gene2])
#   
#   
#   #index_vertex_gene<-match(namePoint,vertex_gene)
#   #print(vertex_gene[index_vertex_gene])
#   
#   #valid_indices_index_edge_gene <- !is.na(index_edge_gene)
#   #valid_indices_index_edge_gene2 <- !is.na(index_edge_gene2)
#  # valid_indices_index_vertex_gene <- !is.na(index_vertex_gene)
#   
#   #edgetest[[1]][index_edge_gene[valid_indices_index_edge_gene]]
#   #edgetest[[2]][index_edge_gene2[valid_indices_index_edge_gene2]]
#  # vertextest[[1]][index_vertex_gene[valid_indices_index_vertex_gene]]
#   
#   #edgetest[[1]][index_edge_gene[valid_indices_index_edge_gene]]<-namechange[valid_indices_index_edge_gene]
#  #edgetest[[2]][index_edge_gene2[valid_indices_index_edge_gene2]]<-namechange[valid_indices_index_edge_gene2]
#   #vertextest[[1]][index_vertex_gene[valid_indices_index_vertex_gene]]<-namechange[valid_indices_index_vertex_gene]
#   
#   #edgetest[[1]][index_edge_gene[valid_indices_index_edge_gene]]
#   #edgetest[[2]][index_edge_gene2[valid_indices_index_edge_gene2]]
#   #vertextest[[1]][index_vertex_gene[valid_indices_index_vertex_gene]]
#   #########################################################################################
#   
#   #print("################################")
#   #print(edgetest)
#   #print("################################")
#   #print(vertextest)
#   #print("################################")
#  # interactions_E<-edgetest
#   #interactions_V<-vertextest
#   
#   #unique_vertices_in_edges <- unique(c(interactions_E$source, interactions_E$target))
#  # missing_vertices <- setdiff(unique_vertices_in_edges, interactions_V$name)
#   #head(rownames(paste0("Raw data/network_merge/edge/", data_name, ".csv")))
#   
#  # print(missing_vertices)
#   
# 
#   g <- graph_from_data_frame(interactions_E, directed = TRUE, vertices = interactions_V)
#   
#   g <- decompose(g, min.vertices = 15)[[1]]
#   g <- set_vertex_attr(g, "degree", index = V(g)$name, degree(g))
#   
#   
#   print(class(degree(g)))
#   print(length(degree(g)))
#   print(class(V(g)$log2FoldChange))
#   print(length(V(g)$log2FoldChange))
#   
#   
#   
#   
#   Deg_max <- sort(degree(g), decreasing = TRUE)[1:2]
#   log2fc_up <- sort(V(g)$log2FoldChange, decreasing = TRUE)[1:2]
#   log2fc_down <- sort(V(g)$log2FoldChange)[1:2]
#   
#   graph <- ggraph(g, layout = "grid") +
#     geom_edge_link(aes(colour = weight), alpha = 0.5) +
#     
#     geom_node_point(aes(size = degree, colour = log2FoldChange, stroke = 1.5)) + 
#     
#     scale_edge_colour_gradientn(name = "weight", colors = rev(brewer.pal(5, "RdYlBu")), 
#                                 limits = c(-1, 1), space = "Lab") +
#     
#     scale_color_gradientn(name = "log2FC", colours = rev(brewer.pal(5, "PiYG")), 
#                           limits = c(-5, 5), space = "Lab")  +
#     
#     scale_size_continuous(name = "degree", range = c(0.01, 2))  +
#     
#     coord_fixed() +
#     
#     ggtitle(label = paste0(data_name,  " nodes:", 
#                            length(V(g)), ' edges:', length(E(g)))) +
#     
#     
#     geom_node_text(aes(filter = (log2FoldChange %in% log2fc_up), label = name),
#                    colour = "red", repel = T) +
#     
#     geom_node_text(aes(filter = (log2FoldChange %in% log2fc_down), label = name),
#                    colour = "blue", repel = T) +
#     
#     geom_node_text(aes(filter = (degree %in% Deg_max), label = name),
#                    colour = "black", repel = T) +
#     
#     theme_void() +
#     
#     theme(plot.title = element_text(size = 20, hjust = 0.5), 
#           legend.title = element_text(size = 18),
#           legend.text = element_text(size = 18))
#   # print(graph)
#   net_res <<- append(net_res, list(graph))
#   return (g)
# }
# 
# 
# #net_res <- list()
# #for(index in 1:length(dataset_TCGA))
# #{
# #network_graph(dataset_TCGA[index])
# #}
# 
# graph_comp<-network_graph(cancer_type)


## vertex_data <- as.data.frame(vertex_attr(graph_comp))
## print(vertex_data)
## edge_data <- as.data.frame(edge_attr(graph_comp))
## print(edge_data)

##class(graph_comp)
##graph_comp.vs

###########################################################




#is_directed(graph_comp)

#print(graph_comp)
#View(graph_comp)
#dir.create("Result", showWarnings = FALSE)
#pdf(paste0("Result/",cancer_type,"_Specific_network_visualization_(Maximum_subnet)02.pdf"), width = 40.0, height = 50)

#print(ggarrange(net_res[[1]], nrow = 6, ncol = 3))
#dev.off()

################################################################################################
### Network community detection
#class(graph_comp[1])
#print(graph_comp[1])
#View(as.data.frame(graph_comp[1]))
#print(row.names(as.data.frame(graph_comp[1])))
#print(rownames(as.data.frame(graph_comp[1])))
#row_names<-row.names(as.data.frame(graph_comp[1]))
#saveRDS(row_names,"./data/fan_network.RData")
#write.csv(row_names,"./data/fan_network.csv")



setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
spg <- list() # List of spinglass simulations
spg_mod <- numeric() # List of modularity simulations
count<-1
# #for (n in c(50, 100, 1000)){
#   #for (k in 1:n){
#   #  print(Sys.time())
#   #  print(count)
#   #  spg[[count]] <- cluster_spinglass(graph_comp,
#                                       spins = 25,
#                                       weights = NULL)
#    # spg_mod[count] <- spg[[count]]$modularity
#     count = count + 1
#     print(Sys.time())
#   }
# }

#count<-1
for (k in 1:50){
  print(Sys.time())
  print(k)
  spg[[k]] <- cluster_spinglass(graph_comp,
                                    spins = 25,
                                    weights = E(graph_comp)$weight,
                                    implementation = "neg")

  spg_mod[k] <- spg[[k]]$modularity
  #count = count + 1
  print(Sys.time())
}





####

########################################################################
# 
# max_index <- which.max(spg_mod)  # 返回最大模块性的索引
# 
# best_spg <- spg[[max_index]]  # 获取对应的社群对象
# print(max_index)
# num_communities <- length(unique(membership(best_spg)))  # 计算社群数量
# colors <- rainbow(num_communities)  # 为每个社群生成一个颜色
# 
# community_ids <- membership(best_spg)  # 获取社群成员身份
# node_colors <- colors[community_ids]  # 根据社群ID分配颜色
# #plot(graph_comp, vertex.color = node_colors, main="Network Visualization with Community Coloring")
# plot(graph_comp, vertex.color = node_colors, vertex.size = 5, vertex.label = NA,
#      main = "Network Visualization with Community Coloring")






########################################################################


# community_ids <- membership(spg[[max_index]])
# 
# 
# print(community_ids)
# 
# 
# 
# num_communities <- max(community_ids)
# 
# print(num_communities)
# colors <- RColorBrewer::brewer.pal(min(num_communities, 9), "Set1")  # 使用颜色方案
# 
# # 将颜色分配给节点
# 
# V(graph_comp)$color <- colors[community_ids]
# 
# # 转换为 visNetwork 需要的格式
# nodes <- data.frame(id = V(graph_comp)$name, label = V(graph_comp)$name, color = V(graph_comp)$color)
# edges <- get.data.frame(graph_comp, what="edges")
# 
# # 如果你的图是有向的，添加 arrows
# edges$arrows <- "to"
# 
# visNetwork(nodes, edges) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
#   visEdges(smooth = FALSE) %>%
#   visNodes(shape = "dot", scaling = list(label = list(enabled = TRUE, min = 8, max = 20))) %>%
#   visLayout(randomSeed = 123)  # 确保布局的一致性

###########################################################################################################
# # 假设已经计算得到了社区 membership
# library(visNetwork)
# 
# nodes <- data.frame(id = names(community_ids), label = names(community_ids), group = community_ids)
# edges <- get.data.frame(graph_comp, what="edges")
# 
# # 生成颜色
# colors <- rainbow(length(unique(community_ids)))
# 
# # 创建可视化
# visNetwork(nodes, edges) %>%
#   visNodes(color = {background = colors[community_ids], border = '#2b2b2b'},
#            shadow = TRUE) %>%
#   visEdges(smooth = FALSE) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
#   visLayout(randomSeed = 42)


###########################################################################
#length(spg[[max_index]])

saveRDS(spg, paste0("./data/spinglass",cancer_type,"_spg.RData"))
saveRDS(spg_mod, paste0("./data/spinglass",cancer_type,"_spg_mod.RData"))

###########################################################################
# Select the fully connected partitions
# full_conn_runs <- numeric() 
# for(ind in order(spg_mod, decreasing = T)){
#   melanet_spg <- spg[[ind]]
#   full_connected <- T
#   for(i in 1:length(melanet_spg)){
#     graph_comp_sub <- induced_subgraph(graph = graph_comp, v = which(melanet_spg$membership == i)) 
#     if(is.connected(graph_comp_sub) == F){
#       full_connected <- F
#     }
#   }
#   if(full_connected){
#     full_conn_runs <- append(full_conn_runs, ind)
#   }
# }
# 
# 
# 
# for(i in 1:length(spg[[max_index]])){
#   graph_comp_sub <- induced_subgraph(graph = graph_comp, v = which(melanet_spg$membership == i)) 
#   if(is.connected(graph_comp_sub) == F){
#     full_connected <- F
#   }
#   else{
#     full_conn_runs <- append(full_conn_runs,)
#     
#   }
# }
#print(max_index)
#full_conn_runs <- append(full_conn_runs,)
# The fully connected partition with the maximum modularity 
#spg_ind <- full_conn_runs[1]


print(spg)
melanet_spg <- spg[[max_index]]
#head(melanet_spg)
saveRDS(melanet_spg, paste0("./data/spinglass/",cancer_type,"_melanet_spg_test.RData"))
melanet_cmt <- as.list(communities(melanet_spg)) # Melanoma network communities
#print(melanet_cmt[[1]])
saveRDS(melanet_cmt, paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))

melanet_cmt_load<-readRDS(paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))

#print(melanet_cmt_load)
length(melanet_cmt_load)
