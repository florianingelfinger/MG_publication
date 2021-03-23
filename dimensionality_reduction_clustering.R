##### Dimensionality reduction and clustering pipeline used in Ingelfinger et al. "Single-cell profiling of myasthenia gravis identifies a pathogenic T cell signature" #####

# load library
library("ggplot2")
library("umap")
library("ggthemes")
library("reshape2")
library("gplots")
library("ComplexHeatmap")
library("ConsensusClusterPlus")


# read data
setwd("/Thymus_aurora_run/2020-06-08_surface/transformation/")
load(file="2020-07-24_thymus_surf_trans.RDa")
data <- data_trans

# subsample cohort for dimensionality reduction
n_samples <- unique(data[,"sample_id"])

data_sub <- data.frame(NULL)
n_sub_i <- 1500
for(i in n_samples){
  data_i <- data[data[,"sample_id"]==i,]
  n_i <- nrow(data_i)
  set.seed(123)
  if(n_sub_i > n_i){
    ix_i <- sample(1:n_i, n_i)
  }
  else{
    ix_i <- sample(1:n_i, n_sub_i)
  }
  data_sub_i <- data_i[ix_i,]
  data_sub <- rbind(data_sub, data_sub_i)
}

# prepare data for UMAP
umap_cols <- colnames(data_sub)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22)]
data_rumap <- data_sub[ , umap_cols]

# calculate UMAP coordinates
data_umap <-  umap(data_rumap, random_state=123, verbose =T)
umap <- as.data.frame(data_umap$layout)
colnames(umap) <- c("UMAP1", "UMAP2")

# plot UMAP black
plot_umap <- ggplot(umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1)
plot_umap

# plot UMAP with expression overlayed
data_plot_u <- cbind(data_rumap, umap)
data_melt_u <- melt(data_plot_u, id.vars = c("UMAP1", "UMAP2"))

color_grad_flow2 <- c("#331820", "#4d4f55", "#55626b", "#5a767e", "#628b8a", "#709b91", "#82aa96", "#98b89a", "#b0c6a2", "#c9d3ab", "#e4e0b6", "#feedc3")

plot_umap_expr <- ggplot(data_melt_u, aes(x = UMAP1, y = UMAP2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = color_grad_flow2, limits = c(0,1)) +
  facet_wrap(~ variable, ncol = 6) +
  theme_few() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
plot_umap_expr

# prepare umap coordinates for export
umap_exp <- cbind(data_sub[,c("gate_source", "cell_id", "sample_id")], data_plot_u[,c("UMAP1", "UMAP2")])
save(umap_exp, file="data_umap_coordinates.RDa")

## FlowSOM clustering
data_flow <- data
data_flow <- data.matrix(data_flow)
clusteringcol_flow <- colnames(data_flow)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22)]

# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(flowFrame(exprs = data_flow, desc = list(FIL = 1)), transform = F, scale = F, compensate = F)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusteringcol_flow)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]

# Metaclustering step
max  <- 20
gate_source <- as.factor(data[,"gate_source"])
meta_results <- data.frame(gate_source)

# do a manual metaclustering for all values up to max
for (i in 3:max) {
  set.seed(123)
  out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = i)
  meta_results <- cbind(meta_results, as.factor(out_meta[labels]))}
meta_results <- meta_results[,2:ncol(meta_results)]
colnames(meta_results) <- paste("k.", 3:max, sep = "")
meta_results <- cbind(meta_results, data[,"cell_id"])
colnames(meta_results)[19] <- "cell_id"

# prepare the joined metaclustering data
data_meta <- umap_exp[,c("cell_id", "UMAP1", "UMAP2")]
meta.melt <- melt(meta_results, id.vars = "cell_id", variable.name = "k.value", value.name = "cluster.assigment")
meta.melt$cluster.assigment <- as.factor(meta.melt$cluster.assigment)
joined.meta <- merge(meta.melt, data_meta, by = "cell_id")


# plot UMAPs with metaclusters overlayed
color.scale <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

plot_umap_meta <- ggplot(joined.meta, aes(x = UMAP1, y = UMAP2, color = cluster.assigment)) +
  geom_point(size = 0.25) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(name = NULL,values = c(color.scale)) +
  facet_wrap(~ k.value, ncol =6) +
  guides(colour = guide_legend(override.aes = list(size=5), title="Cluster"))

plot_umap_meta

# consensus clustering to determine number of clusters
heat_mat <- matrix(NA, nrow = 100, ncol = length(clusteringcol_flow))
for(i in 1:100) {
  temp_mat <- data_flow[labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
rownames(heat_mat) <- 1:100
colnames(heat_mat) <- clusteringcol_flow

set.seed(123)
results <- ConsensusClusterPlus(t(heat_mat), maxK = 20, reps = 1000, pItem = 0.8,
                                pFeature = 1,  clusterAlg = "hc", verbose = F,
                                distance = "euclidean", seed = 123, plot = "png",
                                writeTable = T)


# decide on number of metaclusters and perform metaclustering
choosen_k <- 14
set.seed(123)
out_meta <- metaClustering_consensus(out_fSOM$map$codes, k = choosen_k)
pop_labels <- out_meta[labels]

# Calculate median marker expression for heatmap
heat_mat <- matrix(NA, nrow = choosen_k, ncol = length(clusteringcol_flow))
for(i in 1:choosen_k) {
  temp_mat <- data_flow[pop_labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(heat_mat) <- paste("cluster", 1:choosen_k, sep = "")
colnames(heat_mat) <- clusteringcol_flow

# make a cluster heatmap
heatmap.2(heat_mat,
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks,
          mar=c(6,9))

# manually merge metaclusters
CD4 <- c(8, 7, 14)
CD8 <- c(12, 6, 11)
Thymocytes <- c(1, 2, 9)
NK <- c(5)
NKT <- c(3)
B <- c(13, 10)
myeloid <- c(4)

cluster_labels <- rep(0, length(pop_labels))
cluster_labels[pop_labels %in% CD4] <- 1
cluster_labels[pop_labels %in% CD8] <- 2
cluster_labels[pop_labels %in% Thymocytes] <- 3
cluster_labels[pop_labels %in% NK] <- 4
cluster_labels[pop_labels %in% NKT] <- 5
cluster_labels[pop_labels %in% B] <- 6
cluster_labels[pop_labels %in% myeloid] <- 7

# Heatmap of manually merged clusters
heat_mat <- matrix(NA, nrow = 7, ncol = length(clusteringcol_flow))
for(i in 1:7) {
  temp_mat <- data[cluster_labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
clusternames <- c("CD4 T cells", "CD8 T cells", "Thymocytes", "NK cells", "NKT cells", "B cells", "myeloid cells")
rownames(heat_mat) <- clusternames
colnames(heat_mat) <- clusteringcol_flow
white.black <- colorRampPalette(c("white", "black"))(n = 9)

# make a cluster heatmap
pheatmap(mat = heat_mat, 
         scale= "none",
         color = white.black,
         breaks = breaks,
         cluster_rows=F,
         cluster_cols=F,
         cellwidth = 10, cellheight = 10,
         border_color = NA,
         filename= "main_clusters_heatmap.pdf"
)

# Project clustering on UMAP
data_clust <- cbind(data, cluster_labels)
data_umap_compl <- merge(umap_exp, data_clust[,c("cluster_labels", "cell_id")], by = "cell_id")
data_umap_compl$cluster_labels <- as.factor(data_umap_compl$cluster_labels)

color_qual_flow <- c("#9F5668", "#DAD95E", "#F69165", "#609E74", "#A0C984", "#262626", "#3CA0B2", "#44C6AA")

plot_umap_clust <- ggplot(data_umap_compl, aes(x = UMAP1, y = UMAP2, color = cluster_labels)) +
  geom_point(size = .5) + coord_fixed(ratio = 1) +
  scale_color_manual(values = color_qual_flow, name= "cluster_labels", breaks=1:7,
                     labels= clusternames) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title=element_blank(),
        legend.background = element_rect())
plot_umap_clust

# save clustered data
colnames(data_clust)[ncol(data_clust)] <- "main_clusters"
save(data_clust, file="data_main_clust_k7.RDa")
