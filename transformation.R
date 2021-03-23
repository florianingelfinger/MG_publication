##### Transformation pipeline used in Ingelfinger et al. "Single-cell profiling of myasthenia gravis identifies a pathogenic T cell signature" #####

# load library
library("ncdfFlow")
library("flowCore")
library("gridExtra")
library("grid")
library("ggplot2")


# Load and concatenate different samples 
setwd("~/Thymus_aurora_run/2020-06-08_surface/compensated/2020-07-24/")
file.names <- list.files(path = getwd(), pattern = ".fcs$")
flowset  <- read.ncdfFlowSet(files = file.names,
                              transformation = F,
                              phenoData = ,
                              truncate_max_range = F)


# combine data in a matrix
data_surf <- fsApply(flowset, exprs)
head(data_surf)
dim(data_surf)

# Retrieve markernames and rename
colnames(data_surf)[7:29] <- markernames(fs_surf1)

# generate ids for sample and cell tracing
dim.vec_surf <- fsApply(fs_surf1, dim)
dim.vec_surf <- as.numeric(dim.vec_surf[,1])
gate.source_surf <- as.vector(x = NULL)
for(i in 1:length(dim.vec_surf)) {temp.source <- rep(i, dim.vec_surf[i])
gate.source_surf <- c(gate.source_surf, temp.source)}
cell_id_surf <- 1:dim(data_surf)[1]

# add to master data
data_surf <- cbind(data_surf, gate.source_surf)
data_surf <- cbind(data_surf, cell_id_surf)
colnames(data_surf1)[31:32] <- c("gate_source", "cell_id")


# extract sample id from string in fcs files and replace
sample_id_surf <- as.numeric(gsub("[^0-9]", "", file.names))
sample_id_surf <- as.numeric(gsub("2020060819", "", sample_id_surf))
sample_gate_surf <- cbind(unique(data_surf[,"gate_source"]), sample_id_surf)
colnames(sample_gate_surf1) <- c("gate_source", "sample_id")

data_surf_sid <- merge(data_surf, sample_gate_surf, by= "gate_source")

data_surf <-  data_surf_sid

# save untransformed files
save(data_surf, file="2020-07-24_Thymus_surf_untransf.RDa")

# convert to matrix
data_matrix_surf1 <- data.matrix(data_surf1[,c(markernames(fs_surf), "cell_id", "gate_source", "sample_id")])

### transform data
# note: cofactors have been chosen on na indiividual basis depending on the separation of positive and negative fraction
# as observed in the output pdf file.

asinh_scale <- 150

data.trans_surf <- asinh(data_matrix_surf/ asinh_scale)
data.trans_surf[,"CD3"] <- asinh(150*sinh(data.trans_surf[,"CD3"])/3500)
data.trans_surf[,"CD56"] <- asinh(150*sinh(data.trans_surf[,"CD56"])/1500)
data.trans_surf[,"CD4"] <- asinh(150*sinh(data.trans_surf[,"CD4"])/1500)
data.trans_surf[,"IgG"] <- asinh(150*sinh(data.trans_surf[,"IgG"])/3500)
data.trans_surf[,"CD8"] <- asinh(150*sinh(data.trans_surf[,"CD8"])/2000)
data.trans_surf[,"CD103"] <- asinh(150*sinh(data.trans_surf[,"CD103"])/6000)
data.trans_surf[,"CD19"] <- asinh(150*sinh(data.trans_surf[,"CD19"])/8000)
data.trans_surf[,"CD11c"] <- asinh(150*sinh(data.trans_surf[,"CD11c"])/3000)
data.trans_surf[,"CD14"] <- asinh(150*sinh(data.trans_surf[,"CD14"])/5000)
data.trans_surf[,"CD45R0"] <- asinh(150*sinh(data.trans_surf[,"CD45R0"])/5500)
data.trans_surf[,"CD69"] <- asinh(150*sinh(data.trans_surf[,"CD69"])/12000)
data.trans_surf[,"CD27"] <- asinh(150*sinh(data.trans_surf[,"CD27"])/4500)
data.trans_surf[,"HLA-DR"] <- asinh(150*sinh(data.trans_surf[,"HLA-DR"])/2000)
data.trans_surf[,"CD127"] <- asinh(150*sinh(data.trans_surf[,"CD127"])/2000)
data.trans_surf[,"IgD"] <- asinh(150*sinh(data.trans_surf[,"IgD"])/2500)
data.trans_surf[,"CD38"] <- asinh(150*sinh(data.trans_surf[,"CD38"])/5000)
data.trans_surf[,"CD25"] <- asinh(150*sinh(data.trans_surf[,"CD25"])/12000)
data.trans_surf[,"CD5"] <- asinh(150*sinh(data.trans_surf[,"CD5"])/2000)
data.trans_surf[,"CXCR5"] <- asinh(150*sinh(data.trans_surf[,"CXCR5"])/1000)
data.trans_surf[,"IgM"] <- asinh(150*sinh(data.trans_surf[,"IgM"])/2000)
data.trans_surf[,"ICOS"] <- asinh(150*sinh(data.trans_surf[,"ICOS"])/2000)

data.trans_surf[,"cell_id"] <- data_surf[,"cell_id"]


## shift to 0
q.vector <- apply(data.trans_surf, 2, function(x) quantile(x, 0.001, names = F))
data.shift <- data.trans_surf
data.shift <- sweep(data.shift, 2, q.vector)

data.shift[,c("CD8", "gate_source")] <- sweep(data.trans_surf1[,c("CD8", "gate_source")], 2, quantile(data.trans_surf1[,c("CD8", "gate_source")], 0.01, names = F))


# normalize based on the 99.99th percentile to have everything between 0 and 1
data.trans.new_surf <- data.shift
per.vector_surf <- apply(data.trans.new_surf, 2, function(x) quantile(x, 0.9999, names = F))
per.vector_surf
data.trans.new_surf <- t(t(data.trans.new_surf) / as.numeric(per.vector_surf))

# recreate original ids
data.trans.new_surf[,"cell_id"] <- data_surf[,"cell_id"]
data.trans.new_surf[,"gate_source"] <- data_surf[,"gate_source"]
data.trans.new_surf[,"sample_id"] <- data_surf[,"sample_id"]




## plot all marker vs each other in order to evaluate transformation
# subsample cells for plotting
n_cells <- 50000
n <- nrow(data.trans.new_surf)
set.seed(123)
ix <- sample(1:n, n_cells)
sub_surf <- rbind(data.trans.new_surf[ix,])
s_surf <- data.frame(sub_surf)

color.scale <- c("#1B9E77","gold" , "#7570B3", "#E7298A", "#D95F02" ,  "lightpink1", "#666666",
                 "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
                 "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" ,
                 "#BC80BD", "#80B1D3")

plot_surf_list <- list()
for(i in 1:length(colnames(data.trans_surf))){
  plot_surf <- ggplot(s_surf, aes_string(x = colnames(s_surf)[i], y = "CD3")) +
    geom_hex(bins = 100) +
    xlim(-0.1, 1.5) +
    ylim(-0.1, 1.5) +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colours = color.scale, trans = "sqrt") +
    theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 3), axis.title = element_text(size = 3))
  plot_surf_list[[i]] <-  plot_surf1
}

m1 <- arrangeGrob(grobs = plot_surf_list, ncol = 1, top=textGrob("Experiment 1", gp=gpar(fontsize=6)))
ggsave(file="marker_comp_total.pdf", m1, width = 1, height = 25)

# save transformed dataframe
data_trans <- data.trans.new_surf
save(data_trans, file = "2020-07-24_thymus_surf_trans.RDa")


