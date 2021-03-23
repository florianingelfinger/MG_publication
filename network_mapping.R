##### Network mapping approach to combine blood/thymus and mass/spectral flow cytometry data used in Ingelfinger et al. 
# "Single-cell profiling of myasthenia gravis identifies a pathogenic T cell signature" #####

# load library
library("ggplot2")
library("ggthemes")
library("flowCore")
library("ncdfFlow")
library("panorama")
library("vite")
library("grappolo")

## prepare landmark data
# load untransformed files
setwd("/Myasthenia gravis/Pooled run/")
fs_surf1  <- read.ncdfFlowSet(files = "MG_surf1_conc_untransf.fcs",
                              transformation = F,
                              phenoData = ,
                              truncate_max_range = F)
fs_surf2  <- read.ncdfFlowSet(files = "MG_surf2_conc_untransf.fcs",
                              transformation = F,
                              phenoData = ,
                              truncate_max_range = F)
data_surf1 <- fsApply(fs_surf1, exprs)
data_surf2 <- fsApply(fs_surf2, exprs)
data_blood_untrans <- rbind(data_surf1, data_surf2)


# subset on MG patients
setwd("/Myasthenia gravis/Pooled run/surf analysis/01_FlowSOM clustering k6")
md <- read.csv("surf_meta_2019-08-13.csv")
gate_MG <- subset(md, dis == "MG")$gate_source
data_blood_untrans_MG <- data_blood_untrans[data_blood_untrans[,"gate_source"]%in% gate_MG,]


# subsample cells per MG patient
n_samples <- unique(data_blood_untrans_MG[,"gate_source"])

data_sub <- data.frame(NULL)
n_sub_i <- 1500
for(i in n_samples){
  data_i <- data_blood_untrans_MG[data_blood_untrans_MG[,"gate_source"]==i,]
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

# load clustered master file
setwd("/Myasthenia gravis/Pooled run/surf analysis/01_FlowSOM clustering k6")
fs  <- read.ncdfFlowSet(file= "Surface_merged_clustered_k9_2019-01-07.fcs",
                        transformation = F,
                        phenoData = ,
                        truncate_max_range = F)
data_blood_master <- fsApply(fs, exprs)
data_blood_master <- cbind(data_blood_master, data_blood_master[,"th_labels"])
colnames(data_blood_master)[ncol(data_blood_master)] <- "pool_level"

# modify clustering labels so blood and thymus correspond
data_blood_master[data_blood_master[,"pool_level"]==3,"pool_level"] <- 0
data_blood_master[data_blood_master[,"pool_level"]==4,"pool_level"] <- 1
data_blood_master[data_blood_master[,"pool_level"]==9,"pool_level"] <- 8

# transfer cluster labels to untransformed subset
data_blood_untrans_MG_clust <- merge(data_blood_untrans_MG, data_blood_master[,c("cell_id", "pool_level")], by= "cell_id")

# export landmark nodes
setwd("/Myasthenia gravis/Thymus biopsies/Thymus_aurora_run/2020-07-03_scaffold/gated/")
n_landmark <- unique(data_blood_untrans_MG_clust[,"pool_level"])
clusternames <- c("CD4 T cells", "CD8 T cells", "gd T cells", "Tregs","NK cells", "NKT cells", "B cells",  "Monocytes", "DCs")
data_blood_untrans_MG_clust_m <- data.matrix(data_blood_untrans_MG_clust)

for (i in n_landmark) {
  sample_i <- data_blood_untrans_MG_clust_m[data_blood_untrans_MG_clust_m[,"pool_level"]==i,]
  ff_i <- flowFrame(sample_i)
  print(paste("sample_surface_", clusternames[i], ".fcs", sep=""))
  write.FCS(ff_i, filename =  paste("sample_surface_", clusternames[i], ".fcs", sep=""))
}

## generate blood fcs file for MG
# subsample
data_blood_master_MG <- data_blood_master[data_blood_master[,"gate_source"]%in% gate_MG,]
n_samples <- unique(data_blood_master_MG[,"gate_source"])

data_sub <- data.frame(NULL)
n_sub_i <- 1000
for(i in n_samples){
  data_i <- data_blood_master_MG[data_blood_master_MG[,"gate_source"]==i,]
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

# export fcs file
ff_i <- flowFrame(data.matrix(data_sub))
write.FCS(ff_i, filename =  "MG_blood.fcs")


## generate thymus fcs file for MG and CTRL
setwd("/Myasthenia gravis/Thymus biopsies/Thymus_aurora_run/2020-06-08_surface/1_umap_Flowsom_main/2020_06_25/")
load(file="data_main_clust_k7.RDa")
data_thymus <- data_clustered

setwd("/Myasthenia gravis/Thymus biopsies/Thymus_aurora_run/")
md <- read.xls("md_thymus_aurora.xlsx")
gate_thy_MG <- subset(md, group =="MG")$sample_id
data_thymus_MG <- data_thymus[data_thymus[,"sample_id"]%in%gate_thy_MG,]

n_samples <- unique(data_thymus_MG[,"sample_id"])

data_sub <- data.frame(NULL)
n_sub_i <- 15000
for(i in n_samples){
  data_i <- data_thymus_MG[data_thymus_MG[,"sample_id"]==i,]
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

# export fcs file
setwd("/Myasthenia gravis/Thymus biopsies/Thymus_aurora_run/2020-07-03_scaffold")
ff_i <- flowFrame(data.matrix(data_sub))
write.FCS(ff_i, filename =  "MG_thymus.fcs")


# CTRL
gate_thy_CTRL <- subset(md, group =="CTRL")$sample_id
data_thymus_CTRL <- data_thymus[data_thymus[,"sample_id"]%in%gate_thy_CTRL,]
n_samples <- unique(data_thymus_CTRL[,"sample_id"])

data_sub <- data.frame(NULL)
n_sub_i <- 10000
for(i in n_samples){
  data_i <- data_thymus_CTRL[data_thymus_CTRL[,"sample_id"]==i,]
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

setwd("/Myasthenia gravis/Thymus biopsies/Thymus_aurora_run/2020-07-03_scaffold")
ff_i <- flowFrame(data.matrix(data_sub))
write.FCS(ff_i, filename =  "CTRL_thymus.fcs")


## Cluster files
col.names <- colnames(data_thymus)[c(1, 2, 3, 5, 7, 8, 9, 12, 14, 15, 17, 19, 21)]
files.list <- list.files(path = getwd(), pattern = ".fcs$")
cluster_fcs_files(files.list, num.cores = 4, col.names = col.names, num.clusters = 200,
                  asinh.cofactor = NULL)

input.files <- paste(files.list, "clustered.txt", sep=".")
landmarks.data <- load_landmarks_from_dir("gated/", asinh.cofactor = 5, transform.data = T)

# run scaffold analysis
run_scaffold_analysis(input.files, ref.file = input.files[2], 
                      landmarks.data = landmarks.data, col.names = col.names, process.clusters.data = F)