# sudo apt list -a r-base-core;
# sudo apt-get install r-base-core=4.2.0-1.2204.0;
# Install RStudio via their deb package;
# Installation of bio3D in RStudio: install.packages("bio3d", dependencies=TRUE);
# Installation of packages in terminal: sudo apt-get install gfortran libblas-dev liblapack-dev libnetcdf-dev;
# Installation of shape in RStudio: install.packages("shape");
# Installation of adegenet in RStudio: install.packages("adegenet");
# The majority of this code originates from http://thegrantlab.org/bio3d/articles/online/traj_vignette/Bio3D_md.html; 

#################
#PCA using bio3d#
#################

#Load libraries:
library(bio3d)
library(adegenet)
library(shape)

#Check if libraries are loaded:
sessionInfo()

#Load .pdb file into environment:
mypdbfile <- "/home/wout/AR/7_ARE/1_molecular_dynamics/MMTV/2_third_run/md_nc_la.pdb"

#Read the .pdb file:
pdb <- read.pdb(mypdbfile)

#Check the file:
print(pdb)

#Load .dcd file into environment:
mydcdfile <- "/home/wout/AR/7_ARE/1_molecular_dynamics/MMTV/2_third_run/md_nc_la.dcd"

#Read the dcd file:
dcd <- read.dcd(mydcdfile)

#Check the file: 
print(dcd)

#Select all C-alpha atoms and all atoms for trajectory frame superposition:
ca.inds <- atom.select(pdb, elety="CA")
#ca.inds2 <- atom.select(pdb, "protein")

#Perform superposition on both CA and protein:
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)
#xyz2 <- fit.xyz(fixed = pdb$xyz, mobile = dcd, 
#                  fixed.inds = ca.inds2$xyz, 
#                  mobile.inds = ca.inds2$xyz)
               
#Check if the coordinates match:
dim(xyz) == dim(dcd)
#dim(xyz2) == dim(dcd)

#Check RMSD of CA:
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
plot(rd, typ="l", ylab="RMSD (Å)", xlab="Frame No.", main="AR DBD monomer RMSD")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)

#Check RMSF of CA:
rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF (Å)", xlab="Residue Position", typ="l", main="AR DBD monomer RMSF")

#Principal Component Analysis of CA (PCA = dimension reduction method):
set.seed(1)
pc <- pca.xyz(xyz[,ca.inds$xyz], use.svd = TRUE)
plot(pc, col=bwr.colors(nrow(xyz)))
colorlegend(col= bwr.colors(nrow(xyz)), zlim = range(0:130), main="Time Step(ns)", posx = c(0.49, 0.51), posy = c(0.25, 0.75))

#Principal Component Analysis of protein (PCA = dimension reduction method):
#set.seed(1)
#pc2 <- pca.xyz(xyz2[,ca.inds2$xyz], use.svd = TRUE)
#plot(pc2, col=bwr.colors(nrow(xyz2)))
#colorlegend(col= bwr.colors(nrow(xyz2)), zlim = range(0:150), main="Time Step(ns)", posx = c(0.49, 0.51), posy = c(0.25, 0.75))

#####################################
#Cluster optimization using adegenet#
#####################################

#Cluster identification of CA using Bayesian Information Criterion (BIC = method to optimize choosing between 2 or more alternative models):
clust <- find.clusters(pc$z, n.pca = 100)
print(clust$Kstat[1:50])
plot(clust$Kstat[1:50], type="b", col="blue", xlab = "", ylab = "")
title(xlab="Number of Clusters", ylab= "BIC", main = "AR DBD monomer")

#Cluster identification of protein using Bayesian Information Criterion (BIC = method to optimize choosing between 2 or more alternative models):
#clust2 <- find.clusters(pc2$z, n.pca = 100)
#print(clust2$Kstat[1:20])
#graphics.off()
#plot(clust2$Kstat[1:30], type="b", col="blue", xlab = "", ylab = "")
#title(xlab="Number of Clusters", ylab= "BIC", main = "AR DBD monomer")

#Kmeans clustering of CA:
#kcl <- kmeans(pc$z, 5)
#plot(pc, col = kcl$cluster)

#Kmeans clustering of protein:
#kcl2 <- kmeans(pc2$z, 5)
#plot(pc2, col = kcl$cluster)

#Complete clustering of CA:
#hc <- hclust(dist(pc$z[,1:9]))
#grps <- cutree(hc, k=4)
#table(grps)
#plot(pc, col=grps)

#Complete clustering of protein:
#hc2 <- hclust(dist(pc2$z[,1:7]))
#grps2 <- cutree(hc2, k=7)
#table(grps2)
#plot(pc2, col=grps2) 

#UPGMA clustering of CA:
hc3 <- hclust(dist(pc$z[,1:8]), method = "average")
grps3 <- cutree(hc3, k=25)
table(grps3)
plot(pc, col = grps3)

#UPGMA clustering of protein:
#hc4 <- hclust(dist(pc2$z[,1:7]), method = "average")
#grps4 <- cutree(hc4, k = 7)
#table(grps4)
#plot(pc2, col = grps4)

#Selection of most central conformation of clusters:

cluster_center = aggregate(pc$z,list(cluster=grps3),mean)
pc.z <- cbind(1:nrow(pc$z), pc$z)

index_closest_to_center <- function(n){
  cluster <- data.frame(pc.z[grps3==n,])
  names(cluster_center) <- names(cluster)
  dist = c(0, nrow(cluster))
  angle = c(0, nrow(cluster))
  for (i in 1:nrow(cluster)) {
    dist[i] = dist(rbind(cluster[i,-1], cluster_center[n,-1]))
    angle[i] = sum(as.vector(as.numeric(cluster[i,-1]))*as.vector(as.numeric(cluster_center[n,-1]))) / (sqrt(sum(as.vector(as.numeric(cluster[i,-1])) * as.vector(as.numeric(cluster[i,-1])))) * sqrt(sum(as.vector(as.numeric(cluster_center[n,-1])) * as.vector(as.numeric(cluster_center[n,-1])))))
  }
  cluster$dist <- dist
  cluster$angle <- angle
  closest_dist <- cluster[which.min(cluster$dist), ]$X1
  closest_angle <- cluster[which.max(cluster$angle), ]$X1
  angle_1 <- cluster[which.max(cluster$angle), ]$angle
  return(c(closest_dist, closest_angle, angle_1))
}

cluster1 <- data.frame(pc.z[grps3==1,])
cluster2 <- data.frame(pc.z[grps3==2,])
cluster3 <- data.frame(pc.z[grps3==3,])
cluster4 <- data.frame(pc.z[grps3==4,])

#Change the amount of clusters:
index_closest_to_center(3)
index_closest_to_center(21)
index_closest_to_center(7)

#Chose the first frame if the cosine is is smaller than 0.9:
"
Smallest distance frame		Smallest angle 		frame cosine 
44.0000000 			45.0000000  		0.9344525
"

#Select the frames in vmd and count +1 in the frame:
vmd .dcd chosenstartframe.pdb

#Energy minimize structures after changing CYM back to CYS

#Upload to FTMap

#Plot residue-residue cross correlation:
cij<-dccm(xyz[,ca.inds$xyz])
plot(cij, main="AR DBD monomer residue-residue cross correlation")

#3D PyMol visualisation of residue-residue cross correlation:
pymol(cij, pdb, type="launch",exefile="/opt/pymol/bin/pymol")
