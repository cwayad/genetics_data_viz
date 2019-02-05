# Redirection vers le bon path
setwd("F:/M2-MLDS/Visualisation/projet")


#Chargement des données
pomeroy <- read.table("pomeroy-2002-v2_database.txt",header = T, sep = "\t", fill = T)
#View(pomeroy)
#Transposition de la table pomeroy
pomeroy_trans=t(pomeroy)
#View(pomeroy_trans)
#Renommer les les lignes
colnames(pomeroy_trans)= pomeroy_trans[1,]
#Enlever la premiere ligne
pomeroy_trans = pomeroy_trans[-1,]
#View(pomeroy_trans)
#Résumé des deux tables transposées
#summary(pomeroy_trans)
#Affecter 2 à MPM, et 1 à AD
pomeroy_trans[,1] <- as.factor(pomeroy_trans[,1])
#Transformer les valeurs
pomeroy_trans<-apply(pomeroy_trans, 2,as.numeric)
#Transformer la table pomeroy en data frame
pomeroy_frame <- data.frame(pomeroy_trans)
# Scaling et Mean normalisation
s= apply(pomeroy_trans, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
View(s)

###########################################################################
#                          ACP                                            #
###########################################################################

library(FactoMineR)
#ACP sur les données centrées réduites
#PCA sans enlever les variables corrélées
pomeroy.pca <- PCA(pomeroy_trans[,-1], scale.unit = TRUE, graph = FALSE)
pca.pomeroycoord = cbind(pomeroy.pca$ind$coord[,1],pomeroy.pca$ind$coord[,2],pomeroy_frame[,1])
plot(pca.pomeroycoord[,1],pca.pomeroycoord[,2],pch=21,bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="ACP_pomeroy")
legend("bottomleft", c("MD", "Mglio", "Rhab", "Ncer", "PNET"), pch=20, col=c("blue","yellow","red","magenta","green"))

###########################################################################
#                          LDA                                            #
###########################################################################
library(MASS)
pomeroy.lda <- lda(pomeroy_frame[,1]~., pomeroy_frame[,-1])
# Corrélation variables avec le facteur discriminant
F=as.matrix(pomeroy_frame[,-1]) %*% pomeroy.lda$scaling
plot(F, col = pomeroy_frame[,1], pch=21,bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="LDA_pomeroy")
legend("top", c("MD", "Mglio", "Rhab", "Ncer", "PNET"), pch=20, col=c("blue","yellow","red","magenta","green"))
#cc=cor(as.matrix(gordon_frame[,-1]),F)


###########################################################################
#                          MDS                                            #
###########################################################################

pomeroy.dist <- dist(pomeroy_trans[,-1], method = "euclidean")
pomeroy.cmd <- cmdscale(pomeroy.dist)
plot(pomeroy.cmd[,1],pomeroy.cmd[,2],pch=21,bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="MDS_pomeroy")
legend("topright", c("MD", "Mglio", "Rhab", "Ncer", "PNET"), pch=20, col=c("blue","yellow","red","magenta","green"))


###########################################################################
#                          LLE                                            #
###########################################################################
library(RDRToolbox)
library(MASS)
par(mfrow=c(2,3))

pomeroy_lle = LLE(data=pomeroy_trans[,-1], dim=2, k=3)
plotDR(data=pomeroy_lle, labels=pomeroy_trans[,1]) 
legend("bottomright",c("K=3"))

pomeroy_lle = LLE(data=pomeroy_trans[,-1], dim=2, k=5)
plotDR(data=pomeroy_lle, labels=pomeroy_trans[,1])
legend("topright",c("K=5"))

pomeroy_lle = LLE(data=pomeroy_trans[,-1], dim=2, k=8)
plotDR(data=pomeroy_lle, labels=pomeroy_trans[,1])
legend("topleft",c("K=8"))

pomeroy_lle = LLE(data=pomeroy_trans[,-1], dim=2, k=10)
plotDR(data=pomeroy_lle, labels=pomeroy_trans[,1])
legend("bottomright",c("K=10"))

pomeroy_lle = LLE(data=pomeroy_trans[,-1], dim=2, k=12)
plotDR(data=pomeroy_lle, labels=pomeroy_trans[,1])
legend("bottomright",c("K=12"))

pomeroy_lle = LLE(data=pomeroy_trans[,-1], dim=2, k=15)
plotDR(data=pomeroy_lle, labels=pomeroy_trans[,1])
legend("bottomright",c("K=15"))


###########################################################################
#                          ISOMAP                                         #
###########################################################################
par(mfrow=c(2,3))
pomeroy_isomap = Isomap(data=pomeroy_trans[,-1],dim=2,k=4)
plot(pomeroy_isomap$dim2[,1], pomeroy_isomap$dim2[,2],pch=21, bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="K=4")

pomeroy_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=8)
plot(pomeroy_isomap$dim2[,1], pomeroy_isomap$dim2[,2],pch=21, bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="K=8")

pomeroy_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=10)
plot(pomeroy_isomap$dim2[,1], pomeroy_isomap$dim2[,2],pch=21, bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="K=10")

pomeroy_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=13)
plot(pomeroy_isomap$dim2[,1], pomeroy_isomap$dim2[,2],pch=21, bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="K=13")

pomeroy_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=16)
plot(pomeroy_isomap$dim2[,1], pomeroy_isomap$dim2[,2],pch=21, bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="K=16")

pomeroy_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=20)
plot(pomeroy_isomap$dim2[,1], pomeroy_isomap$dim2[,2],pch=21, bg=c("blue","yellow","red","magenta","green")[as.numeric(pca.pomeroycoord[,3])], main="K=20")

###########################################################################
#                          SOM                                            #
###########################################################################
library(class)
library(kohonen)
x = scale(pomeroy_trans[,-1], center = TRUE, scale = TRUE)
SOM.pomeroy = som(x, grid = somgrid(4,4, "hexagonal"))
colour1 = tricolor(SOM.pomeroy$grid)
par(mfrow=c(1,3))
plot(SOM.pomeroy, type="mapping", bgcol = rgb(colour1))
plot(SOM.pomeroy, type="mapping", labels=pomeroy_trans[,1], pchs=2)
plot(SOM.pomeroy,type="quality") # qualité de représentation

