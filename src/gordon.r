# Redirection vers le bon path
setwd("F:/M2-MLDS/Visualisation/projet")

#Chargement des données
gordon <- read.table("gordon-2002_database.txt", header = T,sep="\t")
#View(gordon)
#Transposition de la table GORDON 
gordon_trans=t(gordon)
#View(gordon_trans)
#Renommer les lignes
colnames(gordon_trans)= gordon_trans[1,]
#Enlever la premiere ligne
gordon_trans = gordon_trans[-1,]
#View(gordon_trans)
#Résumé de la table transposée
summary(gordon_trans)
#Affecter 2 à MPM, et 1 à AD
gordon_trans[,1] <- as.factor(gordon_trans[,1])
#Transformer les valeurs
gordon_trans<-apply(gordon_trans, 2,as.numeric)
#Transformer la table gordon en data frame
gordon_frame <- data.frame(gordon_trans)
View(gordon_trans)
s= apply(gordon_trans, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
View(s)

###########################################################################
#                          ACP                                            #
###########################################################################

library(FactoMineR)
#ACP sur les données centrées réduites
#PCA sans enlever les variables corrélées
gordon.pca <- PCA(gordon_trans[,-1], scale.unit = TRUE, graph = FALSE)
pca.Gordoncoord = cbind(gordon.pca$ind$coord[,1],gordon.pca$ind$coord[,2],gordon_frame[,1])
plot(pca.Gordoncoord[,1],pca.Gordoncoord[,2],pch=21,bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="ACP_GORDON")
legend("top", c("MPM","AD"), pch=20, col=c("blue","yellow"))

###########################################################################
#                          LDA                                            #
###########################################################################
library(MASS)
gordon.lda <- lda(gordon_frame[,1]~., gordon_frame[,-1])
# Corrélation variables avec le facteur discriminant
F=as.matrix(gordon_frame[,-1]) %*% gordon.lda$scaling
plot(F, col = gordon_frame[,1], pch=21,bg=c("blue","magenta")[as.numeric(pca.Gordoncoord[,3])], main="LDA_GORDON")
legend("top", c("MPM","AD"), pch=20, col=c("blue","magenta"))
#cc=cor(as.matrix(gordon_frame[,-1]),F)

###########################################################################
#                          MDS                                            #
###########################################################################

gordon.dist <- dist(gordon_trans[,-1], method = "euclidean")
gordon.cmd <- cmdscale(gordon.dist)
plot(gordon.cmd[,1],gordon.cmd[,2],pch=21,bg=c("magenta","yellow")[as.numeric(pca.Gordoncoord[,3])], main="MDS_GORDON")
legend("topright", c("MPM","AD"), pch=20, col=c("magenta","yellow"))


###########################################################################
#                          LLE                                            #
###########################################################################

library(RDRToolbox)
library(MASS)
par(mfrow=c(2,3))

gordon_lle = LLE(data=gordon_trans[,-1], dim=2, k=3)
plotDR(data=gordon_lle, labels=gordon_trans[,1]) 
legend("bottomright",c("K=3"))

gordon_lle = LLE(data=gordon_trans[,-1], dim=2, k=5)
plotDR(data=gordon_lle, labels=gordon_trans[,1])
legend("bottomright",c("K=5"))

gordon_lle = LLE(data=gordon_trans[,-1], dim=2, k=8)
plotDR(data=gordon_lle, labels=gordon_trans[,1])
legend("topright",c("K=8"))

gordon_lle = LLE(data=gordon_trans[,-1], dim=2, k=10)
plotDR(data=gordon_lle, labels=gordon_trans[,1])
legend("topright",c("K=10"))

gordon_lle = LLE(data=gordon_trans[,-1], dim=2, k=12)
plotDR(data=gordon_lle, labels=gordon_trans[,1])
legend("bottomright",c("K=12"))

gordon_lle = LLE(data=gordon_trans[,-1], dim=2, k=15)
plotDR(data=gordon_lle, labels=gordon_trans[,1])
legend("topright",c("K=15"))



###########################################################################
#                          ISOMAP                                         #
###########################################################################
par(mfrow=c(2,3))
gordon_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=4)
plot(gordon_isomap$dim2[,1], gordon_isomap$dim2[,2],pch=21, bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="K=4")

gordon_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=8)
plot(gordon_isomap$dim2[,1], gordon_isomap$dim2[,2],pch=21, bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="K=8")

gordon_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=10)
plot(gordon_isomap$dim2[,1], gordon_isomap$dim2[,2],pch=21, bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="K=10")

gordon_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=13)
plot(gordon_isomap$dim2[,1], gordon_isomap$dim2[,2],pch=21, bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="K=13")

gordon_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=16)
plot(gordon_isomap$dim2[,1], gordon_isomap$dim2[,2],pch=21, bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="K=16")

gordon_isomap = Isomap(data=gordon_trans[,-1],dim=2,k=20)
plot(gordon_isomap$dim2[,1], gordon_isomap$dim2[,2],pch=21, bg=c("blue","yellow")[as.numeric(pca.Gordoncoord[,3])], main="K=20")

###########################################################################
#                             SOM                                         #
###########################################################################

library(class)
library(kohonen)
par(mfrow=c(3,3))
x = scale(gordon_trans[,-1], center = TRUE, scale = TRUE)
SOM.GORDON = som(x, grid = somgrid(8,8, "hexagonal"))
#SOM.GORDON$codes # coordonnees des centres (espace d'entrée)
#SOM.GORDON$distances # distances entres neuronnes adjacents
#SOM.GORDON$changes # erreur entre les donnees et les centres
colour1 = tricolor(SOM.GORDON$grid)
#par(mfrow=c(2,3))
plot(SOM.GORDON, type="mapping", bgcol = rgb(colour1))
plot(SOM.GORDON, type="mapping", labels=gordon_trans[,1], pchs=2)
plot(SOM.GORDON,type="quality") 
plot(SOM.GORDON,type="counts") 
plot(SOM.GORDON,type="codes")

# Modification de la palette de couleurs (bleu -> rouge)
coolBlueHotRed = function(n, alpha = 1)
{
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}
# U-matrix (matrice de voisinage)
plot(SOM.GORDON, type="dist.neighbours",palette.name=coolBlueHotRed)
# Component planes
par(mfrow=c(2,2))
for (i in 2:10)
{
  plot(SOM.GORDON, type="property", palette.name=coolBlueHotRed,
       property=SOM.GORDON$codes[,i], main=colnames(SOM.GORDON$codes)[i])
}

