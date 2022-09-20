#create fasta file from protein sequences
library(seqinr)
setwd("c:/users/mudoe/desktop/hetexp working")
seq = read.csv("sequences.csv")
gpos = read.csv("workflow1/bleedproc/6-21gpos.csv", header = FALSE)

colnames(gpos) = gpos[1,]
gpos=gpos[-1,]
df = read.csv("connectors.csv")
df[,1] = NULL
for (i in 1:nrow(seq)) {
  seq[i,1] = gsub(".*_", "", seq[i,1])
  seq[i,7] = gsub("\\*","",seq[i,7])
}
seq$Original.Protein.Sequence = gsub("*","",as.character(seq$Original.Protein.Sequence))
for (i in 1:nrow(seq)) {
  seq[i,1] = unique(df[which(df[,1] %in% seq[i,1]),3])
}

write.fasta(as.list(seq$Original.Protein.Sequence),as.list(seq$JGI.name),
            "c:/users/mudoe/desktop/locus.fasta",as.string=TRUE)

#visualize and compare trees
library("ape")
library("phytools")
library("ggtree")
setwd("c:/users/mudoe/desktop/hetexp working/trees")
tre = read.tree("clap dendro.tre")
tre2 = read.tree("clustal2.tre")
comp = comparePhylo(tre,tre2)
cophylo = cophylo(tre, tre2)
plot(cophylo)

#create heat map
setwd("c:/users/mudoe/desktop/hetexp working")
df = read.csv("connectors.csv")
df[,1] = NULL
seq = read.csv("sequences.csv")
gpos = read.csv("c:/users/mudoe/desktop/full desktop/orderedforfigure.csv", header = FALSE)
###not dynamic, must count number of hydrogenated feature rows###
gpos = gpos[1:22,]
colnames(gpos) = gpos[1,]
gpos=gpos[-1,]
rownames(gpos) = gpos[,1]
gpos[,1] = NULL
gpos = gpos+1
gpos = log10(gpos)
relg = gpos[,1:366]

library("ape")
library("phytools")
#library("ggtree")
setwd("c:/users/mudoe/desktop/hetexp working/trees")
tre = read.tree("clap dendro.tre")
tre2 = read.tree("c:/users/mudoe/desktop/full desktop/matchplasmid.dnd.tre")
order = tre2[[4]]
filler <- matrix(nrow=nrow(relg),
                 ncol=2*length(order))
aces = 1
for (i in seq(1,ncol(filler),by=2)) {
  filler[1,i] = order[aces]
  filler[1,i+1] = order[aces]
  aces = aces + 1
}
filler = data.frame(filler)
colnames(filler) = filler[1,]
x = colnames(filler)
y = colnames(relg)

df2 = df
#for (i in 1:3) {
for (i in 1:ncol(filler)) {
  #for (j in 1:ncol(relg)) {
  zep = which(df2[,1] %in% x[i])
  run = df2[zep,2][1]
  df2[zep[1],1] = "pulled"
  cub = y[which(y %in% run)]
  filler[,i] = relg[cub]
  #}
}
rownames(filler) = rownames(relg)

library(ggplot2)
library(gplots)
library(reshape)
col <- colorRampPalette(c("royalblue4","royalblue4","steelblue1","beige","red2","brown"))(100)
heatmap.2(as.matrix(filler), Rowv = FALSE, Colv = FALSE, density.info = "none", col = col,trace = "none", margins = c(7, 15), labCol = NA )

png("c:/users/mudoe/desktop/guy3.png",width=2000,height=800)
heatmap.2(as.matrix(filler), Rowv = FALSE, Colv = FALSE, density.info = "none", col = col,trace = "none", margins = c(2.2, 4.8),labCol = NA, keysize = 1)
dev.off()
plot(tre2)
tre2
#dendrogram and heatmap must be combined in another program

####mass spec fig
data = read.csv("c:/users/mudoe/desktop/chromat data.csv")
altdata = data
altdata[,2] = log10(altdata[,2] + 1)
altdata[,5] = log10(altdata[,5] + 1)
altdata[,8] = log10(altdata[,8] + 1)
first = altdata[,c(1,2)]
sec = altdata[,c(4,5)]
thir = altdata[,c(7,8)]

library(ggplot2)
library(reshape)


png(file="C:/users/mudoe/desktop/take3.png",
    width=750, height=504)
plot(altdata$Retention.time.2,altdata$M1_C1.intensity,type="l",col="green",
     xlim = c(2.5,4.5),ylim=c(2.5,5.5),lwd = 2,xlab = "RT",ylab = "log10 intensity")
lines(altdata$Retention.time,altdata$X3OHC10.intensity,type="l",col="red",
      lwd = 2)
lines(altdata$Retention.time.1,altdata$AWP201.intensity,type="l",col="blue",
      lwd = 2)
legend("topleft", legend = c("Heterologous Expression", "3OHC10","AWP201"),
       col = c("green", "red","blue"), lwd = 2, cex = .8)
dev.off()

plot(data$Retention.time.2,data$M1_C1.intensity,type="l",col="green",
     lwd = 2,xlab = "RT",ylab = "intensity")
lines(data$Retention.time,data$X3OHC10.intensity,type="l",col="red",
      lwd = 2)
lines(data$Retention.time.1,data$AWP201.intensity,type="l",col="blue",
      lwd = 2)
legend("topleft", legend = c("Heterologous Expression", "3OHC10","AWP201"),
       col = c("green", "red","blue"), lwd = 2, cex = .8)
