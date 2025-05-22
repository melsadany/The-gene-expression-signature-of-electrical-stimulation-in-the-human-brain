#### used manuscript color palettes


## Muhammad's color palettes
boxplot.colors <- c("#aaf0d1", "#b39eb5")
redblu.col <-  c("#ff6961", "#89cff0")
redblack.col <- c("#800000", "black")
six.colors <- c("#800000", "#cc7277", "#4f6162", "#e65236", "#56483a", "#73937e")
ten.colors <- c("#800000", "#cc7277", "#4f6162", "#e65236", "#56483a", 
                "#73937e", "#06241b", "#b8860b", "#e07c4c", "#9a81b0")
antique.colors <- c("#855C75","#D9AF6B","#AF6458","#736F4C","#526A83",
                    "#625377","#68855C","#9C9C5E","#A06177","#8C785D",
                    "#467378","#7C7C7C")
redblu.ni.col <- c("#E43339", "#0033FF")
warm.cold.col <- c("blue", "deepskyblue", "white", "yellow", "orange", "red")
abstract.colors <- c("#CC0A7D", "#EAA91E", "#A46CE1", "#7B8CB2", "#798632", "#4A9A3B", "#FFD9B8", "#781285", "#C746FF", "#F9FF32")
abstract.colors.2 <- c("#86324A", "#333286", "#798632", "#4A9A3B", "#643264", "#3B39D2", "#DCD813", "#00FF00", "#DD9FDD", "#2D8B57", "#86CDEB")


## celltype colors
cell.colors <- c(antique.colors[c(1,5:7,3,4,2,8)])
names(cell.colors) <- c("Astro", "Micro", "Oligo", "OPC", "Exc", "Inh", "Endo", "VLMC")


## participant labels colors
sample.colors <- abstract.colors[c(1,3:7,11)]
names(sample.colors) <- paste0("P", c(4,5,6,18,19,20,22))
sample.colors.2 <- abstract.colors[c(1,3:7,11)]
names(sample.colors.2) <- c(paste0("E", c(3,1,2,4)),paste0("G", c(1,3,2)))


## stimulation colors
stim.colors <- c(six.colors[1], antique.colors[12])
names(stim.colors) <- c("stimulation", "baseline")

## dx colors
dx.colors <- c(antique.colors[11:12])
names(dx.colors) <- c("epilepsy", "glioma")

## rna and atac colors
mod.colors <- c("#f2aa84", "#a6caec")
names(mod.colors) <- c("RNA-Seq", "ATAC-Seq")

## significant genes
sig.colors <- c(abstract.colors.2[1:3],"grey")
names(sig.colors) <- c("DEG in both", "edgeR DEG ONLY", "lmmSeq DEG ONLY", "ns")

## 