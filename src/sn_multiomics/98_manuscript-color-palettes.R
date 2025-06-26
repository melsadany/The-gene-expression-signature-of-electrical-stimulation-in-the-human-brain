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
# cell.colors.2 <- c("Astro" = "#ef7cef", "Micro" = "#e99633", "Oligo"= "#ff8eb9",
#                    "OPC" = "#beb83b", "Exc"= "#42cf8e", "Inh" = "#a3afff")
cell.colors.2 <- c("Astro" = "#C77CFF", "Micro" = "#FF7F00", "Oligo"= "#FF61CC",
                   "OPC" = "#7CAE00", "Exc"= "#00BA38", "Inh" = "#619CFF")
# 
# c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#FF61CC","#00BA38","#619CFF","#F564E3")
# [1]  "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00"
# [6]  "#FFFF33" "#A65628" "#F781BF" "#999999" "#66C2A5"
# [11] "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"
# [16] "#E5C494" "#B3B3B3" "#8DD3C7" "#FFFFB3" "#BEBADA"
# [21] "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5"
# [26] "#D9D9D9"


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