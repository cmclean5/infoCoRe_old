##reset
rm(list=ls())

library(infoCoRe)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## 
adj=get.adjacency(gg)

infoCoRe::test(Adj=as(adj,"generalMatrix"))
