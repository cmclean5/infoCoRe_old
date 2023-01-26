##reset
rm(list=ls())

library(infoCoRe)
library(igraph)

## load graph
##gg = read.graph("PPI_Presynaptic.gml",format="gml")

## 
##adj=get.adjacency(gg)

##infoCoRe::test(Adj=as(adj,"generalMatrix"))

## Directed Graph Example
ba = sample_pa(n=100, power=1, m=1,  directed=T)

adj=get.adjacency(ba)

infoCoRe::test2(Adj=as(adj,"generalMatrix"))
