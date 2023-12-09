##reset
rm(list=ls())

library(infoCoRe)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## 
adj=get.adjacency(gg)

infoCoRe::driver(Adj=as(adj,"generalMatrix"),
                 weighted=0,
                 directed=0,
                 val_only=1,
                 norm=1)

## Directed Graph Example
## http://www-personal.umich.edu/~mejn/netdata/
##gg = sample_pa(n=100, power=1, m=1,  directed=T)
gg = read.graph("polblogs.gml",format="gml")

adj=get.adjacency(gg)

infoCoRe::driver(Adj=as(adj,"generalMatrix"),
                 weighted=0,
                 directed=1,
                 val_only=1,
                 norm=1)

## use igraph's c libraries to decompose a graph
infoCoRe::test(Adj=as(adj,"generalMatrix"),
               directed=is.directed(gg),
               loops=0
               )

## use igraph's R decompose function
cat("> No: of components (igraph R) =", length(decompose(gg)), "\n")
