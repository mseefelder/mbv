# Minimum number of Branch Vertices spanning tree

## Compilation

`g++ -std=c++14 backtrack-unionfind.cpp -o backtrack`

## How to use

### On bash:

```
./backtrack < your-text-file-input.txt
```

## Input format:

The program reads a stream fom standard input and expects the data to be on the following format:

* On the first line three value `n m d`, where `n` is the number of vertices, `m` is the number of edges and `d` is a dummy value that is a requirement for compatiblity reasons (with the input I'm using);
* On the following `m` lines three values on each `u v d`, where `u` and `v` are the vertices connected by the edge and `d` is another dummy value.

## Input recommendations:

I've tested the code with the dataset that [Carrabs](http://www.dipmat2.unisa.it/people/carrabs/www/) et al. used in [Lower and upper bounds for the spanning tree with minimum branch vertices](http://www.dipmat2.unisa.it/people/carrabs/www/pdf/2013_SPT_MinimumBranchVertices.pdf):

[> Dataset](www.dipmat2.unisa.it/people/carrabs/www/DataSet/MBV_Instances.zip)

## About code versions

* `backtrack.cpp`: Works. Worst implementation, uses a structure for disjoint sets that results in inefficiently deep trees and doesn't try to do any bounding or pruning by inviability;
* `backtrack-unionfind.cpp`: Works. Better than previous implementation, uses a efficient unionfind structure with rank optimization. 
* `branch-and-bound.cpp`: Works. Does bounding with a modified kruskal for every branch that prohibits a vertex, giving us pruning by inviability and bound. Improves the performance over previous implementations but has to recalculate and re-sort all edge weights every time the bound is calculated
* `branch-and-bound-no-reorder.cpp`: Works. Faster than previous implementation, also runs a kruskal for every branch that prohibits a vertex, but calculates edge weights and sorts them only once
* `branch-and-bound-no-reorder-heuristic.cpp`: Doesn't work yet. Implementing GRASP-like heuristic for MBV problem to be used as starting pont for branch-and-bound-no-reorder.
