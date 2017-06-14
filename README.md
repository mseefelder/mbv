# Minimum number of Branch Vertices spanning tree

## Compilation

`g++ -std=c++14 backtrack-unionfind.cpp -o backtrack`

## How to use

### On bash:

```
./backtrack < your-text-file-input.txt
```

## Input format:

The program reads s tream fom standard input and expects the data to be on the following format:

* On the first line three value `n m d`, where `n` is the number of vertices, `m` is the number of edges and `d` is a dummy value that is a requirement for compatiblity reasons (with the input I'm using);
* On the following `m` lines three values on each `u v d`, where `u` and `v` are the vertices connected by the edge and `d` is another dummy value.

## Input recommendations:

I've tested the code with the dataset that [Carrabs](http://www.dipmat2.unisa.it/people/carrabs/www/) et al. used in [Lower and upper bounds for the spanning tree with minimum branch vertices](http://www.dipmat2.unisa.it/people/carrabs/www/pdf/2013_SPT_MinimumBranchVertices.pdf):

[> Dataset](www.dipmat2.unisa.it/people/carrabs/www/DataSet/MBV_Instances.zip)
