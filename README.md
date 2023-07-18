# Efficient $k$-clique Listing: An Edge-Oriented Branching Strategy

This repository contains the codes and a technical report for $k$-clique listing with the edge-oriented branching strategy. 

## Setup

To run the codes correctly, please ensure you have installed `cmake (minimum version 3.16)` on your computer or server since we need to use cmake to generate `Makefile` to compile our codes. 
There are two steps to run our codes. The first step is a pre-processing procedure, which is to ensure that the graph contains no self-loops and no duplicate edges. The second step is the main listing procedure, which outputs all $k$-cliques. 

## Step 0 - Compile the codes

When you already download the codes, run the following commands to compile our codes. 

```bash
cd src/
mkdir build && cd build/
cmake ..
make
```

After running the codes, there will be an executable file called `BBkC`, which means you have already compiled our codes. 

## Step 1 - Pre-processing Procedure

When the graph is downloaded from the [Graph Repository](https://networkrepository.com/), the graph may have some comments at the very beginning, please ensure those comments are deleted before running the codes. The graph file may contain a number of edges, with each edge occupying a line (e.g., [nasasrb.edges](./dataset/nasasrb.edges)). Besides, the edges may also have some other attributes (e.g., [facebook.edges](./dataset/facebook.edges)). Our codes can handle both cases by only focusing on the first two elements of each edge, which represent the endpoints of the edge. 

To pre-process the graph, run the following command: 

```bash
./BBkC p /PATH_TO_EDGE_DATA
```

For example, 

```bash
./BBkC p ../../dataset/facebook.edges     # pre-process the `facebook` dataset
./BBkC p ../../dataset/nasasrb.edges      # pre-process the `nasasrb` dataset
```

For `facebook` dataset, since it contains many duplicate edges, it outputs

```
Pre-process ../../dataset/facebook.edges
... duplicates.
Pre-processed in 3917.07 ms
```

For `nasasrb`, since it contains no duplicate edges and no self-loops, it outputs

```
Pre-process ../../dataset/nasasrb.edges
Pre-processed in 4535.87 ms
```

The pre-processing procedure will output two files with the prefix same as the original file but end with `.clean` and `.index`, respectively. We will use the `.index` file for the listing procedure. 

## Step 2 - Listing Procedure

To list $k$-cliques in the graph, we need to specify the value of $k$. In addition, we also need to specify a value of $t$ to control when we should conduct the early termination (details can be found in our paper). Therefore, the running command is: 

```bash
./BBkC e /PATH_TO_INDEX_DATA k_val t_val
```

For example, 

```bash
./BBkC e ../../dataset/facebook.index 20 3    # list 20-clique in `facebook` with early-termination in 3-plex 
./BBkC e ../../dataset/nasasrb.index 12 2     # list 12-clique in `nasasrb` with early-termination in 2-plex
```

For `facebook` dataset, it outputs

```
Reading edges from ../../dataset/facebook.index ...
|V| = 2543, |E| = 62844
Truss number = 35
Building necessary data structure ...
Iterate over all cliques
Number of 20-cliques: 1163924709
EBBkC+ET (t = 3) runtime 1391.73 ms
```

For `nasasrb` dataset, it outputs

```
Reading edges from ../../dataset/nasasrb.index ...
|V| = 54567, |E| = 1299405
Truss number = 23
Building necessary data structure ...
Iterate over all cliques
Number of 12-cliques: 7905736209
EBBkC+ET (t = 2) runtime 12683.40 ms
```

