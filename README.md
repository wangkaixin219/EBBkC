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

After running the codes, there will be two executable files: one is called `BBkC` and the other is called `truss`, which means you have already compiled our codes. 

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

The pre-processing procedure will output a file with the prefix same as the original file but end with `.clean`. 

We will use the `.index` file for the listing procedure. 

## Step 2 - Serial Listing Procedure

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

## Step 3 - Parallel Listing Procedure

To list $k$-cliques in the graph with multiple threads, we need to specify an additional parameter to indicate the number of threads to be used. Therefore, the running command is: 

```bash
./BBkC ep /PATH_TO_INDEX_DATA k_val t_val p_val
```

For example, 

```bash
./BBkC ep ../../dataset/nasasrb.index 12 3 16   # list 12-clique in `nasasrb` with early-termination in 3-plex with 16 threads
```

It should outputs

```
Reading edges from ../../dataset/nasasrb.index ...
|V| = 54567, |E| = 1299405
Truss number = 23
Building necessary data structure ...
Iterate over all cliques
Thread: 9, Runtime = 881.88 ms, handled 80552 edges
Thread: 2, Runtime = 881.33 ms, handled 82149 edges
Thread: 10, Runtime = 881.87 ms, handled 83261 edges
Thread: 12, Runtime = 881.31 ms, handled 80916 edges
Thread: 13, Runtime = 881.61 ms, handled 82643 edges
Thread: 0, Runtime = 881.96 ms, handled 87243 edges
Thread: 15, Runtime = 881.17 ms, handled 73071 edges
Thread: 1, Runtime = 881.66 ms, handled 80179 edges
Thread: 5, Runtime = 881.73 ms, handled 76805 edges
Thread: 7, Runtime = 881.75 ms, handled 80162 edges
Thread: 8, Runtime = 882.12 ms, handled 79822 edges
Thread: 14, Runtime = 881.33 ms, handled 81962 edges
Thread: 3, Runtime = 881.80 ms, handled 82818 edges
Thread: 11, Runtime = 881.95 ms, handled 81079 edges
Thread: 6, Runtime = 881.88 ms, handled 83730 edges
Thread: 4, Runtime = 882.28 ms, handled 83013 edges
Number of 12-cliques: 7905736209
EBBkC+ET (t = 3) runtime 882.28 ms
```
