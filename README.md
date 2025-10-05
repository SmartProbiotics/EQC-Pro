# Maximum Edge-Based Quasi-Clique Computation

This project implements maximum edge-based quasi-clique algorithm——`EQC-Pro`

## Compilation

```sh
$ make clean
$ make
```
It generates an executable `EQC`(i.e. `EQC-Pro`)、`EQC-TD`、`EQC-NH`.

## Execution

To run the program, use the following command:
```sh
$ ./EQC {path_to_graph} {gamma}
```

```sh
$ ./EQC ./data/bin/soc-delicious.bin 0.95
```

### Different variants of the algorithm can be run as follows.

#### 1. Run the algorithm EQC-TD

```sh
$ ./EQC-TD {path_to_graph} {gamma}
```

#### 2. Run the algorithm EQC-NH

```sh
$ ./EQC-NH {path_to_graph} {gamma}
```
#### 3. Run the algorithm EQC-Pro
```sh
$ ./EQC  {path_to_graph} {gamma}
```

## Data format
### Format : `*.bin`
binary graph:
- first $4$ Bytes: **sizeof(uint32_t)**, which should be $4$
- then $4$ Bytes: $n$
- then $4$ Bytes: $2\times m$
- then $4\times n$ Bytes: the degree $d_G(\cdot)$ of $n$ vertices
- then: $n$ parts ($2m\times 4$ Bytes in total), each part has $d_G(u)$ integers which are the neighbors of $u$ ***in ascending order***

We provide some examples of binary graph file  in [**data**](.data/bin).
### Transform other formats to `*.bin`
We provide [**some C++ code**](./data/change_graph_format) to convert four different formats into binary files.\
More details about it can be found in [**data/change_graph_format/readme.txt**](./data/change_graph_format/readme.txt)

## Get datasets
Real-world graphs collection: http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz

Facebook graphs collection: https://networkrepository.com/socfb.php
