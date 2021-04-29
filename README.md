# Privacy-preserving DBSCAN and TRACLUS

By Beyza Bozdemir¹, Sébastien Canard², Orhan Ermis¹, Helen Möllering³, Melek Önen¹, and Thomas Schneider³. (¹[EURECOM](http://www.eurecom.fr); ²[Applied Cryptography Group](https://crypto.orange-labs.fr/), Orange Labs; ³[ENCRYPTO](http://www.encrypto.de), TU Darmstadt) <br>In [16th ACM ASIA Conference on Computer and Communications Security (ACM ASIACCS 2021)](https://asiaccs2021.comp.polyu.edu.hk/). <!--%[Paper available here.](https://encrypto.de/papers/) -->

This work is based on the [ABY framework](https://github.com/encryptogroup/ABY/) for efficient mixed-protocol secure two-party computation.

For simplicity, we will only describe the short version of the build process here, including the differences and additional steps required compared to the ABY build procedure.
Please also refer to the detailed build instructions of ABY.

This code is provided as an experimental implementation for testing purposes and should not be used in a productive environment. We cannot guarantee security and correctness.

### (Short) Requirements and Installation

This code requires all [ABY requirements](https://github.com/encryptogroup/ABY#requirements) and can be set up and executed like any other [ABY example](https://github.com/encryptogroup/ABY#aby-applications).

### Parameters Distance Calculation
```Usage: ./asv_test
 -r [Role: 0/1, required]
 -d address of dataset file
 -e input dataset size
```
More (optional) parameters can be inspected by simply running `./ppDBSCAN_distance`.

### Parameters Grouping
```Usage: ./asv_test
 -r [Role: 0/1, required]
 -d address of dataset file
 -e input dataset size
 -i index of input record that is currently checked if it creates a cluster (including the cluster expansion by neighbors). 
```
More (optional) parameters can be inspected by simply running `./ppDBSCAN_grouping`.

### Running ppDBSCAN on Lsun
For an easy test of our ppDBSCAN implementation, we included two scripts ([1](https://github.com/encryptogroup/ppDBSCAN/blob/main/client_script.sh) and [2](https://github.com/encryptogroup/ppDBSCAN/blob/main/server_script.sh)), each for one of the two computing parties, that run ppDBSCAN on the Lsun dataset [[1]](#1) scaled to integer values. 

As visible in the scripts, ppDBSCAN is split into distance calculation and grouping. Furthermore, the grouping is interrupted after checking each element of the input dataset to reset the ABY circuit and free memory occupied from ABY. This is realized by writing the secret shares at each of the two computing parties into a local file at the respective party. Thereby, nothing about the result is leaked as all data is kept as secret-shares.

The [Lsun data file](https://github.com/encryptogroup/ppDBSCAN/blob/main/data/Lsun_prepared.txt) must be placed in `ABY/build/bin/data` to be automatically found when executing the two scripts.

### Clustering other datasets
To cluster different datasets, the values of ppDBSCAN's 'eps' and 'minPts' as well as the dataset's dimension (indicated by 'dim') should be adapted in this [file](https://github.com/encryptogroup/ppDBSCAN/blob/main/src/examples/ppdbscan/common/ppdbscan_coordination.cpp). Additionally, we recommend to re-use the scripts for clustering Lsun by simply adapting the parameters for the dataset's size and changing the file's address. The dataset size must also be used for the number of grouping iterations (i.e., the number of calls to ppDBSCAN_grouping) as each iteration checks for one element if it creates a cluster including the recursive expansion through potential neighbors. Further, if the new dataset is not given in the same format, the function `loadData(..)` in this [file](https://github.com/encryptogroup/ppDBSCAN/blob/main/src/examples/ppdbscan/common/ppdbscan_coordination.cpp) should be adapted.

### Issues & Notes
ABY has a known problem where it randomly returns invalid results in about 1 out of 10 executions. This problem is still an open and also reported in several issues, e.g., [here](https://github.com/MPC-SoK/frameworks/issues/19) or [here](https://github.com/encryptogroup/ABY/issues/114). Given that we have to do many runs of an ABY program in the splitted grouping phase, it follows that our prototype can only give valid runtimes but not valid clustering results. To get valid clustering results, the split of the grouping phase must be removed.

### References
<a id="1">[1]</a> 
Ultsch, A. (2005).
Clustering with SOM: U*C. 
Proc. Workshop on Self-Organizing Maps, Paris, France, pp. 75-82.
