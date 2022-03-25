
## Version 14 is for memory optimization

Command to find peak memory usages
$ /usr/bin/time -f "%M" ./sqlattice 1024 5
or,
$ time -f "%M" ./sqlattice 1024 5
for L=1024 and 5 runs

### Memory Usage by V13
L=2000
Mem=2.68 GB
Time=8.5 sec per iteration



### Memory Usage by V14
L=2000
Mem=1.53 GB
Time=6.7 sec per iteration


### Memory Usages by Lattice of v14
L=1000 Mem=390340 KB = 381 MB
L=2000 Mem=1550680 KB = 1514 MB


### Log 2022.03.24
Just creating SitePercolationL0_v13 object takes 1.02 GB RAM. then another ~500 MB for running the rest of the code.
I think there are multiple copies of same things in the program, need to clean it.

commented "first_nn_count" and "second_directional_nn_count" in Site.

L=2000 Mem=1.50 GB RAM by "run_v13_rsbd_L0" method using V14 Site/Bond. 6.3 sec per iteration.
L=2000 Mem=2.67 GB RAM by "run_v13_rsbd_L0" method using V13 Site/Bond. 8.5 sec per iteration.

L=2000 Mem=1.48 GB RAM by "run_v13_rsbd_L1" method using V14 Site/Bond. 11.8 sec per iteration.
L=2000 Mem=2.66 GB RAM by "run_v13_rsbd_L1" method using V13 Site/Bond. 30.0 sec per iteration.

So, V14 is both memory efficient and fast. Reason: it takes non-zero finite amount of time to get data from RAM and if you do it
millions of times then the total time is huge and does not improve performance. If it took zero time to get data from RAM, then, I guess,
v13 would have performed better.


### log 2022.03.24
Using "first_nn_count" and "second_directional_nn_count" in Site. and V14 Site/Bond we get
L=2000 Mem=1.51 GB RAM by "run_v13_rsbd_L1" method using V14 Site/Bond. 9.76 sec per iteration. While, L1 is almost unchanged.
So, using it is better. Check for for L2. TODO
