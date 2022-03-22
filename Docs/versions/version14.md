
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

