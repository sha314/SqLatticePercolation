
## Version 14 is for memory optimization

Command to find peak memory usages
$ /usr/bin/time -f "%M" ./sqlattice 1024 5
or,
$ time -f "%M" ./sqlattice 1024 5
for L=1024 and 5 runs

### Memory Usage by V13
Length PeakMemory
100
500
1000
2000



### Memory Usage by V14
Length  PeakMemory(KB)
100     8480
500     111508
1000    432592
2000    1734244


### Memory Usages by Lattice of v14
L=1000 Mem=390340 KB = 381 MB
L=2000 Mem=1550680 KB = 1514 MB

