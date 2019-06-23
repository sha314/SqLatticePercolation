# Percolation on Square Lattice
## Version 9

#### Previously
In previous versions we would merge two clusters by copying all
 elements of one cluster to another then delete the smaller 
 cluster. After deleting one cluster we would have to keep track 
 of id and index of all clusters separately using binary search tree
  (BST) of C++, i.e. std::map, and a linked list. We called that InverseArray.
   Which would take most of the time.
   
#### Paradigm Shift
We do not delete cluster only make it empty so that index and id remains the same.
This way Searching for index of the cluster using id can be omitted.
Performance improved by a considerable amount.
       


#### General Rule
 * First it randomizes the site index list then use it.

 

#### Feature :
1. Can turn on and off both horizontal and boundary condition

2. Uses class "Cluster_v3" for storing clusters

3. Uses Group_ID for Bonds and Sites to identify that they are in the same cluster

4. Occupation probability is calculated by site_index_sequence,
  i.e. number of active site_index_sequence divided by total number of site_index_sequence, in site percolation
 
5. Spanning is calculated by number of bonds in a spanning clusters with periodicity turned off,
       i.e., number of bonds in the spanning clusters divided by total number of bonds
 
6. Weighted and Unweighted relabeling is available

7. Runtime is significantly improved. For example, if L=1000 program will take ~5 sec to place all site_index_sequence.
