## 1# mean cluster size
It is the ratio of the second and the first moment of n_s. That is, 
\sum_s s^2 n_s/sum_s s n_s 
excluding the largest or the spanning cluster.

\sum_s s^2 n(s) / \sum_s s n(s)

basically, (s_1^2+s_2^2+...+s_k^2)/(s_1+s_2+---+s_k) but you exclude the largest or the spanning cluster.
do it for each value of p.


## 2# Cluster Size
In bond percolation, the number of sites in a cluster is the cluster size. 
In New definition of site percolation(Our Definition), the number os bonds in a cluster is the cluster size.


## 3a# Site Index (r,c)
We place the sites in a grid and use cartesian coordinate system here, except for the fact that the coordinate system is rotated by 90 degree clockwise. so that 
horizontal axis -> y axis and y value increasing toward right
vertical axis -> x axis and x value increasies downward

This convention helps with relative indexing.

## 3b# Bond Index
Since we have twice as many bonds...


## 4# Relative Index
To check a wrapping we make use of relative index.

The root site(The site of a cluster that was the first to be occupied) is labeled <0,0>
the site on it's right will have index <0,1> if it belongs to the same cluster. In the following 4x4 grid, the site (1,2) is the root site


Upper Left corner is <0,0>
<x,y> means relative index
90 degree Clockwise rotated Coordinate system
y value increases as we go rightward. Like columns
x value increases as we go downward . Like rows
Format 'gid<x,y>'


-- y -->
  | 0     | 1     | 2     | 3     |
---------------------------------
0 |       |       | <-1,0> | <-1,1>
1 |       |       | <0,0> | <0,1>
2 |
4 |

### 4a. Issue. BondPercolation_v13
    length=5;
    seed=2;

After 4 Iteration I've got
----------------------------------------------------------------------
    |     0     |     1     |     2     |     3     |     4     |
----------------------------------------------------------------------
  0 |21<  0, -1>|21< -4,  0>| 2<  0,  0>| 3<  0,  0>| 4<  0,  0>|
  1 | 5<  0,  0>|11< -1,  0>| 7<  0,  0>| 8<  0,  0>|14< -1,  0>|
  2 |10<  0,  0>|11<  0,  0>|12<  0,  0>|13<  0,  0>|14<  0,  0>|
  3 |15<  0,  0>|16<  0,  0>|17<  0,  0>|18<  0,  0>|19<  0,  0>|
  4 |20<  0,  0>|21<  0,  0>|22<  0,  0>|23<  0,  0>|24<  0,  0>|
----------------------------------------------------------------------

While I should get (for cluster with id 21)
----------------------------------------------------------------------
    |     0     |     1     |     2     |     3     |     4     |
----------------------------------------------------------------------
  0 |21<  1, -1>|21<  1,  0>| 2<  0,  0>| 3<  0,  0>| 4<  0,  0>|
  1 | 5<  0,  0>|11< -1,  0>| 7<  0,  0>| 8<  0,  0>|14< -1,  0>|
  2 |10<  0,  0>|11<  0,  0>|12<  0,  0>|13<  0,  0>|14<  0,  0>|
  3 |15<  0,  0>|16<  0,  0>|17<  0,  0>|18<  0,  0>|19<  0,  0>|
  4 |20<  0,  0>|21<  0,  0>|22<  0,  0>|23<  0,  0>|24<  0,  0>|
----------------------------------------------------------------------


## 5# Merging Clusters. Bond Percolation
A bond connects two sites.
1. Fint the two connected sites by the selected bond.
2. Find the group id of those sites to figure out which cluster they belong to.
3. Relabel site and bond gids of smaller cluster.
4. Relabel relative index of sites of smaller cluster.
5. Add bonds and sites ids to larger cluster.
6. Clear the smaller cluster.


