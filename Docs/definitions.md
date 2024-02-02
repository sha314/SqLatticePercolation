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



  | 0     | 1     | 2     | 3     |
---------------------------------
0 |       |       | <0,1> | <1,1>
1 |       |       | <0,0> |
2 |
4 |
