
# Sometimes we want to run the program for multiple lengths. This bash script can make things easier. It will generate 
# a unique log file for each lengths and each time you run it.
En=2000

for L in 128 256 512
do
	./sqlattice $L $En >> log$L$(date +"%Y%m%d-%H%M%S") &
	echo $L
done

