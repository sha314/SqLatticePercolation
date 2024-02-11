Percolation on a Square Lattice


Method to run the program:
1. using commands as described in cmd.args.txt file
2. using a "source form" file. You take the form file and fill it up.
    Then pass the name of the source file as command line argument.
    The program reads the file and run the program as described in the
    form. If no form is supplied as argument then the program creates
    "input.form.txt" file which contains all information which to be filled.


## v13 vs v14
Looks like v13 is faster in Mac OS than v14. While v14 is specifically built using memory optimization in mind at the cost of little CPU usages.


# How to Build (Terminal)
1. Pull the repository using command
$ git clone https://github.com/sha314/SqLatticePercolation.git

2. Go to the SqLatticePercolation directory and create a new folder where the executable will be built.
$ cd SqLatticePercolation
$ mkdir build
$ cd build

3. Run cmake and make in build folder. I am using 6 threads to build it.
$ cmake ..
$ make -j 6

4. Run the executable
$ ./sqlattice



# Installing Required Tools to properly build this program (g++, gdb, cmake, make)
1. Ubuntu OS

$ apt install git cmake make

2. Fedora OS

$ dnf install cmake g++
$ dnf groupinfo "Development Tools"
$ dnf groupinstall "Development Tools"

2. Mac OS



3. Manjaro OS




