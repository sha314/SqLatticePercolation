cmake_minimum_required(VERSION 3.0)
project(SqLatticePercolation)

set(CMAKE_CXX_STANDARD 11)

#set (CMAKE_C_COMPILER               /usr/bin/gcc)
#set (CMAKE_CXX_COMPILER             /home/shahnoor/software/pgi/linux86-64-llvm/2018/bin/pgc++)
#set (CMAKE_MAKE_PROGRAM             /usr/bin/make)
SET( CMAKE_CXX_FLAGS  "-pthread -fopenmp")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

set(SOURCE_FILES
        src/cmdLineArgs/cmd_args.cpp
        src/cmdLineArgs/cmd_args.h
        src/cmdLineArgs/cmd_args_json.cpp
        src/cmdLineArgs/cmd_args_json.h
        src/cmdLineArgs/percolation_fractal.cpp
        src/cmdLineArgs/percolation_fractal.h
        src/cmdLineArgs/sq_lattice_methods.h
        src/ensemble/site_percolation.cpp
        src/ensemble/site_percolation.h
        src/exception/exceptions.cpp
        src/exception/exceptions.h
        src/ext_libs/inverse_array.h
        src/fractal/fractal_tool.cpp
        src/fractal/fractal_tool.h
        src/index/delta.cpp
        src/index/delta.h
        src/index/id.cpp
        src/index/id.h
        src/index/index.cpp
        src/index/index.h
        src/lattice/bond.cpp
        src/lattice/bond.h
        src/lattice/bond_type.h
        src/lattice/connector.cpp
        src/lattice/connector.h
        src/lattice/lattice.cpp
        src/lattice/lattice.h
        src/lattice/site.cpp
        src/lattice/site.h
        src/percolation/boundary.cpp
        src/percolation/boundary.h
        src/percolation/calculation.cpp
        src/percolation/calculation.h
        src/percolation/cluster.cpp
        src/percolation/cluster.h
        src/percolation/cluster_v2.cpp
        src/percolation/id_index.cpp
        src/percolation/id_index.h
        src/percolation/percolation.cpp
        src/percolation/percolation.h
        src/percolation/percolation_bond_v0.cpp
        src/percolation/percolation_site_ballistic_deps.cpp
        src/percolation/percolation_site_explosive.cpp
        src/percolation/percolation_site_v8.cpp
        src/percolation/percolation_type.h
        src/percolation/site_position.cpp
        src/percolation/site_position.h
        src/tests/test_bonds.cpp
        src/tests/test_bonds.h
        src/tests/test_lattice.cpp
        src/tests/test_lattice.h
        src/tests/test_percolation.cpp
        src/tests/test_percolation.h
        src/util/printer.h
        src/util/time_tracking.cpp
        src/util/time_tracking.h
        src/flags.h
        src/helper.cpp
        src/helper.h
        src/latticeConfiguration.cpp
        src/latticeConfiguration.h
        src/logs.cpp
        src/logs.h
        src/types.h
        src/ziff_algorithm.cpp
        src/ziff_algorithm.h
        src/util/printer.cpp
        src/percolation/percolation_site_v9.cpp
        src/util/checking.cpp
        src/util/checking.h
        src/percolation/percolation_bond_v1.cpp
        main.cpp
        )

add_executable(SqLatticePercolation ${SOURCE_FILES})