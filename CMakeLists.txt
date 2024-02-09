cmake_minimum_required(VERSION 3.0)
project(SqLatticePercolation)

set(CMAKE_CXX_STANDARD 17)

#set (CMAKE_C_COMPILER               /usr/bin/gcc)
#set (CMAKE_CXX_COMPILER             /home/shahnoor/software/pgi/linux86-64-llvm/2018/bin/pgc++)
#set (CMAKE_MAKE_PROGRAM             /usr/bin/make)
#SET( CMAKE_CXX_FLAGS  "-Wall -D DEBUG_FLAG")
#SET( CMAKE_CXX_FLAGS  "-Wall -D UNIT_TEST")
#SET( CMAKE_CXX_FLAGS  "-Wall -D DEBUG_FLAG -D UNIT_TEST")
#SET( CMAKE_CXX_FLAGS  "-Wall")
SET( CMAKE_CXX_FLAGS  "-Wall -O3 -DCMAKE_BUILD_TYPE=RELEASE")
#SET( CMAKE_CXX_FLAGS  "-Wall -O3 -pthread -fopenmp")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

set(SOURCE_FILES
        src/flags.h
        src/types.h
        src/main.cpp
        src/cmdLineArgs/cmd_args.cpp
        src/cmdLineArgs/cmd_args.h
        src/cmdLineArgs/cmd_args_json.cpp
        src/cmdLineArgs/cmd_args_json.h
        src/cmdLineArgs/sq_lattice_methods.h
        src/ensemble/site_percolation.cpp
        src/ensemble/site_percolation.h
        src/exception/exceptions.cpp
        src/exception/exceptions.h
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
        src/percolation/cluster/cluster.cpp
        src/percolation/cluster/cluster.h
        src/percolation/id_index.cpp
        src/percolation/id_index.h
        src/percolation/percolation.cpp
        src/percolation/percolation.h
        src/percolation/site/percolation_site_ballistic_deps_v2.cpp
        src/percolation/site/percolation_site_ballistic_deps_v3.cpp
        src/percolation/site/percolation_site_ballistic_deps_v4.cpp
        src/percolation/site/percolation_site_explosive.cpp
        src/percolation/percolation_type.h
        src/percolation/site_position.cpp
        src/percolation/site_position.h
        src/util/printer.h
        src/util/time_tracking.cpp
        src/util/time_tracking.h
        src/util/printer.cpp
        src/percolation/site/percolation_site_v9.h
        src/percolation/site/percolation_site_v9.cpp
        src/percolation/site/percolation_site_v10.h
        src/percolation/site/percolation_site_v10.cpp
        src/percolation/site/percolation_site_v11.h
        src/percolation/site/percolation_site_v11.cpp
        src/util/checking.cpp
        src/util/checking.h
        src/percolation/bond/percolation_bond_v1.h
        src/percolation/bond/percolation_bond_v2.h
        src/percolation/bond/percolation_bond_v3.h
        src/percolation/bond/percolation_bond_v1.cpp
        src/percolation/bond/percolation_bond_v2.cpp
        src/percolation/bond/percolation_bond_v3.cpp
        src/tests/test_percolation.cpp
        src/tests/test_percolation.h
        src/percolation/percolation_v2.cpp
        src/percolation/percolation_v2.h
        src/lattice/lattice_v2.cpp
        src/tests/test_v11.cpp
        src/tests/test_v11.h
        src/percolation/site/percolation_site_explosive.h
        src/percolation/cluster/cluster_v4.cpp
        src/percolation/site/percolation_site_ballistic.h
        src/tests/test_v9.cpp
        src/tests/test_v9.h
        src/tests/test_v10.cpp
        src/tests/test_v10.h
        src/percolation/percolation_v12.cpp
        src/percolation/percolation_v12.h
        src/v13/bond_v13.cpp
        src/v13/bond_v13.h
        src/v13/cluster.cpp
        src/v13/cluster.h
        src/v13/element_v13.cpp
        src/v13/element_v13.h
        src/v13/index_v13.cpp
        src/v13/index_v13.h
        src/v13/percolation_v13.cpp
        src/v13/percolation_v13.h
        src/v13/site_v13.cpp
        src/v13/site_v13.h
        src/tests/test_v13.cpp
        src/tests/test_v13.h
        src/v13/lattice_v13.cpp
        src/v13/lattice_v13.h
        src/v13/run_v13.cpp
        src/v13/run_v13.h
        src/v13/percolation_rsbd_v13.h
        src/v13/status.h
        src/v13/percolation_rsbd_v13.cpp
        src/v13/percolation_v13_bond.h
        src/v13/percolation_v13_bond.cpp
        )

add_executable(sqlattice ${SOURCE_FILES})
