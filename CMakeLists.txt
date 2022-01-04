cmake_minimum_required(VERSION 3.11)
project(ACPC_Preparation)


#Put this in your CMakeLists.txt file (change any options from OFF to ON if you want):



include_directories("C:/local/boost_1_78_0")
find_package(Boost)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_STACK_SIZE /STACK:4000000000)


add_executable(ACPC_Preparation linear_algebra.h modular_arithmetic.h number_theory.h
        ring_extension.h polynomial.h fft.h optimisation.h string.h data_structures.h binary_operation.h
        topology.h statistic_tree.h main.cpp )

        
target_include_directories(ACPC_Preparation PUBLIC "C:/local/boost_1_78_0/boost)")
