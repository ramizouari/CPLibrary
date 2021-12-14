cmake_minimum_required(VERSION 3.0)
project(ACPC_Preparation)

set(CMAKE_CXX_STANDARD 20)

add_executable(ACPC_Preparation linear_algebra.h modular_arithmetic.h number_theory.h
        ring_extension.h polynomial.h fft.h optimisation.h string.h data_structures.h binary_operation.h
        topology.h statistic_tree.h main.cpp)
