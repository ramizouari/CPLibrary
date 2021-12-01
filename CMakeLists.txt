cmake_minimum_required(VERSION 3.0)
project(ACPC_Preparation)

set(CMAKE_CXX_STANDARD 17)

add_executable(ACPC_Preparation number_theory.h polynomial.h abstract_algebra.h ring_extension.cpp fft.h)
