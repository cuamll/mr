#!/bin/bash

# i am too fucking lazy to sort out the makefile for the moment, so

# comment out for linux
# gfortran-mp-5 -g -c common.f90
# gfortran-mp-5 -g -c io.f90
# gfortran-mp-5 -g -c linear_solver.f90
# gfortran-mp-5 common.o io.o linear_solver.o main.f90 -o mr_test

# comment out for mac
gfortran -g -c common.f90
gfortran -g -c io.f90
gfortran -g -c linear_solver.f90
gfortran common.o io.o linear_solver.o main.f90 -o mr_test
