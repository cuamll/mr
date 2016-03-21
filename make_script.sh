#!/bin/bash

# i am too fucking lazy to sort out the makefile for the moment, so

gfortran-mp-5 -c common.f90
gfortran-mp-5 -c linear_solver.f90
gfortran-mp-5 common.o linear_solver.o main.f90 -o linsol
