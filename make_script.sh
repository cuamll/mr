#!/bin/bash

# i am too fucking lazy to sort out the makefile for the moment, so

gfortran-mp-5 -g -c common.f90
gfortran-mp-5 -g -c io.f90
gfortran-mp-5 -g -c linear_solver.f90
gfortran-mp-5 common.o io.o linear_solver.o main.f90 -o linsol
