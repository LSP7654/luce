#####################################################
###                                               ###
### Makefile for luce                             ###
###                                               ###
### Duncan H. Forgan 01/02/2019          	  ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables:
CC     = g++

VERSION = $(shell git describe --tags --abbrev=0)

# Git version
GIT_HASH = $(shell git describe --abbrev=9 --dirty --always)

# Define compiler flags
#CFLAGS = -O3 -DVERSION=\"$(VERSION)\" -DGIT_HASH=\"${GIT_HASH}\" -std=c++11
CFLAGS = -DVERSION=\"$(VERSION)\" -DGIT_HASH=\"${GIT_HASH}\" -std=c++11 -ggdb -g3


# Create object files:
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<
%.o: %.f90
	$(CC) $(CFLAGS) -c $<

# Source files (.F90)
SOURCESA = main.cpp Vector3D.cpp Body.cpp Star.cpp Planet.cpp PlanetSurface.cpp System.cpp parFile.cpp
OBJECTSA    = $(SOURCESA:.cpp=.o)

# Create executable files:
build: luce

luce:  $(OBJECTSA)
	$(CC) $(CFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o luce

# run:
# 	\mkdir -p ../data;
# 	\cd ../data;
# 	\../src/luce ../src/kepler.test
#fix later

# End Makefile
