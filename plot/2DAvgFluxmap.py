# Written by D Forgan, 17/12/2014
# Reads in output dumps from nbody_2DFlux (C++ code)
# Produces 2D maps of each timestep (PNG)

import numpy as np
import matplotlib.pyplot as plt
import infofile
#from string import split
from os import system
#from sys import exit

pi = 3.1415926585

longcol = 0
latcol = 1
fluxcol = 2
darkcol = 3


# Read in input parameters

#prefix = input("What is the file prefix? ")
#planetname = input("What is the planet name? ")

prefix = "test"
planetname = "Kepler-16b"



nfiles, nstars, starname, starradius, startemp, starcolor, fluxmax, fluxmin, Avgfluxmax, Pmax, Pmin = infofile.read_infofile(prefix)

inputfile = prefix +'_'+planetname+'.avg'
avgfluxfile = 'avgflux_'+prefix+'_'+planetname+'.png'

f = open(inputfile, 'r')

line = f.readline()

numbers=line.split() 

nlat = int(numbers[0])
nlong = int(numbers[1])
        
f.close()

data = np.genfromtxt(inputfile, skip_header=1)

latitude = data[:,latcol].reshape(nlat,nlong)*180.0/pi
longitude= data[:,longcol].reshape(nlat,nlong)*180.0/pi
avgflux = data[:,fluxcol].reshape(nlat,nlong)
darkness = data[:,darkcol].reshape(nlat,nlong)

fig1 = plt.figure(1)
ax = fig1.add_subplot(111)
ax.set_xlabel('Longitude (degrees)')
ax.set_ylabel('Latitude (degrees)')

import matplotlib.colors as colors

#print("max =", fluxmax)
#print("min =", fluxmin)
pcm = ax.pcolor(longitude,latitude,avgflux,norm=colors.Normalize(vmin= fluxmin, vmax= Avgfluxmax), cmap='nipy_spectral')
fig1.colorbar(pcm)

plt.savefig(avgfluxfile, format= 'png')
plt.clf()
