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
planetname = "Earth"

nfiles, nstars, starname, starradius, startemp, starcolor, fluxmax, fluxmin = infofile.read_infofile(prefix)

# print(row.split(',')[0] row in infofile.read_infofile)
# first_col = (int(row.split(',')[0]) for row in infofile.read_infofile(prefix))
# fluxmin = fluxmax = next(first_col, None)
# for i in first_col:
#     if i < fluxmin:
#         fluxmin = i
#     elif i > fluxmax:
#         fluxmax = i

moviechoice = input("Make an animated gif at end? (y/n) ")
deletechoice = 'n'
if(moviechoice=='y'):
    deletechoice = input("Delete .png files? (y/n) ")

nzeros = int(np.log10(nfiles))+1

# Loop over files
print(range(nfiles))

for i in range(nfiles):

    # Create filename - how many zeros needed?
    num = str(i+1)
    k=np.log10(i+1)

       
    ##For Avg Flux Plots 
    #Written 1/25/24 by LSP7654
    #inputfile = prefix +'_'+planetname+'_'+num+'.avg'
    inputfile = prefix +'_'+planetname+'.avg'
    avgfluxfile = 'avgflux_'+prefix+'_'+planetname+'_'+num+'.png'

    # Read in header - time, position data etc

    f = open(inputfile, 'r')

    line = f.readline()

    numbers=line.split() 

    # time = 1
    # print("warning: time is hardcoded T=1")
    nlat = int(numbers[0])
    nlong = int(numbers[1])
        
    f.close()
    
    print('File', str(i+1),'Time', time, "yr")
    
    # Read in rest of file
    
    data = np.genfromtxt(inputfile, skip_header=1)
        
    # Reshape to fit 2D array
    
    latitude = data[:,latcol].reshape(nlat,nlong)*180.0/pi
    longitude= data[:,longcol].reshape(nlat,nlong)*180.0/pi
    avgflux = data[:,fluxcol].reshape(nlat,nlong)
    darkness = data[:,darkcol].reshape(nlat,nlong)

    # Plot 2D maps of this timestep    
    
    # Flux
    
    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    #vmax was set to fluxmax but I changed it
    import matplotlib.colors as colors


    print("max =", fluxmax)
    print("min =", fluxmin)
    pcm = ax.pcolor(longitude,latitude,avgflux,norm=colors.Normalize(vmin= fluxmin, vmax= fluxmax), cmap='Spectral')
    fig1.colorbar(pcm)


    plt.savefig(avgfluxfile, format= 'png')
    plt.clf()

# end of loop


# Command for converting images into gifs - machine dependent

convertcommand = '/opt/homebrew/Cellar/imagemagick/7.1.1-22/bin/convert '

# Create movie if requested
if(moviechoice=='y'):
    print('Creating animated gif of avg flux pattern, filename avgfluxmovie.gif')
    system(convertcommand +'-delay 10 avgflux_'+prefix+'*.png avgfluxmovie.gif')

    if(deletechoice=='y'):
        print('Deleting png files')
        system('rm '+prefix+'_'+planetname+'*.png')                        
        

print('Complete')
