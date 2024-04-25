import numpy as numpy
import matplotlib.pyplot as plt

# https://lucydot.github.io/python_novice/

a = 1.0e-5
b = 1.0e-3
g = 2.0
P = 0

I = numpy.linspace(0,5000,100)
P = I/(a*I*I + b*I + g) - 0.3*P

plt.plot(I,P)
plt.xlabel("Intensity")
plt.ylabel("Photosynthesis Rate")
plt.show() 