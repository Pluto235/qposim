import numpy as np
import matplotlib.pyplot as plt

d=np.loadtxt("weight.txt")

plt.plot(d[0:,0]/1e9,np.exp(d[:,1]))
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.ylim(1e16,1e56)
plt.show()
