import numpy as np
import matplotlib.pyplot as plt
G = 6.674e-8
MSUN = 1.989e33
C = 2.998e10
KPC = 3.086e21

MBH=6.5e9
d=16.8e3

th=90
dth=20
solid=2*np.pi*(np.cos((th-dth/2)*np.pi/180)-np.cos((th+dth/2)*np.pi/180))
print(f'solid angle={solid} degree')

dx=10
dy=60
df=(G*MBH*MSUN/C**2)/(d*KPC)
fov=dx*df*dy*df
print(f'df={df} FOV={fov}')

d0=np.loadtxt('bhac_data/spec_1_mpi.dat')
d1=np.loadtxt('bhac_data/spec_1_mp.dat')
d2=np.loadtxt('bhac_data/spec_1_84.dat')
d3=np.loadtxt('spec_1.dat')

plt.figure(dpi=150)
plt.step(d0[:,0]/1e9,d0[:,10],lw=2,alpha=0.5,label=r'MPI')
plt.step(d0[:,0]/1e9,d1[:,10],lw=2,alpha=0.5,label=r'OPENMP')
plt.step(d0[:,0]/1e9,d2[:,10],lw=2,alpha=0.5,label=r'8c*4t')
plt.step(d0[:,0]/1e9,d3[:,10],lw=2,alpha=0.5,label=r'test')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlim(1e-1,1e10)
plt.ylim(1e-8,10)
plt.xlabel('GHz')
plt.ylabel('Flux[Jy]')
plt.grid()
plt.minorticks_on()
plt.show()
