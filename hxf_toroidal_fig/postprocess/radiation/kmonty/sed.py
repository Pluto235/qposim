import numpy as np
import matplotlib.pyplot as plt

G = 6.674e-8
MSUN = 1.989e33
C = 2.998e10
KPC = 3.086e21

MBH=6.5e9
d=1.68e3*KPC
#print("{:e}".format(G*MBH*MSUN/C**2))
#print(d*KPC)

th=90
dth=20
solid=2*np.pi*(np.cos((th-dth/2)*np.pi/180)-np.cos((th+dth/2)*np.pi/180))
print(f'solid angle={solid}')

dx=10
dy=60
df=(G*MBH*MSUN/C**2)/(d*KPC)
fov=dx*df*dy*df
print(f'df={df} FOV={fov}')
print(fov*1e23)
#d0=np.loadtxt('spec_96_1e2.dat')
d1=np.loadtxt('spec_96.dat')
d3=np.loadtxt('th/spec_96s.dat')
d2=np.loadtxt('th/spectrum_40_90.00.dat')

plt.figure(dpi=150)
'''
#plt.step(d0[:,0]/1e9,d0[:,10],lw=2,alpha=0.5,label=r'1e2')
plt.step(d1[:,0]/1e9,d1[:,10],lw=2,alpha=0.5,label=r'IC')
plt.step(d2[:,0]/1e9,d2[:,1],lw=2,alpha=0.5,label=r'Raptor')
plt.step(d3[:,0]/1e9,d3[:,10],lw=2,alpha=0.5,label=r'Syn')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlim(1,1e15)
plt.ylim(1e-8,1e2)
plt.title(r'Thermal $90^\circ$')
plt.xlabel('GHz')
plt.ylabel('Flux [Jy]')
plt.grid()
plt.minorticks_on()
plt.show()


'''
plt.step(d2[:,0]/1e9,d2[:,5]/(4*np.pi*d**2),lw=2,alpha=0.5,label=r'Raptor')
plt.step(d1[:,0]/1e9,d1[:,9]/(4*np.pi*d**2),lw=2,alpha=0.5,label=r'IC')
plt.step(d3[:,0]/1e9,d3[:,9]/(4*np.pi*d**2),lw=2,alpha=0.5,label=r'Syn')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlim(1,1e15)
plt.ylim(1e-14,1e-7)
plt.title('Thermal $90^\circ$')
plt.xlabel('GHz')
plt.ylabel(r'$\nu F_{\nu}$ [$erg\cdot cm^{-2}\cdot s^{-1}$]')
plt.grid()
plt.minorticks_on()
plt.show()

