# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 11:42:53 2021

@author: frottenberg
"""
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')


def array_pattern(theta,M):
  a=np.exp(-1j*np.pi*np.cos(theta)*np.arange(M) )
  return a

M=32 # Number of antennas
theta=3*np.pi/4
theta=80*np.pi/180


h=array_pattern(theta,M)

w_MRT=h.conj()#./np.sqrt(M)
w_MRT=w_MRT/np.sqrt(np.sum(abs(w_MRT)**2)/M)
M_s=1
w_pure_Z3RO=h.conj()
w_pure_Z3RO[0:M_s]=-(np.sum(abs(h[M_s:])**4)/np.sum(abs(h[0:M_s])**4))**(1/3)*w_pure_Z3RO[0:M_s]
w_pure_Z3RO=w_pure_Z3RO/np.sqrt(np.sum(abs(w_pure_Z3RO)**2)/M)

#w=w_MRT.copy() # Uncomment for MRT precoder
w=w_pure_Z3RO.copy() # Uncomment for Z3RO precoder


M_pt=2000
theta_range=np.linspace(0,np.pi,M_pt)
phi_range=np.pi*np.cos(theta_range)

P_sig_w=np.zeros(M_pt)
P_dist_w=np.zeros(M_pt)


for index_theta in range(0, M_pt):    
    temp=array_pattern(theta_range[index_theta],M)
    P_sig_w[index_theta]=P_sig_w[index_theta]+abs(np.sum(w*temp))**2
    P_dist_w[index_theta]=P_dist_w[index_theta]+abs(np.sum(w*abs(w)**(3-1)*temp))**2  
    

P_tot=2*np.sum(P_sig_w)*(theta_range[1]-theta_range[0])
P_sig_w=P_sig_w/(P_tot/2/np.pi)

P_tot=2*np.sum(P_dist_w)*(theta_range[1]-theta_range[0])
P_dist_w=P_dist_w/(P_tot/2/np.pi)


fig = plt.figure()

ax = plt.subplot(projection='polar')

ax.plot(theta_range,10*np.log10(P_sig_w) ,label="Signal")
ax.plot(theta_range,10*np.log10(P_dist_w),'r--', label="Third-order distortion")
ax.set_rmax(np.max(10*np.log10(P_sig_w)))
ax.set_rmin(np.max(10*np.log10(P_sig_w))-30)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rticks(np.round(np.max(10*np.log10(P_sig_w)))+[-30, -25, -20, -15, -10, -5, 0])  # Less radial ticks
ax.grid(True)
plt.xlabel('[dB]')
ax.set_title(f"Directivity pattern with Z3RO precoding (M={M}, Ms=1)")
plt.legend()
plt.tight_layout()
plt.show()

