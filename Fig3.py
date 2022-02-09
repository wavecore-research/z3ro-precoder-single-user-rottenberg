# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 11:42:53 2021

@author: frottenberg
"""
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')





def array_pattern(theta,M):
  a=np.exp(-1j*np.pi*np.cos(theta)*np.arange(M) )#*np.sin(theta)
  return a


def compute_radiation_pattern(w):
    P_sig_w=np.zeros(N_pt)
    P_dist_w=np.zeros(N_pt)
        
    for index_theta in range(0, N_pt):    
        temp=array_pattern(theta_range[index_theta],M)
        P_sig_w[index_theta]=P_sig_w[index_theta]+abs(np.sum(w*temp))**2
        P_dist_w[index_theta]=P_dist_w[index_theta]+abs(np.sum(w*abs(w)**(3-1)*temp))**2  
    return P_sig_w, P_dist_w



M=32
theta=3*np.pi/4
theta=80*np.pi/180



w_MRT=array_pattern(theta,M).conj()#./np.sqrt(M)

M_s=1
M_2=M-M_s
w_pure_3rd=w_MRT.copy()
w_pure_3rd[0:M_s]=-(M_2/M_s)**(1/3)*w_pure_3rd[0:M_s]
w_pure_3rd=w_pure_3rd/np.sqrt(np.sum(abs(w_pure_3rd)**2)/M)

M_s=2
M_2=M-M_s
w_pure_3rd_2=w_MRT.copy()
w_pure_3rd_2[0:M_s]=-(M_2/M_s)**(1/3)*w_pure_3rd_2[0:M_s]
w_pure_3rd_2=w_pure_3rd_2/np.sqrt(np.sum(abs(w_pure_3rd_2)**2)/M)

M_s=4
M_2=M-M_s
w_pure_3rd_4=w_MRT.copy()
w_pure_3rd_4[0:M_s]=-(M_2/M_s)**(1/3)*w_pure_3rd_4[0:M_s]
#w_pure_3rd_4=np.random.normal(0,1/2,M)+1j*np.random.normal(0,1/2,M)
w_pure_3rd_4=w_pure_3rd_4/np.sqrt(np.sum(abs(w_pure_3rd_4)**2)/M)



N_pt=2000
theta_range=np.linspace(0,np.pi,N_pt)
phi_range=np.pi*np.cos(theta_range)


P_sig_w_MRT, P_dist_w_MRT=compute_radiation_pattern(w_MRT)
P_sig_w_3rd, P_dist_w_3rd=compute_radiation_pattern(w_pure_3rd)
P_sig_w_3rd_2, P_dist_w_3rd_2=compute_radiation_pattern(w_pure_3rd_2)
P_sig_w_3rd_4, P_dist_w_3rd_4=compute_radiation_pattern(w_pure_3rd_4)


fig = plt.figure()

ax = plt.subplot(projection='polar')

ax.plot(theta_range,10*np.log10(P_dist_w_MRT),'k',label='MRT')
ax.plot(theta_range,10*np.log10(P_dist_w_3rd),'b--',label='$M_{s}=1$')
ax.plot(theta_range,10*np.log10(P_dist_w_3rd_2),'r--',label='$M_{s}=2$')
ax.plot(theta_range,10*np.log10(P_dist_w_3rd_4),'g--',label='$M_{s}=4$')
ax.set_rmax(np.max(10*np.log10(P_dist_w_MRT)))
ax.set_rmin(np.max(10*np.log10(P_dist_w_MRT))-30)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rticks(np.round(np.max(10*np.log10(P_dist_w_MRT)))+[-30, -25, -20, -15, -10, -5, 0])  # Less radial ticks
ax.grid(True)
#plt.xlabel('Power [dB]')
plt.legend(bbox_to_anchor=(-0.1, 0.05),loc="lower left", ncol = len(ax.lines) )


fig = plt.figure()

ax = plt.subplot(projection='polar')

ax.plot(theta_range,10*np.log10(P_sig_w_MRT),'k',label='MRT')
ax.plot(theta_range,10*np.log10(P_sig_w_3rd),'b--',label='$M_{s}=1$')
ax.plot(theta_range,10*np.log10(P_sig_w_3rd_2),'r--',label='$M_{s}=2$')
ax.plot(theta_range,10*np.log10(P_sig_w_3rd_4),'g--',label='$M_{s}=4$')
ax.set_rmax(np.max(10*np.log10(P_sig_w_MRT)))
ax.set_rmin(np.max(10*np.log10(P_sig_w_MRT))-30)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rticks(np.round(np.max(10*np.log10(P_sig_w_MRT)))+[-30, -25, -20, -15, -10, -5, 0])  # Less radial ticks
ax.grid(True)
plt.legend(bbox_to_anchor=(-0.1, 0.05),loc="lower left", ncol = len(ax.lines) )
plt.xlabel('Power [dB]')
plt.tight_layout()
plt.show()
