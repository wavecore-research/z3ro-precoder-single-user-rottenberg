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

def rapp_model(r_in,r_sat,S):
    r_out=r_in/((1+(r_in/r_sat)**(2*S))**(1/2/S))
    return r_out

def compute_signal_distortion_power_simu(w,p,r_sat,h):
    N_data=10000
    s=np.sqrt(p)/np.sqrt(2)*(np.random.randn(N_data)+1j*np.random.randn(N_data))
#    s=np.sqrt(p)/np.sqrt(2)*(np.random.randint(2,size=N_data)*2-1+1j*(np.random.randint(2,size=N_data)*2-1))
    x=np.zeros((M,N_data),dtype=complex)
    y=np.zeros((M,N_data),dtype=complex)    
    for index_M in range(0,M):
        x[index_M,:]=w[index_M] * s
        S=2
        y[index_M,:]=rapp_model(abs(x[index_M,:]),r_sat,S)*np.exp(1j*np.angle(x[index_M,:]))
    r=h @ y   
    C_rs=np.average(s*np.conj(r))
    C_ss=np.average(s*np.conj(s))
    C_rr=np.average(r*np.conj(r))    
    B=C_rs/C_ss
    sigma2_s=abs(B)**2*C_ss
    sigma2_d=C_rr-abs(B)**2*C_ss  
    return sigma2_s, sigma2_d

M_p=20
back_off_set=10**(np.linspace(-15,0,M_p)/10)
p=1
p_sat_set=p/back_off_set





M=64

theta=3*np.pi/4
theta=60*np.pi/180

#h=array_pattern(theta,M)
h=np.random.normal(0,1/2,M)+1j*np.random.normal(0,1/2,M)
h=h/np.sqrt(np.sum(abs(h)**2)/M)

w_MRT=h.conj()#./np.sqrt(M)
w_MRT=w_MRT/np.sqrt(np.sum(abs(w_MRT)**2)/M)

M_s=2
w_Z3RO=h.conj()
w_Z3RO[0:M_s]=-(np.sum(abs(h[M_s:])**4)/np.sum(abs(h[0:M_s])**4))**(1/3)*w_Z3RO[0:M_s]
w_Z3RO=w_Z3RO/np.sqrt(np.sum(abs(w_Z3RO)**2)/M)


#sigma2_nu=1e0
sigma2_nu=1e1

sigma2_s_MRT=np.zeros(M_p)
sigma2_s_Z3RO=np.zeros(M_p)
sigma2_d_MRT=np.zeros(M_p)
sigma2_d_Z3RO=np.zeros(M_p)
w_pure_SNDR_m=np.zeros((M,M_p))
for index_p in range(0,M_p):
    p_sat=p_sat_set[index_p]    
    sigma2_s_MRT[index_p], sigma2_d_MRT[index_p] = np.real(compute_signal_distortion_power_simu(w_MRT,p,np.sqrt(p_sat),h))    
    sigma2_s_Z3RO[index_p], sigma2_d_Z3RO[index_p] = np.real(compute_signal_distortion_power_simu(w_Z3RO,p,np.sqrt(p_sat),h))    
    


SDR_MRT=sigma2_s_MRT/sigma2_d_MRT
SDR_Z3RO=sigma2_s_Z3RO/sigma2_d_Z3RO

SNR_MRT=sigma2_s_MRT/sigma2_nu
SNR_Z3RO=sigma2_s_Z3RO/sigma2_nu

SNDR_MRT=1/(1/SDR_MRT+1/SNR_MRT)
SNDR_Z3RO=1/(1/SDR_Z3RO+1/SNR_Z3RO)



MFB_MRT=(sigma2_s_MRT+sigma2_d_MRT)/sigma2_nu
MFB_Z3RO=(sigma2_s_Z3RO+sigma2_d_Z3RO)/sigma2_nu




fig = plt.figure()
plt.plot(10*np.log10(back_off_set),10*np.log10(SDR_MRT),':r', linewidth=2,label='SDR - MRT')
plt.plot(10*np.log10(back_off_set),10*np.log10(SNR_MRT),'--r', linewidth=2,label='SNR - MRT')
plt.plot(10*np.log10(back_off_set),10*np.log10(SNDR_MRT),'-r', linewidth=2,label='SNDR - MRT')
plt.plot(10*np.log10(back_off_set),10*np.log10(SDR_Z3RO),':b', linewidth=2,label='SDR - Z3RO')
plt.plot(10*np.log10(back_off_set),10*np.log10(SNR_Z3RO),'--b', linewidth=2,label='SNR - Z3RO')
plt.plot(10*np.log10(back_off_set),10*np.log10(SNDR_Z3RO),'-b', linewidth=2,label='SNDR - Z3RO')


#plt.plot(10*np.log10(back_off_set),10*np.log10(MFB_MRT),'-xr', linewidth=2,label='MFB - MRT')
#plt.plot(10*np.log10(back_off_set),10*np.log10(MFB_Z3RO),'-xb', linewidth=2,label='MFB - Z3RO')

#plt.plot(np.array((-6,-6)),np.array((-100,100)),'-.k', linewidth=0.5)
plt.text(-4,27.5,'Saturation regime')
plt.text(-9.5,27.5,'Linear regime')

plt.xlim((-15,0))
plt.ylim((10,30))
plt.legend()
plt.xlabel('$p/p_{sat}$ [dB]')
plt.ylabel('Magnitude [dB]')
plt.show()