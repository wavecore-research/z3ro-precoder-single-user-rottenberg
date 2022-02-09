# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 11:42:53 2021

@author: frottenberg
"""
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')


M_set=np.arange(1,500)
M_M=np.size(M_set)



array_gaiM_MRT=np.zeros(M_M)
array_gaiM_pure_3rd=np.zeros(M_M)
array_gaiM_pure_3rd_2=np.zeros(M_M)
array_gaiM_pure_3rd_4=np.zeros(M_M)

for index_M in range(0,M_M):
    M=M_set[index_M]
    w_MRT=np.ones(M)#./np.sqrt(M)
    array_gaiM_MRT[index_M]=abs(np.sum(w_MRT))**2/M
    
    M_1=1
    M_2=M-M_1
    w_pure_3rd=w_MRT.copy()
    w_pure_3rd[0:M_1]=-(M_2/M_1)**(1/3)*w_pure_3rd[0]
    w_pure_3rd=w_pure_3rd/np.sqrt(np.sum(abs(w_pure_3rd)**2)/M)
    if M_1 > M/2:
        array_gaiM_pure_3rd[index_M]=0
    else:
        array_gaiM_pure_3rd[index_M]=abs(np.sum(w_pure_3rd))**2/M    
            
    M_1=2
    M_2=M-M_1
    w_pure_3rd_2=w_MRT.copy()
    w_pure_3rd_2[0:M_1]=-(M_2/M_1)**(1/3)*w_pure_3rd_2[0]
    w_pure_3rd_2=w_pure_3rd_2/np.sqrt(np.sum(abs(w_pure_3rd_2)**2)/M)
    if M_1 > M/2:
        array_gaiM_pure_3rd_2[index_M]=1e-6
    else:
        array_gaiM_pure_3rd_2[index_M]=abs(np.sum(w_pure_3rd_2))**2/M
    
    M_1=4
    M_2=M-M_1
    w_pure_3rd_4=w_MRT.copy()
    w_pure_3rd_4[0:M_1]=-(M_2/M_1)**(1/3)*w_pure_3rd_4[0]
    w_pure_3rd_4=w_pure_3rd_2/np.sqrt(np.sum(abs(w_pure_3rd_4)**2)/M)
    if M_1 > M/2:
        array_gaiM_pure_3rd_4[index_M]=1e-6
    else:
        array_gaiM_pure_3rd_4[index_M]=abs(np.sum(w_pure_3rd_4))**2/M
    
    

fig = plt.figure()
plt.plot(M_set,10*np.log10(np.ones(M_M)),'--k')
plt.plot(M_set,10*np.log10(array_gaiM_pure_3rd/array_gaiM_MRT), '-b', label='$M_{s}=1$')
plt.plot(M_set,10*np.log10(array_gaiM_pure_3rd_2/array_gaiM_MRT), '--r', label='$M_{s}=2$')
plt.plot(M_set,10*np.log10(array_gaiM_pure_3rd_4/array_gaiM_MRT), 'g:', label='$M_{s}=4$')
#plt.plot(M_set,array_gaiM_pure_3rd)
plt.legend()
plt.xlim(0, np.max(M_set))  
plt.ylim(-10, 1)  
plt.xlabel('Mumber of antennas $M$')
plt.ylabel('Array gain penalty versus MRT [dB]')
plt.tight_layout()
plt.show()

