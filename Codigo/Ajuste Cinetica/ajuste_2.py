# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 19:49:37 2020

@author: rozog
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import pandas as pd
from scipy.optimize import minimize
from scipy import optimize
from scipy.optimize import fmin
datos=pd.read_excel(r'./Datos/ajuste.xlsx')



# Metodo Cinetica
def f(t,u):
    global k1
    k=k1
    du_dt=np.zeros((len(u)))
    # Siendo u[0]=POOH, u[1]=POO, u[2]=P u[3]=PH;
    #Reacciones k0: 2POOH = POO +P; k1: P+O2= POO ; 
    PO2=1e5;
    so2=4e-8;
    O2=so2*PO2;
    du_dt[0]=-2*k[0]*u[0]**2+k[2]*u[1]*u[3];
    du_dt[1]=k[0]*u[0]**2+k[1]*u[2]*O2-k[2]*u[1]*u[3]-2*k[3]*u[1]**2;
    du_dt[2]=k[0]*u[0]**2-k[1]*u[2]*O2+k[2]*u[1]*u[3];
    du_dt[3]=-k[2]*u[1]*u[3];
    return du_dt;


def modelo_minimizar(k):
   global data_t, data_u, T 
   u=solve_ode(T,k,data_t);
   rho=930; #g/L
   PO2=1e5;
   so2=4e-8;
   O2=so2*PO2;
   dm_dt=(32/rho)*(k[1]*(u[2,:].T)*O2);
   m=[];                 
   for i in range(len(data_t)):
       m.append(np.trapz(dm_dt[0:i], data_t[0:i]));
       
   return  np.sum((np.asarray(m)-np.asarray(data_u))**2); 
   
def solve_ode(T,uo,t,val=0):
    global k1
    tf=t[-1]
    sol=solve_ivp(f,[0,tf],uo,method='BDF',dense_output=True)
    u=sol.sol(t);
    if val==1: return sol;
    return u


datos_2=pd.read_excel(r'./80.xlsx',sheet_name=0)
data_t=datos_2['Tiempo (s)'].values.tolist();
data_u=datos_2['Masa (m/mo)'].values.tolist();


fig2 = plt.figure()
ax2=plt.subplot(1, 1, 1)



Ti=[80, 90 ,100];
for j in range(len(Ti)):
    T=Ti[j]+273; 
    rho=930; #g/L
    PO2=101325*0.21;
    so2=4e-8;
    O2=so2*PO2;
    global k1
    k1=np.asarray([4.09e10*np.exp(-101.5e3/(8.314*T)),1e7,3.52e6*np.exp(-31e3/(8.314*T)),3.05e15*np.exp(-48.4e3/(8.314*T))])
    t=data_t
    t=np.linspace(0,24*3600,100)
    uo=[0.5e-2,0,0,6.8]
    sol=solve_ode(T,uo,t,val=1);
    u=sol.sol(t)
    dm_dt=(32.0/rho)*(k1[1]*(u[2,:].T)*O2);
    m=[];                 
    for i in range(len(t)):
        m.append(np.trapz(dm_dt[0:i], t[0:i]));
    plt.plot(t,m, label=np.str(Ti[j]))
    plt.xlabel('Tiempo [s]')
    plt.ylabel('$O_2$ absorbido (g/$g_{oil}$)')
    plt.legend()
#    
#    POOH=u[0,:].T
#    POOH*=1.0e6/(3*rho+0.0)
#    ax2.plot(t,POOH,label='Model '+np.str(Ti[j]))
#    ax2.ticklabel_format(axis="y", style='sci', scilimits=(0,0))  
#    ax2.set_xlabel('Tiempo (h)')
#    ax2.set_ylabel('Peroxide Value [mmol/kg]')
#    datos_3=pd.read_excel(r'./Peroxide_value.xlsx',sheet_name=j)
#    data_t_PV=datos_3['t[s]'].values.tolist();
#    data_u_PV=datos_3['PV'].values.tolist();
#    ax2.plot(data_t_PV, data_u_PV, label='Data ' +np.str(Ti[j]) )
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    
    
#plt.plot(data_t,data_u,'-',label='data')
#plt.legend()
#plt.show()



#fig3 = plt.figure()
#ax3=plt.subplot(1, 1, 1)
#
#fig4 = plt.figure()
#ax4=plt.subplot(1, 1, 1)
#
#fig5= plt.figure()
#ax5=plt.subplot(1, 1, 1)
#
#ax3.plot(t,(u[1,:].T),label='POO')
#ax3.ticklabel_format(axis="y", style='sci', scilimits=(0,0))  
#ax4.plot(t,u[2,:].T,label='P')
#ax4.ticklabel_format(axis="y", style='sci', scilimits=(0,0))     
#ax5.plot(t,u[3,:].T,label='PH')
#ax5.ticklabel_format(axis="y", style='sci', scilimits=(0,0)) 

   




