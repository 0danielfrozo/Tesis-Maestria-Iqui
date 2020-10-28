import numpy as np
import scipy.optimize as spo
import pandas as pd
import matplotlib
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix, csr_matrix, tril,triu ,diags
from scipy.sparse.linalg import cg
from scipy.sparse.linalg import spsolve 
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import MaxNLocator



#-------------------------------------------------------------------------------
    
def Propiedades_Material(L):
    global materiales,R,T,x,D,S
    n=len(x);
    
    Ea=36500 #J/mol 
    Do=0.58 #cm2/s PP
    #D=Do*np.exp(-Ea/(R*T))#cm2/s
    D=1e-11*1e4    
    Delta_H=2775.7 #J/mol
    So=3.14E-11 #mol/cm3*Pa
    S=So*np.exp((-Delta_H/R)*(1.0/T -1.0/298.0))
    
    resp_D=D*np.ones((len(x)));   
    resp_S=S*np.ones((len(x)));  
    #resp_S=np.ones((len(x)));
    
    D_nuevo=np.zeros((5*len(x)));
    D_nuevo[4*n:]=resp_D;
    D=diags(D_nuevo).tolil();
    S=diags(resp_S*1e6*R*T).tolil();
    #S=diags(resp_S).tolil();


def Reacc(u):
    k=np.asarray([4.09e13*np.exp(-105.5e3/(R*T)),1e10,3.52e9*np.exp(-31e3/(R*T)),3e14,1e13,1e11])
    dx=x[1]-x[0]
    especies=np.reshape(u,(5,len(x)));
    du_dt=np.zeros(np.shape(especies));
    #-----------------Especies-------------------------
    #u[0],=POOH, u[1]=POO°, u[2]= P°, u[3]=PH, u[4]=O_2
    #----------------Reacciones------------------------
    #k0: 2POOH=POO° +P°,  k1: P°+O_2=POO°
    #k2: POO°+PH=POOH+P°, k3: 2P°= ...
    #k4: P° +POO°=...,    k5: POO° +POO° = ...
    #----------Calculo de la derivada temporal -------
    du_dt[0,:]=-2*k[0]*especies[0,:]**2+k[2]*especies[1,:]*especies[3,:];
    du_dt[1,:]=k[0]*especies[0,:]**2+k[1]*especies[2,:]*especies[4,:]-k[2]*especies[1,:]*especies[3,:]-k[4]*especies[2,:]*especies[1,:]-2*k[5]*especies[1,:]**2;
    du_dt[2,:]=k[0]*especies[0,:]**2-k[1]*especies[2,:]*especies[4,:]+k[2]*especies[1,:]*especies[3,:]-2*k[3]*especies[2,:]**2-k[4]*especies[2,:]*especies[1,:];
    du_dt[3,:]=-k[2]*especies[1,:]*especies[3,:]-k[0]*especies[0,:]**2;
    du_dt[4,:]=-k[1]*especies[2,:]*especies[4,:]+k[5]*especies[1,:]**2;
    du_dt=np.reshape(du_dt,5*len(x)); 
    return du_dt   

def f(t,u):
    global x
    Dif=D[4*n:,4*n:];
    dx=x[1];
    nabla_2=diags([1/(dx**2),-2/(dx**2), 1/(dx**2)], [-1, 0, 1],shape=(n,n)).tolil();
    diffusivo=(Dif*nabla_2);
    resp=Reacc(u)
    nx=len(x)
    do2_dx=diffusivo.dot(u[4*nx:]);
    resp[4*nx:]+= do2_dx;
    resp[4*nx]=0;
    resp[-1]=0;
    return resp

def ic(L,OS_0,co,ci):
    global S, T, R, x 
    #---Concentracion inicial de  Espcies en Aceite de Linaza---
    #Lin_oil_Co=np.asarray([0.05/1000.0,0,0,6.8/1000.0]); #mol/cm3 Lin_oil
    rho_oil=0.960 #g/cm3
    rho_ipp=0.936 #g/cm3
    Lin_oil_Co=np.asarray([0.5/1000.0,0,0,6e-2/(rho_oil)]);
    
    #--------------Vector de concentrracion inicial ------------
    uo=np.zeros((5,len(x)));
    for k  in  range(4): uo[k,:]=OS_0[0]*(rho_ipp/rho_oil)*Lin_oil_Co[k];   
    #-------------Concentracion Inicial de O2 -------------------
    uo[4,0]=co*S[0,0];
    uo[4,-1]=ci*S[-1,-1];
    uo=np.reshape(uo,5*len(x));
    return uo

L=0.015;
n=51;
R=8.314;#J/mol K
T=60+273.0;# K
A=25 #cm2
V=np.inf #cm3
x=np.linspace(0,L,n); 
dx=x[1]
Propiedades_Material(L)
co=(0.21*101325/(R*T))*1e-6; #mol/cm3
ci=(0.03*101325/(R*T))*1e-6; #mol/cm3
OS_0=[0.01];
uo=ic(L,OS_0,co,co)

tf=7200*12*3;
t=np.linspace(0,tf,tf/3600 +1);
sol = solve_ivp(f, [0,tf],uo,method='BDF',t_eval=t,rtol=1e-11,atol=1e-15);




#k0: 2POOH=POO° +P°,  k1: P°+O_2=POO°
#k2: POO°+PH=POOH+P°, k3: 2P°= ...
#k4: P° +POO°=...,    k5: POO° +POO° = ...

#-----------------Especies-------------------------
#u[0],=POOH, u[1]=POO°, u[2]= P°, u[3]=PH, u[4]=O_2
k=np.asarray([4.09e13*np.exp(-105.5e3/(R*T)),1e10,3.52e9*np.exp(-31e3/(R*T)),3e14,1e13,3.42e9])
dO2=k[1]*sol.y[2*n:3*n,:]*sol.y[4*n:,:]-k[5]*sol.y[n:2*n,:]**2

FFO2=0.5*(dO2[:,1:]+dO2[:,:-1])*(t[1:]-t[:-1]);

O2=np.zeros((n,len(t)))
for i in range(1,len(t)-1):
    O2[:,i]=O2[:,i-1]+FFO2[:,i]
O2[:,-1]=O2[:,-2]+FFO2[:,-1]    

f1=plt.figure();
a = f1.gca(projection='3d')
X, Y = np.meshgrid(x, sol.t/3600.0)
surf = a.plot_surface(X, Y, O2.T, cmap=cm.jet,linewidth=0,edgecolor='none', antialiased=False)
c1=f1.colorbar(surf,ax = a,aspect=20, shrink=0.5,format='%.2e')
c1.set_label(r'$O_2$ [mol/$cm^3$]')
a.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
a.view_init(azim=-75, elev=30)
a.set_xlabel('Posicion [cm]')
a.set_ylabel('tiempo [h]')
a.set_zlabel('$O_2$ [mol/$cm^3$]')


f1=plt.figure();
a = f1.gca(projection='3d')
X, Y = np.meshgrid(x, sol.t/3600.0)
surf = a.plot_surface(X, Y, sol.y[4*n:,:].T, cmap=cm.jet,linewidth=0,edgecolor='none', antialiased=False)
c1=f1.colorbar(surf,ax = a,aspect=20, shrink=0.5,format='%.2e')
c1.set_label(r'$O_2$ [mol/$cm^3$]')
a.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
a.view_init(azim=-75, elev=30)
a.set_xlabel('Posicion [cm]')
a.set_ylabel('tiempo [h]')
a.set_zlabel('$O_2$ [mol/$cm^3$]')


f2=plt.figure()
plt.plot(sol.t/(3600.0*24),sol.y[-1,:])
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}'))
plt.ylabel('$O_2$ [mol/$cm^3$]')
plt.xlabel('tiempo [d]')
