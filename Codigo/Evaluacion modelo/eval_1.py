
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
#Metodo que calcula la difusividad y solublidad del O2 en la pelicula
#Este metodo no devuelve nada pero genera dos variables globales D y S
# D: matriz diagonal con la difusividad del oxigeno en cada posicion de la pelicula [cm2/s]
# S: matriz diagonal con la solubilidad del oxigeno en cada posicion de la pelicula 
# Para el calculo de la difusividad se utiliza: D(T)=Do exp(-Ea/RT)
# Para el calculo de la difusividad se utiliza: S(T)=So exp(-Delta_H/R [1/T-1/298K])

def Propiedades_Material():
    global D,S
    
#    #Calculo de la Difusividad de O2 en PP por Ley de Arrehnius     
#    Ea=36500 #J/mol 
#    Do=0.58 #cm2/s 
#    D=Do*np.exp(-Ea/(R*T))#cm2/s 
#    
#    #Calculo de la Solubilidad de O2 en PP por Ley de Arrehnius  
#    Delta_H=2775.7 #J/mol
#    So=3.14E-11 #mol/cm3*Pa
#    S=So*np.exp((-Delta_H/R)*(1.0/T -1.0/298.0)) #mol/cm3*Pa
    
    D=2.30E-07;
    S=5.17E-12;
    
    
    
    #Almacenamiento de la diusividad y la solubilidad en vectores.   
    resp_D=D*np.ones((len(x)));   
    resp_S=S*np.ones((len(x)));  
    
    #Conversion de la difusividad y la solubilidad en matrices.
    #Los unicos valores de la matriz diagonal diferentes de cero son las posiciones correspondientes
    #a la concentracion de oxigeno
    # la solubilidad esta dada en unidades [mol/cm3_PP* cm3/mol_aire]
    D_nuevo=np.zeros((5*len(x)));
    D_nuevo[4*n:]=resp_D;
    D=diags(D_nuevo).tolil();
    S=diags(resp_S*1e6*R*T).tolil();
    
# Metodo que devuelve la condicon iniciai de la pelicula 
# Este metodo tiene como  parametros:
#   OS_0: Concentracion masica de OS en pelicula [g_os/g_PP]
#   co: Concentracion de oxigeno en el aire en x=0 [mol/cm3]
#   ci: Concentracion de oxigeno en el aire en x=L [mol/cm3]
#   Constante:  inidica si la concentracion del oxigeno externo es constante a lo largo del tiempo. Booleano
#   c_ext: inidica el valor incial de la concentracion de oxigeno exterior que va a cambiar con el tiempo. [mol/cm3]  
#   Este metodo regresa un vector de concentraciones iniciales a lo largo de la pelicula. 

def ic(OS_0,co,ci,constante, c_ext):
    #---Concentracion inicial de  Espcies en Aceite de Linaza---
    Lin_oil_Co=np.asarray([0.005/1000.0,0,0,6.8/1000.0]); #mol/cm3 Lin_oil
    rho_oil=0.960 #g/cm3
    rho_ipp=0.936 #g/cm3    
    #--------------Vector de concentrracion inicial ------------
    uo=np.zeros((5,len(x)));
    for k  in  range(4): uo[k,:]=OS_0[0]*(rho_ipp/rho_oil)*Lin_oil_Co[k];   
    #-------------Concentracion Inicial de O2 -------------------
    uo[4,0]=co*S[0,0];
    uo[4,-1]=ci*S[-1,-1];
    uo=np.reshape(uo,5*len(x));
    if constante==True:
        return uo 
    else:
        u_final=np.zeros(5*len(x)+1);
        u_final[:-1]=uo;
        u_final[-1]=c_ext;
        return u_final;    
    
    
def Reacc(u):
    if constante==True:
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
        du_dt[3,:]=-k[2]*especies[1,:]*especies[3,:];
        du_dt[4,:]=-k[1]*especies[2,:]*especies[4,:];
        du_dt=np.reshape(du_dt,5*len(x)); 
        return du_dt   
    else:
        especies=np.reshape(u[:-1],(5,len(x)));
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
        du_dt[3,:]=-k[2]*especies[1,:]*especies[3,:];
        du_dt[4,:]=-k[1]*especies[2,:]*especies[4,:];
        du_dt=np.reshape(du_dt,5*len(x)); 
        du_dt_def=np.zeros(5*len(x)+1);
        du_dt_def[:-1]=du_dt[:];
        return du_dt_def  



def f(t,u):
    Dif=D[4*n:,4*n:];
    dx=x[1];
    nabla_2=diags([1/(dx**2),-2/(dx**2), 1/(dx**2)], [-1, 0, 1],shape=(n,n)).tolil();
    diffusivo=(Dif*nabla_2);
    resp=Reacc(u)
    nx=len(x)
    
    if constante==True:
        do2_dx=diffusivo.dot(u[4*nx:]);
        resp[4*nx:]+= do2_dx;   
        if bc[0]==0:
            resp[4*nx]=0;
        elif bc[0]==1:
            resp[4*nx]=-2*u[4*nx]/(dx**2)+2*u[4*nx+1]/(dx**2);
        if bc[1]==0:
            resp[-1]=0;
        elif bc[1]==1:
            resp[-1]=-2*u[-1]/(dx**2)+2*u[-2]/(dx**2);
            
    else:
       do2_dx=diffusivo.dot(u[4*nx:-1]);
       resp[4*nx:-1]+= do2_dx;
       resp[-1]=(A/(V+0.0))*(D[4*n,4*n]/(2.0*dx))*(-3*u[4*nx]+4*u[4*nx+1]-u[4*nx+2]-3*u[-2]+4*u[-3]-u[-4])
       if bc[0]==0:
           resp[4*nx]=S[0,0]*resp[-1];
       elif bc[0]==1:   
           resp[4*nx]=-2*u[4*nx]/(dx**2)+2*u[4*nx+1]/(dx**2);
       if bc[1]==0:   
           resp[-2]=S[-1,-1]*resp[-1];
       elif bc[1]==1:
           resp[-2]=-2*u[-2]/(dx**2)+2*u[-3]/(dx**2);     
    return resp

def f_4(t,u):
    global x, A, V
    Dif=D[4*n:,4*n:];
    dx=x[1];
    nabla_2=diags([1/(dx**2),-2/(dx**2), 1/(dx**2)], [-1, 0, 1],shape=(n,n)).tolil();
    diffusivo=(Dif*nabla_2);
    resp=Reacc_2(u)
    nx=len(x)
    do2_dx=diffusivo.dot(u[4*nx:-1]);
    resp[4*nx:-1]+= do2_dx;
    resp[-1]=(A/(V+0.0))*(D[4*n,4*n]/(2.0*dx))*(-3*u[-2]+4*u[-3]-u[-4])
    resp[4*nx]=0;
    resp[-2]=S[-1,-1]*resp[-1];
    return resp


OS_0_lista=[0,0.125,0.25,0.5,1];

L=0.06;
n=51;
R=8.314;
T=80+273.0;
A=10 #cm2
V=100 #cm3
x=np.linspace(0,L,n); 
dx=x[1]
Propiedades_Material()
co=(0.21*101325/(R*T))*1e-6; #mol/cm3
ci=(0.05*101325/(R*T))*1e-6; #mol/cm3
bc=[0,1];
constante=False;
k=np.asarray([4.09e13*np.exp(-101.5e3/(R*T)),1e10,3.52e9*np.exp(-31e3/(R*T)),3e14,1e13,3.05e18*np.exp(-48.4e3/(R*T))])
tf=7200*12*32;
t=np.linspace(0,tf,tf/3600 +1);

Longitudes=[0.06]
for i in range (len(Longitudes)):
    L=Longitudes[i];
    n=51;
    R=8.314;
    T=80+273.0;
    A=10 #cm2
    V=100 #cm3
    x=np.linspace(0,L,n); 
    dx=x[1]
    Propiedades_Material()
    co=(0.21*101325/(R*T))*1e-6; #mol/cm3
    ci=(0.05*101325/(R*T))*1e-6; #mol/cm3
    bc=[0,1];
    constante=False;
    k=np.asarray([4.09e13*np.exp(-101.5e3/(R*T)),1e10,3.52e9*np.exp(-31e3/(R*T)),3e14,1e13,3.05e18*np.exp(-48.4e3/(R*T))])
    tf=7200*50;
    t=np.linspace(0,tf,tf/3600 +1);
    OS_0=[0.5];
    uo=ic(OS_0,co,0,constante,co)
    sol = solve_ivp(f, [0,tf],uo,method='BDF',t_eval=t,rtol=1e-11,atol=1e-15);
    

##k0: 2POOH=POO° +P°,  k1: P°+O_2=POO°
##k2: POO°+PH=POOH+P°, k3: 2P°= ...
##k4: P° +POO°=...,    k5: POO° +POO° = ...
#
##-----------------Especies-------------------------
##u[0],=POOH, u[1]=POO°, u[2]= P°, u[3]=PH, u[4]=O_2
#
#    if constante==True:
#        dO2=k[1]*sol.y[2*n:3*n,:]*sol.y[4*n:,:];
#    else:
#        dO2=k[1]*sol.y[2*n:3*n,:]*sol.y[4*n:-1,:];
#
#    FFO2=0.5*(dO2[:,1:]+dO2[:,:-1])*(t[1:]-t[:-1]);
#
#    O2=np.zeros((n,len(t)))
#    for i in range(1,len(t)-1):
#        O2[:,i]=O2[:,i-1]+FFO2[:,i]
#        O2[:,-1]=O2[:,-2]+FFO2[:,-1]    

#    f1=plt.figure();
#    a = f1.gca(projection='3d')
#    X, Y = np.meshgrid(x, sol.t/3600.0)
#    surf = a.plot_surface(X, Y, O2.T, cmap=cm.jet,linewidth=0,edgecolor='none', antialiased=False)
#    c1=f1.colorbar(surf,ax = a,aspect=20, shrink=0.5,format='%.2e')
#    c1.set_label(r'$O_2$ [mol/$cm^3$]')
#    a.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
#    a.view_init(azim=-75, elev=30)
#    a.set_xlabel('Posicion [cm]')
#    a.set_ylabel('tiempo [h]')
#    a.set_zlabel('$O_2$ [mol/$cm^3$]')

    if constante==False:
        #f2=plt.figure()
        plt.plot(sol.t/(3600.0),sol.y[-1,:]/co,label=str(L) +'cm')
        I=sol.y[-1,:]/co;
        #plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}'))
        plt.ylabel('% $O_2$ residual')
        plt.xlabel('tiempo [h]')
        leg = plt.legend();


#f3=plt.figure();
#a = f3.gca(projection='3d')
#X, Y = np.meshgrid(x, sol.t/3600.0)
#if constante==True:
#    surf = a.plot_surface(X, Y, sol.y[4*n:,:].T, cmap=cm.jet,linewidth=0,edgecolor='none', antialiased=False)
#else:
#    surf = a.plot_surface(X, Y, sol.y[4*n:-1,:].T, cmap=cm.jet,linewidth=0,edgecolor='none', antialiased=False)
#    c1=f3.colorbar(surf,ax = a,aspect=20, shrink=0.5,format='%.2e')
#    c1.set_label(r'$O_2$ [mol/$cm^3$]')
#    a.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
#    a.view_init(azim=-75, elev=30)
#    a.set_xlabel('Posicion [cm]')
#    a.set_ylabel('tiempo [h]')
#    a.set_zlabel('$O_2$ [mol/$cm^3$]')


