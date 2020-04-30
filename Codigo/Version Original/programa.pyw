# -*- coding: utf-8 -*-
"""
@author: rozog
"""

from tkinter import *
from tkinter import ttk
from PIL.ImageTk import PhotoImage
from PIL import Image
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.integrate import solve_ivp


raiz=Tk()

#------------Vetana Principal --------------
raiz.title("Herramienta de Diseno V1.0")
raiz.resizable(1,1) #width,height
raiz.iconbitmap("logo.ico")
#raiz.geometry("650x350")


#----------Panel Derecha-------------------
Frame_Der=Frame(raiz,width=2400,height=600)
Frame_Der.config(bd=10,relief="sunken")
Frame_Der.pack(side="right",anchor="n")
titulo_param=Label(Frame_Der,text="Parametros")
titulo_param.grid(row=0,column=0,columnspan=3)


##Label Panel Derecha
label_tipo=Label(Frame_Der,text="Tipo de Pelicula");
label_tipo.grid(row=1,column=0)
label_grosor=Label(Frame_Der,text="Grosor Pelicula");
label_grosor.grid(row=2,column=0)
label_unid_grosor=Label(Frame_Der,text="cm");
label_unid_grosor.grid(row=2,column=2)
label_Area=Label(Frame_Der,text="Area Pelicula");
label_Area.grid(row=3,column=0)
label_unid_Area=Label(Frame_Der,text="cm^2");
label_unid_Area.grid(row=3,column=2)
label_vol_finito=Label(Frame_Der,text="Volumen Finito: ");
label_vol_finito.grid(row=4,column=0)
label_vol=Label(Frame_Der,text="Volumen: ");
label_vol.grid(row=5,column=0)
label_unid_vol=Label(Frame_Der,text='cm^3');
label_unid_vol.grid(row=5,column=2)
label_O2_ini=Label(Frame_Der,text='Concentracion O_2:');
label_O2_ini.grid(row=6,column=0)
label_O2_ini_unid=Label(Frame_Der,text='mol/cm^3');
label_O2_ini_unid.grid(row=6,column=2)
label_OS_ini=Label(Frame_Der,text='Concentracion OS:');
label_OS_ini.grid(row=7,column=0)
label_OS_ini_unid=Label(Frame_Der,text='mol/cm^3');
label_OS_ini_unid.grid(row=7,column=2)
label_Dr=Label(Frame_Der,text='Difusividad Pel.React:');
label_Dr.grid(row=8,column=0)
label_Dr_unid=Label(Frame_Der,text='cm^2/s');
label_Dr_unid.grid(row=8,column=2)
label_tiempo_total=Label(Frame_Der,text='Tiempo:');
label_tiempo_total.grid(row=9,column=0)
label_tiempo_total_un=Label(Frame_Der,text='h');
label_tiempo_total_un.grid(row=9,column=2)

#EntradasPanel Derecha
longitud=StringVar()
Area=StringVar();
Volumen=StringVar();
O2_inicial=StringVar();
OS_inicial=StringVar();
Difusividad_reac=StringVar();
tiempo_tot=StringVar();

entry_grosor=Entry(Frame_Der,width=12,textvariable=longitud,justify="right");
entry_grosor.grid(row=2,column=1)
entry_area=Entry(Frame_Der,width=12,textvariable=Area,justify="right");
entry_area.grid(row=3,column=1)
entry_volumen=Entry(Frame_Der,width=12,textvariable=Volumen,justify="right");
entry_volumen.grid(row=5,column=1)
entry_O2_inicial=Entry(Frame_Der,width=12,textvariable=O2_inicial,justify="right");
entry_O2_inicial.grid(row=6,column=1)
entry_OS_inicial=Entry(Frame_Der,width=12,textvariable=OS_inicial,justify="right");
entry_OS_inicial.grid(row=7,column=1)
entry_Dr=Entry(Frame_Der,width=12,textvariable=Difusividad_reac,justify="right");
entry_Dr.grid(row=8,column=1)
entry_t=Entry(Frame_Der,width=12,textvariable=tiempo_tot,justify="right");
entry_t.grid(row=9,column=1)

finito=IntVar()
finito.set(1);
def cambiar_vol():
    if finito.get()==2:
        Volumen.set("");
        entry_volumen.config(state='disabled')
    else:
        entry_volumen.config(state='normal')

radio_si=Radiobutton(Frame_Der,text='Si',variable=finito,value=1,command=cambiar_vol)
radio_si.grid(row=4,column=1)
radio_no=Radiobutton(Frame_Der,text='No',variable=finito,value=2,command=cambiar_vol)
radio_no.grid(row=4,column=2)



##Codigo Boton Correr
Boton_correr=Button(Frame_Der,text='Calcular Desempeno',width=15,command=lambda:correr());
Boton_correr.grid(row=10,column=0,columnspan=3,pady=2)




##--------------Subpanel Multicapa--------------
def tipo_pelicula(event=None):
    if tipo.get()=="Homogenea":
        Frame_multicapa.grid_forget();
        Boton_correr.grid(row=10,column=0,columnspan=3,pady=2) 
    else:
        Frame_multicapa.grid(row=10,column=0,columnspan=3) 
        Boton_correr.grid(row=11,column=0,columnspan=3,pady=2) 

tipo=StringVar()
tipo.set("Homogenea")
lista_tipo=ttk.Combobox(Frame_Der,width=15,textvariable=tipo,state="readonly")
lista_tipo.bind("<<ComboboxSelected>>", tipo_pelicula)
lista_tipo['values']=("Homogenea","Multicapa ABA")
lista_tipo.config(justify='center')
lista_tipo.grid(row=1,column=1,columnspan=2)
lista_tipo.current()

Frame_multicapa=Frame(Frame_Der, highlightbackground="black",highlightthickness=1);
load = Image.open("multicapa.png")
load=load.resize((250, 250), Image.ANTIALIAS)
render = PhotoImage(load,master=Frame_multicapa)
Imagen=Label(Frame_multicapa,image=render).grid(row=0,column=0,columnspan=3,pady=2)
label_L1=Label(Frame_multicapa,text='L1:');
label_L1.grid(row=1,column=0)
label_L1_un=Label(Frame_multicapa,text='cm');
label_L1_un.grid(row=1,column=2)
label_L2=Label(Frame_multicapa,text='L2:');
label_L2.grid(row=2,column=0)
label_L2_un=Label(Frame_multicapa,text='cm');
label_L2_un.grid(row=2,column=2)
label_Di=Label(Frame_multicapa,text='Difusividad Pel. Inerte:');
label_Di.grid(row=3,column=0)
label_Di_unid=Label(Frame_multicapa,text='cm^2/s');
label_Di_unid.grid(row=3,column=2)

#Entradas  subPanel multicapa
Difusividad_inert=StringVar();
longitud_L1=StringVar();
longitud_L2=StringVar();
entry_l1=Entry(Frame_multicapa,width=12,textvariable=longitud_L1,justify="right");
entry_l1.grid(row=1,column=1)
entry_l2=Entry(Frame_multicapa,width=12,textvariable=longitud_L2,justify="right");
entry_l2.grid(row=2,column=1)
entry_Di=Entry(Frame_multicapa,width=12,textvariable=Difusividad_inert,justify="right");
entry_Di.grid(row=3,column=1)

##------------------Panel Izquierda-----$$
Frame_izq=Frame(raiz,width=2400,height=600)
Frame_izq.config(bd=10,relief="sunken")
Frame_izq.pack(side="left",anchor="n")
titulo_result=Label(Frame_izq,text="Resultados")
titulo_result.grid(row=0,column=0,columnspan=3)

cuaderno=ttk.Notebook(Frame_izq)
cuaderno.grid(row=1,column=0,columnspan=3)


tiempo_result=StringVar();
O2_result=StringVar();
Label_tiempo_result=Label(Frame_izq, text='Tiempo S.S:');
Label_tiempo_result.grid(row=2,column=0)
entry_tiempo_result=Entry(Frame_izq,width=12,textvariable=tiempo_result,justify="right");
entry_tiempo_result.grid(row=2,column=1)
entry_tiempo_result.config(state='readonly')
Label_tiempo_result_un=Label(Frame_izq, text='h');
Label_tiempo_result_un.grid(row=2,column=2)
Label_O2_result=Label(Frame_izq, text='O_2 Total Absorbido:');
Label_O2_result.grid(row=3,column=0)
entry_O2_result=Entry(Frame_izq,width=12,textvariable=O2_result,justify="right");
entry_O2_result.grid(row=3,column=1)
entry_O2_result.config(state='readonly')
Label_O2_result_un=Label(Frame_izq, text='g');
Label_O2_result_un.grid(row=3,column=2)




grafico_1=Frame(cuaderno)
grafico_12=Frame(cuaderno)
grafico_2=Frame(cuaderno)
grafico_22=Frame(cuaderno)
grafico_3=Frame(cuaderno)
cuaderno.add(grafico_1,text='O2 Pel.')
cuaderno.add(grafico_12,text='O2 Pel. Contorno')
cuaderno.add(grafico_2,text='OS Pel.')
cuaderno.add(grafico_22,text='OS Pel. Contorno')
cuaderno.add(grafico_3,text='O2 HeadSpace')
    
f1 = Figure(figsize=(6,5), dpi=100)
f12 = Figure(figsize=(6,5), dpi=100)
f2 = Figure(figsize=(6,5), dpi=100)
f22 = Figure(figsize=(6,5), dpi=100)
f3 = Figure(figsize=(6,5), dpi=100)

##Figura 1-----------------     
canvas = FigureCanvasTkAgg(f1, grafico_1)
canvas.draw()
canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=False)
a = f1.gca(projection='3d')
toolbar = NavigationToolbar2Tk(canvas,grafico_1)
toolbar.update()
canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=False)

##Figura 1.2-----------------     
canvas12 = FigureCanvasTkAgg(f12, grafico_12)
canvas12.draw()
canvas12.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=False)
a12 = f12.add_subplot(111)
toolbar = NavigationToolbar2Tk(canvas12,grafico_12)
toolbar.update()
canvas12._tkcanvas.pack(side=TOP, fill=BOTH, expand=False)




##Figura 2-----------------   
canvas2 = FigureCanvasTkAgg(f2, grafico_2)
canvas2.draw()
canvas2.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=False)
a2 = f2.gca(projection='3d')
toolbar2 = NavigationToolbar2Tk(canvas2,grafico_2)
toolbar2.update()
canvas2._tkcanvas.pack(side=TOP, fill=BOTH, expand=False)


##Figura 1.2-----------------     
canvas22 = FigureCanvasTkAgg(f22, grafico_22)
canvas22.draw()
canvas22.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=False)
a22 = f22.add_subplot(111)
toolbar = NavigationToolbar2Tk(canvas22,grafico_22)
toolbar.update()
canvas22._tkcanvas.pack(side=TOP, fill=BOTH, expand=False)



##Figura 3-----------------
canvas3 = FigureCanvasTkAgg(f3, grafico_3)
canvas3.draw()
canvas3.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=False)
a3 = f3.add_subplot(111)

toolbar3 = NavigationToolbar2Tk(canvas3,grafico_3)
toolbar3.update()
canvas3._tkcanvas.pack(side=TOP, fill=BOTH, expand=False)


##PArametros Iniciales
Difusividad_reac.set(np.str(2.4e-9))
longitud.set(np.str(35e-4))
longitud_L1.set(np.str(13e-4))
longitud_L2.set(np.str(9e-4))
Difusividad_inert.set(np.str(4.81e-9))
Area.set(np.str(25))
Volumen.set(np.str(70));
OS_inicial.set(np.str(1.186e-5))
O2_inicial.set(np.str(8.561e-6))
tiempo_tot.set(2)








#----------------Funcionamiento diseno

def correr():
 global Di, Dr,L,L1,L2,nu,k,A,V,no,s,co,n,x,to,tf
 Dr=np.float(Difusividad_reac.get());
 Di=Dr;
 L=np.float(longitud.get());
 if tipo.get()=="Homogenea": 
     L1=0;
     L2=L;
 else:
     L1=np.float(longitud_L1.get());
     L2=np.float(longitud_L2.get());
     Di=np.float(Difusividad_inert.get());
 nu=0.5;
 k=1e4;
 A=np.float(Area.get());
 V=np.float(Volumen.get());
 no=np.float(OS_inicial.get());
 s=7.16e-2;
 co=np.float(O2_inicial.get());
 n=31;
 x=np.linspace(0,L,n);
 to=0
 tf=np.float(tiempo_tot.get())*3600;
 pde(x,to,tf)


def D(x):
    global Di, Dr,L,L1,L2,nu,k,A,V,no,s,co,n,to,tf
    resp=np.zeros((len(x)));
    for i in range(len(x)):
        if x[i]<L1 or x[i]>(L1+L2):
            resp[i]=Di;
        else:
            resp[i]=Dr;
    return resp;    
    
def f(t,u):
    #Parametros del sistema a estudiar;
    global Di, Dr,L,L1,L2,nu,k,A,V,no,s,co,n,x,to,tf
    dx=x[1]-x[0]
    #Calculo de la derivada temporal 
    du_dt=np.zeros((2*n+1));
    du_dt[1:n-1]=D(x[1:n-1])*((-2*u[1:n-1]+u[2:n]+u[0:n-2])/(dx**2))-k*u[1:n-1]*u[n+1:(2*n-1)]
    du_dt[n:(2*n)]=-nu*k*u[:n]*u[n:(2*n)];
    du_dt[2*n]=2*(A/V)*Di*(-u[2]+4*u[1]-3*u[0])/(2*dx);
    du_dt[0]=du_dt[2*n]*s
    du_dt[n-1]=du_dt[2*n]*s
    return du_dt

def ic(x):
    global Di, Dr,L,L1,L2,nu,k,A,V,no,s,co,n,to,tf
    vo=np.zeros((len(x)));
    for i in range(len(x)):
        if x[i]<= (L1+L2) and x[i]>=L1:
            vo[i]=no;
    uo=np.zeros((len(x)));
    uo[0]=co*s;
    uo[-1]=uo[0];
    uo=np.append(uo,vo);
    uo=np.append(uo,co);
    return uo
    


def pde(x,to,tf):
    dx=x[1];
    uo=ic(x);
    
    u_3=solve_ivp(f,[to,tf],uo,method='BDF');
    t=u_3.t;
    u=u_3.y.T;
     
    O2=u[:,:n];
    OS=u[:,n:2*n];
    HS=u[:,-1];
    
    
    f1.clear()
    a = f1.gca(projection='3d')
    X, Y = np.meshgrid(x, t)
    surf = a.plot_surface(X, Y, O2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    f1.colorbar(surf, shrink=0.5, aspect=5)
    a.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
    a.view_init(azim=-75, elev=30)
    a.set_xlabel('Posicion [cm]')
    a.set_ylabel('tiempo [s]')
    a.set_zlabel('$O_2$ [mol/$cm^3$]')
    canvas.draw()
    
    f12.clear()
    a12=f12.add_subplot(111)
    cp = a12.contourf(X, Y, O2)
    f12.colorbar(cp) # Add a colorbar to a plot
    a12.set_title('Contorno Concentracion $O_2$ [mol/$cm^3$]')
    a12.set_xlabel('x (cm)')
    a12.set_ylabel('t (s)')
    canvas12.draw()
    

    f2.clear()
    a2 = f2.gca(projection='3d')
    a2.clear()
    surf = a2.plot_surface(X, Y, OS, cmap=cm.rainbow,linewidth=0, antialiased=False)
    a2.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
    f2.colorbar(surf, shrink=0.5, aspect=5)
    a2.set_xlabel('Posicion [cm]')
    a2.set_ylabel('tiempo [s]')
    a2.set_zlabel('concentracion [mol/$cm^3$]')
    canvas2.draw()
    
    
    f22.clear()
    a22=f22.add_subplot(111)
    cp = a22.contourf(X, Y, OS)
    f22.colorbar(cp) # Add a colorbar to a plot
    a22.set_title('Contorno Concentracion OS [mol/$cm^3$]')
    a22.set_xlabel('x (cm)')
    a22.set_ylabel('t (s)')
    canvas22.draw()
    
    
    a3.clear()
    a3.plot(t,HS)
    a3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    a3.set_xlabel('t (s)')
    a3.set_ylabel('$O_2$ Head Space [mol/$cm^3$]')
    canvas3.draw()
    
raiz.mainloop()
