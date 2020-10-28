# -*- coding: utf-8 -*-
"""
@author: rozog
"""

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter.filedialog import asksaveasfile
from PIL.ImageTk import PhotoImage
from PIL import Image
import numpy as np
import scipy.optimize as spo
import pandas as pd
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
from scipy.sparse import lil_matrix, csr_matrix, tril,triu ,diags
from scipy.sparse.linalg import cg
from scipy.sparse.linalg import spsolve 


raiz=Tk()

#------------Vetana Principal --------------
raiz.title("Herramienta de Diseno V2.2")
raiz.resizable(1,1) #width,height
raiz.iconbitmap("logo.ico")
#raiz.geometry("650x350")

#------------Barra de Menu---------------------

menubar = Menu(raiz);
raiz.config(menu=menubar) 
filemenu = Menu(menubar, tearoff=0)
ayudamenu= Menu(menubar, tearoff=0)
filemenu.add_command(label="Modo Desempeno", command=raiz.destroy)
filemenu.add_command(label="Modo Analisis", command=raiz.destroy)
filemenu.add_separator()
filemenu.add_command(label="Salir", command=raiz.destroy)
ayudamenu.add_command(label="Acerca de la Herramienta", command=lambda:crearVentanaAyuda())

menubar.add_cascade(label="Archivo", menu=filemenu)
menubar.add_cascade(label="Ayuda", menu=ayudamenu)

#----Informacion de Materiales------------------
data=pd.read_excel(r'./Datos/Data.xlsx')
materiales=data['Material'].values.tolist()
#----------Panel Derecha-------------------
Frame_Der=Frame(raiz,width=2400,height=600)
Frame_Der.config(bd=2,relief="raised")
Frame_Der.pack(side="right",anchor="n")

#---------------------Label Panel Derecha
Label(Frame_Der,text="Parametros").grid(row=0,column=0,columnspan=3)
Label(Frame_Der,text="Tipo de Pelicula:").grid(row=1,column=0);
Label(Frame_Der,text="Area Pelicula [cm^2]:").grid(row=2,column=0)
Label(Frame_Der,text="Volumen Finito: ").grid(row=3,column=0)
Label(Frame_Der,text="Volumen [cm^3]: ").grid(row=4,column=0);
Label(Frame_Der,text='Concentracion Oxigeno [mol/cm^3]:').grid(row=5,column=0);
Label(Frame_Der,text='Tiempo [h]: ').grid(row=6,column=0);
Label(Frame_Der,text='Temperatura [K]').grid(row=7,column=0);

#-------------------- Variables Panel Derecha
Area=StringVar();
Volumen=StringVar();
O2_inicial=StringVar();
tiempo_tot=StringVar();
temp=StringVar();
Entry(Frame_Der,width=12,textvariable=Area,justify="right").grid(row=2,column=1);
Entry(Frame_Der,width=12,textvariable=Volumen,justify="right").grid(row=4,column=1);
Entry(Frame_Der,width=12,textvariable=O2_inicial,justify="right").grid(row=5,column=1);
Entry(Frame_Der,width=12,textvariable=tiempo_tot,justify="right").grid(row=6,column=1);
Entry(Frame_Der,width=12,textvariable=temp,justify="right").grid(row=7,column=1);

finito=IntVar()
finito.set(1);
def cambiar_vol():
    if finito.get()==2:
        Volumen.set("");
        entry_volumen.config(state='disabled')
    else:
        entry_volumen.config(state='normal')

radio_si=Radiobutton(Frame_Der,text='No',variable=finito,value=1,command=cambiar_vol)
radio_si.grid(row=3,column=1)
radio_no=Radiobutton(Frame_Der,text='Si',variable=finito,value=2,command=cambiar_vol)
radio_no.grid(row=3,column=2)

#----------------------------------Combobox Panel Derecha
def tipo_pelicula(event=None):
    global Frame_Homogeneo, Frame_Multicapa, Frame_Der
    if tipo.get()=="Homogenea": 
        try:
            Frame_Multicapa.grid_forget();
        except:
            None
        Panel_homogeneo();
    else:
        try:
            Frame_Homogeneo.grid_forget();
        except:
            None    
        Panel_Multicapa();
    Frame_Der.update();

tipo=StringVar()
tipo.set("Homogenea")
lista_tipo=ttk.Combobox(Frame_Der,width=15,textvariable=tipo,state="readonly")
lista_tipo.bind("<<ComboboxSelected>>", tipo_pelicula)
lista_tipo['values']=("Homogenea","Multicapa")
lista_tipo.config(justify='center')
lista_tipo.grid(row=1,column=1,columnspan=2)
lista_tipo.current()

#----------------------Boton Correr
Button(Frame_Der,text='Calcular Desempeno',width=15,command=lambda:correr()).grid(row=10,column=0,columnspan=3,pady=2);

#////////////////////////////////////////////////////////////////////////////////////////

#------------------------Panel homogeneo
def Panel_homogeneo():
    global Frame_Homogeneo, Frame_Der
    global longitud, lista_variables_material,inicial_OS_n
    Frame_Homogeneo=Frame(Frame_Der);
            
    Label(Frame_Homogeneo,text="Grosor Pelicula [cm]:").grid(row=0,column=0);
    Label(Frame_Homogeneo,text='Concentracion OS [mol/cm^3]:').grid(row=1,column=0);
    Label(Frame_Homogeneo,text='Material:').grid(row=2,column=0);    

    #Entradas Panel Homogeneo
    longitud=[]
    lon=StringVar()
    inicial_OS_n=[]
    OS_inicial=StringVar();
    Entry(Frame_Homogeneo,width=12,textvariable=lon,justify="right").grid(row=0,column=1);
    Entry(Frame_Homogeneo,width=12,textvariable=OS_inicial,justify="right").grid(row=1,column=1);
    longitud.append(lon);
    inicial_OS_n.append(OS_inicial);   
    
    lista_variables_material=[]
    mat=StringVar();
    combo=ttk.Combobox(Frame_Homogeneo,width=5,textvariable=mat,state="readonly",values=materiales) 
    combo.grid(row=2,column=1)
    combo.config(justify='center')
    lista_variables_material.append(mat);
    
    Frame_Homogeneo.grid(row=8,column=0,columnspan=3)   

Panel_homogeneo();
#//////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------Subpanel Multicapa--------------              
def Panel_Multicapa():
    global capas, Frame_Multicapa, Frame_Der 
    Frame_Multicapa=Frame(Frame_Der)
    
     
    Label(Frame_Multicapa,text="Cantidad de capas:").grid(row=0,column=0)
    capas=StringVar();
        
    combo_capas=ttk.Combobox(Frame_Multicapa,width=5,textvariable=capas,state="readonly")
    combo_capas.bind("<<ComboboxSelected>>",numero_capas);
    combo_capas['values']=("2","3","4","5")
    combo_capas.config(justify='center')
    combo_capas.grid(row=0,column=1,columnspan=2)
    combo_capas.current()
    
    Frame_Multicapa.grid(row=8,column=0,columnspan=3)   
    
    
class Checkbar(Frame):
   def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
      Frame.__init__(self, parent)
      global inicial_OS_n
      inicial_OS_n=[]
      self.vars = []
      for pick in picks:
         var = IntVar()
         var2 = StringVar()
         chk = Checkbutton(self, text=pick, variable=var, command=imprimir)
         chk.pack(side=side, anchor=anchor, expand=YES)
         self.vars.append(var)
         inicial_OS_n.append(var2);
      
   def state(self):
      return map((lambda var: var.get()), self.vars)
      
       

  
def numero_capas(event=None):
    global capas, Frame_Multicapa, reactiva,Frame_longitud, Frame_material,Frame_multicapa
    global longitud, lista_variables_material 
    numero=np.int(capas.get());
    lista=np.linspace(1,numero,numero)
    try:
        Frame_longitud.destroy()
        Frame_material.destroy()
        reactiva.destroy()
        Frame_multicapa.destroy()
    except:
        None
    Frame_longitud=Frame(Frame_Multicapa)
    Frame_material=Frame(Frame_Multicapa)
    reactiva = Checkbar(Frame_Multicapa, lista);
    Label(Frame_Multicapa,text='Seleccione las capas que sean reactivas:').grid(row=1,columnspan=3)
    reactiva.grid(row=2,column=0, columnspan=3)
    longitud=[]
    lista_variables_material=[]
    for i in range(0,numero):
        lon=StringVar();
        mat=StringVar();
        Label(Frame_longitud,text='L'+str(i+1)+':').grid(row=i,column=0)
        Entry(Frame_longitud,textvariable=lon,width=12,justify="right").grid(row=i,column=1) 
        Label(Frame_material,text='Material Capa '+str(i+1)+':').grid(row=i,column=0)
        combo=ttk.Combobox(Frame_material,width=5,textvariable=mat,state="readonly",values=materiales) 
        combo.grid(row=i,column=1)
        combo.config(justify='center')
        longitud.append(lon)
        lista_variables_material.append(mat)
    Frame_longitud.grid(row=4,column=0,columnspan=3);
    Frame_material.grid(row=5,column=0,columnspan=3);
    

def imprimir(event=None):
    global reactiva,Frame_Multicapa,inicial_OS_n, Frame_multicapa
    estado=list(reactiva.state()); 
    numero=len(estado)
    try:
       Frame_multicapa.destroy();
       Frame_Multicapa.update()
    except:
       None
    Frame_multicapa=Frame(Frame_Multicapa, highlightbackground="black",highlightthickness=1);
    fila_or=0;
    for i in range(0,numero):
        if estado[i]==1:
            Label(Frame_multicapa,text='Concentracion OS Capa ' + np.str(i+1)+' :').grid(row=fila_or,column=0);
            Entry(Frame_multicapa,width=12,textvariable=inicial_OS_n[i],justify="right").grid(row=fila_or,column=1);
            fila_or+=1;            
    Frame_multicapa.grid(row=3,column=0,columnspan=3)

#//////////////////////////////////////////////////////////////////////////////////////////////////

##------------------Panel Izquierda-----$$
Frame_izq=Frame(raiz,width=2400,height=600)
Frame_izq.config(bd=10,relief="sunken")
Frame_izq.pack(side="left",anchor="n")
titulo_result=Label(Frame_izq,text="Resultados")
titulo_result.grid(row=0,column=0,columnspan=2)

cuaderno=ttk.Notebook(Frame_izq)
cuaderno.grid(row=1,column=0,columnspan=2)

#----------Label Frames Izquierdo-------------------
Label(Frame_izq, text='Tiempo S.S. [h]:').grid(row=2,column=0);
Label(Frame_izq, text=r'$O_2$ Total Absorbido [g]:').grid(row=3,column=0);
#------------Resultados Frame Izquierdo------------
tiempo_result=StringVar();
O2_result=StringVar();
Entry(Frame_izq,width=12,textvariable=tiempo_result,justify="right",state='readonly').grid(row=2,column=1);
Entry(Frame_izq,width=12,textvariable=O2_result,justify="right",state='readonly').grid(row=3,column=1);

#--------------Boton de almacenamiento de resultados------------------
Button(Frame_izq,text='Guardar Resultados',width=15,command=lambda:guardar()).grid(row=4,column=0,columnspan=2,pady=2);

#---------------Inclusion de Pesetañas Frame Izquierdo ------------------
grafico_1=Frame(cuaderno)
grafico_12=Frame(cuaderno)
grafico_2=Frame(cuaderno)
grafico_22=Frame(cuaderno)
grafico_3=Frame(cuaderno)
grafico_4=Frame(cuaderno)
cuaderno.add(grafico_1,text='O2 Pel.')
cuaderno.add(grafico_12,text='O2 Pel. Contorno')
cuaderno.add(grafico_2,text='OS Pel.')
cuaderno.add(grafico_22,text='OS Pel. Contorno')
cuaderno.add(grafico_3,text='O2 HeadSpace')
cuaderno.add(grafico_4,text='O2 Masa Absorbida')

f1 = Figure(figsize=(6,5), dpi=100)
f12 = Figure(figsize=(6,5), dpi=100)
f2 = Figure(figsize=(6,5), dpi=100)
f22 = Figure(figsize=(6,5), dpi=100)
f3 = Figure(figsize=(6,5), dpi=100)
f4 = Figure(figsize=(6,5), dpi=100)

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

##Figura 4-----------------
canvas4 = FigureCanvasTkAgg(f4, grafico_4)
canvas4.draw()
canvas4.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=False)
a4 = f4.add_subplot(111)

toolbar4 = NavigationToolbar2Tk(canvas4,grafico_4)
toolbar4.update()
canvas4._tkcanvas.pack(side=TOP, fill=BOTH, expand=False)


#///////////////////////////////////////////////////////////////////////////////////////////////////////

##PArametros Iniciales
#Difusividad_reac.set(np.str(2.4e-9))
#longitud_L1.set(np.str(13e-4))
#longitud_L2.set(np.str(9e-4))
#Difusividad_inert.set(np.str(4.81e-9))
Area.set(np.str(25))
Volumen.set(np.str(70));
#OS_inicial.set(np.str(1.186e-5))
O2_inicial.set(np.str(8.561e-6))
tiempo_tot.set(2)

#///////////////////////////////////////////////////////////////////////////////////////////////

#----------------Funcionamiento diseno
def Propiedades_Material(L, mat):
    global materiales,R,T,x,D,S
    actual=materiales.index(mat[0].get())                         
   
    Ea=np.float(data.iloc[actual,2]);
    Do=np.float(data.iloc[actual,1]);
    D=Do*np.exp(-Ea/(R*T))
        
    Delta_H=np.float(data.iloc[actual,4]);
    So=np.float(data.iloc[actual,3]);
    S=So*np.exp((-Delta_H/R)*(1.0/T -1.0/298.0))
    
    resp_D=D*np.ones((len(x)));   
    resp_S=S*np.ones((len(x)));  
    for i in range(len(L)-1):
        actual=materiales.index(mat[i+1].get())
                
        Ea=np.float(data.iloc[actual,2]);
        Do=np.float(data.iloc[actual,1]);
        D=Do*np.exp(-Ea/(R*T));
        
        Delta_H=np.float(data.iloc[actual,4]);
        So=np.float(data.iloc[actual,3]);
        S=So*np.exp((-Delta_H/R)*(1.0/T -1.0/298));
        
        for j in range(len(x)):
            if x[j]>=L[i]:
                resp_D[j]=D;
                resp_S[j]=S;
                
    D=diags(resp_D).tolil();
    S=diags(resp_S*1e6*R*T).tolil();
    
def A_matriz(dt):
    global D ,x
    dx=x[1];
    n=len(x);
    #----------Calculo de las Matrices de Difussion-----------------
    nabla_2=diags([1/(dx**2),-2/(dx**2), 1/(dx**2)], [-1, 0, 1],shape=(n,n)).tolil();
    identity=diags([1], [0],shape=(n,n)).tolil();
    A_izq=identity-((dt)/2.0)*D*nabla_2;
    A_der=identity+((dt)/2.0)*D*nabla_2;
    
    #------Condiciones de FRonterea de Dirichlet-------------------
    A_izq[0,1]=0; A_izq[0,0]=1; A_der[0,1]=0; A_der[0,0]=1;
    A_izq[-1,-2]=0;A_izq[-1,-1]=1;A_der[-1,-2]=0;A_der[-1,-1]=1;
      
    return A_izq, A_der

def Reacc(u):
    global x,T,R
    k=np.asarray([4.09e10*np.exp(-101.5e3/(R*T)),1e7,3.52e6*np.exp(-31e3/(R*T)),0,0,3.05e15*np.exp(-48.4e3/(R*T))])
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
    du_dt[3,:]=-k[2]*especies[1,:]*especies[3,:];
    du_dt[4,:]=-k[1]*especies[2,:]*especies[4,:];
    du_dt=np.reshape(du_dt,5*len(x)); 
    return du_dt


def dCext_dt(u):
    global D,x
    # ----------------Ecuacion------------------------
    # dCext_dt=A/V ([D*dO_2/dx]_x=0 - [D*dO_2/dx]_x=L)
    #----------Calculo de la derivada temporal -------
    return (A/(V+0.0))*(D[0,0]*(-u[2]+4*u[1]-3*u[0])/(2.0*x[1]) +D[-1,-1]*(-u[-3]+4*u[-2]-3*u[-1])/(2.0*x[1])) ;


def ic(L,OS_0,co):
    global S, T, R, x 
    #---Concentracion inicial de  Espcies en Aceite de Linaza---
    Lin_oil_Co=np.asarray([0.05,0,0,6.8]);
    #--------------Vector de concentrracion inicial ------------
    uo=np.zeros((5,len(x)));
    for k  in  range(4): uo[k,:]=OS_0[0]*Lin_oil_Co[k];   
    
    for j in range(len(L)-1):
        for i in range(len(x)):
            if x[i]>=L[j]: 
                for k  in  range(4): uo[k,i]=OS_0[j+1]*Lin_oil_Co[k]; 
    #-------------Concentracion Inicial de O2 -------------------
    uo[4,0]=co*S[0,0];
    uo[4,-1]=co*S[-1,-1];
    uo=np.reshape(uo,5*len(x));
    return uo 

def frac_met(uo,co,t):
    dx=x[1];
    dt=t[1];
    n=len(x);
    # --------------- Barra de Progreso--------------
    top = Toplevel();
    top.title("Progreso")
    barra= ttk.Progressbar(top, orient = HORIZONTAL, 
              length = 100,mode = 'determinate');
    barra.pack(pady = 10) 
    barra.update()
    
    #------- Inicializacion y matrices a solucinar---
    u=np.zeros((5*len(x), len(t)));
    c_ext=np.zeros((len(t)));
    u[:,0]=uo;
    c_ext[0]=co;
    A_izq, A_der=A_matriz(dt);
        
    for i in range(0,len(t)-1):
        u_actual=u[:,i].copy();
        u_medio=u_actual.copy();
        #---------- Solucion Difusion ---------------
        O2_actual=u[4*n:,i];
        O2_actual[0]=S[0,0]*c_ext[i];O2_actual[-1]=S[-1,-1]*c_ext[i];
        b=A_der.dot(O2_actual);
        u_medio[4*n:]=spsolve(A_izq,b);       
        #---------- Solucion Reaccion ---------------
        def g(u): return u-0.5*dt*(Reacc(u_medio)+Reacc(u))-u_medio;
        u[:,i+1]=spo.fsolve(g,u_medio)
        
        #---------- Actualizacion Headspace ---------
        c_ext[i+1]=dt*dCext_dt(u_actual[4*n:])+c_ext[i];
                
         #-------Actualizacion Barra de Progreso -----
        barra['value'] = ((i+1.0)/len(t))*100;
        barra.update()
    
    top.destroy()
    return c_ext,u

def frac_met_paso_varible(uo,co,tf):
    dx=x[1];
    n=len(x);
    #-------Parametros Metodo Adaptativo----------
    dt=0.001;
    nu_min=1e-7;
    nu_max=0.01;
    dt_max=1000;
    dt_min=1e-3;
    rho=5;
    sigma=0.2;
    # --------------- Barra de Progreso--------------
    top = Toplevel();
    top.title("Progreso")
    barra= ttk.Progressbar(top, orient = HORIZONTAL, 
              length = 100,mode = 'determinate');
    barra.pack(pady = 10) 
    barra.update()
    
    #-------------- Inicializacion-----------
    u=[]
    c_ext=[]
    u.append(uo)
    c_ext.append(co);
    u_actual=uo;
    c_ext_actual=co;
    
    tiempo=[]
    tiempo.append(0)
    t=0;
    while t<tf: 
        u_medio=u_actual.copy();
        #---------- Solucion Difusion ---------------
        O2_actual=u_actual[4*n:];
        O2_actual[0]=S[0,0]*c_ext_actual;O2_actual[-1]=S[-1,-1]*c_ext_actual;
        A_izq, A_der=A_matriz(dt); 
        b=A_der.dot(O2_actual);
        u_medio[4*n:]=spsolve(A_izq,b);       
        #---------- Solucion Reaccion ---------------
        def g(u): return u-0.5*dt*(Reacc(u_medio)+Reacc(u))-u_medio;
        u_nuevo=spo.fsolve(g,u_medio)  
        #---------- Actualizacion Headspace ---------
        c_ext_nuevo=dt*dCext_dt(u_actual[4*n:])+c_ext_actual;
        nu=np.sum(np.abs(u_nuevo-u_actual))/(np.sum(np.abs(u_actual))+1)
            
        if nu>=nu_min and nu<=nu_max:            
            u.append(u_nuevo);
            c_ext.append(c_ext_nuevo);
            u_actual=u_nuevo.copy();
            c_ext_actual=c_ext_nuevo.copy();
            t+=dt;
            if t> tf: 
                t=tf;
                dt=tf-tiempo[-1];
            tiempo.append(t);
            
        elif nu<nu_min:
            u.append(u_nuevo);
            c_ext.append(c_ext_nuevo);
            u_actual=u_nuevo.copy();
            c_ext_actual=c_ext_nuevo.copy();
            t+=dt;
            dt*=rho;
            if t> tf: 
                t=tf;
                dt=tf-tiempo[-1];    
            if dt>dt_max: dt=dt_max;
            tiempo.append(t);
        else:
            dt*=sigma;
            if dt <dt_min: dt=dt_min;        
                
        #-------Actualizacion Barra de Progreso -----
        barra['value'] = ((t+1.0)/tf)*100;
        barra.update()
    
    top.destroy()
    return tiempo,c_ext,np.array(u).T



def masa_absorbida(u,t):
    global x,D,A
    dx=x[1]
    O2=u[4*len(x):,:];
    
    Flux=-D[0,0]*(-O2[2,:]+4*O2[1,:]-3*O2[0,:])/(2.0*dx)-D[-1,-1]*(-O2[-3,:]+4*O2[-2,:]-3*O2[-1,:])/(2.0*dx)
        
    masa=[];
    for i in range(len(t)): 
        masa.append(64*A*np.trapz(Flux[:i], t[:i]));
    
    return masa

def tiempo_lag(masa,t):
    t_ss='-';
    duplicados = [y for n, y in enumerate(masa) if y in masa[:n]]
    if len(duplicados)!=0:t_ss=t[masa.index(duplicados[0])]/3600.0;
    return t_ss;

def variables_respuesta(masa,t_ss):
    O2_result.set(format(masa[-1], "10.2E"))
    tiempo_result.set(t_ss);
    if type(t_ss) is float: tiempo_result.set(format(t_ss, "10.2f"));

def graficador(c_ext,u,masa,t_ss,t):
    global x
    
    nx=len(x);
    HS=c_ext; 
    O2=u[4*nx:,:];
    OS=u[:nx,:];

    #--------------Gráfica O_2 en pelicula 3d------------------
    f1.clear()
    a = f1.gca(projection='3d')
    X, Y = np.meshgrid(x, t)
    surf = a.plot_surface(X, Y, O2.T, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    f1.colorbar(surf, shrink=0.5, aspect=5)
    a.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
    a.view_init(azim=-75, elev=30)
    a.set_xlabel('Posicion [cm]')
    a.set_ylabel('tiempo [s]')
    a.set_zlabel('$O_2$ [mol/$cm^3$]')
    canvas.draw()
    
    #--------------Gráfica O_2 en pelicula en contorno------------------
    f12.clear()
    a12=f12.add_subplot(111)
    cp = a12.contourf(X, Y, O2.T)
    f12.colorbar(cp) # Add a colorbar to a plot
    a12.set_title('Contorno Concentracion $O_2$ [mol/$cm^3$]')
    a12.set_xlabel('x (cm)')
    a12.set_ylabel('t (s)')
    canvas12.draw()
    
    #--------------Gráfica POOH en pelicula 3d------------------
    f2.clear()
    a2 = f2.gca(projection='3d')
    a2.clear()
    surf = a2.plot_surface(X, Y, OS.T, cmap=cm.rainbow,linewidth=0, antialiased=False)
    a2.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
    f2.colorbar(surf, shrink=0.5, aspect=5)
    a2.set_xlabel('Posicion [cm]')
    a2.set_ylabel('tiempo [s]')
    a2.set_zlabel('concentracion [mol/$cm^3$]')
    canvas2.draw()
    
    #--------------Gráfica POOH en pelicula en contorno------------------
    f22.clear()
    a22=f22.add_subplot(111)
    cp = a22.contourf(X, Y, OS.T)
    f22.colorbar(cp) # Add a colorbar to a plot
    a22.set_title('Contorno Concentracion OS [mol/$cm^3$]')
    a22.set_xlabel('x (cm)')
    a22.set_ylabel('t (s)')
    canvas22.draw()
    
    #---------------Grafica de Concentracion O_2 en Headspace ------------
    a3.clear()
    a3.plot(t,HS,'-')
    a3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    a3.set_xlabel('t (s)')
    a3.set_ylabel('$O_2$ Head Space [mol/$cm^3$]')
    canvas3.draw()
    
    #---------------Grafica Absoricon de masa de O_2 en Pelicula ------------
    a4.clear()
    a4.plot(t,masa,'-')
    a4.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    a4.set_xlabel('t (s)')
    a4.set_ylabel('$O_2$ Total Mass [g]')
    canvas4.draw()


def almacenamiento_parametros(L,OS_0,co,t, material):
    global datos_parametros_entrada 
    global A, V,T
    
    variables=['Temperatura [K]:','Tiempo Total [h]:','Area Pelicula [cm^2]:','Concentracion O_2 [mol/cm^3]','Volumen [cm^3]','numero de capas :','Longitud Capa 1 [cm]:','Material Capa 1','Carga Aceite Linaza Capa 1 [goil/g_pelicula]']
    valores=[T,t[-1]/3600.0,A,co,V,len(L),L[0] , material[0].get(), OS_0[0]]
    for i  in range(1,len(L)):
        variables.append('Longitud Capa '+ i +' [cm]:','Material Capa '+ i,'Carga Aceite Linaza Capa '+ i +' [goil/g_pelicula]');
        valores.append(L[i],material[i].get(), OS_0[i])
    
    datos_parametros_entrada=pd.DataFrame({'Variable':variables, 'Valor':valores})
    
def almacenamiento_resultados(u,c_ext,masa,t):
    global datos_O2,datos_POOH,datos_PH, datos_HS, datos_O2_absorbido
    global x
    
    n=len(x)
    O2=u[4*n:,:];
    POOH=u[:n,:];
    PH=u[2*n:3*n,:];
    
    tabla_O2=np.zeros((len(t)+1,len(x)+1));
    tabla_O2[1:,0]=t;
    tabla_O2[0,1:]=x;
    tabla_O2[1:,1:]=O2.T;
    tabla_O2=tabla_O2.tolist();
    tabla_O2[0][0]='Tiempo [s]\ Posicion [m]'
    
    tabla_POOH=np.zeros((len(t)+1,len(x)+1));
    tabla_POOH[1:,0]=t;
    tabla_POOH[0,1:]=x;
    tabla_POOH[1:,1:]=POOH.T;
    tabla_POOH=tabla_POOH.tolist();
    tabla_POOH[0][0]='Tiempo [s]\ Posicion [m]'
    
    tabla_PH=np.zeros((len(t)+1,len(x)+1));
    tabla_PH[1:,0]=t;
    tabla_PH[0,1:]=x;
    tabla_PH[1:,1:]=PH.T;
    tabla_PH=tabla_PH.tolist();
    tabla_PH[0][0]='Tiempo [s]\ Posicion [m]'

    datos_O2=pd.DataFrame(tabla_O2)
    datos_POOH=pd.DataFrame(tabla_POOH)
    datos_PH=pd.DataFrame(tabla_PH)
    datos_HS=pd.DataFrame({'Tiempo [s]':t, 'O2 HeadSpace [mol/cm3]':c_ext})
    datos_O2_absorbido=pd.DataFrame({'Tiempo [s]':t, 'O2 Absorbido [g/g_oil]':masa})
    
    
def almacenamiento_variables_desempeno(masa,t_ss):
    global datos_var_desempeno
    variables=['Masa O_2 Total Absorbida [g/g_oil]', 'Tiempo Lag [h]']
    valores=[masa[-1],t_ss]
    datos_var_desempeno=pd.DataFrame({'Variable':variables, 'Valor':valores})

def guardar():
    global datos_parametros_entrada 
    global datos_O2,datos_POOH,datos_PH, datos_HS, datos_O2_absorbido
    global datos_var_desempeno
    
    files = [('Excel', '*.xlsx')] 
    file = filedialog.asksaveasfilename(filetypes = files, defaultextension = files)  
    
    writer = pd.ExcelWriter(file, engine='xlsxwriter')
    datos_parametros_entrada.to_excel(writer, sheet_name='Parametros de Entrada', index=False)
    datos_O2.to_excel(writer, sheet_name='Concentracion O2', index=False)
    datos_POOH.to_excel(writer, sheet_name='Concentracion POOH', index=False)
    datos_PH.to_excel(writer, sheet_name='Concentracion PH', index=False)
    datos_HS.to_excel(writer, sheet_name='Headspace O2 Concentration', index=False)
    datos_O2_absorbido.to_excel(writer, sheet_name='Masa O2 absorbido', index=False)
    datos_var_desempeno.to_excel(writer, sheet_name='Variables Desempeno', index=False)
    writer.save()

def correr():
    global temp,reactiva,longitud,lista_variables_material
    global T,x,A,V,R
    
    #--------------- Definicion Temperatura -------------
    T=np.float(temp.get());
    #--------------Cte de Gases-----------
    R=8.314;     
    #---------------Definicion mallas ------
    L=[];
    L.append(np.float(longitud[0].get()));
    n=21;
    x=np.linspace(0,L[-1],n);   
    # ----------------MAterial y Area de la Pelicula-------
    Propiedades_Material(L, lista_variables_material)
    A=np.float(Area.get());   

    #------- --Definicion de Capas Reactivas ---------
    capas=len(longitud)
    OS_0=np.zeros(capas);
    try:
        estado=list(reactiva.state());
    except:
         estado=np.zeros((1));     
    if capas==1: estado[0]=1; 
    for i in range(capas):
        if inicial_OS_n[i].get()!='' and estado[i]!=0:
            OS_0[i]=np.float(inicial_OS_n[i].get());
    for i in range(1,capas):
        L.append(L[i-1]+np.float(longitud[i].get()));
    
    #--------Parametros HeadSpace--------------
    V=np.float(Volumen.get());
    co=np.float(O2_inicial.get());
    
    #----------Tiempo Total Simulado -------------
    tf=np.float(tiempo_tot.get())*3600;
    t=np.linspace(0,tf,1001)
    #----------- Condicion inicial------------------
    uo=ic(L,OS_0,co)
    
    #------------Resolucion del PDE-----------------
    t,c_ext,u=frac_met_paso_varible(uo,co,tf)
    
    #--------- Calculo de masa de O_2 absorbida----
    O2_absorbido=masa_absorbida(u,t);
    t_ss=tiempo_lag(O2_absorbido,t);
  
    #----------Actualizacion de Interfaz----------
    variables_respuesta(O2_absorbido,t_ss)
    graficador(c_ext,u,O2_absorbido,t_ss,t)
    
    #----------Almacenamiento de Resultados-------
    almacenamiento_resultados(u,c_ext,O2_absorbido,t)
    almacenamiento_variables_desempeno(O2_absorbido,t_ss)
    almacenamiento_parametros(L,OS_0,co,t, lista_variables_material)

#-------------------------------------------------------------------------------
    
    
def imex(x,T,t):
    global co, s, D, L
    R=8.314;
    dx=x[1];
    dt=t[1];
    n=len(x);
    u=np.zeros((5*len(x), len(t)));
    D=Difusividad(T,x,L);
    A=A_izq(dx,dt,len(x),D);
    pseudo_b=A_der(dx,dt,len(x),D);
    c_ext=np.zeros((len(t)));
    c_ext[0]=co;
    u[:,0]=ic(x);    
 
    top = Toplevel();
    top.title("Progreso")
    barra= ttk.Progressbar(top, orient = HORIZONTAL, 
              length = 100,mode = 'determinate');
    barra.pack(pady = 10) 
    barra.update()
    
    u_actual=u[:,0];
    c_ext[1]=c_ext[0]+(dt)*(f(u_actual[4*n:])) 
    A[4*n,:]=0;
    A[5*n-1,:]=0;
    A[4*n,4*n]=1;
    A[5*n-1,5*n-1]=1;
    b=pseudo_b.dot(u_actual)+(dt)*(Reacc(u_actual));
    b[4*n]=c_ext[1]*s[0]*1e6*R*T;
    b[5*n-1]=c_ext[1]*s[-1]*1e6*R*T;
    u_nuevo=spsolve(A,b);
    u[:,1]=u_nuevo.copy(); 
    
    barra['value'] =1;
    barra.update()

#-----------------Menu de Ayuda---------------------
    
def crearVentanaAyuda():
    global tree, texto_label
    newWindow = Toplevel(raiz)
    newWindow.iconbitmap("logo.ico")
    tree=ttk.Treeview(newWindow)
    Intro=tree.insert("", END, text="Introduccion")
    Tutorial=tree.insert("", END, text="Tutorial")
    tree.heading("#0", text="Categoria")
    tree.insert(Tutorial, END, text="Modo Desempeno")
    tree.bind("<Double-1>",item_seleccionado)
    tree.pack(side=LEFT)
    texto_Frame=Frame(newWindow);
    texto_Frame.config(bd=2,relief="groove")
    texto_Frame.pack(side=RIGHT, fill=Y)
    
     
    barra=Scrollbar(texto_Frame);
    posicion_instr="""HAMLET: To be, or not to be--that is the question:
Whether 'tis nobler in the mind to suffer
The slings and arrows of outrageous fortune
Or to take arms against a sea of troubles
And by opposing end them. To die, to sleep--
No more--and by a sleep to say we end
The heartache, and the thousand natural shocks
That flesh is heir to. 'Tis a consummation
Devoutly to be wished."""
    barra=Scrollbar(texto_Frame);
    texto_label=Text(texto_Frame, height=4, width=50);
    barra.pack(side=RIGHT, fill=Y)
    texto_label.pack(side=LEFT, fill=Y)
    barra.config(command=texto_label.yview);
    texto_label.config(yscrollcommand=barra.set)
    texto_label.insert("insert",posicion_instr)
    texto_label.config(state=DISABLED,font=("Helvetica", 10, "italic"))
    
    

def item_seleccionado(event=None):
    global tree, texto_label
    seleccionado=tree.item(tree.selection());
    valor=seleccionado['text'];
    
    texto_label.config(state=NORMAL)
    texto_label.delete(1.0, END)
    if valor=="Tutorial":
        texto_label.insert("insert",valor);
        
    texto_label.config(state=DISABLED)      
    
    

    


raiz.mainloop()
