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


raiz=Tk()

#------------Vetana Principal --------------
raiz.title("Herramienta de Diseno V2.1")
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
    global longitud, material,inicial_OS_n
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
    
    material=[]
    mat=StringVar();
    combo=ttk.Combobox(Frame_Homogeneo,width=5,textvariable=mat,state="readonly",values=materiales) 
    combo.grid(row=2,column=1)
    combo.config(justify='center')
    material.append(mat);
    
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
    global longitud, material 
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
    material=[]
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
        material.append(mat)
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

def correr():
 global nu,D,k,A,V,s,co,n,x,tf,temp,T,inicial_OS_n,reactiva,longitud,material,OS_0,L
 T=np.float(temp.get());
 L=[];
 L.append(np.float(longitud[0].get()));
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
 n=101;
 x=np.linspace(0,L[-1],n);
 D=Difusividad(T,x,L,material);
 s=Solubilidad(T,x,L,material);

 A=np.float(Area.get());
 V=np.float(Volumen.get());

 co=np.float(O2_inicial.get());

 to=0
 tf=np.float(tiempo_tot.get())*3600;
 pde(x,to,tf)



def Solubilidad(T,x,L,mat):    
    R=8.314;
    actual=materiales.index(mat[0].get()) 
    Delta_H=np.float(data.iloc[actual,4]);
    So=np.float(data.iloc[actual,3]);
    S=So*np.exp((-Delta_H/R)*(1.0/T -1.0/298.0))
    resp=S*np.ones((len(x)));    
    for i in range(len(L)-1):
        actual=materiales.index(mat[i+1].get()) 
        Delta_H=np.float(data.iloc[actual,4]);
        So=np.float(data.iloc[actual,3]);
        S=So*np.exp((-Delta_H/R)*(1.0/T -1.0/298))
        for j in range(len(x)):
            if x[j]>=L[i]:
                resp[j]=S;   
    return resp;

def Difusividad(T,x,L,mat):
    global materiales 
    R=8.314;
    actual=materiales.index(mat[0].get())                         
    Ea=np.float(data.iloc[actual,2]);
    Do=np.float(data.iloc[actual,1]);
    D=Do*np.exp(-Ea/(R*T))
    resp=D*np.ones((len(x)));    
    for i in range(len(L)-1):
        actual=materiales.index(mat[i+1].get())
        Ea=np.float(data.iloc[actual,2]);
        Do=np.float(data.iloc[actual,1]);
        D=Do*np.exp(-Ea/(R*T));
        for j in range(len(x)):
            if x[j]>=L[i]:
                resp[j]=D;   
    return resp;
    
def f(t,u):
    #Parametros del sistema a estudiar;
    global D,L,L1,L2,nu,k,A,V,no,s,co,n,x,T
    R=8.314;
    k=np.asarray([4.09e10*np.exp(-101.5e3/(R*T)),1e7,3.52e6*np.exp(-31e3/(R*T)),3e11,1e10,3.05e15*np.exp(-48.4e3/(R*T))])
    dx=x[1]-x[0]
    #Calculo de la derivada temporal 
    especies=np.reshape(u[:-1],(5,len(x)));
    du_dt=np.zeros(np.shape(especies));  
    du_dt[0,:]=-2*k[0]*especies[0,:]**2+k[2]*especies[1,:]*especies[3,:];
    du_dt[1,:]=k[0]*especies[0,:]**2+k[1]*especies[2,:]*especies[4,:]-k[2]*especies[1,:]*especies[3,:]-2*k[5]*especies[1,:]**2-k[4]*especies[2,:]*especies[1,:];
    du_dt[2,:]=k[0]*especies[0,:]**2-k[1]*especies[2,:]*especies[4,:]+k[2]*especies[1,:]*especies[3,:]-2*k[3]*especies[2,:]**2-k[4]*especies[2,:]*especies[1,:];
    du_dt[3,:]=-k[2]*especies[1,:]*especies[3,:];
    du_dt[4,1:-1]=D[1:-1]*((-2*especies[4,1:-1]+especies[4,2:]+especies[4,:-2])/(dx**2))-k[1]*especies[2,1:-1]*especies[4,1:-1];
    dc_ext_dt=(A/V)*(D[0]*(-especies[4,2]+4*especies[4,1]-3*especies[4,0])/(2.0*dx) +D[-1]*(-especies[4,-3]+4*especies[4,-2]-3*especies[4,-1])/(2.0*dx)) ;
    du_dt[4,0]=dc_ext_dt*s[0]*1e6*R*T;
    du_dt[4,-1]=dc_ext_dt*s[-1]*1e6*R*T;
    
    du_dt=np.reshape(du_dt,(1,-1))
    du_dt=np.append(du_dt,dc_ext_dt);
    return du_dt

def ic(x):
    global OS_0, s, L, T,co
    c_o=[0.05,0,0,6.8];
    c_o=np.asarray(c_o);
    uo=np.ones((5,len(x)));
    for o  in  range(4): uo[o,:]=OS_0[0]*c_o[o]*uo[o,:];   
    R=8.314;
    for j in range(len(L)-1):
        for i in range(len(x)):
            if x[i]>=L[j]: 
                for o  in  range(4): uo[o,i]=OS_0[j+1]*c_o[o]*uo[o,i]; 
    uo[4,:]=0;
    uo[4,0]=co*s[0]*1e6*R*T;
    uo[4,-1]=co*s[-1]*1e6*R*T;
    
   
    uo=np.reshape(uo,5*len(x));
    uo=np.append(uo,co);
    return uo 

  
def pde(x,to,tf):
    global D, A,O2_result
    global t, O2, OS, HS , masa 
    dx=x[1];
    uo=ic(x);
    
    u_3=solve_ivp(f,[to,tf],uo,method='BDF');
    t=u_3.t;
    u=u_3.y.T;
    HS=u[:,-1]; 
    O2=u[:,4*n:5*n];
    OS=u[:,:n];
#    u=np.reshape(u[:,:-1],(5,len(t),len(x)));
     
    
    
    Flux=-D[0]*(-O2[:,2]+4*O2[:,1]-3*O2[:,0])/(2.0*dx) -D[-1]*(-O2[:,-3]+4*O2[:,-2]-3*O2[:,-1])/(2.0*dx)
    
    masa=[];
    rep=[];
    masa_tot=64*A*np.trapz(Flux, t)
    for i in range(1,len(t)+1): 
        masa.append(64*A*np.trapz(Flux[:i], t[:i]))
    for i in range(len(t)):   
        if np.str(format(masa[i], "10.2E"))==np.str(format(masa_tot, "10.2E")):rep.append(i); 
    
    O2_result.set(format(masa[-1], "10.2E"))
    tiempo_SS=masa.index(masa[-1]);
    tiempo_result.set('-');
    if len(rep)>1:tiempo_result.set(format(t[rep[0]]/3600.0, "10.2f"))
     
    
    
    
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
    a3.plot(t,HS,'-')
    a3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    a3.set_xlabel('t (s)')
    a3.set_ylabel('$O_2$ Head Space [mol/$cm^3$]')
    canvas3.draw()
    
    
    a4.clear()
    a4.plot(t,masa,'-')
    a4.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    a4.set_xlabel('t (s)')
    a4.set_ylabel('$O_2$ Total Mass [g]')
    canvas4.draw()

def guardar():
    global x, t, O2, OS, HS
    global nu,D,k,A,V,s,co,n,x,tf,temp,T,inicial_OS_n,reactiva,longitud,material,OS_0,L
    files = [('Excel', '*.xlsx')] 
    file = filedialog.asksaveasfilename(filetypes = files, defaultextension = files)
    tabla =np.zeros((len(t)+1,len(x)+1));
    tabla[1:,0]=t;
    tabla[0,1:]=x;
    tabla[1:,1:]=O2;
    tabla=tabla.tolist();
    tabla[0][0]='Tiempo [s]\ Posicion [m]'
    
    
    tabla_2 =np.zeros((len(t)+1,len(x)+1));
    tabla_2[1:,0]=t;
    tabla_2[0,1:]=x;
    tabla_2[1:,1:]=OS;
    tabla_2=tabla_2.tolist();
    tabla_2[0][0]='Tiempo [s]\ Posicion [m]'
    
    
    writer = pd.ExcelWriter(file, engine='xlsxwriter')
    df=pd.DataFrame(tabla)
    df_2=pd.DataFrame(tabla_2)
    df_3=pd.DataFrame({'Tiempo [s]':t, 'O2 HeadSpace [mol/cm3]':HS})
    df_4=pd.DataFrame({'Tiempo [s]':t,'Masa O2 absorbida [g]':masa })
    df.to_excel(writer, sheet_name='Concentracion O2', index=False)
    df_2.to_excel(writer, sheet_name='Concentracion OS', index=False)
    df_3.to_excel(writer, sheet_name='Headspace O2 Concentration', index=False)
    df_4.to_excel(writer, sheet_name='Masa O2 absorbido', index=False)
    writer.save()



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
