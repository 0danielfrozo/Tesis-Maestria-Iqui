from matplotlib import animation, rc
from IPython.display import HTML
import matplotlib.pyplot as plt
from fipy import * 
import numpy as np
import matplotlib.animation as animation

#Condiciones calculo camino libre medio
P=101325 #Pa
T=298 #K

#Constantes 
R=8.314 #J/mol K
Na=6.022e23 
#Propiedades oxigeno
d_m=346e-12 #m diametro kinetico O_2
c=0.21*P/(R*T) #mol/m3 concentracion de oxigeno en el aire

#Camino libre medio 
lamnda=1/(np.sqrt(2)*np.pi*d_m**2*c*Na) #m

#Diametro de poro
sigma=20e-10 #m Tomado de Tesis Reynel Gomez (2018)

#numero de Knudsen
K_n=lamnda/sigma

#Dimensiones Capsula
r_ext=7.5e-8 #[m]
r_in=4.182e-8 #[m]

#Vector de radio discretizado
n=200
r=np.linspace(0,r_ext,n)
delta_r=r[1];

#Funcion Escalon 
H=(r_in<r).astype(int);
frontera=np.ravel(np.asarray(H==1).nonzero())[0]


#Propiedades oxigeno
M=32e-3 #Kg/mol

#Propiedades Membrana 
vol_capsula=(4/3)*np.pi*(350e-7)**3 #cm3
peso_capsula=2.5996e-15 #g
vol_total=vol_capsula/peso_capsula
vol_poro=0.002 #cm3/g

E=vol_poro/vol_total #Porosidad 
t_2=E**(1-1.5) #Tortuosidad 

#Calculo difusividad
def dif_knudsen(T=298):
    Ko=1/3 *d_m*E/t_2;
    v=np.sqrt(8*R*T/(np.pi*M))
    D_k=Ko*v
    return D_k

#Difusividad y solubilidad de oxigeno en aceite
def Difusividad(T=298):
    Dac=0.4e-9 #m2/s
    D_k=dif_knudsen(T) #m2/s
    return D_k,Dac

def Solubilidad(T=298):
    R=8.314;
    Sac=6.858e6*np.exp(-1.478e4/(R*T)); #Pa m3/mol
    return 1/(R*T),1/Sac

def cons_cinetica(T=298):
    R=8.314;
    return np.asarray([52691244*np.exp(-101.5e3/(R*T)),6581,572.68*np.exp(-31e3/(R*T)),2.2914e+08,10362358.9,96814085.92*np.exp(-48.4e3/(R*T))])


D_si,D_ac=Difusividad()
S_si,S_ac=Solubilidad()

n=200
grid_x,grid_y = np.meshgrid(np.linspace(-150e-9,150e-9,n),np.linspace(-150e-9,150e-9,n))
D=(grid_x**2+grid_y**2<r_in**2)*D_ac+(grid_x**2+grid_y**2<r_ext**2)*D_si+0.219e-4*(grid_x**2+grid_y**2>r_ext**2)
D=D-(grid_x**2+grid_y**2<r_in**2)*D_si

S=(grid_x**2+grid_y**2<r_in**2)*S_ac+S_si
S=S-(grid_x**2+grid_y**2<r_in**2)*S_si

k=cons_cinetica()
dx=np.linspace(0,2*150e-9,n)[1]
mesh = Grid2D(nx=n, ny=n, dx=dx, dy=dx)

O2 = CellVariable(name=r"$O_2$", mesh=mesh,hasOld=True)
ROOH = CellVariable(name=r"$ROOH$", mesh=mesh,hasOld=True)
RH = CellVariable(name=r"$RH$", mesh=mesh,hasOld=True)
ROO = CellVariable(name=r"$ROO\cdot$", mesh=mesh,hasOld=True)
R_reac = CellVariable(name=r"$R\cdot$", mesh=mesh,hasOld=True)

O2o=0.21*P;
ROOH_o=19.1789709;
RH_o=1559.37052;

O2_inicial=(grid_x**2+grid_y**2>r_ext**2)*O2o
ROOH_inicial=(grid_x**2+grid_y**2<r_in**2)*ROOH_o
RH_inicial=(grid_x**2+grid_y**2<r_in**2)*RH_o

O2.value=np.reshape(O2_inicial,n*n);
ROOH.value=np.reshape(ROOH_inicial,n*n);
RH.value=np.reshape(RH_inicial,n*n);


D1 = CellVariable(mesh=mesh, value=np.reshape(D,n*n))
S1=CellVariable(mesh=mesh, value=np.reshape(S,n*n))

elapsed = 0
fig=plt.figure()
fig.suptitle('t = %.2g' %(elapsed), fontsize=16)
ax1=fig.add_subplot(111)
#viewer_def=MultiViewer(viewers=(MatplotlibViewer(vars=O2,datamin=0,datamax=O2o,axes=fig.add_subplot(121)),(MatplotlibViewer(vars=O2*S1,axes=fig.add_subplot(122),datamin=0,datamax=O2o*S_si,))) )

viewer_def =MultiViewer(viewers=(MatplotlibViewer(vars=O2,datamin=0,datamax=O2o,axes=ax1)))
def inicial():
    viewer_def.plot()
    fig.suptitle('t = %.2g' %(elapsed), fontsize=16)
    ax1.update(props)
    fig.canvas.draw_idle()


dexp = -5
elapsed = 0
props = {'autoscalex_on': False}

reacc_1=k[0]*ROOH
reacc_2=k[1]*O2*S1
reacc_2_alt=k[1]*R_reac
reacc_3=k[2]*RH
reacc_4=k[3]*R_reac
reacc_5=k[4]*R_reac
reacc_6=k[5]*ROO


eqn0= (TransientTerm(var=O2) == DiffusionTerm(coeff=D1,var=O2)-ImplicitSourceTerm(coeff=reacc_2_alt,var=O2))
eqn1= (TransientTerm(var=ROOH) == -2*ImplicitSourceTerm(coeff=reacc_1,var=ROOH)+ImplicitSourceTerm(coeff=reacc_3,var=ROO))
eqn2= (TransientTerm(var=ROO) == ImplicitSourceTerm(coeff=reacc_1,var=ROOH)+ImplicitSourceTerm(coeff=reacc_2, var=R_reac)-ImplicitSourceTerm(coeff=reacc_3,var=ROO)-ImplicitSourceTerm(coeff=reacc_5, var=ROO)-2*ImplicitSourceTerm(coeff=reacc_6, var=ROO))
eqn3= (TransientTerm(var=R_reac) == ImplicitSourceTerm(coeff=reacc_1,var=ROOH)-ImplicitSourceTerm(coeff=reacc_2, var=R_reac)+ImplicitSourceTerm(coeff=reacc_3,var=ROO)-2*ImplicitSourceTerm(coeff=reacc_4, var=R_reac)-ImplicitSourceTerm(coeff=reacc_5, var=ROO))
eqn4= (TransientTerm(var=RH) == -ImplicitSourceTerm(coeff=reacc_3,var=ROO))

eq = eqn0 & eqn1 & eqn2 &  eqn3 &  eqn4

def update_plot(parametro):
    global elapsed, O2,ROOH,ROO,R_reac,RH,eq,viewer_def,fig,dexp,ax1

    if elapsed < 1.5:
        dt = min(100, np.exp(dexp))
        elapsed += dt
        dexp += 0.01
   
        elapsed += dt
    
        O2.updateOld()
        ROOH.updateOld()
        ROO.updateOld()
        R_reac.updateOld()
        RH.updateOld()
    
        eq.solve(dt=dt)
        viewer_def.plot()
        fig.suptitle('t = %.2g' %(elapsed), fontsize=16)
        ax1.update(props)
        fig.canvas.draw_idle()

fps=5;
frn=50;
ani = animation.FuncAnimation(fig, update_plot,init_func=inicial)
fn = 'plot_surface_animation_funcanimation_3'
Writer = animation.writers['pillow']
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
ani.save(fn+'.gif',writer=writer)
