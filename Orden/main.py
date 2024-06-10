import pfsim
import pandas as pd
import numpy as np
import time
import gurobipy as gp
from gurobipy import *

sim = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#sim = pfsim.PowerFactorySim('Ejemplo Clase')
#sim = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')



ptdf_dataframe, indices_bus = sim.export_csv('AC')
#ptdf_dataframe = pd.read_csv(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.csv',skiprows = 0,delimiter=';')
#ptdf_dataframe = pd.read_csv(r'C:\Users\lldie\Desktop\SF_datos.csv',skiprows = 0,delimiter=';')
SF, indices_obj = pfsim.ShiftFactors(ptdf_dataframe)

(dict_barras, dict_cargas, dict_lineas, dict_trafos, dict_gen, dict_genstat) = sim.get_data(indices_bus)



# demanda por barra
dda_barra = np.zeros((len(indices_bus)))
for load in dict_cargas:
    dda_barra[dict_cargas[load][0]] += dict_cargas[load][1]/sim.Sb
Flujo_dda = SF @ dda_barra   # Flujo en elemntos debido a la demanda en barra




# GENERADORES EFECTIVOS - EN SERVICIO
dict_gen_eff = dict()
for gen in dict_gen:
    if dict_gen[gen][2] == 0:
        dict_gen_eff[gen] = dict_gen[gen]
ngen_eff=len(dict_gen_eff)

pos_gen = np.zeros(ngen_eff).astype(int)
ngen_par = np.zeros(ngen_eff)
Cvar_gen = np.zeros(ngen_eff)
Ccte_gen = np.zeros(ngen_eff)
Pmin_gen = np.zeros(ngen_eff)
Pmax_gen =np.zeros(ngen_eff)
cont=0
for gen in dict_gen_eff:
    pos_gen[cont] = dict_gen_eff[gen][0]
    ngen_par[cont] = dict_gen_eff[gen][1]
    Pmin_gen[cont] = dict_gen_eff[gen][3]/sim.Sb
    Pmax_gen[cont] = dict_gen_eff[gen][4]/sim.Sb
    Cvar_gen[cont] = dict_gen_eff[gen][6]
    Ccte_gen[cont] = dict_gen_eff[gen][7]
    cont+=1


# GENERADORES ESTATICOS EFECTIVOS - EN SERVICIO
dict_genstat_eff = dict()
for gen in dict_genstat:
    if dict_genstat[gen][2] == 0:
        dict_genstat_eff[gen] = dict_genstat[gen]
ngenstat_eff=len(dict_genstat_eff)

pos_genstat = np.zeros(ngenstat_eff).astype(int)
ngenstat_par =np.zeros(ngenstat_eff)
Pmin_genstat = np.zeros(ngenstat_eff)
Pmax_genstat =np.zeros(ngenstat_eff)
Cvar_genstat = np.zeros(ngenstat_eff)
Ccte_genstat = np.zeros(ngenstat_eff)
cont=0
for gen in dict_genstat_eff:
    pos_genstat[cont] = dict_genstat_eff[gen][0]
    ngenstat_par[cont] = dict_genstat_eff[gen][1]
    Pmin_genstat[cont] = dict_genstat_eff[gen][3]/sim.Sb
    Pmax_genstat[cont] = dict_genstat_eff[gen][4]/sim.Sb
    Cvar_genstat[cont] = dict_genstat_eff[gen][6]
    Ccte_genstat[cont] = dict_genstat_eff[gen][7]
    cont+=1

n_elem = len(indices_obj)   # N° de líneas + 'trf2'
# Obtención de fmax, R y X
dict_full = dict()
FMax = np.zeros(n_elem)
R = np.zeros(n_elem)
X = np.zeros(n_elem)
cont = -1
for i in indices_obj:
    cont += 1
    if i in dict_lineas:
        dict_full[i] = dict_lineas[i]
        Fmax_i = dict_lineas[i][2]
        R_i = dict_lineas[i][0]
        X_i = dict_lineas[i][1]
    elif i in dict_trafos:
        dict_full[i] = dict_trafos[i]
        Fmax_i = dict_trafos[i][3]
        R_i = dict_trafos[i][1]
        X_i = dict_trafos[i][0]
    else:
        print('Elemento no encontrado en diccs: ' + i)
    FMax[cont] = Fmax_i/sim.Sb
    R[cont] = R_i
    X[cont] = X_i


val_min = 1e-24 #Valor mínimo para evitar valores nulos en matriz X
yprim = np.zeros(n_elem).astype(complex)
for i in range(n_elem):
    if X[i] == 0:
        X[i] = val_min
    yprim[i] = 1/(complex(R[i],X[i]))
G = np.real(yprim)
B = np.imag(yprim)

# Definir Modelo
m = gp.Model('Modelo 1')
m.Params.MIPGap = 1e-5
m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

p_g = m.addMVar(ngen_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg')
#pdis_g = m.addMVar(ngen_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_dis')
p_statg = m.addMVar(ngenstat_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Psg')
#pdis_statg = m.addMVar(ngenstat_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Psg_dis')

f = m.addMVar(n_elem,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas

# Función objetivo
f_obj = 0
#costo_gen = 0
#costo_genstat = 0

costo_gen = p_g* ngen_par @ Cvar_gen*sim.Sb + Ccte_gen.sum()
costo_genstat = p_statg* ngenstat_par @ Cvar_genstat*sim.Sb + Ccte_genstat.sum()

f_obj = costo_gen + costo_genstat

m.setObjective(f_obj, GRB.MINIMIZE)
m.getObjective()

###  Restricciones
## Balance Nodal
m.addConstr(p_g @ ngen_par + p_statg @ ngenstat_par == dda_barra.sum(), name='Balance')

## Restricciones de generación
# Limitaciones de Potencia
m.addConstr(p_g >= Pmin_gen, name='Pmin_gen')
m.addConstr(-p_g >= -Pmax_gen, name='Pmax_gen')
m.addConstr(p_statg >= Pmin_genstat, name='Pmin_genstat')
m.addConstr(-p_statg >= -Pmax_genstat, name='Pmin_genstat')

## Sistema de transmisión
m.addConstr(-f == SF[:,pos_gen] @ (p_g*ngen_par) + SF[:,pos_genstat] @ (p_statg*ngenstat_par)- Flujo_dda, name='f')
m.addConstr(-f >= -FMax, name = 'fmax+')
m.addConstr(f >= -FMax, name = 'fmax-')

# %%
m.write('OPF.lp')

t2 = time.time()
m.optimize()
t3=time.time()

print('-----------------------------')
print('La demanda total del sistema es: %.2f (MW)' % (dda_barra.sum()*sim.Sb))    
status = m.Status
if status == GRB.Status.OPTIMAL:
    print ('Costo => %.2f ($/h)' % m.objVal) 
    #print ('Las perdidas son %.2f (MW)' % sum(pk_loss[i].X for i in range(len(pk_loss)))) 
elif status == GRB.Status.INF_OR_UNBD or \
    status == GRB.Status.INFEASIBLE  or \
    status == GRB.Status.UNBOUNDED:
    print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
    m.computeIIS() 
    m.write("GTCEP.ilp")

cont=0
for gen in dict_gen_eff:
    print('%s => %.2f (MW)' % (gen,ngen_par[cont]*p_g[cont].X*sim.Sb))
    cont += 1

cont=0
for genstat in dict_genstat_eff:
    print('%s => %.2f (MW)' % (genstat,ngenstat_par[cont]*p_statg[cont].X*sim.Sb))
    cont += 1

cont=0
for elem in indices_obj:
    print('Línea: %s => %.2f (MW)' % (elem,f[cont].X*sim.Sb))
    cont += 1

print('finish')
#  %% Comprobación de resultados



#indices = pd.read_csv(r'C:\Users\lldie\Desktop\SF_Indices.csv',skiprows = 0,delimiter=';', usecols=[2], encoding= 'ISO-8859-1')

#array_indices = indices.to_numpy()[:,0]
#array_indices_corregidos = list()
#for i in array_indices:
#    if '\\' in i:
#        array_indices_corregidos.append(i.split('\\')[1])
#    else:
#        array_indices_corregidos.append(i)
#