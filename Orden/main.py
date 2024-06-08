import pfsim
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import *

sim = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#sim = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')



ptdf_dataframe, indices_bus = sim.export_csv()
#ptdf_dataframe = pd.read_csv(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.csv',skiprows = 0,delimiter=';')
#ptdf_dataframe = pd.read_csv(r'C:\Users\lldie\Desktop\SF_datos.csv',skiprows = 0,delimiter=';')
shift_factors, indices_obj = pfsim.ShiftFactors(ptdf_dataframe)

(dict_barras, dict_cargas, dict_lineas, dict_trafos, dict_gen, dict_genstat) = sim.get_data(indices_bus)

# demanda por barra
dda_barra = np.zeros((len(indices_bus)))
for load in dict_cargas:
    dda_barra[dict_cargas[load][0]] += dict_cargas[load][1]

Flujo_dda = shift_factors @ dda_barra   # Flujo en elemntos debido a la demanda en barra




# GENERADORES EFECTIVOS - EN SERVICIO
dict_gen_eff = dict()
for gen in dict_gen:
    if dict_gen[gen][2] == 0:
        dict_gen_eff[gen] = dict_gen[gen]
ngen_eff=len(dict_gen_eff)

# GENERADORES ESTATICOS EFECTIVOS - EN SERVICIO
dict_genstat_eff = dict()
for gen in dict_genstat:
    if dict_genstat[gen][2] == 0:
        dict_genstat_eff[gen] = dict_genstat[gen]
ngenstat_eff=len(dict_genstat_eff)


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
    FMax[cont] = Fmax_i
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
model = gp.Model('Modelo 1')
model.Params.MIPGap = 1e-5
model.Params.OutputFlag = 1 # eliminar mensajes adicioneales Gurobi

p_g = model.addMVar(ngen_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg')
p_sg = model.addMVar(ngenstat_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Psg')

f = model.addVars(n_elem,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas


#indices = pd.read_csv(r'C:\Users\lldie\Desktop\SF_Indices.csv',skiprows = 0,delimiter=';', usecols=[2], encoding= 'ISO-8859-1')

#array_indices = indices.to_numpy()[:,0]
#array_indices_corregidos = list()
#for i in array_indices:
#    if '\\' in i:
#        array_indices_corregidos.append(i.split('\\')[1])
#    else:
#        array_indices_corregidos.append(i)
#