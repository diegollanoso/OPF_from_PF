# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 00:20:53 2022

@author: danie
"""
# %% Preambulo
# Importar Gurobi
import gurobipy as gp  
from gurobipy import GRB,quicksum
import pandas as pd
import numpy as np
from numpy import array, dot as mult, zeros , imag , real

# Importar tiempo
import time
t0 = time.time()

# Importar DIgSILENT
import sys
sys.path.append(r'C:\Program Files\DIgSILENT\PowerFactory 2023 SP3A\Python\3.11')
import powerfactory

# Activar proyecto
app = powerfactory.GetApplication()
PROJECT = 'Ejemplo Clase'
app.ActivateProject(PROJECT)
app.Show() 

app.ResetCalculation() 

# Comandos
ldf = app.GetFromStudyCase('ComLdf')

factores = app.GetFromStudyCase('ComVstab') # SF


# %% Datos del sistema 

#### Listas de elementos de red
generadores = app.GetCalcRelevantObjects('*.ElmSym') 
cargas = app.GetCalcRelevantObjects('*.ElmLod')
lineas = app.GetCalcRelevantObjects('*.ElmLne')
barras = app.GetCalcRelevantObjects('*.ElmTerm')

#### Cantidad de elementos 
ng = len(generadores) # N° de generadores
nc = len(cargas) # N° de cargas
nl = len(lineas) # N° de lineas
nb = len(barras) # N° de barras

# activar generadores
for g in generadores:
    g.outserv = 0

#### Demanda del sistema
demanda_total = 0
for c in cargas:
    demanda_total += c.plini

fact_perd = 1.00725
#### Parámetros de generadores
nom_gen = list()
p_min = list()
p_max = list()
Cv = list()
Cf = list()

# valores base 
Sb = 500
Vb = 110
raiz3 = 1.73205080757
Ib = Sb/(raiz3*Vb)
Zb = Vb*Vb/Sb

for g in generadores:
    nom_gen.append(g.loc_name)
    p_min.append(g.Pmin_uc)
    p_max.append(g.Pmax_uc)
    Cv.append(g.penaltyCosts)
    Cf.append(g.fixedCosts)
    

    
    
    
# %% Obtención de Shift-Factors mediante DIgSILENT
factores.iopt_method = 2 # En DC
factores.pResult.loc_name = "shift_factors"
factores.frmLimitsBrc = 1
factores.cdpflim = 0 
factores.Execute()


archivo_resultados = app.GetFromStudyCase('shift_factors.ElmRes') 
exportar = app.GetFromStudyCase('ComRes')
exportar.pResult = archivo_resultados
exportar.iopt_vars = 0
exportar.iopt_exp = 6 #tipo de archivo a exportar
exportar.f_name = r'C:\Users\danie\OneDrive\Escritorio\Clase OPF\SF_datos.csv'
exportar.iopt_csel = 0 
exportar.locn = 1
exportar.iopt_sep=0
exportar.ciopt_head = 0
exportar.dec_Sep = '.'
exportar.Execute()

# Extraer Datos 
ptdf_dataframe = pd.read_csv('SF_datos.csv',skiprows = 0,delimiter=';')

ejes = ptdf_dataframe.axes
eliminar_col = list([0])
cont = 0 
for i in ejes[1]:
    if '.' in i:        
        eliminar_col.append(cont)
    cont += 1
    
ptdf_array_T = np.asarray(ptdf_dataframe)
ptdf_array_T = np.delete(ptdf_array_T ,eliminar_col, axis = 1)

for i in range(len(ptdf_array_T)):
    for u in range(len(ptdf_array_T[i])):
        if ptdf_array_T[i,u] == '   ----':
            ptdf_array_T[i,u] = 0
            
sf = np.transpose(ptdf_array_T).astype(float) 

exportar.f_name = r'C:\Users\danie\OneDrive\Escritorio\Clase OPF\indices.csv'
exportar.iopt_vars = 2
exportar.Execute()

indices = pd.read_csv('indices.csv',skiprows = 0,delimiter=';', usecols=[2])

array_indices= indices.to_numpy()[:,0]
                       
ind_lin = array_indices[0:nl]
ind_bus = array_indices[nl:nl+nb]

print("Shift Factors")
print(sf)


# %% Listas de datos del sistema
dict_barras = {}
cont = 0 
for b in ind_bus:
    dict_barras[str(b)]= cont
    cont += 1
    
dict_gen = {}
cont = 0 
for g in generadores:
    dict_gen[str(g.loc_name)]= cont
    cont += 1    
    
dict_lineas = {}
cont = 0 
for l in ind_lin:
    dict_lineas[str(l)]= cont
    cont += 1  

ldf.Execute() # Ejecutar flujo de potencia

# Lista de cargas
list_cargas = np.zeros((nc, 2))
cont = 0
for c in cargas:
    list_cargas[cont] = [dict_barras[c.bus1.cterm.loc_name], c.GetAttribute('m:Psum:bus1')]
    cont += 1

# Lista de lineas
# "from","to","R", "X", "Lmin", "Lmáx"
list_lineas = np.zeros((nl, 7))
cont1 = 0
for line in lineas:
    list_lineas[cont1] = [dict_lineas[line.loc_name],dict_barras[line.bus1.cterm.loc_name], dict_barras[line.bus2.cterm.loc_name], round(line.R1/Zb,6), round(line.X1/Zb,6),round(-line.Inom_a*Vb*raiz3,4),round(line.Inom_a*Vb*raiz3,4)]
    cont1 += 1

list_lineas = np.array(sorted(list_lineas, key=lambda x:x[0]))
fm = list_lineas[:,6]  # flujo máximo

# Lista generadores
list_gen = np.zeros((ng, 6))
pos_g = []
cont = 0
for g in generadores:
    pos_g.append(dict_barras[g.bus1.cterm.loc_name])
    potencia = 0 
    if g.outserv == 0:
        potencia = round(g.GetAttribute('m:Psum:bus1'),3)
    list_gen[cont] = [dict_gen[g.loc_name],dict_barras[g.bus1.cterm.loc_name], g.ip_ctrl ,potencia,g.Pmin_uc ,g.P_max]
    cont += 1
    

# #### Lista de barras
# # N° de barra, balance nodal
# list_barras = np.zeros((nb, 2))
# cont = 0
# for b in ind_bus:
#     Pd = 0
#     Pg = 0
#     for g in range(ng):
#         if list_gen[g,1] == dict_barras[b]:  
#             Pg += float(list_gen[g,3])
#     for c in range(nc):
#         if list_cargas[c,0] == dict_barras[b]:           
#             Pd += float(list_cargas[c,1])    
#     list_barras[cont] = [dict_barras[b], round(Pd - Pg,2)]
#     cont += 1
    
demanda_por_barra = np.zeros((nb, 2))
cont = 0
for b in ind_bus:
    Pd = 0
    for c in range(nc):
        if list_cargas[c,0] == dict_barras[b]:           
            Pd += float(list_cargas[c,1])    
    demanda_por_barra[cont] = [dict_barras[b], round(Pd,2)]
    cont += 1

    
total_G0 = sum(list_gen[:,3])
total_L0 = sum(list_cargas[:,1])
perd_0 = total_G0 - total_L0

PF = mult(sf, demanda_por_barra[:,1])

R = list_lineas[:,3]
X = list_lineas[:,4]
yprim = zeros(nl).astype(complex)
    
for l in range(nl):
    yprim[l] = 1/(complex(R[l],X[l]))
    
G = real(yprim)
B = imag(yprim)
# %% Optmización
# Definir modelo
model = gp.Model('Modelo 1')
#model.Params.MIPGap = 1e-6
model.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi
#### Definición de variables
p = model.addVars(ng,vtype=GRB.CONTINUOUS,ub=GRB.INFINITY,lb=0, name='Pg') # Potencia generadores
f = model.addVars(nl,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas

# Variables para pérdidas
L=6
pk_mas = model.addVars(nl, vtype=GRB.CONTINUOUS, lb=0, name='Pk_mas')
pk_menos = model.addVars(nl, vtype=GRB.CONTINUOUS, lb=0, name='Pk_menos')
dpk = model.addVars(nl,L, vtype=GRB.CONTINUOUS, lb=0, name='dpk')
pk_loss = model.addVars(nl, vtype=GRB.CONTINUOUS, lb=0 , name= 'pk_loss')

#### Función objetivo
f_obj = 0   
for i in range(ng):
    f_obj += Cv[i]*p[i] + Cf[i]
    
model.setObjective(f_obj,GRB.MINIMIZE)  

#### Restricciones     
# Flujo máximo por las lineas
for i in range(nl):#
    for u in range(L):#
        model.addConstr( -dpk[i,u] >= -fm[i]/L) #

for i in range(nl):
    kl = zeros(L)
    for l in range(0,L):
        kl[l] = (2*(l+1)-1)*(fm[i]/Sb)/L#
    model.addConstr(pk_loss[i] == (G[i]/(B[i])**2) * quicksum(kl*dpk.select(i,'*')))
    
    f_g = quicksum(p.select('*')*sf[i,pos_g])
    model.addConstr(f[i] == PF[i] - f_g)
    model.addConstr(f[i] == pk_mas[i] - pk_menos[i])
    model.addConstr(sum(dpk.select(i,'*')) == pk_mas[i] + pk_menos[i])
    model.addConstr(float(fm[i]) >= pk_mas[i] + pk_menos[i] + pk_loss[i])
    
    
# Potencia máxima y mínima
for i in range(ng):
    model.addConstr(p[i] >= p_min[i], name = 'Pmin'+ str(i))
    model.addConstr(-p[i] >= -p_max[i], name = 'Pmax' + str(i))
    
# Balance de potencia 
model.addConstr(p.sum('*') == demanda_total + pk_loss.sum('*'), name ='Balance')    

t2 = time.time()
model.optimize()
model.write('ejemplo1.lp')    
t3 = time.time()
  

# salida de información
print('-----------------------------')

print('La demanda total del sistema es: %.2f (MW)' % (demanda_total*fact_perd))    
status = model.Status
if status == GRB.Status.OPTIMAL:
    print ('Costo => %.2f ($/h)' % model.objVal) 
    # for v in model.getConstrs():
    #     if v.pi > 1e-2:
    #         print('%s = %g ($/MWh)' % (v.ConstrName,v.pi))
elif status == GRB.Status.INF_OR_UNBD or \
    status == GRB.Status.INFEASIBLE  or \
    status == GRB.Status.UNBOUNDED:
    print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
    model.computeIIS() 
    model.write("GTCEP.ilp")


cont = 0
for g in generadores:
    g.pgini = p[cont].X
    print('%s => %.2f (MW)' % (g.loc_name,p[cont].X))
    cont += 1

ldf.Execute()

# %% Corroboración de Resultados

flujos_AC = zeros(nl)
for l in lineas:
    flujos_AC[dict_lineas[l.loc_name]] = l.GetAttribute('m:Psum:bus2')
    
for i in range(nl):
    print('%s ; Gurobi (DC) = %.3f [MW]; DIgSILENT (AC) = %.3f [MW]' % (ind_lin[i],f[i].X,flujos_AC[i]))

print('Tiempo de total: %.4f (s)'% (t3-t0))
print('Tiempo de optimización : %.4f (s)' % (t3-t2))
