import gurobipy as gb
from gurobipy import *
import time
import numpy as np
from numpy import array, dot as mult, zeros , imag , real
import pandas as pd

# Importar DIgSILENT
import sys
sys.path.append(r'C:\Program Files\DIgSILENT\PowerFactory 2023 SP3A\Python\3.11')
import powerfactory

# valores base 
Sb = 500
Vb = 110
raiz3 = 1.73205080757
Ib = Sb/(raiz3*Vb)
Zb = Vb*Vb/Sb


t0 = time.time() #formulation time
# Activar proyecto
app = powerfactory.GetApplication()
app.Show() 
PROJECT = 'IEEE 14-bus con PE y PV'
#PROJECT = 'Ejemplo Clase'
app.ActivateProject(PROJECT)

app.ResetCalculation() 

# Comando flujo de carga
ldf = app.GetFromStudyCase('ComLdf')
# Factores Shift-factors
sf = app.GetFromStudyCase('ComVstab') 


# %% Datos del sistema 

#### Listas de elementos de red
generadores = app.GetCalcRelevantObjects('*.ElmSym') 
PV_plants = app.GetCalcRelevantObjects('*.ElmPvsys')
PE_plants = app.GetCalcRelevantObjects('*.ElmGenstat')
cargas = app.GetCalcRelevantObjects('*.ElmLod')
lineas = app.GetCalcRelevantObjects('*.ElmLne')
barras = app.GetCalcRelevantObjects('*.ElmTerm')

#### Cantidad de elementos 
ng = len(generadores) # N° de generadores
n_pe = len(PE_plants) # N° de generadores
n_pv = len(PV_plants) # N° de generadores
nc = len(cargas) # N° de cargas
nl = len(lineas) # N° de lineas
nb = len(barras) # N° de barras

# Activar genradores
for g in generadores:
    g.outserv = 0

#### Demanda del sistema
demanda_total = 0
for c in cargas:
    demanda_total += c.plini

#### Parámetros de generadores
nom_gen = list()
p_min = list()
p_max = list()
Cv = list()
Cf = list()

pv_nom_gen = list()
pv_p_min = list()
pv_p_max = list()
pv_Cv = list()
pv_Cf = list()

pe_nom_gen = list()
pe_p_min = list()
pe_p_max = list()
pe_Cv = list()
pe_Cf = list()
pe_parallel = list()

for g in generadores:
    nom_gen.append(g.loc_name)
    p_min.append(g.Pmin_uc)
    p_max.append(g.Pmax_uc)
    Cv.append(g.penaltyCosts)
    Cf.append(g.fixedCosts)

for pv in PV_plants:
    pv_nom_gen.append(pv.loc_name)
    pv_p_min.append(pv.Pmin_uc)
    pv_p_max.append(pv.Pmax_uc)
    pv_Cv.append(pv.penaltyCosts)
    pv_Cf.append(pv.fixedCosts)

for pe in PE_plants:
    pe_nom_gen.append(pe.loc_name)
    pe_p_min.append(pe.Pmin_uc)
    pe_p_max.append(pe.Pmax_uc)
    pe_Cv.append(pe.penaltyCosts)
    pe_Cf.append(pe.fixedCosts)
    pe_parallel.append(pe.ngnum)


# %% Shift-Factors
#proj = app.GetProjectFolder('study')
#studycase = proj.GetContents('*.IntCase')[0]

#if len(studycase.GetContents('shift_factors.ElmRes')) == 0:
#    Results = studycase.CreateObject('ElmRes','shift_factors')
#    Results.calTp = 31
#    Line_Results = Results.CreateObject('IntMon','Line')
#    Line_Results.classnm = 'ElmLne'
#    Line_Results.AddVar('dPdP:bus1')
#    Line_Results.AddVar('dPdP:bus2')

    #Results.Delete()



sf.iopt_method = 0 # En AC
sf.pResult.loc_name = 'shift_factors'
Results = sf.pResult
for i in Results.GetContents():
    i.Delete()
Line_Results = Results.CreateObject('IntMon','Line')
Line_Results.classnm = 'ElmLne'
Line_Results.AddVar('m:dPdP:bus1')
Line_Results.AddVar('m:dPdP:bus2')
sf.frmElmFilt4Res = 0
sf.frmLimitsBrc = 1
sf.cdpflim = 0 
sf.frmLimitsBus = 0
#print(sf.pResult)
sf.Execute()



Results = app.GetFromStudyCase('shift_factors.ElmRes') # Archivo Resultados
#print(Results.GetContents())
export = app.GetFromStudyCase('ComRes') # Método para exportar 
export.pResult = Results
export.iopt_vars = 0
export.iopt_exp = 6 #tipo de archivo a exportar
export.f_name = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.csv'
export.iopt_csel = 0 
export.locn = 1
export.iopt_sep=0
export.ciopt_head = 0
export.dec_Sep = '.'
export.Execute()

ptdf_dataframe = pd.read_csv('SF_datos.csv',skiprows = 0,delimiter=';')
ejes = ptdf_dataframe.axes

list_bus = list()
for i in sf.p_bus.All():
    list_bus.append(i.loc_name)

list_lin = list()
cont=0
for i in ejes[1]:
    if not '.' in i and cont>2:        
        list_lin.append(i)
    cont+=1


ptdf_array_T = np.asarray(ptdf_dataframe)
ptdf_array_T = np.delete(ptdf_array_T ,[0,1,2], axis = 1)

sf = np.transpose(ptdf_array_T).astype(float) 

#print("Shift Factors")
#print(sf)

# %% asdf

dict_barras = {}
cont = 0 
for b in list_bus:
    dict_barras[str(b)]= cont
    cont += 1
    
dict_gen = {}
cont = 0 
for g in generadores:
    dict_gen[str(g.loc_name)]= cont
    cont += 1    
    
dict_lineas = {}
cont = 0 
for l in list_lin:
    dict_lineas[str(l)]= cont
    cont += 1  

#ldf.Execute() # Ejecutar flujo de potencia

# Lista de cargas
list_cargas = zeros((nc, 2))
cont = 0
for c in cargas:
    #list_cargas[cont] = [dict_barras[c.bus1.cterm.loc_name], c.GetAttribute('m:Psum:bus1')]
    list_cargas[cont] = [dict_barras[c.bus1.cterm.loc_name], c.plini]
    cont += 1

# Lista de lineas
# "Name","from","to","R", "X", "Lmin", "Lmáx"
list_lineas = zeros((nl, 7))
cont1 = 0
for line in lineas:
    list_lineas[cont1] = [dict_lineas[line.loc_name],dict_barras[line.bus1.cterm.loc_name], dict_barras[line.bus2.cterm.loc_name], round(line.R1/Zb,6), round(line.X1/Zb,6),round(-line.Inom_a*Vb*raiz3,4),round(line.Inom_a*Vb*raiz3,4)]
    cont1 += 1

list_lineas = np.array(sorted(list_lineas, key=lambda x:x[0]))
fm = list_lineas[:,6]  # flujo máximo

list_gen = np.zeros((ng, 6))
pos_g = []
cont = 0
for g in generadores:
    pos_g.append(dict_barras[g.bus1.cterm.loc_name])
    potencia = 0 
    if g.outserv == 0:
        potencia = round(g.pgini,3)
        #potencia = round(g.GetAttribute('m:Psum:bus1'),3)
    list_gen[cont] = [dict_gen[g.loc_name],dict_barras[g.bus1.cterm.loc_name], g.ip_ctrl ,potencia,g.Pmin_uc ,g.P_max]
    cont += 1


# #### Lista de barras
demanda_por_barra = np.zeros((nb, 2))
cont = 0
for b in list_bus:
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
t1 = time.time()
# Definir modelo
model = gb.Model('Modelo 1')
#model.Params.MIPGap = 1e-6
model.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

#### Definición de variables
p = model.addVars(ng,vtype=GRB.CONTINUOUS,ub=GRB.INFINITY,lb=0, name='Pg')
#e = model.addVars(ng,vtype=GRB.BINARY, name='e') # Estado del generador
p_pv = model.addVars(n_pv,vtype=GRB.CONTINUOUS,ub=GRB.INFINITY,lb=0, name='Pg_pv')
p_pe = model.addVars(n_pe,vtype=GRB.CONTINUOUS,ub=GRB.INFINITY,lb=0, name='Pg_pe')

#### Función objetivo
f_obj = 0   
for i in range(ng):
#    f_obj += e[i]*(Cv[i]*p[i] + Cf[i])
    f_obj += Cv[i]*p[i] + Cf[i]
for i in range(n_pv):
    f_obj += pv_Cv[i]*p_pv[i] + pv_Cf[i]
for i in range(n_pe):
    f_obj += pe_Cv[i]*p_pe[i]*pe_parallel[i] + pe_Cf[i]    
    

model.setObjective(f_obj,GRB.MINIMIZE)  

#### Restricciones     
# Potencia máxima y mínima
for i in range(ng):
    #model.addConstr(p[i] >= p_min[i]*e[i],'Pmin')
    #model.addConstr(-p[i] >= -p_max[i]*e[i], 'Pmax')
    model.addConstr(p[i] >= p_min[i],'Pmin')
    model.addConstr(-p[i] >= -p_max[i], 'Pmax')
for i in range(n_pv):
    model.addConstr(p_pv[i] >= pv_p_min[i],'Pmin_pv')
    model.addConstr(-p_pv[i] >= -pv_p_max[i], 'Pmax_pv')
for i in range(n_pe):
    model.addConstr(p_pe[i] >= pe_p_min[i],'Pmin_pe')
    model.addConstr(-p_pe[i] >= -pe_p_max[i], 'Pmax_pe')

# Balance de potencia 
#model.addConstr(p.sum('*') + p_pv.sum('*') + p_pe.sum('*') == demanda_total*fact_perd,'Balance')   
model.addConstr(p.sum('*') + p_pv.sum('*') + quicksum(pe_parallel[i]*p_pe[i] for i in range(n_pe)) == demanda_total,'Balance')    


# %%

t2 = time.time()
model.optimize()
model.write('opf.lp')    
t3 = time.time()

# Salida de información
print('-----------------------------')
print('La demanda total del sistema es: %.2f (MW)' % demanda_total)    
status = model.Status
if status == GRB.Status.OPTIMAL:
    print ('Costo => %.2f ($/h)' % model.objVal) 
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
    #if e[cont].X == 0.:
    #    g.outserv = 1
    cont += 1

cont = 0
for pv in PV_plants:
    pv.pgini = p_pv[cont].X
    print('%s => %.2f (MW)' % (pv.loc_name,p_pv[cont].X))
    cont += 1

cont = 0
for pe in PE_plants:
    pe.pgini = p_pe[cont].X
    print('%s => %.2f (MW)' % (pe.loc_name,pe_parallel[cont]*p_pe[cont].X))
    cont += 1
    
print('Tiempo de total: %.4f (s)'% (t3-t0))
print('Tiempo de optimización : %.4f (s)' % (t3-t2))

ldf.Execute()


# %% Corroboración de Resultados

flujos_AC = zeros(nl)
for l in lineas:
    flujos_AC[dict_lineas[l.loc_name]] = l.GetAttribute('m:Psum:bus2')
    
for i in range(nl):
    print('%s ; Gurobi (DC) = %.3f [MW]; DIgSILENT (AC) = %.3f [MW]' % (list_lin[i],f[i].X,flujos_AC[i]))

print('Tiempo de total: %.4f (s)'% (t3-t0))
print('Tiempo de optimización : %.4f (s)' % (t3-t2))


# %%