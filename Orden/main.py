import pfsim
import pandas as pd
import numpy as np
import time
from scipy.sparse import csr_matrix as sparse, identity as sparseI
import gurobipy as gp
from gurobipy import *

#sim = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#sim = pfsim.PowerFactorySim('Ejemplo Clase')
#sim = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')
sim = pfsim.PowerFactorySim('Articulo2')

Gen_AGC = ['G4_3', 'G6_3', 'G2_3', 'G4_1', 'G7_1', 'G2_1', 'G3_1', 'G4_4', 'G6_4', 'G2_4', 'G3_4', 'G4_2', 'G3_2']

before=False
L=6

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
Gen_on_AGC = np.zeros(ngen_eff)
pos_gen_agc_list = list()
Pgen_pre = np.zeros(ngen_eff)
cont=0
for gen in dict_gen_eff:
    pos_gen[cont] = dict_gen_eff[gen][0]
    ngen_par[cont] = dict_gen_eff[gen][1]
    Pmin_gen[cont] = dict_gen_eff[gen][3]/sim.Sb
    Pmax_gen[cont] = dict_gen_eff[gen][4]/sim.Sb
    Pgen_pre[cont] = dict_gen_eff[gen][5]/sim.Sb
    Cvar_gen[cont] = dict_gen_eff[gen][6]
    Ccte_gen[cont] = dict_gen_eff[gen][7]
    if gen in Gen_AGC:
        Gen_on_AGC[cont] = 1
        pos_gen_agc_list.append(cont)
    cont+=1

n_gen_agc = len(pos_gen_agc_list)

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
Genstat_on_AGC = np.zeros(ngenstat_eff)
pos_genstat_agc_list = list()
Pgenstat_pre = np.zeros(ngenstat_eff)
cont=0
for gen in dict_genstat_eff:
    pos_genstat[cont] = dict_genstat_eff[gen][0]
    ngenstat_par[cont] = dict_genstat_eff[gen][1]
    Pmin_genstat[cont] = dict_genstat_eff[gen][3]/sim.Sb
    Pmax_genstat[cont] = dict_genstat_eff[gen][4]/sim.Sb
    Pgenstat_pre[cont] = dict_genstat_eff[gen][5]/sim.Sb
    Cvar_genstat[cont] = dict_genstat_eff[gen][6]
    Ccte_genstat[cont] = dict_genstat_eff[gen][7]
    if gen in Gen_AGC:
        Genstat_on_AGC[cont] = 1
        pos_genstat_agc_list.append(cont)
    cont+=1

n_genstat_agc = len(pos_genstat_agc_list)

P_pre = Pgen_pre.sum() + Pgenstat_pre.sum()


n_elem = len(indices_obj)   # N° de líneas + 'trf2'
# Obtención de fmax, R y X
dict_full = dict()
pos_line = list()
FMax = np.zeros(n_elem)
R = np.zeros(n_elem)
X = np.zeros(n_elem)
i_buses = np.zeros(n_elem).astype(int)
j_buses = np.zeros(n_elem).astype(int)
cont = -1
for i in indices_obj:
    cont += 1
    if i in dict_lineas:
        dict_full[i] = dict_lineas[i]
        Fmax_i = dict_lineas[i][2]
        R_i = dict_lineas[i][0]
        X_i = dict_lineas[i][1]
        i_bus = dict_barras[dict_lineas[i][4]]
        j_bus = dict_barras[dict_lineas[i][5]]
        pos_line.append(cont)
    elif i in dict_trafos:
        dict_full[i] = dict_trafos[i]
        Fmax_i = dict_trafos[i][3]
        R_i = dict_trafos[i][1]
        X_i = dict_trafos[i][0]
        i_bus = dict_barras[dict_trafos[i][4]]
        j_bus = dict_barras[dict_trafos[i][5]]
    else:
        print('Elemento no encontrado en diccs: ' + i)
    FMax[cont] = Fmax_i/sim.Sb
    R[cont] = R_i
    X[cont] = X_i
    i_buses[cont] = i_bus
    j_buses[cont] = j_bus

I = np.r_[range(n_elem), range(n_elem)]
S = sparse((np.r_[np.ones(n_elem), -np.ones(n_elem)], (I, np.r_[i_buses, j_buses])), (n_elem, len(indices_bus)))
A = np.array(S.todense())


val_min = 1e-24 #Valor mínimo para evitar valores nulos en matriz X
yprim = np.zeros(n_elem).astype(complex)
for i in range(n_elem):
    if X[i] == 0:
        X[i] = val_min
    yprim[i] = 1/(complex(R[i],X[i]))
G = np.real(yprim)
B = np.imag(yprim)

# Switch event
# Gen Out
name_events = list()
events_folder = sim.IntEvt.GetContents()
for e in events_folder:
    name_events.append(e.loc_name)
    e.outserv = 1   


#df = pd.DataFrame()

Ns = len(dict_gen_eff)
# Dynamic Simulation until 50 seg - Finishing CPF
t_initial = 0.5
t_final = 300
tstop_cpf = 49.9
Pgen_pfc = np.zeros((ngen_eff, Ns))
Pgenstat_pfc = np.zeros((ngenstat_eff, Ns))
Gen_on_AGC = np.tile(Gen_on_AGC,(Ns,1)).T
Genstat_on_AGC = np.tile(Genstat_on_AGC,(Ns,1)).T
P_out = np.zeros(Ns)
Pos_gen_out = np.zeros(Ns)
cont = 0
if False:
    for gen_out in dict_gen_eff:
        P_out[cont] = dict_gen_eff[gen_out][5]/sim.Sb
        Pos_gen_out[cont] = dict_gen_eff[gen_out][0]
        evt = events_folder[name_events.index('Salida Gen')]
        evt.outserv = 0
        evt.time = t_initial
        evt.p_target = sim.generadores[list(dict_gen).index(gen_out)]
    
        #sim.prepare_dynamic_sim({'*.ElmSym' : ['m:Psum:bus1']}, 'rms', end_time=300)
        sim.prepare_dynamic_sim({}, 'rms', end_time=t_final)
        sim.run_dynamic_sim(end_sim=tstop_cpf)

        Pgen_pfc[:,cont] = np.array(list(map(lambda x: x.GetAttribute('m:Psum:bus1'), sim.generadores)))  # Potencia generadores al final del CPF
        Pgenstat_pfc[:,cont] = np.array(list(map(lambda x: x.GetAttribute('m:Psum:bus1'), sim.genstate)))  # Potencia generadores al final del CPF

        if (gen_out in Gen_AGC) and (gen_out in dict_gen_eff):
            Gen_on_AGC[list(dict_gen_eff.keys()).index(gen_out),cont] = 0

        elif (gen_out in Gen_AGC) and (gen_out in dict_genstat_eff):
            Genstat_on_AGC[list(dict_genstat_eff.keys()).index(gen_out),cont] = 0

        cont += 1

    df = pd.DataFrame(Pgen_pfc).T
    df.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenPFC.csv", index=False)
    df1 = pd.DataFrame(Pgenstat_pfc).T
    df1.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenstatPFC.csv", index=False)
    df2 = pd.DataFrame(Gen_on_AGC).T
    df2.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenAGC.csv", index=False)
    df3 = pd.DataFrame(Genstat_on_AGC).T
    df3.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenStatAGC.csv", index=False)
    df4 = pd.DataFrame(P_out).T
    df4.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\Pout.csv", index=False)
    df4 = pd.DataFrame(Pos_gen_out).T
    df4.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PosPout.csv", index=False)


if True:
    Pgen_pfc = pd.read_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenPFC.csv").to_numpy().T
    try:
        Pgenstat_pfc = pd.read_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenstatPFC.csv").to_numpy().T
    except:
        Pgenstat_pfc = np.zeros((ngenstat_eff, Ns))
    Gen_on_AGC = pd.read_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenAGC.csv").to_numpy().T
    try:
        Genstat_on_AGC = pd.read_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenStatAGC.csv").to_numpy().T
    except:
        Genstat_on_AGC = np.tile(Genstat_on_AGC,(Ns,1)).T
    P_out = pd.read_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\Pout.csv").to_numpy().T
    Pps_gen_out = pd.read_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PosPout.csv").to_numpy().T



Neg_Gen_on_AGC = np.where((Gen_on_AGC==0)|(Gen_on_AGC==1), Gen_on_AGC.astype(int)^1, Gen_on_AGC.astype(int))
Neg_Genstat_on_AGC = np.where((Genstat_on_AGC==0)|(Genstat_on_AGC==1), Genstat_on_AGC.astype(int)^1, Genstat_on_AGC.astype(int))



VOLL = 300

m = gp.Model('Modelo AGC')
#m.setParam('DualReductions', 0)
m.Params.MIPGap = 1e-5
m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

pg_inc = m.addMVar((ngen_eff, Ns), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_inc')
vg_inc = m.addMVar(n_gen_agc, vtype=GRB.BINARY, name='vg_inc')

pstatg_inc = m.addMVar((ngenstat_eff, Ns), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pstatg_inc')
vstatg_inc = m.addMVar(n_genstat_agc, vtype=GRB.BINARY, name='vstatg_inc')


p_ens = m.addMVar((len(indices_bus),Ns), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='P_ens')

f = m.addMVar((n_elem, Ns),vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas
fp = m.addMVar((n_elem, Ns),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='fp') # Flujo-p por las líneas
fn = m.addMVar((n_elem, Ns),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='fn') # Flujo-n por las líneas


#perdidas
ploss = m.addMVar((n_elem, Ns),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='losses') # Flujo por las líneas
dpk = m.addMVar((n_elem,L, Ns), vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo
n_l = m.addMVar((n_elem, Ns), vtype=GRB.BINARY, name='n_l')                                          # variable binaria complentaridad
n_a = m.addMVar((n_elem,L, Ns), vtype=GRB.BINARY, name='n_a')                                          # variable binaria adyacencia



f_obj = 0
#costo_gen = 0
#costo_genstat = 0
Csfc = Ccte_gen[pos_gen_agc_list] @ vg_inc
Csfc += Ccte_genstat[pos_genstat_agc_list] @ vstatg_inc

Cop = 0
Cue = 0
Cpl = 0
for s in range(Ns):
    Cop += Gen_on_AGC[:,s] * Cvar_gen*sim.Sb @ pg_inc[:,s]
    #Cop += Neg_Gen_on_AGC[:,s] * VOLL*sim.Sb @ pg_inc[:,s]
    if np.size(Cvar_genstat) != 0:
        Cop += Genstat_on_AGC[:,s] * Cvar_genstat*sim.Sb @ pstatg_inc[:,s]
        #Cop += Neg_Genstat_on_AGC[:,s] * VOLL *sim.Sb @ pstatg_inc[:,s]

    Cue += (VOLL * p_ens[:,s]).sum() * sim.Sb
    Cpl += (VOLL * ploss[:,s]).sum() * sim.Sb

f_obj = Csfc + Cop + Cue + Cpl
m.setObjective(f_obj,GRB.MINIMIZE)

#Flujo1
f_gen =  SF[:,pos_gen] @ Pgen_pre
f_genstat =  SF[:,pos_genstat] @ Pgenstat_pre

#Balance
for s in range(Ns):
    m.addConstr(P_pre - P_out[s] + pg_inc[pos_gen_agc_list,s].sum() + pstatg_inc[pos_genstat_agc_list,s].sum() + p_ens[:,s].sum() == dda_barra.sum() + ploss[:,s].sum() , name ='Balance')
#Flujo 2
    f_ens =  SF @ p_ens[:,s]
    f_gen_agc = SF[:,pos_gen[pos_gen_agc_list]] @ pg_inc[pos_gen_agc_list,s]
    f_genstat_agc = 0
    if np.size(Cvar_genstat) != 0:
        f_genstat_agc = Genstat_on_AGC[:,s] * SF[:,pos_genstat] @ pstatg_inc[:,s]
    f_gen_out = SF[:,int(Pos_gen_out[s])]*P_out[s]
    f_loss = 0.5 * SF @ abs(A.T) @ ploss[:,s]
    m.addConstr(f[:,s] == f_gen + f_genstat + f_ens + f_gen_agc + f_genstat_agc - Flujo_dda - f_gen_out - f_loss, name='flujos') 
    m.addConstr(f[:,s] == fp[:,s] - fn[:,s], name = 'f')
    m.addConstr(fp[:,s] + fn[:,s] == dpk[:,:,s].sum(1), name='SumaDpk')
    kl = np.zeros((n_elem, L))
    for l in range(L):
        kl[:,l] = (2*(l+1)-1)*(FMax)/L
    m.addConstr(ploss[:,s] == G/(B**2)*(quicksum(kl[:,i]*dpk[:,i,s] for i in range(L))), name = 'Ploss')  
    m.addConstr(-f[:,s] - 0.5*ploss[:,s] >= -FMax, name = 'fp')
    m.addConstr(f[:,s] - 0.5*ploss[:,s] >= -FMax, name = 'fn')
    m.addConstr(-fp[:,s] >= -FMax, name = 'fp+')
    m.addConstr(-fn[:,s] >= -FMax, name = 'fn+')
    for l in range(L):
        m.addConstr(-dpk[:,l,s] >= -FMax/L, name = 'LimiteDpk')
    
    m.addConstr(-fp[:,s]>=-n_l[:,s]*FMax, name='fp')   #flujo positvo-restriccion de complementaridad
    m.addConstr(-fn[:,s]>=(1-n_l[:,s])*(-FMax), name='fn') #flujo nefativo-restriccion de complementaridad

    for l in range(L): 
        if l==0:
            m.addConstr(-dpk[:,l,s]>=-FMax/(L), name='d_f_Res_max_A_l')
            m.addConstr(dpk[:,l,s]>=n_a[:,l,s]*(FMax/(L)), name='d_f_Res_min_A_l')
        elif l==L-1:
            m.addConstr(-dpk[:,l,s]>=-n_a[:,l-1,s]*FMax/(L), name='d_f_Res_max_A_L')
            m.addConstr(dpk[:,l,s]>=0, name='d_f_Res_min_A_L')
        else:
            m.addConstr(-dpk[:,l,s]>=-n_a[:,l-1,s]*(FMax/(L)), name='d_f_Res_max_A_L-1')
            m.addConstr(dpk[:,l,s]>=n_a[:,l,s]*(FMax/(L)), name='d_f_Res_min_A_L-1')


# Generadores
    m.addConstr(-pg_inc[pos_gen_agc_list,s] >= -vg_inc * (Pmax_gen[pos_gen_agc_list] * ngen_par[pos_gen_agc_list] - Pgen_pre[pos_gen_agc_list]))
    if np.size(Cvar_genstat) != 0:
        m.addConstr(-pg_inc[pos_genstat_agc_list,s] >= -vstatg_inc * (Pmax_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list] - Pgen_pre[pos_genstat_agc_list]))

    m.addConstr(-p_ens[:,s] >= -dda_barra)

# %%
m.write('OPF.lp')

t2 = time.time()
m.optimize()
t3=time.time()

print('-----------------------------')
print('La demanda total del sistema es: %.2f (MW)' % (dda_barra.sum()*sim.Sb))    
print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($) + Cpl = %.2f ($)' % (m.objVal,Csfc.getValue(),Cop.getValue(),Cue.getValue(),Cpl.getValue()))
print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (m.NumVars,m.NumConstrs,m.DNumNZs)) #print('num_Vars =  %d / num_Const =  %d' % (len(m.getVars()), len(m.getConstrs())))      
print ('Total P_loss = %.2f [MW]'%(sim.Sb*ploss.sum().getValue()))
#print('=> Formulation time: %.4f (s)'% (t1-t0))
print('=> Solution time: %.4f (s)' % (t3-t2))
print('=> Solver time: %.4f (s)' % (m.Runtime))
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

part_factors = np.zeros((ngen_eff,Ns))
for gen, s in np.ndindex(pg_inc.x.shape):
    if pg_inc.x[gen,s] != 0:
        part_factors[gen,s] = pg_inc.x[gen,s]/P_out[s]

#df = pd.DataFrame(part_factors).T
#df.to_csv(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PartFactors.csv", index=False)




#results_line = np.vstack((np.array(indices_obj),f.x))[:,pos_line]

#sim.check_results(dict_gen_eff, dict_genstat_eff, p_g, p_statg, results_line)

#cont=0
#for gen in dict_gen_eff:
#    print('%s => %.2f (MW)' % (gen,ngen_par[cont]*p_g[cont].X*sim.Sb))
#    cont += 1
#
#cont=0
#for genstat in dict_genstat_eff:
#    print('%s => %.2f (MW)' % (genstat,ngenstat_par[cont]*p_statg[cont].X*sim.Sb))
#    cont += 1
#
#cont=0
#for elem in indices_obj:
#    print('Línea: %s => %.2f (MW)' % (elem,f[cont].X*sim.Sb))
#    #print("f_p = %.3f // f_n = %.3f" % (fp.X[cont],fn.X[cont]))
#    #print("d_f[0]=%.3f     -     d_f[1]=%.3f " % (dpk.X[cont,0],dpk.X[cont,1]))    
#    cont += 1


print('finish!')
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

#### OLD CODE FOR LOSSES

#if before:
#    # Definir Modelo
#    m = gp.Model('Modelo 1')
#    m.Params.MIPGap = 1e-5
#    m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi
#
#    p_g = m.addMVar(ngen_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg')
#    pdis_g = m.addMVar(ngen_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_dis')
#    p_statg = m.addMVar(ngenstat_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Psg')
#    pdis_statg = m.addMVar(ngenstat_eff, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Psg_dis')
#
#    f = m.addMVar(n_elem,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas
#    fp = m.addMVar(n_elem,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='fp') # Flujo-p por las líneas
#    fn = m.addMVar(n_elem,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='fn') # Flujo-n por las líneas
#    #perdidas
#    ploss = m.addMVar(n_elem,vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='losses') # Flujo por las líneas
#    dpk = m.addMVar((n_elem,L), vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo
#    n_l = m.addMVar(n_elem, vtype=GRB.BINARY, name='n_l')                                          # variable binaria complentaridad
#    n_a = m.addMVar((n_elem,L), vtype=GRB.BINARY, name='n_a')                                          # variable binaria adyacencia
#
#
#
#
#    # Función objetivo
#    f_obj = 0
#    #costo_gen = 0
#    #costo_genstat = 0
#
#    costo_gen = p_g* ngen_par @ Cvar_gen*sim.Sb + Ccte_gen.sum()
#    costo_genstat = p_statg* ngenstat_par @ Cvar_genstat*sim.Sb + Ccte_genstat.sum()
#
#    f_obj = costo_gen + costo_genstat
#
#    m.setObjective(f_obj, GRB.MINIMIZE)
#    m.getObjective()
#
#    ###  Restricciones
#    ## Balance Nodal
#    m.addConstr(p_g @ ngen_par + p_statg @ ngenstat_par == dda_barra.sum(), name='Balance')
#    m.addConstr(pdis_g @ ngen_par + p_statg @ ngenstat_par >= 1.1*dda_barra.sum(), name="Reserva")
#
#
#    ## Restricciones de generación
#    # Limitaciones de Potencia
#    m.addConstr(p_g >= Pmin_gen, name='Pmin_gen')
#    m.addConstr(-p_g >= -Pmax_gen, name='Pmax_gen')
#    m.addConstr(pdis_g >= 0, name='Pmin_gen_disp')
#    m.addConstr(-pdis_g >= -Pmax_gen, name='Pmax_gen_dip')
#
#    m.addConstr(p_statg >= Pmin_genstat, name='Pmin_genstat')
#    m.addConstr(-p_statg >= -Pmax_genstat, name='Pmin_genstat')
#    m.addConstr(pdis_statg >= 0, name='Pmin_genstat_disp')
#    m.addConstr(-pdis_statg >= -Pmax_genstat, name='Pmax_genstat_dip')
#
#
#    ## Sistema de transmisión
#    m.addConstr(f == SF[:,pos_gen] @ (p_g*ngen_par) + SF[:,pos_genstat] @ (p_statg*ngenstat_par)- Flujo_dda, name='sf')
#    m.addConstr(f == fp - fn, name = 'f')
#    m.addConstr(fp + fn == dpk.sum(1), name = 'SumDpk')
#    kl = np.zeros((n_elem,L))
#    for l in range(L):
#        kl[:,l] = (2*(l+1)-1)*FMax/L
#        m.addConstr(-dpk[:,l] >= -FMax/L, name = 'LimiteDpk')
#    m.addConstr(ploss == G/(B**2)*(quicksum(kl[:,i]*dpk[:,i] for i in range(L))), name = 'Ploss')
#    m.addConstr(-f - 0.5*ploss >= -FMax, name = 'fp')
#    m.addConstr(f - 0.5*ploss >= -FMax, name = 'fn')
#    m.addConstr(-fp >= -FMax, name = 'fp+')
#    m.addConstr(-fn >= -FMax, name = 'fn+')
#
#    m.addConstr(-fp>=-n_l*FMax, name='Cfp')   #flujo positvo-restriccion de complementaridad
#    m.addConstr(-fn>=(-1+n_l)*FMax, name='Cfn') #flujo negativo-restriccion de complementaridad
#
#    for l in range(L): 
#        if l==0:
#            m.addConstr(-dpk[:,l]>=-FMax/(L), name='d_f_Res_max_A_l')
#            m.addConstr(dpk[:,l]>=n_a[:,l]*(FMax/(L)), name='d_f_Res_min_A_l')
#        elif l==L-1:
#            m.addConstr(-dpk[:,l]>=-n_a[:,l-1]*FMax/(L), name='d_f_Res_max_A_L')
#            m.addConstr(dpk[:,l]>=0, name='d_f_Res_min_A_L')
#        else:
#            m.addConstr(-dpk[:,l]>=-n_a[:,l-1]*(FMax/(L)), name='d_f_Res_max_A_L-1')
##            m.addConstr(dpk[:,l]>=n_a[:,l]*(FMax/(L)), name='d_f_Res_min_A_L-1')
