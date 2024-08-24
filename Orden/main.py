import sys
import pfsim
import pandas as pd
import numpy as np
import time
from scipy.sparse import csr_matrix as sparse, identity as sparseI
import gurobipy as gp
from gurobipy import *
import optm

#sim = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#sim = pfsim.PowerFactorySim('Ejemplo Clase')
#sim = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')

Estudio = 1
if Estudio == 1:
    CASES = ['A'] # Sin sistema de transmisión y VOLL = 150000
    VOLL = 150000
elif Estudio == 2:
    CASES = ['B']
    VOLL = 300
    Costo_ts = 50
    M = 1e6
elif Estudio == 3:
    CASES = ['C']









t0 = time.time()
sim = pfsim.PowerFactorySim('Articulo2')

Gen_AGC = ['G2_1', 'G3_1', 'G4_1', 'G7_1', 'G2_3', 'G4_3', 'G6_3', 'G3_2', 'G4_2', 'G2_4', 'G3_4', 'G4_4', 'G6_4']
#Gen_outages = ['G5_1']
#Gen_outages = ['G4_1', 'G5_1', 'G6_1', 'G2_2', 'G5_2']
Gen_outages = ['G4_1', 'G5_1', 'G6_1', 'G2_2', 'G5_2', 'G5_3', 'G3_4', 'G5_4', 'G7_4']
Ns = len(Gen_outages)
escenarios = ['peak', 'media', 'valle']

TS_lines = ['']
nl_ts=len(TS_lines)
TS = False



Nt = 3 # N° de escenarios
before=False
L=11
T_Sfc = 15 # 15 Min, <-- Tiempo de operación AGC

ptdf_dataframe, indices_bus = sim.export_csv('AC')
#ptdf_dataframe = pd.read_excel(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.xlsx',skiprows = 0,delimiter=';')
#ptdf_dataframe = pd.read_excel(r'C:\Users\lldie\Desktop\SF_datos.xlsx',skiprows = 0,delimiter=';')
SF, indices_obj = pfsim.ShiftFactors(ptdf_dataframe)

(dict_barras, dict_cargas, dict_lineas, dict_trafos, dict_gen, dict_genstat, bus_slack) = sim.get_data(indices_bus, Nt)



# demanda por barra
Nb = len(indices_bus)
dda_barra = np.zeros((Nb, Nt))
for load in dict_cargas:
    for ti in range(Nt):
        dda_barra[dict_cargas[load][0], ti] += dict_cargas[load][ti+2]/sim.Sb

#Flujo_dda = SF @ dda_barra   # Flujo en elemntos debido a la demanda en barra





# GENERADORES EFECTIVOS - EN SERVICIO
#dict_gen = dict()
#for gen in dict_gen:
#    if dict_gen[gen][2] == 0:
#        dict_gen[gen] = dict_gen[gen]
ngen=len(dict_gen)

pos_gen = np.zeros(ngen).astype(int)
ngen_par = np.zeros(ngen)
Cvar_gen = np.zeros(ngen)
Ccte_gen = np.zeros(ngen)
Pmin_gen = np.zeros(ngen)
Pmax_gen =np.zeros(ngen)
Gen_on_AGC = np.zeros(ngen)
pos_gen_agc_list = list()
name_gen_agc_list = list()
Pgen_pre = np.zeros((ngen,Nt))
Ramp_gen = np.zeros(ngen)
all_gen = list()
cont=0
for gen in dict_gen:
    all_gen.append(gen)
    pos_gen[cont] = dict_gen[gen][0]
    ngen_par[cont] = dict_gen[gen][1]
    Pmin_gen[cont] = dict_gen[gen][3]/sim.Sb
    Pmax_gen[cont] = dict_gen[gen][4]/sim.Sb
    Cvar_gen[cont] = dict_gen[gen][5]
    Ccte_gen[cont] = dict_gen[gen][6]
    Ramp_gen[cont] = dict_gen[gen][7]/sim.Sb
    for ti in range(Nt):
        Pgen_pre[cont,ti] = dict_gen[gen][8+ti]/sim.Sb
    if gen in Gen_AGC:
        Gen_on_AGC[cont] = 1
        pos_gen_agc_list.append(cont)
        name_gen_agc_list.append(gen)
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
Ramp_genstat = np.zeros(ngenstat_eff)
cont=0
for gen in dict_genstat_eff:
    pos_genstat[cont] = dict_genstat_eff[gen][0]
    ngenstat_par[cont] = dict_genstat_eff[gen][1]
    Pmin_genstat[cont] = dict_genstat_eff[gen][3]/sim.Sb
    Pmax_genstat[cont] = dict_genstat_eff[gen][4]/sim.Sb
    Pgenstat_pre[cont] = dict_genstat_eff[gen][5]/sim.Sb
    Cvar_genstat[cont] = dict_genstat_eff[gen][6]
    Ccte_genstat[cont] = dict_genstat_eff[gen][7]
    Ramp_genstat[cont] = dict_genstat_eff[gen][8]
    if gen in Gen_AGC:
        Genstat_on_AGC[cont] = 1
        pos_genstat_agc_list.append(cont)
    cont+=1

n_genstat_agc = len(pos_genstat_agc_list)

#P_pre = Pgen_pre.sum() + Pgenstat_pre.sum()


n_elem = len(indices_obj)   # N° de líneas + 'trf2'
# Obtención de fmax, R y X
dict_full = dict()
pos_line = list()
all_line = list()
FMax = np.zeros(n_elem)
R = np.zeros(n_elem)
X = np.zeros(n_elem)
i_buses = np.zeros(n_elem).astype(int)
j_buses = np.zeros(n_elem).astype(int)
pos_ts = list()
cont = -1
for i in indices_obj:
    cont += 1
    if i in dict_lineas:
        all_line.append(i)
        dict_full[i] = dict_lineas[i]
        Fmax_i = dict_lineas[i][2]
        R_i = dict_lineas[i][0]
        X_i = dict_lineas[i][1]
        i_bus = dict_barras[dict_lineas[i][4]]
        j_bus = dict_barras[dict_lineas[i][5]]
        pos_line.append(cont)
        if i in TS_lines:
            pos_ts.append(cont)

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

pos_nots = [pos for pos in pos_line if pos not in pos_ts]
# Lineas Candidatas TS
i_buses_ts = i_buses[pos_ts]
j_buses_ts = j_buses[pos_ts]

I = np.r_[range(n_elem), range(n_elem)]
S = sparse((np.r_[np.ones(n_elem), -np.ones(n_elem)], (I, np.r_[i_buses, j_buses])), (n_elem, Nb))
A = np.array(S.todense())


val_min = 1e-24 #Valor mínimo para evitar valores nulos en matriz X
yprim = np.zeros(n_elem).astype(complex)
for i in range(n_elem):
    if X[i] == 0:
        X[i] = val_min
    yprim[i] = 1/(complex(R[i],X[i]))
G = np.real(yprim)
B = np.imag(yprim)

## CALCULO DE SHIFT FACTORS

n_bus_slack = list(dict_barras.keys()).index(bus_slack)
noslack = np.flatnonzero(np.arange(Nb) != n_bus_slack)


BfR = sparse((np.r_[np.imag(yprim), -np.imag(yprim)], (I, np.r_[i_buses, j_buses])), (n_elem,Nb))
BbusR = S.T * BfR
SFR = np.zeros((n_elem,Nb))
SFR[:,noslack] = BfR[:, noslack].todense()*np.linalg.inv(BbusR[np.ix_(noslack, noslack)].todense())    

difference = (SFR - SF)/SFR

SF_PF = pd.DataFrame(np.vstack((list(dict_barras.keys()),SF)).T, columns=np.insert(indices_obj, 0, 'x', axis=0))
SF_DC_Lossy = pd.DataFrame(np.vstack((list(dict_barras.keys()),SFR)).T, columns=np.insert(indices_obj, 0, 'x', axis=0))
diff = pd.DataFrame(np.vstack((list(dict_barras.keys()),difference)).T, columns=np.insert(indices_obj, 0, 'x', axis=0))

SF= SFR


with pd.ExcelWriter('3ShiftFactors.xlsx') as writer:
    SF_PF.to_excel(writer, sheet_name='SF_PF', index=False)
    SF_DC_Lossy.to_excel(writer, sheet_name='SF_DC', index=False)
    diff.to_excel(writer, sheet_name='Diff', index=False)





# Switch event
# Gen Out
name_events = list()
events_folder = sim.IntEvt.GetContents()
for e in events_folder:
    name_events.append(e.loc_name)
    e.outserv = 1   


#df = pd.DataFrame()

class Simulacion(object):
    def __init__(self):
        self.pos_gen_agc_list = pos_gen_agc_list


    
    # Dynamic Simulation until 50 seg - Finishing CPF
    t_initial = 0.5 # Tiempo Perturbación
    t_final = 420 # Tiempo total
    tstop_cpf = 29.9 # Tiempo hasta el CPF
    D_pfc = np.zeros((Nb, Ns, Nt)) # Demanda al finaliza el CPF
    PL_pre_line = np.zeros(Nt) # Perdidas al finalizar el CPF
    #PL_pre_trafo2 = np.zeros((len(dict_trafos),Ns))
    Pgen_pfc = np.zeros((ngen, Ns, Nt))
    Pgenstat_pfc = np.zeros((ngenstat_eff, Ns, Nt))
    Gen_on_AGC = np.tile(Gen_on_AGC,(Ns,1)).T
    Gen_on_AGC = np.stack([Gen_on_AGC]*3,axis=0)
    Genstat_on_AGC = np.tile(Genstat_on_AGC,(Ns,1)).T
    P_out = np.zeros((Ns, Nt))
    Barra_gen_out = np.zeros((Ns))
    for ti in range(Nt):
        sim.scens[ti].Activate()
        print('Escenario: ' + str(ti+1))
        cont = 0
        for gen_out in Gen_outages:
            sys.stdout.write('Progreso: ' + str(cont+1+Ns*ti) + '/' + str(Ns*Nt))
            sys.stdout.write('Generador Out: ' + gen_out)
            sys.stdout.flush()
            P_out[cont,ti] = dict_gen[gen_out][8+ti]/sim.Sb
            Barra_gen_out[cont] = dict_gen[gen_out][0]
            evt = events_folder[name_events.index('Salida Gen')]
            evt.outserv = 0
            evt.time = t_initial
            evt.p_target = sim.generadores[list(dict_gen).index(gen_out)]

            #sim.prepare_dynamic_sim({'*.ElmSym' : ['m:Psum:bus1']}, 'rms', end_time=300)
            sim.prepare_dynamic_sim({}, 'rms', end_time=t_final)
            if PL_pre_line[ti] == 0: # Obtener perdidas solo una vez para cada contingencia y cada operación demanda
                PL_pre_line[ti] = np.array(list(map(lambda x: x.GetAttribute('c:Losses'), sim.lineas))).sum()
            sim.run_dynamic_sim(end_sim=tstop_cpf)

            for load in list(map(lambda x: (x.GetAttribute('m:Psum:bus1'),dict_barras[x.bus1.cterm.loc_name]), sim.cargas)):  # Demanda al final del CPF por carga
                D_pfc[load[1], cont, ti] += load[0]/sim.Sb

            #PL_pre_trafo2[:,cont] = np.array(list(map(lambda x: x.GetAttribute('c:Ploss'), sim.trafos)))  


            for gen in sim.generadores:
                if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                    potencia = gen.GetAttribute('m:Psum:bus1')/sim.Sb
                else:
                    potencia = 0
                Pgen_pfc[list(all_gen).index(gen.loc_name),cont, ti] = potencia

            for gen in sim.genstate:
                Pgenstat_pfc[list(all_gen).index(gen.loc_name),cont, ti] = gen.GetAttribute('m:Psum:bus1')/sim.Sb

            #Pgen_pfc[:,cont] = np.array(list(map(lambda x: x.GetAttribute('m:Psum:bus1'), sim.generadores)))/sim.Sb  # Potencia generadores al final del CPF
            #Pgenstat_pfc[:,cont,ti] = np.array(list(map(lambda x: x.GetAttribute('m:Psum:bus1'), sim.genstate)))/sim.Sb  # Potencia generadores al final del CPF

            if (gen_out in Gen_AGC) and (gen_out in dict_gen):
                Gen_on_AGC[ti, list(all_gen).index(gen_out),cont] = 0

            elif (gen_out in Gen_AGC) and (gen_out in dict_genstat_eff):
                Genstat_on_AGC[ti, list(dict_genstat_eff.keys()).index(gen_out),cont] = 0

            cont += 1

#df0 = pd.DataFrame(np.array(D_pfc)).T
#df0.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\DPFC.xlsx", index=False)
#df = pd.DataFrame(Pgen_pfc).T
#df.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenPFC.xlsx", index=False)
#df1 = pd.DataFrame(Pgenstat_pfc).T
#df1.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenstatPFC.xlsx", index=False)
#df2 = pd.DataFrame(Gen_on_AGC).T
#df2.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenAGC.xlsx", index=False)
#df3 = pd.DataFrame(Genstat_on_AGC).T
#df3.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenStatAGC.xlsx", index=False)
#df4 = pd.DataFrame(P_out).T
#df4.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\Pout.xlsx", index=False)
#df5 = pd.DataFrame(Barra_gen_out).T
#df5.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\BarraGenOut.xlsx", index=False)
#df6 = pd.DataFrame(PL_pre_line).T
#df6.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PLLine.xlsx", index=False)



if False:
    PL_pre_line = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PLLine.xlsx").to_numpy().T
    D_pfc = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\DPFC.xlsx").to_numpy().T
    Pgen_pfc = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenPFC.xlsx").to_numpy().T
    try:
        Pgenstat_pfc = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PgenstatPFC.xlsx").to_numpy().T
    except:
        Pgenstat_pfc = np.zeros((ngenstat_eff, Ns))
    Gen_on_AGC = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenAGC.xlsx").to_numpy().T
    try:
        Genstat_on_AGC = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\GenStatAGC.xlsx").to_numpy().T
    except:
        Genstat_on_AGC = np.tile(Genstat_on_AGC,(Ns,1)).T
    P_out = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\Pout.xlsx").to_numpy()[0]
    Barra_gen_out = pd.read_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\BarraGenOut.xlsx").to_numpy()[0]



#Neg_Gen_on_AGC = np.where((Gen_on_AGC==0)|(Gen_on_AGC==1), Gen_on_AGC.astype(int)^1, Gen_on_AGC.astype(int))
#Neg_Genstat_on_AGC = np.where((Genstat_on_AGC==0)|(Genstat_on_AGC==1), Genstat_on_AGC.astype(int)^1, Genstat_on_AGC.astype(int))

t1 = time.time()







for Case in CASES:
    m = optm.Modelo()

    (pg_inc,pg_dec,vg_inc,vg_dec,p_ens) = optm.VarPotencia(m,n_gen_agc,Ns,Nt,Nb)

    (f, fp, fn) = optm.Flujos(m,n_elem,Ns,Nt)

    (ploss,dpk) = optm.Perdidas(m,n_elem,Ns,Nt,L)


    m = gp.Model('Modelo AGC')
    #m.setParam('DualReductions', 0)
    m.Params.MIPGap = 1e-5
    m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi
    
    pg_inc = m.addMVar((n_gen_agc, Ns, Nt), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_inc')
    pg_dec = m.addMVar((n_gen_agc, Ns, Nt), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_dec')
    vg_inc = m.addMVar(n_gen_agc, vtype=GRB.BINARY, name='vg_inc')
    vg_dec = m.addMVar(n_gen_agc, vtype=GRB.BINARY, name='vg_dec')
    
    #pstatg_inc = m.addMVar((ngenstat_eff, Ns), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pstatg_inc')
    #pstatg_dec = m.addMVar((ngenstat_eff, Ns), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pstatg_dec')
    #vstatg_inc = m.addMVar(n_genstat_agc, vtype=GRB.BINARY, name='vstatg_inc')
    #vstatg_dec = m.addMVar(n_genstat_agc, vtype=GRB.BINARY, name='vstatg_dec')
    
    
    p_ens = m.addMVar((Nb,Ns, Nt), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='P_ens')
    
    f = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas
    fp = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='fp') # Flujo-p por las líneas
    fn = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='fn') # Flujo-n por las líneas
    
    
    #perdidas
    ploss = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='losses') # Flujo por las líneas
    dpk = m.addMVar((n_elem, L, Ns, Nt), vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo
    
    # transmission switching
    s_ts = m.addMVar((nl_ts, Ns, Nt), vtype=GRB.BINARY, name='s_ts') 
    
    f_obj = 0
    f_obj_1 = 0
    #costo_gen = 0
    #costo_genstat = 0
    Csfc = Ccte_gen[pos_gen_agc_list] @ (vg_inc + vg_dec)
    #Csfc += Ccte_genstat[pos_genstat_agc_list] @ (vstatg_inc + vstatg_dec)
    Cop = 0
    Cue = 0
    Cpl = 0
    Cts = 0
    for ti in range(Nt):
        Cop_s = 0
        Cue_s = 0
        Cpl_s = 0
        Cts_s = 0
        for s in range(Ns):
            Cop_s += Gen_on_AGC[ti,:,s][pos_gen_agc_list] * Cvar_gen[pos_gen_agc_list]*sim.Sb @ (pg_inc[:,s,ti] + pg_dec[:,s,ti])
            if np.size(Cvar_genstat) != 0:
                Cop_s += Genstat_on_AGC[:,s] * Cvar_genstat*sim.Sb @ pstatg_inc[:,s]
    
            Cue_s += (VOLL * p_ens[:,s,ti]).sum() * sim.Sb
            Cpl_s += (VOLL * ploss[:,s,ti]).sum() * sim.Sb
            Cts_s += Costo_ts*(1-s_ts[:,s,ti]).sum()
    
        Cop += Cop_s
        Cue += Cue_s
        Cpl += Cpl_s
        Cts += Cts_s
    
    f_obj = Csfc + Cop + Cue + Cpl + Cts
    m.setObjective(f_obj,GRB.MINIMIZE)
    
    
    for ti in range(Nt):
        for s in range(Ns):
            PL_pre = PL_pre_line[ti]/sim.Sb/1000 #+ PL_pre_trafo2[:,s].sum()
            #Balance
            #m.addConstr(P_pre - P_out[s] + pg_inc[pos_gen_agc_list,s].sum() - pg_dec[pos_gen_agc_list,s].sum() + pstatg_inc[pos_genstat_agc_list,s].sum() - pstatg_dec[pos_genstat_agc_list,s].sum() + p_ens[:,s].sum() == dda_barra.sum() + ploss[:,s].sum() , name ='Balance')
            #m.addConstr(P_pre - P_out[s] + pg_inc[pos_gen_agc_list,s].sum() - pg_dec[pos_gen_agc_list,s].sum() + pstatg_inc[pos_genstat_agc_list,s].sum() - pstatg_dec[pos_genstat_agc_list,s].sum() + p_ens[:,s].sum() == D_pfc[:,s].sum() + ploss[:,s].sum() , name ='Balance')
            m.addConstr(pg_inc[:,s,ti].sum() - pg_dec[:,s,ti].sum() + pstatg_inc[:,s].sum() - pstatg_dec[:,s].sum() + p_ens[:,s].sum() == P_out[s,ti] - dda_barra[:,ti].sum() + D_pfc[:,s,ti].sum() - PL_pre + ploss[:,s,ti].sum(), name ='Balance')
            
            #Flujo Basic
            if TS == False:
                f_gen =  SF[:,pos_gen] @ Pgen_pre[:,ti]
                f_genstat =  SF[:,pos_genstat] @ Pgenstat_pre
                f_ens =  SF @ p_ens[:,s,ti]
                f_gen_agc = SF[:,pos_gen[pos_gen_agc_list]] @ (pg_inc[:,s,ti] - pg_dec[:,s,ti])
                f_genstat_agc = 0
                #if np.size(Cvar_genstat) != 0:
                #    f_genstat_agc = Genstat_on_AGC[:,s] * SF[:,pos_genstat] @ (pstatg_inc[:,s] - pstatg_dec[:,s])
                f_gen_out = SF[:,int(Barra_gen_out[s])]*P_out[s,ti]
                f_loss = 0.5 * SF @ abs(A.T) @ ploss[:,s,ti]
                Flujo_dda = SF @ D_pfc[:,s,ti]
                m.addConstr(f[:,s,ti] == f_gen + f_genstat + f_ens + f_gen_agc + f_genstat_agc - Flujo_dda - f_gen_out - f_loss, name='flujos') 
                m.addConstr(f[:,s,ti] == fp[:,s,ti] - fn[:,s,ti], name = 'f')
                m.addConstr(fp[:,s,ti] + fn[:,s,ti] == dpk[:,:,s,ti].sum(1), name='SumaDpk')
                kl = np.zeros((n_elem, L))
                for l in range(L):
                    kl[:,l] = (2*(l+1)-1)*(FMax)/L
                m.addConstr(ploss[:,s,ti] == G/(B**2)*(quicksum(kl[:,i]*dpk[:,i,s,ti] for i in range(L))), name = 'Ploss')  
                m.addConstr(-f[:,s,ti] - 0.5*ploss[:,s,ti] >= -FMax, name = 'fp')
                m.addConstr(f[:,s,ti] - 0.5*ploss[:,s,ti] >= -FMax, name = 'fn')
                m.addConstr(-fp[:,s,ti] >= -FMax, name = 'fp+')
                m.addConstr(-fn[:,s,ti] >= -FMax, name = 'fn+')
                for l in range(L):
                    m.addConstr(-dpk[:,l,s,ti] >= -FMax/L, name = 'LimiteDpk')
            #Flujo con transmission switching
            elif TS == True:
                #Líneas no candidatas
                fe_gen =  SF[pos_nots,:][:,pos_gen] @ Pgen_pre[:,ti]
                fe_ens =  SF[pos_nots,:] @ p_ens[:,s,ti]
                fe_gen_agc = SF[pos_nots,:][:,pos_gen[pos_gen_agc_list]] @ (pg_inc[:,s,ti] - pg_dec[:,s,ti])
                fe_gen_out = SF[pos_nots,:][:,int(Barra_gen_out[s])]*P_out[s,ti]
                fe_loss = 0.5 * SF[pos_nots,:] @ abs(A.T) @ ploss[:,s,ti]
                Flujoe_dda = SF[pos_nots,:] @ D_pfc[:,s,ti]
                fv = (SF[pos_nots,:] @ A[pos_ts,:].T) @ f[:,s,ti]
                fe = fe_gen + fe_ens + fe_gen_agc - Flujoe_dda - fe_gen_out - fe_loss
                m.addConstr(-(fe+fv) >= -FMax[pos_nots], name = 'fe_p')
                m.addConstr(fe+fv >= -FMax[pos_nots], name = 'fe_n')
    
    
            # Generadores
            m.addConstr(-pg_inc[:,s,ti] >= -vg_inc * (Pmax_gen[pos_gen_agc_list] * ngen_par[pos_gen_agc_list] - Pgen_pre[pos_gen_agc_list,ti]), name = 'PMax')
            m.addConstr(-pg_dec[:,s,ti] >= -vg_dec * (Pgen_pre[pos_gen_agc_list,ti] - Pmin_gen[pos_gen_agc_list] * ngen_par[pos_gen_agc_list]), name = 'PMin')
            #Nueva Lista que solo tiene unidades participantes en el AGC / se quita de unidades participantes la unidad que sale de servicio
            pos_part_gen_agc_list = pos_gen_agc_list[:]
            x = list(range(n_gen_agc))
            if Gen_outages[s] in name_gen_agc_list:
                pos_part_gen_agc_list.remove(all_gen.index(Gen_outages[s]))
                x.pop(name_gen_agc_list.index(Gen_outages[s]))
                m.addConstr(pg_inc[name_gen_agc_list.index(Gen_outages[s]), s, ti] == 0, name = 'GenOut+')
                m.addConstr(pg_dec[name_gen_agc_list.index(Gen_outages[s]), s, ti] == 0, name = 'GenOut-')
            m.addConstr(-pg_inc[x, s, ti] >= -T_Sfc * Ramp_gen[pos_part_gen_agc_list] - (Pgen_pfc[pos_part_gen_agc_list,s,ti] - Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + Gen_outages[s] +'_Ramp+')
            m.addConstr(-pg_dec[x, s, ti] >= -T_Sfc * Ramp_gen[pos_part_gen_agc_list] + (Pgen_pfc[pos_part_gen_agc_list,s,ti] - Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + Gen_outages[s] +'_Ramp-')
            if np.size(Cvar_genstat) != 0:
                m.addConstr(-pg_inc[pos_genstat_agc_list, s, ti] >= -vstatg_inc * (Pmax_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list] - Pgen_pre[pos_genstat_agc_list,ti]))
                m.addConstr(-pg_dec[pos_genstat_agc_list, s, ti] >= -vstatg_dec * (Pgen_pre[pos_genstat_agc_list,ti] - Pmin_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list]))
                #m.addConstr(-pg_inc[pos_genstat_agc_list,s] >= -T_Sfc * Ramp_genstat[pos_genstat_agc_list] - (Pgenstat_pfc[pos_genstat_agc_list, s] - Pgenstat_pre[pos_genstat_agc_list]))
                #m.addConstr(-pg_dec[pos_genstat_agc_list,s] >= -T_Sfc * Ramp_genstat[pos_genstat_agc_list] + (Pgenstat_pfc[pos_genstat_agc_list, s] - Pgenstat_pre[pos_genstat_agc_list]))
    
    
            m.addConstr(-p_ens[:,s,ti] >= -D_pfc[:,s,ti], 'LimENS')
    
    # %%
    m.write('OPF.lp')
    
    t2 = time.time()
    m.optimize()
    t3=time.time()
    
    status = m.Status
    if status == GRB.Status.OPTIMAL:
        print('-----------------------------')
        #print('La demanda total del sistema es: %.2f (MW)' % (dda_barra[:,ti].sum()*sim.Sb))    
        print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($) + Cpl = %.2f ($)' % (m.objVal,Csfc.getValue(),Cop.getValue(),Cue.getValue(),Cpl.getValue()))
        print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (m.NumVars,m.NumConstrs,m.DNumNZs)) #print('num_Vars =  %d / num_Const =  %d' % (len(m.getVars()), len(m.getConstrs())))      
        #print ('Total P_loss = %.2f [MW]'%(sim.Sb*ploss.sum().getValue()))
        print('=> Extract data time: %.4f (s)'% (t1-t0))
        print('=> Formulation time: %.4f (s)' % (t2-t1))
        print('=> Solution time: %.4f (s)' % (t3-t2))
        print('=> Solver time: %.4f (s)' % (m.Runtime))
        print ('Costo => %.2f ($/h)' % m.objVal) 
        #print ('Las perdidas son %.2f (MW)' % sum(pk_loss[i].X for i in range(len(pk_loss)))) 
    
        Post_gen_out = np.zeros((Ns,Nt))
    
        part_factors = np.zeros((n_gen_agc,Ns,Nt))
        for gen, s, ti in np.ndindex(pg_inc.x.shape):
            Post_gen_out[s, ti] = P_out[s,ti] - dda_barra[:,ti].sum() + D_pfc[:,s,ti].sum() - PL_pre_line[ti]/sim.Sb/1000 +  ploss[:,s,ti].x.sum()
            if pg_inc.x[gen,s,ti] != 0 or pg_dec.x[gen,s,ti] != 0:
                part_factors[gen,s,ti] = (pg_inc.x[gen,s,ti]-pg_dec.x[gen,s,ti])/Post_gen_out[s,ti]
    
    
    elif status == GRB.Status.INF_OR_UNBD or \
        status == GRB.Status.INFEASIBLE  or \
        status == GRB.Status.UNBOUNDED:
        print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
        m.computeIIS() 
        m.write("GTCEP.ilp")
    
    t4 = time.time()
    ### Chekeos
    #1
    all_gen = np.array(all_gen)
    Datos_gen = [all_gen[pos_gen_agc_list], Cvar_gen[pos_gen_agc_list], Ccte_gen[pos_gen_agc_list], Pmin_gen[pos_gen_agc_list]*sim.Sb,Pmax_gen[pos_gen_agc_list]*sim.Sb,Ramp_gen[pos_gen_agc_list]*sim.Sb, Pgen_pre[pos_gen_agc_list,2]*sim.Sb, Pgen_pre[pos_gen_agc_list,1]*sim.Sb, Pgen_pre[pos_gen_agc_list,0]*sim.Sb, vg_inc.x, vg_dec.x]
    
    datos_gen = pd.DataFrame(np.array(Datos_gen).T, columns=['Gen', 'C_fuel', 'C_AGC', 'P min', 'P max', 'Rampa', 'P valle', 'P media', 'P peak', 'v_inc', 'v_dec'])
    datos_gen = datos_gen.sort_values(by = 'Gen', key = lambda col:col.str[3]+col.str[1])
    
    #2  
    all_P_out = np.zeros(Nt*Ns)
    for ti in range(Nt):
        all_P_out[ti*Ns:ti*Ns+Ns] = P_out[:,ti]
    dddf = np.vstack((range(1,Ns*Nt+1),np.tile(Gen_outages,Nt),all_P_out,np.tile(Barra_gen_out,Nt)))
    Res_escenarios = pd.DataFrame(dddf.T,index=range(1,Ns*Nt+1), columns=['Escenario', 'Gen', 'Potencia ', 'Posición'])
    
    
    #df_genOnAGC = pd.DataFrame(np.vstack((Gen_outages, Gen_on_AGC[ti,:,:])).T, columns=np.insert(all_gen,0,'Gen', axis=0))
    
    # %%
    with pd.ExcelWriter("0Resultados.xlsx") as writer:
        Res_escenarios.to_excel(writer, sheet_name='Escenarios - Gen OUT', index=False)
        datos_gen.to_excel(writer,sheet_name='DatosGenAGC', index = False)
    
        for ti in range(Nt):
            name_column=['Gen', 'Potencia']
            all_list= [all_gen]
            all_list.append(Pgen_pre[:,ti])
            for scen in range(Ns):
                all_list.append(Pgen_pfc[:,scen,ti])
                name_column.append('E' + str(scen + 1))
    
                delta_inc = np.zeros(ngen)
                delta_dec = np.zeros(ngen)
                delta_inc[pos_gen_agc_list] = pg_inc[:,scen,ti].x
                delta_dec[pos_gen_agc_list] = pg_dec[:,scen,ti].x
    
                scenarios = pd.DataFrame(np.vstack((all_gen,Pgen_pre[:,ti],Pgen_pfc[:,scen,ti], delta_inc - delta_dec,Pgen_pre[:,ti] + delta_inc ,Gen_on_AGC[ti,:,scen], T_Sfc * Ramp_gen)).T, 
                                     columns=['Gen', 'Potencia Previa', 'Potencia PFC', 'Aumento', 'Potencia SFC', 'Gen en AGC', 'Rampa'])
                name_hoja = 'D' + str(ti) + '_E ' + str(scen) + '--' + Gen_outages[scen]
                scenarios.to_excel(writer, sheet_name= name_hoja,index=False)
    
    
            powers_gen = pd.DataFrame(np.array(all_list).T,columns=name_column)
            powers_gen.to_excel(writer, sheet_name='Gen D' + str(ti+1), index=False)
        
        #df_genOnAGC.to_excel(writer, sheet_name='Generador en AGC', index=False)
    
    
    
    #PowerFactory
    results_line = np.zeros((len(dict_lineas),Ns,Nt))
    result_gen = np.zeros((ngen,Ns,Nt))
    for ti in range(Nt):
        sim.scens[ti].Activate()
        print('Escenario: ' + str(ti+1))
        for scen in range(Ns):
            for e in events_folder:
                e.outserv = 1   
            print('Progreso: ' + str(scen+1+Ns*ti) + '/' + str(Ns*Nt))
            print('Generador Out: ' + Gen_outages[scen])
            evt = events_folder[name_events.index('Salida Gen')]
            evt.outserv = 0
            evt.time = t_initial
            evt.p_target = sim.generadores[list(dict_gen).index(Gen_outages[scen])]
            count=0
            for gen in Gen_AGC:
                evt2 = events_folder[name_events.index('Evento Gamma ' + str(count+1))]
                evt2.outserv = 0
                evt2.value = str(part_factors[:,scen,ti][name_gen_agc_list.index(gen)])
                count += 1
            sim.prepare_dynamic_sim({}, 'rms', end_time=t_final)
            sim.run_dynamic_sim(end_sim=t_final)
    
            for gen in sim.generadores:
                if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                    result_gen[list(all_gen).index(gen.loc_name),scen,ti] = gen.GetAttribute('m:Psum:bus1')
    
            for line in sim.lineas:
                results_line[all_line.index(line.loc_name),scen,ti] = (line.GetAttribute('m:Psum:bus1') - line.GetAttribute('m:Psum:bus2'))/2
            #results_line[:,scen] = np.array(list(map(lambda x: (x.GetAttribute('m:Psum:bus1') + x.GetAttribute('m:Psum:bus2'))/2 , sim.lineas)))  
    
    with pd.ExcelWriter('1Comparaciones.xlsx') as writer:
        for ti in range(Nt):    
            for scen in range(Ns):
                delta_inc = np.zeros(ngen)
                delta_dec = np.zeros(ngen)
                delta_inc[pos_gen_agc_list] = pg_inc[:,scen,ti].x
                delta_dec[pos_gen_agc_list] = pg_dec[:,scen,ti].x
    
                P_final = Pgen_pre[:,ti] + delta_inc - delta_dec
                P_final[list(all_gen).index(Gen_outages[scen])] = 0
                Res_gen_export = pd.DataFrame(np.vstack((all_gen,result_gen[:,scen,ti],P_final*sim.Sb)).T, columns=['Gen', 'P_PF', 'P_Py'])
                name_hoja = 'D' + str(ti) + '_E' + str(scen) + ' ' + Gen_outages[scen]
                Res_gen_export.to_excel(writer, sheet_name= name_hoja,index=False)
    
    
                Res_line_export = pd.DataFrame(np.vstack((all_line,results_line[:,scen,ti],f[pos_line,scen,ti].x*sim.Sb, FMax[pos_line]*sim.Sb)).T, columns=['Linea', 'P_PF', 'P_Py', 'FMax'])
                Res_line_export.to_excel(writer, sheet_name= name_hoja,index=False, startcol = 6)
    
    escenariosPy_col = ['Py_peak', 'Py_media', 'Py_valle']
    escenariosPF_col = ['PF_peak', 'PF_media', 'PF_valle']
    
    with pd.ExcelWriter('3Casos.xlsx') as writer:
        for ti in range(Nt):
            PL_pre = PL_pre_line[ti]/sim.Sb/1000 #+ PL_pre_trafo2[:,s].sum()
            P_agc=np.zeros(Ns)
            for scen in range(Ns):
                P_agc[scen] = P_out[scen,ti] - dda_barra[:,ti].sum() + D_pfc[:,scen,ti].sum() - PL_pre + ploss[:,scen,ti].x.sum()
            diff_power = np.array(pg_inc.x[:,:,ti]-pg_dec.x[:,:,ti])*sim.Sb
            diff_power = pd.DataFrame(np.vstack((Gen_outages,diff_power,P_agc*sim.Sb,P_out[:,ti]*sim.Sb)).T, columns = [escenariosPy_col[ti]] + name_gen_agc_list + ['P_agc', 'P_out'])
            diff_power.to_excel(writer, index = False, startrow= ti*12)
    
    
            diff_gen_pf = result_gen[pos_gen_agc_list,:,ti] - np.tile(Pgen_pre[pos_gen_agc_list,ti]*sim.Sb,(Ns,1)).T
            results_pf = pd.DataFrame(np.vstack((Gen_outages,diff_gen_pf)).T, columns = [escenariosPF_col[ti]] + name_gen_agc_list)
            results_pf.to_excel(writer, index = False, startrow= 12*(3+ti))
    
    
    t5 = time.time()
    #Excel
    #dict_scen = dict()
    
    #for s in range(Ns):
    #    pos_part_gen_agc_list = pos_gen_agc_list[:]
    #    if Barra_gen_out[s] in pos_gen_agc_list:
    #        pos_part_gen_agc_list.remove(Barra_gen_out[s])
    #    y = -pg_inc[pos_part_gen_agc_list,s].x >= -T_Sfc * Ramp_gen[pos_part_gen_agc_list] - (Pgen_pfc[pos_part_gen_agc_list,s] - Pgen_pre[pos_part_gen_agc_list])
    #    if False in y:
    #        print('Escenario n°:' + str(s))
    #        print(y)
    
    
    
    #df = pd.DataFrame(part_factors).T
    #df.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\PartFactors.xlsx", index=False)




#results_line = np.vstack((np.array(indices_obj),f.x))[:,pos_line]

#sim.check_results(dict_gen, dict_genstat_eff, p_g, p_statg, results_line)

#cont=0
#for gen in dict_gen:
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

print('=> Check time: %.4f (s)' % (t5-t4))


print('finish!')
#  %% Comprobación de resultados



