import gurobipy as gp
from gurobipy import *

def Modelo():
    m = gp.Model('Modelo AGC')
    #m.setParam('DualReductions', 0)
    m.Params.MIPGap = 1e-5
    m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

    return m

def VarPotencia(m,n_gen_agc, Ns, Nt,Nb):
    pg_inc = m.addMVar((n_gen_agc, Ns, Nt), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_inc')
    pg_dec = m.addMVar((n_gen_agc, Ns, Nt), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg_dec')
    vg_inc = m.addMVar(n_gen_agc, vtype=GRB.BINARY, name='vg_inc')
    vg_dec = m.addMVar(n_gen_agc, vtype=GRB.BINARY, name='vg_dec')
    p_ens = m.addMVar((Nb,Ns, Nt), vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='P_ens')

    return (pg_inc,pg_dec,vg_inc,vg_dec,p_ens)

def Flujos(m,n_elem,Ns,Nt):
    f = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY, name='flujo') # Flujo por las líneas
    fp = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='fp') # Flujo-p por las líneas
    fn = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='fn') # Flujo-n por las líneas

    return (f, fp, fn)

def Perdidas(m,n_elem,Ns,Nt,L):
    ploss = m.addMVar((n_elem, Ns, Nt),vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='losses') # Flujo por las líneas
    dpk = m.addMVar((n_elem, L, Ns, Nt), vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo

    return (ploss,dpk)

def ObjFunction(m,):
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

    return f_obj