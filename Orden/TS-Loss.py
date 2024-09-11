import pfsimGeneral as pfsim
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import *

Sb = 100
Vb = 110
Zb = Vb**2/Sb

if False:
    z12 = complex(0.02, 0.1)
    z13 = complex(0.05, 0.1)
    z23 = complex(0.09, 0.3)
    z24 = complex(0.06, 0.15)
    z34 = complex(0.06, 0.15)

    z12 = complex(0.09, 0.3)
    z13 = complex(0.04, 0.2)
    z23 = complex(0.02, 0.05)


    lineas = [z12, z13, z23, z24, z34]
    name_lineas = ['z12', 'z13', 'z23', 'z24', 'z34']

    cont=0
    for i in lineas:
        val = i*Zb
        complex_val = round(val.real, 3) + round(val.imag, 3) * 1j
        print(name_lineas[cont] + ' = ' + str(complex_val))
        cont+=1


#pf = pfsim.PowerFactorySim('TS-Loss', False)
pf = pfsim.PowerFactorySim('3Bus_TS', False)


bin_comp = True
bin_ady = True


pf.flujo_dc = 0
pf.TS = 1
TsLoss = 1
costo_ts = 0
LOSS = True


M = 1e5
L = 2
pf.Nt = 1
pf.export_csv('AC')
pf.use_sf = 1


dict_linea = pf.get_data()

m = gp.Model('Modelo TS-Loss')
#m.setParam('DualReductions', 0)
m.Params.MIPGap = 1e-10
m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi
m.Params.IntegralityFocus = 1


pg = m.addMVar((pf.ngen),vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='Pg')
f = m.addMVar((pf.n_elem),vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=-GRB.INFINITY, name='f')
fp = m.addMVar((pf.n_elem),vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='fp')
fn = m.addMVar((pf.n_elem),vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=0, name='fn')
ploss = m.addMVar((pf.n_elem), vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='ploss')   # Perdidas por cada línea
dpk = m.addMVar((pf.n_elem,L), vtype=GRB.CONTINUOUS, lb=0, ub=GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo
n_l = m.addMVar(pf.n_elem, vtype=GRB.BINARY, name='n_l')                                          # variable binaria complentaridad
n_a = m.addMVar((pf.n_elem,L), vtype=GRB.BINARY, name='n_a')                                          # variable binaria adyacencia
s_ts = m.addMVar(len(pf.pos_ts), vtype=GRB.BINARY, name='s_ts')                                    # variable binaria de TS
f_ts = m.addMVar(len(pf.pos_ts),vtype=GRB.CONTINUOUS, ub=GRB.INFINITY, lb=-GRB.INFINITY, name='f')



Cop = pf.Ccte_gen @ pg * pf.Sb
Cts = costo_ts * (1 - s_ts).sum()
f_obj = Cop + Cts

m.setObjective(f_obj,gp.GRB.MINIMIZE)

m.addConstr(pg.sum() == pf.dda_barra.sum() + ploss.sum(), name = 'Balance')

m.addConstr(-pg >= -pf.Pmax_gen, name = 'PMax')
m.addConstr(pg >= pf.Pmin_gen, name = 'PMin')


if pf.TS and not TsLoss:
    #Líneas NO Candidatas TS
    fe = pf.SF[pf.pos_nots,:][:,pf.pos_gen] @ pg - pf.SF[pf.pos_nots,:]@pf.dda_barra[:,0]
    fv = (pf.SF[pf.pos_nots,:] @ pf.A[pf.pos_ts,:].T) @ f_ts

    m.addConstr(-(fe+fv) >= -pf.FMax[pf.pos_nots], name = 'fe_p')
    m.addConstr(fe+fv >= -pf.FMax[pf.pos_nots], name = 'fe_n')

    #Líneas candidatas
    f1 = pf.SF[pf.pos_ts,:][:,pf.pos_gen] @ pg - pf.SF[pf.pos_ts,:]@pf.dda_barra[:,0]
    f2 = f_ts - (pf.SF[pf.pos_ts,:] @ pf.A[pf.pos_ts,:].T) @ f_ts

    m.addConstr(f1-f2 <= (pf.FMax[pf.pos_ts]) * s_ts, name = 'fs1_p') # 1
    m.addConstr(f1-f2 >= -(pf.FMax[pf.pos_ts]) * s_ts, name = 'fs1_n')

    m.addConstr(f_ts <= M*(1 - s_ts), name = 'fs2_p') # 2
    m.addConstr(f_ts >= -M*(1 - s_ts), name = 'fs2_n') # 2

elif TsLoss:

    kl = np.zeros((pf.n_elem, L))
    for l in range(L):
        kl[:,l] = (2*(l+1)-1)*(pf.FMax)/L
    
    m.addConstr(ploss == pf.G/(pf.B**2)*(gp.quicksum(kl[:,i]*dpk[:,i] for i in range(L))), name = 'Ploss')  


    f_loss_nc = 0.5 * pf.SF[pf.pos_nots,:] @ abs(pf.A[pf.pos_nots,:].T) @ ploss[pf.pos_nots]
    f_loss_c = 0.5 * pf.SF[pf.pos_ts,:] @ abs(pf.A[pf.pos_ts,:].T) @ ploss[pf.pos_ts]

    #Líneas NO Candidatas TS
    fe = pf.SF[pf.pos_nots,:][:,pf.pos_gen] @ pg - pf.SF[pf.pos_nots,:]@pf.dda_barra[:,0]
    fv = (pf.SF[pf.pos_nots,:] @ pf.A[pf.pos_ts,:].T) @ f_ts

    m.addConstr(-(fe-f_loss_nc+fv) >= -pf.FMax[pf.pos_nots], name = 'fe_p')
    m.addConstr(fe-f_loss_nc+fv>= -pf.FMax[pf.pos_nots], name = 'fe_n')

    #Líneas candidatas
    f1 = pf.SF[pf.pos_ts,:][:,pf.pos_gen] @ pg - pf.SF[pf.pos_ts,:]@pf.dda_barra[:,0]
    f2 = f_ts - (pf.SF[pf.pos_ts,:] @ pf.A[pf.pos_ts,:].T) @ f_ts

    m.addConstr(f1-f_loss_c-f2 <= (pf.FMax[pf.pos_ts]) * s_ts, name = 'fs1_p') # 1
    m.addConstr(f1-f_loss_c-f2 >= -(pf.FMax[pf.pos_ts]) * s_ts, name = 'fs1_n')

    m.addConstr(f_ts <= M*(1 - s_ts), name = 'fs2_p') # 2
    m.addConstr(f_ts >= -M*(1 - s_ts), name = 'fs2_n') # 2


    m.addConstr(dpk.sum(1) == fp + fn, 'dpk')

    m.addConstr(f == fp - fn, name = 'flujo')

    m.addConstr(f[pf.pos_ts] == f1-f_loss_c-f2)
    m.addConstr(f[pf.pos_nots] == fe-f_loss_nc+fv)

    m.addConstr(-fp - fn- 0.5*ploss >= -pf.FMax, name = 'suma_f')

    if bin_comp:
        m.addConstr(-fp >= -n_l*pf.FMax, name='fp_n')   #flujo positvo-restriccion de complementaridad
        m.addConstr(-fn >= (1-n_l)*(-pf.FMax), name='fn_n') #flujo nefativo-restriccion de complementaridad

    if bin_ady:
        for l in range(L): 
            if l==0:
                m.addConstr(-dpk[:,l]>= -pf.FMax/L, name='d_f_Res_max_A_l')
                m.addConstr(dpk[:,l] >= n_a[:,l]*(pf.FMax/L), name='d_f_Res_min_A_l')
            elif l==L-1:
                m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*pf.FMax/L, name='d_f_Res_max_A_L')
                m.addConstr(dpk[:,l] >=0, name='d_f_Res_min_A_L')
            else:
                m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*pf.FMax/L, name='d_f_Res_max_A_L-1')
                m.addConstr(dpk[:,l] >=n_a[:,l]*pf.FMax/L, name='d_f_Res_min_A_L-1')


    #m.addConstr(-f >= -pf.FMax, name = 'f+')
    #m.addConstr(f >= -pf.FMax, name = 'f-')






else:
    f_gen = pf.SF[:,pf.pos_gen] @ pg
    f_loss = 0.5 * pf.SF @ abs(pf.A.T) @ ploss
    f_dda = pf.SF @ pf.dda_barra

    if LOSS:
        m.addConstr(f == f_gen - f_dda[:,0] - f_loss, name = 'flujoSF')
    else:
        m.addConstr(f == f_gen - f_dda[:,0], name = 'flujoSF')

    m.addConstr(f == fp - fn, name = 'flujo')


    m.addConstr(-f >= -pf.FMax, name = 'f+')
    m.addConstr(f >= -pf.FMax, name = 'f-')

    m.addConstr(dpk.sum(1) == fp + fn, 'dpk')

    if LOSS:
        kl = np.zeros((pf.n_elem, L))
        for l in range(L):
            kl[:,l] = (2*(l+1)-1)*(pf.FMax)/L
        m.addConstr(ploss == pf.G/(pf.B**2)*(gp.quicksum(kl[:,i]*dpk[:,i] for i in range(L))), name = 'Ploss')  

    m.addConstr(-fp - fn- 0.5*ploss >= -pf.FMax, name = 'suma_f')
    #m.addConstr(-f - 0.5*ploss >= -pf.FMax)
    #m.addConstr(f - 0.5*ploss >= -pf.FMax)


    m.addConstr(-fp >= -pf.FMax, name = 'fp+')
    m.addConstr(-fn >= -pf.FMax, name = 'fn+')
    for l in range(L):
        m.addConstr(-dpk[:,l] >= -pf.FMax/L, name = 'LimiteDpk')

    if bin_comp:
        m.addConstr(-fp >= -n_l*pf.FMax, name='fp_n')   #flujo positvo-restriccion de complementaridad
        m.addConstr(-fn >= (1-n_l)*(-pf.FMax), name='fn_n') #flujo nefativo-restriccion de complementaridad

    if bin_ady:
        for l in range(L): 
            if l==0:
                m.addConstr(-dpk[:,l]>= -pf.FMax/L, name='d_f_Res_max_A_l')
                m.addConstr(dpk[:,l] >= n_a[:,l]*(pf.FMax/L), name='d_f_Res_min_A_l')
            elif l==L-1:
                m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*pf.FMax/L, name='d_f_Res_max_A_L')
                m.addConstr(dpk[:,l] >=0, name='d_f_Res_min_A_L')
            else:
                m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*pf.FMax/L, name='d_f_Res_max_A_L-1')
                m.addConstr(dpk[:,l] >=n_a[:,l]*pf.FMax/L, name='d_f_Res_min_A_L-1')



    


m.write('OPF-TSLoss.lp')


m.optimize()

status = m.Status
if status == gp.GRB.Status.OPTIMAL:
    print('-----------------------------')
    #print('La demanda total del sistema es: %.2f (MW)' % (pf.dda_barra[:,ti].sum()*pf.Sb))
            
    print ('Cost = %.2f ($)' % (m.objVal))
            
    print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (m.NumVars,m.NumConstrs,m.DNumNZs)) #print('num_Vars =  %d / num_Const =  %d' % (len(m.getVars()), len(m.getConstrs())))      
    #print ('Total P_loss = %.2f [MW]'%(pf.Sb*ploss.sum().getValue()))
    print('=> Solver time: %.4f (s)' % (m.Runtime))
    for flujo, name_flujo in zip(f.x, pf.all_line):
        print(name_flujo + ' = ' + str(flujo))
    
    for val_gen, name_gen in zip(pg.x, pf.all_gen):
        print(name_gen + ' = ' + str(val_gen))


    #print ('Costo => %.2f ($/h)' % m.objVal) 
    #print ('Las perdidas son %.4f (MW)' % ploss.x.sum()*100) 

elif status == gp.GRB.Status.INF_OR_UNBD or \
    status == gp.GRB.Status.INFEASIBLE  or \
    status == gp.GRB.Status.UNBOUNDED:
    print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
    m.computeIIS() 
    m.write("GTCEP-TSLoss.ilp")


print('end!')