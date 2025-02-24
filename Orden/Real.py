import sys
import pfsim
import pandas as pd
import numpy as np
import time
import gurobipy as gp
from gurobipy import *
import optm

#pf = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#pf = pfsim.PowerFactorySim('Ejemplo Clase')
#pf = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')

Sim_PF = False
flujo_dc = 0

Estudio = 6
Nt = 1 # N° de escenarios
TS = False
Costo_ts = 0

if Estudio == 0: # Flujo DC
    VOLL = 2500
    Pot_Down = 0
    Flujos = 1
    Perdidas = 0
    PerdidasPlus = 0
    flujo_dc = 1
    Sim_PF = True



if Estudio == 1: # Sin sistema de transmisión y VOLL = 150000
    VOLL = 2500
    Pot_Down = 0
    Flujos = 0
    Perdidas = 0
    PerdidasPlus = 0

elif Estudio == 2:# Sistema de transmisión sin perdidas y con efecto en Voltage
    VOLL = 2500
    Pot_Down = 0
    Flujos = 1
    Perdidas = 0
    PerdidasPlus = 0

elif Estudio == 3:# Sistema de transmisión con perdidas
    VOLL = 300
    Pot_Down = 0
    Flujos = 1
    Perdidas = 1
    PerdidasPlus = 1

elif Estudio == 4: # No se utilizan las condiciones de complementaridad / exclusividad
    VOLL = 2500
    Pot_Down = 0
    Flujos = 1
    Perdidas = 1
    PerdidasPlus = 0

elif Estudio == 5: # Se puede bajar la potencia de los generadores
    VOLL = 300
    Pot_Down = 1
    Flujos = 1
    Perdidas = 1
    PerdidasPlus = 1

elif Estudio == 6: # Transmission Switching con perdidas
    VOLL = 300
    Pot_Down = 1
    Flujos = 1
    Perdidas = 1
    PerdidasPlus = 1

    Costo_ts = 10
    TS = 1


t0 = time.time()
#project_file = '3Bus_TS'
project_file = 'Articulo2'
pf = pfsim.PowerFactorySim(project_file)

#Gen_AGC = ['G2_1', 'G3_1', 'G4_1', 'G7_1', 'G2_3', 'G4_3', 'G6_3', 'G3_2', 'G4_2', 'G2_4', 'G3_4', 'G4_4', 'G6_4']

#Gen_outages = ['G4_1', 'G5_1', 'G6_1', 'G2_2', 'G5_2']
#Gen_outages = ['G4_1', 'G5_1', 'G6_1', 'G2_2', 'G5_2', 'G5_3', 'G3_4', 'G5_4', 'G7_4']

#lines = ['Line3_1', 'Line3_3', 'Line3_4']
#val = [600, 725, 590.2]

#pf.ChangeMaxLine(lines, val)


escenarios = ['peak', 'media', 'valle']

pf.TS = TS

L=11
#T_Sfc = 15 # 15 Min, <-- Tiempo de operación AGC

pf.Nt = Nt # N° de escenarios
if Estudio == 0:
    pf.export_csv('DC')
    pf.flujo_dc = flujo_dc
else:
    pf.export_csv('AC')
#ptdf_dataframe = pd.read_excel(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.xlsx',skiprows = 0,delimiter=';')
#ptdf_dataframe = pd.read_excel(r'C:\Users\lldie\Desktop\SF_datos.xlsx',skiprows = 0,delimiter=';')
#SF, indices_obj = pfsim.ShiftFactors()



pf.use_sf = 1 # 1 = SF calculados con matrices;  0 = SF de PF

dict_lineas = pf.get_data()


if False:
    # demanda por barra
    Nb = len(indices_bus)
    pf.dda_barra = np.zeros((Nb, Nt))
    for load in dict_cargas:
        for ti in range(Nt):
            pf.dda_barra[dict_cargas[load][0], ti] += dict_cargas[load][ti+2]/pf.Sb

    #Flujo_dda = SF @ pf.dda_barra   # Flujo en elemntos debido a la demanda en barra





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
        Pmin_gen[cont] = dict_gen[gen][3]/pf.Sb
        Pmax_gen[cont] = dict_gen[gen][4]/pf.Sb
        Cvar_gen[cont] = dict_gen[gen][5]
        Ccte_gen[cont] = dict_gen[gen][6]
        Ramp_gen[cont] = dict_gen[gen][7]/pf.Sb
        for ti in range(Nt):
            Pgen_pre[cont,ti] = dict_gen[gen][8+ti]/pf.Sb
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
        Pmin_genstat[cont] = dict_genstat_eff[gen][3]/pf.Sb
        Pmax_genstat[cont] = dict_genstat_eff[gen][4]/pf.Sb
        Pgenstat_pre[cont] = dict_genstat_eff[gen][5]/pf.Sb
        Cvar_genstat[cont] = dict_genstat_eff[gen][6]
        Ccte_genstat[cont] = dict_genstat_eff[gen][7]
        Ramp_genstat[cont] = dict_genstat_eff[gen][8]
        if gen in Gen_AGC:
            Genstat_on_AGC[cont] = 1
            pos_genstat_agc_list.append(cont)
        cont+=1

    n_genstat_agc = len(pos_genstat_agc_list)

    #P_pre = Pgen_pre.sum() + Pgenstat_pre.sum()

    n_elem = len(pf.indices_obj)   # N° de líneas + 'trf2'
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
    for i in pf.indices_obj:
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
        FMax[cont] = Fmax_i/pf.Sb
        R[cont] = R_i
        X[cont] = X_i
        i_buses[cont] = i_bus
        j_buses[cont] = j_bus

    pos_nots = [pos for pos in pos_line if pos not in pos_ts]
    # Lineas Candidatas TS
    i_buses_ts = i_buses[pos_ts]
    j_buses_ts = j_buses[pos_ts]
    
    I = np.r_[range(n_elem), range(n_elem)]
    S = sparse((np.r_[np.ones(n_elem), -np.ones(n_elem)], (I, np.r_[i_buses, j_buses])), (n_elem, pf.Nb))
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
    noslack = np.flatnonzero(np.arange(pf.Nb) != n_bus_slack)


    BfR = sparse((np.r_[np.imag(yprim), -np.imag(yprim)], (I, np.r_[i_buses, j_buses])), (n_elem,Nb))
    BbusR = S.T * BfR
    SFR = np.zeros((n_elem,Nb))
    SFR[:,noslack] = BfR[:, noslack].todense()*np.linalg.inv(BbusR[np.ix_(noslack, noslack)].todense())    

    difference = (SFR - SF)/SFR

    SF_PF = pd.DataFrame(np.vstack((list(dict_barras.keys()),SF)).T, columns=np.insert(indices_obj, 0, 'x', axis=0))
    SF_DC_Lossy = pd.DataFrame(np.vstack((list(dict_barras.keys()),SFR)).T, columns=np.insert(indices_obj, 0, 'x', axis=0))
    diff = pd.DataFrame(np.vstack((list(dict_barras.keys()),difference)).T, columns=np.insert(indices_obj, 0, 'x', axis=0))

    SF= SFR

    #Excel con Shift Factors
    with pd.ExcelWriter('3ShiftFactors.xlsx') as writer:
        SF_PF.to_excel(writer, sheet_name='SF_PF', index=False)
        SF_DC_Lossy.to_excel(writer, sheet_name='SF_DC', index=False)
        diff.to_excel(writer, sheet_name='Diff', index=False)


#df = pd.DataFrame()

# Dynamic Simulation until 50 seg - Finishing CPF
t_initial = 0 # Tiempo Perturbación
#t_final = 420 # Tiempo total
t_final = 200 + t_initial # Tiempo total

tstop_cpf = 30 + t_initial# Tiempo hasta el CPF


sim = pfsim.Simulacion(pf, t_initial, tstop_cpf)

t1 = time.time()

m_optm = optm.Modelo()
m_optm.pot_down = Pot_Down # 1 = Aumento y disminución de potencia / 0 = Solo aumento de potencia
m_optm.flujos = Flujos # 1 = OPF / 0 = ED
m_optm.losses = Perdidas # 1 = Perdidas / 0 = Sin perdidas
m_optm.losses_plus = PerdidasPlus # 1 = Añadir Complementaridad / Exclusividad
m_optm.Voll = VOLL # VOLL
m_optm.TS = TS
m_optm.costo_ts = Costo_ts

m_optm(pf, sim)

t2 = time.time()

m_optm.run()

t3=time.time()

m_optm.Results(pf, sim)

print('=> Extract data time: %.4f (s)'% (t1-t0))
print('=> Formulation time: %.4f (s)' % (t2-t1))
print('=> Solution time: %.4f (s)' % (t3-t2))


pf.app.ResetCalculation()

# Tiempo entre cada simulación
t_step = 4
n_step = int((t_final-tstop_cpf)/t_step  +1)


part_factors = np.zeros((pf.n_gen_agc,n_step,pf.Nt))

# Variabilidad
for ti in range(pf.Nt): 
    pf.scens[ti].Activate()
    scen = -1
    for gen_out in pf.Gen_Outages:
        scen += 1
        
        #Desactivar todos los eventos
        events_folder = pf.IntEvt.GetContents()
        for e in events_folder:
                e.outserv = 1
        
        ### NUEVO SF Con TRANSMISSION SWITCHING
        new_SF = pfsim.new_SF(pf, m_optm, ti, scen)
                
        ### Eventos de desconexiones de líneas al inicio del CSF
        pfsim.CreateEvents_line(pf, m_optm, scen, ti)

        events_folder = pf.IntEvt.GetContents()

        evt = events_folder[sim.name_events.index('Salida Gen')]
        evt.outserv = 0
        evt.time = t_initial
        evt.p_target = pf.generadores[list(pf.dict_gen).index(gen_out)]


        P_rec = 0
        for t_int in range(tstop_cpf, t_final, t_step):
            if t_int == tstop_cpf:
                P_before = sim.Pgen_pfc[:,scen,ti].sum() - sim.D_pfc[:,scen,ti].sum()
            
            ShortSim = pfsim.ShortSim(pf, m_optm, sim, gen_out, t_int)
            
            P_rec += ShortSim.p_gen.sum() - ShortSim.D_t.sum() - P_before
            P_before = ShortSim.p_gen.sum() - ShortSim.D_t.sum()
            
            PartialModel = optm.PartialModel(pf,m_optm,sim,ShortSim,new_SF, scen, ti)

            part_factors[:,int((t_int-tstop_cpf)/t_step),ti] = PartialModel.part_factors


            cont=0
            for gen in pf.gen_agc:
                name_event =  'Evt Gamma ' + str(gen.loc_name) + ' - ' + str(t_int)
                if not name_event in pf.IntEvt.GetContents() and float(PartialModel.part_factors[cont]) != 0.0:
                    evento = pf.IntEvt.CreateObject('EvtParam', 'Evt Gamma ' + str(gen.loc_name) + ' - ' + str(t_int))    
                    evento.variable = pf.signal_list[cont]
                    evento.time = t_int
                    evento.p_target = pf.dsl_agc_bloques
                    evento.value = str(PartialModel.part_factors[cont])
                cont+=1

                #PartialModel.part_factors[gen] 
            


if True:
    t4 = time.time()
    ### Chekeos
    #1
    all_gen = np.array(pf.all_gen)
    PgenPre_all = list()
    NamePPre_all =  ['P valle', 'P media', 'P peak']
    NamePPre_all = NamePPre_all[0:Nt]
    for ti in range(Nt):
        PgenPre_all.append(pf.Pgen_pre[pf.pos_gen_agc_list,ti]*pf.Sb)
    
    if m_optm.pot_down:
        Datos_gen = [all_gen[pf.pos_gen_agc_list], pf.Cvar_gen[pf.pos_gen_agc_list], pf.Ccte_gen[pf.pos_gen_agc_list], pf.Pmin_gen[pf.pos_gen_agc_list]*pf.Sb,pf.Pmax_gen[pf.pos_gen_agc_list]*pf.Sb,pf.Ramp_gen[pf.pos_gen_agc_list]*pf.Sb] + PgenPre_all + [m_optm.vg_inc.x, m_optm.vg_dec.x]
        datos_gen = pd.DataFrame(np.array(Datos_gen).T, columns=['Gen', 'C_fuel', 'C_AGC', 'P min', 'P max', 'Rampa'] + NamePPre_all + ['v_inc', 'v_dec'])
        datos_gen = datos_gen.sort_values(by = 'Gen', key = lambda col:col.str[3]+col.str[1])

    else:
        Datos_gen = [all_gen[pf.pos_gen_agc_list], pf.Cvar_gen[pf.pos_gen_agc_list], pf.Ccte_gen[pf.pos_gen_agc_list], pf.Pmin_gen[pf.pos_gen_agc_list]*pf.Sb,pf.Pmax_gen[pf.pos_gen_agc_list]*pf.Sb,pf.Ramp_gen[pf.pos_gen_agc_list]*pf.Sb] + PgenPre_all + [m_optm.vg_inc.x]
        datos_gen = pd.DataFrame(np.array(Datos_gen).T, columns=['Gen', 'C_fuel', 'C_AGC', 'P min', 'P max', 'Rampa'] + NamePPre_all + ['v_inc'])
        datos_gen = datos_gen.sort_values(by = 'Gen', key = lambda col:col.str[3]+col.str[1])

    
    
    #2  
    all_P_out = np.zeros(Nt*pf.Ns)
    for ti in range(Nt):
        all_P_out[ti*pf.Ns:ti*pf.Ns+pf.Ns] = sim.P_out[:,ti]
    dddf = np.vstack((range(1,pf.Ns*Nt+1),np.tile(pf.Gen_Outages,Nt),all_P_out,np.tile(sim.Barra_gen_out,Nt)))
    Res_escenarios = pd.DataFrame(dddf.T,index=range(1,pf.Ns*Nt+1), columns=['Escenario', 'Gen', 'Potencia ', 'Posición'])
    
    
    #df_genOnAGC = pd.DataFrame(np.vstack((Gen_outages, Gen_on_AGC[ti,:,:])).T, columns=np.insert(all_gen,0,'Gen', axis=0))
    
    # %%
    with pd.ExcelWriter('00Resultados'+ project_file+ '.xlsx') as writer:
        Res_escenarios.to_excel(writer, sheet_name='Escenarios - Gen OUT', index=False)
        datos_gen.to_excel(writer,sheet_name='DatosGenAGC', index = False)
    
        for ti in range(Nt):
            name_column=['Gen', 'Potencia']
            all_list= [all_gen]
            all_list.append(pf.Pgen_pre[:,ti])
            for scen in range(pf.Ns):
                all_list.append(sim.Pgen_pfc[:,scen,ti])
                name_column.append('E' + str(scen + 1))
    
                delta_inc = np.zeros(pf.ngen)
                delta_inc[pf.pos_gen_agc_list] = m_optm.pg_inc[:,scen,ti].x
                
                if m_optm.pot_down:
                    delta_dec = np.zeros(pf.ngen)
                    delta_dec[pf.pos_gen_agc_list] = m_optm.pg_dec[:,scen,ti].x
                    delta = delta_inc - delta_dec
                
                else:
                    delta = delta_inc
    
                scenarios = pd.DataFrame(np.vstack((all_gen,pf.Pgen_pre[:,ti],sim.Pgen_pfc[:,scen,ti], delta,pf.Pgen_pre[:,ti] + delta_inc ,sim.Gen_on_AGC[ti,:,scen], m_optm.T_Sfc * pf.Ramp_gen)).T, 
                                     columns=['Gen', 'Potencia Previa', 'Potencia PFC', 'Aumento', 'Potencia SFC', 'Gen en AGC', 'Rampa'])
                name_hoja = 'D' + str(ti) + '_E ' + str(scen) + '--' + pf.Gen_Outages[scen]
                scenarios.to_excel(writer, sheet_name= name_hoja,index=False)
    
    
            powers_gen = pd.DataFrame(np.array(all_list).T,columns=name_column)
            powers_gen.to_excel(writer, sheet_name='Gen D' + str(ti+1), index=False)
        
        #df_genOnAGC.to_excel(writer, sheet_name='Generador en AGC', index=False)
    
    


    
    #PowerFactor
    results_line = np.zeros((len(dict_lineas),pf.Ns,Nt))
    result_gen = np.zeros((pf.ngen,pf.Ns,Nt))
    if Sim_PF:
        if flujo_dc:
            for ti in range(Nt):
                pf.scens[ti].Activate()
                print('Escenario: ' + str(ti+1))
                s=0
                for gen_out in pf.Gen_Outages:
                    for gen in pf.generadores:
                        if gen_out == gen.loc_name:
                            print('Gen Out : ' + gen.loc_name)
                            s_gen_out = gen
                            s_gen_out.outserv = 1
                            original_power = np.zeros((len(pf.Gen_AGC)))
                            cont= 0
                            for gen_agc in pf.name_gen_agc_list:
                                for gen_change in pf.generadores:
                                    if gen_agc == gen_change.loc_name:
                                        original_power[cont] = gen_change.pgini
                                        gen_change.pgini += pf.Sb*m_optm.pg_inc[cont, s,ti].x
                                        cont += 1
                            
                            pf.ldf.Execute()



                            for gen in pf.generadores:
                                if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                                    result_gen[list(all_gen).index(gen.loc_name),s,ti] = gen.GetAttribute('m:P:bus1')

                            for line in pf.lineas:
                                results_line[pf.all_line.index(line.loc_name),s,ti] = (line.GetAttribute('m:P:bus1') - line.GetAttribute('m:P:bus2'))/2
                            s += 1
                            print(original_power)
                            s_gen_out.outserv = 0
                            cont=0
                            for gen_agc in pf.Gen_AGC:
                                for gen_change in pf.generadores:
                                    if gen_agc == gen_change.loc_name:
                                        gen_change.pgini = original_power[cont]
                                        cont+=1    


        else:
            
            signal_list = list()
            for gen in pf.gen_agc:
                comp = gen.c_pmod
                for slot in comp.pblk:
                    ## Slot pcu
                    if slot.loc_name[:9] == 'pcu Slot' or slot.loc_name[:9] == 'Gov Slot':
                        val = comp.pelm[comp.pblk.index(slot)].signal[3][-2:]
                        if not val[0].isdigit():
                            val = val[1]
                        signal_list.append('Evento Gamma ' + val)


            for ti in range(Nt):
                pf.scens[ti].Activate()
                print('Escenario: ' + str(ti+1))
                for scen in range(pf.Ns):
                    for e in pf.events_folder:
                        e.outserv = 1   
                    print('Progreso: ' + str(scen+1+pf.Ns*ti) + '/' + str(pf.Ns*Nt))
                    print('Generador Out: ' + pf.Gen_Outages[scen])
                    evt = pf.events_folder[sim.name_events.index('Salida Gen')]
                    evt.outserv = 0
                    evt.time = t_initial
                    evt.p_target = pf.generadores[list(pf.dict_gen).index(pf.Gen_Outages[scen])]
                    count=0
                    for gen in pf.Gen_AGC:
                        evt2 = pf.events_folder[sim.name_events.index(signal_list[count])]
                        evt2.outserv = 0
                        evt2.value = str(m_optm.part_factors[:,scen,ti][pf.name_gen_agc_list.index(gen)])
                        count += 1
                    cont = 0
                    line_sim_out = list()
                    for line in pf.lineas:
                        if line.loc_name in pf.TS_lines:
                            if m_optm.s_ts.x[cont,scen,ti] == 0:
                                line.outserv = 1
                                line_sim_out.append(line)
                        cont += 1

                    pf.prepare_dynamic_sim({}, 'rms')
                    pf.run_dynamic_sim(end_sim=t_final)

                    for gen in pf.generadores:
                        if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                            result_gen[list(all_gen).index(gen.loc_name),scen,ti] = gen.GetAttribute('m:Psum:bus1')

                    for line in pf.lineas:
                        if not line in line_sim_out:
                            results_line[pf.all_line.index(line.loc_name),scen,ti] = (line.GetAttribute('m:Psum:bus1') - line.GetAttribute('m:Psum:bus2'))/2
                    #results_line[:,scen] = np.array(list(map(lambda x: (x.GetAttribute('m:Psum:bus1') + x.GetAttribute('m:Psum:bus2'))/2 , pf.lineas)))  
    

                    for line in line_sim_out:
                        line.outserv = 0

    with pd.ExcelWriter('01Comparaciones'+ project_file +'.xlsx') as writer:
        for ti in range(Nt):    
            for scen in range(pf.Ns):
                delta_inc = np.zeros(pf.ngen)
                delta_inc[pf.pos_gen_agc_list] = m_optm.pg_inc[:,scen,ti].x
                if m_optm.pot_down:
                    delta_dec = np.zeros(pf.ngen)
                    delta_dec[pf.pos_gen_agc_list] = m_optm.pg_dec[:,scen,ti].x
                    delta = delta_inc - delta_dec
                else:
                    delta = delta_inc
    
                P_final = pf.Pgen_pre[:,ti] + delta
                P_final[list(all_gen).index(pf.Gen_Outages[scen])] = 0
                Res_gen_export = pd.DataFrame(np.vstack((all_gen,result_gen[:,scen,ti],P_final*pf.Sb)).T, columns=['Gen', 'P_PF', 'P_Py'])
                name_hoja = 'D' + str(ti) + '_E' + str(scen) + ' ' + pf.Gen_Outages[scen]
                Res_gen_export.to_excel(writer, sheet_name= name_hoja,index=False)
    
                if m_optm.flujos:
                    Res_line_export = pd.DataFrame(np.vstack((pf.all_line,results_line[:,scen,ti],m_optm.f[pf.pos_line,scen,ti].x*pf.Sb, pf.FMax[pf.pos_line]*pf.Sb)).T, columns=['Linea', 'P_PF', 'P_Py', 'FMax'])
                    Res_line_export.to_excel(writer, sheet_name= name_hoja,index=False, startcol = 6)
    
    escenariosPy_col = ['Py_peak', 'Py_media', 'Py_valle']
    escenariosPF_col = ['PF_peak', 'PF_media', 'PF_valle']
    
    with pd.ExcelWriter('03Casos_DC' + project_file +'.xlsx') as writer:
        for ti in range(Nt):
            P_agc = np.zeros(pf.Ns)
            costo = np.zeros(pf.Ns)
            for scen in range(pf.Ns):
                #j=0
                #for gen in pf.name_gen_agc_list:
                #    print(gen + ' = ' + str(m_optm.pg_inc.x[j,scen,ti]*100))
                #    j+=1
                
                #j=0
                #for line in pf.all_line:
                #    if line in pf.TS_lines:
                #        if  m_optm.s_ts.x[j,scen,ti] != 1 and m_optm.s_ts.x[j,scen,ti] != 0:
                #            print(line + ' = ' +str(m_optm.s_ts.x[j,scen,ti]))
                #        j+=1

                if m_optm.pot_down:
                    costo[scen] = pf.Cvar_gen[pf.pos_gen_agc_list]*pf.Sb @ (m_optm.pg_inc[:,scen,ti].x + m_optm.pg_dec[:,scen,ti].x)
                else:
                    costo[scen] = pf.Cvar_gen[pf.pos_gen_agc_list]*pf.Sb @ (m_optm.pg_inc[:,scen,ti].x)
                
                if m_optm.losses:
                    P_agc[scen] = sim.P_out[scen,ti] - pf.dda_barra[:,ti].sum() + sim.D_pfc[:,scen,ti].sum() - sim.PL_pre_line[ti] + m_optm.ploss[:,scen,ti].x.sum()
                elif m_optm.flujos and not pf.flujo_dc:
                    P_agc[scen] = sim.P_out[scen,ti] - pf.dda_barra[:,ti].sum() + sim.D_pfc[:,scen,ti].sum()
                else:
                    P_agc[scen] = sim.P_out[scen,ti] 

            if m_optm.pot_down:
                diff_power = np.array(m_optm.pg_inc.x[:,:,ti]-m_optm.pg_dec.x[:,:,ti])*pf.Sb
            else:
                diff_power = np.array(m_optm.pg_inc.x[:,:,ti])*pf.Sb
            
            columns_diffpower = [escenariosPy_col[ti]] + pf.name_gen_agc_list + ['P_agc', 'P_out', 'Costo']
            diff_power = pd.DataFrame(np.vstack((pf.Gen_Outages,diff_power,P_agc*pf.Sb,sim.P_out[:,ti]*pf.Sb, costo)).T, columns = columns_diffpower)
            diff_power.to_excel(writer, index = False, sheet_name='Case', startrow= ti*12)
    
    
            diff_gen_pf = result_gen[pf.pos_gen_agc_list,:,ti] - np.tile(pf.Pgen_pre[pf.pos_gen_agc_list,ti]*pf.Sb,(pf.Ns,1)).T
            results_pf = pd.DataFrame(np.vstack((pf.Gen_Outages,diff_gen_pf)).T, columns = [escenariosPF_col[ti]] + pf.name_gen_agc_list)
            results_pf.to_excel(writer, index = False, sheet_name='Case', startrow= 12*(3+ti))
    
            diff_line = pd.DataFrame(np.vstack((m_optm.f.x[:,:,ti])).T, columns = list(pf.indices_obj))
            diff_line.to_excel(writer, index = False, sheet_name='Case', startrow= ti*12, startcol=len(columns_diffpower)+2)

    
    t5 = time.time()



print('=> Check time: %.4f (s)' % (t5-t4))

if m_optm.m.SolCount > 1:
    m_optm.m.setParam(m_optm.gp.GRB.Param.SolutionNumber,1)
    print('Another Solution: ' + str(m_optm.m.PoolObjVal))

    #for e in range(m_optm.m.SolCount):
    #    m_optm.m.setParam(m_optm.gp.GRB.Param.SolutionNumber,e)
        


print('finish!')
#  %% Comprobación de resultados