import sys
import pfsim
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import optm
import logging

logging.basicConfig(filename='Main.log',
                    filemode='a',
level=logging.INFO,
format='%(asctime)s - %(levelname)s - %(message)s')

def setup_simulation(estudio):
    logging.info("Setting up simulation with estudio: %d", estudio)
    params = {
        'VOLL': 2500,
        'Pot_Down': 0,
        'Flujos': 1,
        'Perdidas': 0,
        'PerdidasPlus': 0,
        'Sim_PF': True,
        'Costo_ts': 0,
        'TS': False
    }

    if estudio == 1:  # Sin sistema de tranmisión
        params.update({'VOLL': 300, 
                        'Flujos': 0})
    elif estudio == 2: # Sistema de transmisión sin perdidas y con efecto de voltage
        params.update({'VOLL': 300})
    elif estudio == 3:  # Sistema de transmisión con perdidas
        params.update({'VOLL': 300, 
                       'Perdidas': 1, 
                       'PerdidasPlus': 1})
    elif estudio == 5:  # Se puede bajar la potencia de los generadores
        params.update({'VOLL': 300, 
                       'Pot_Down': 1, 
                       'Perdidas': 1, 
                       'PerdidasPlus': 1})
    elif estudio == 6:  # Transmission Switching con perdidas
        params.update({'VOLL': 300, 
                       'Pot_Down': 1, 
                       'Perdidas': 1, 
                       'Costo_ts': 10, 
                       'TS': 1})

    logging.info("Simulation parameters: %s", params)
    return params

def change_max_line_values(pf, lines, values):
    logging.info("Changing max line values for lines: %s with values: %s", lines, values)
    logging.info("Changing max line values")
    pf.ChangeMaxLine(lines, values)

def export_simulation_data(pf, estudio):
    logging.info("Exporting simulation data for estudio: %d", estudio)
    logging.info("Exporting CSV")
    #if estudio == 0:
    #    pf.export_csv('DC')
        #else:
    #    pf.export_csv('AC')
    pf.export_csv('AC')


def run_optimization_model(pf, sim, params):
    logging.info("Running optimization model with params: %s", params)
    logging.info("Running optimization model")
    m_optm = optm.Modelo()
    m_optm.pot_down = params['Pot_Down']
    m_optm.flujos = params['Flujos']
    m_optm.losses = params['Perdidas']
    m_optm.losses_plus = params['PerdidasPlus']
    m_optm.Voll = params['VOLL']
    m_optm.TS = params['TS']
    m_optm.costo_ts = params['Costo_ts']

    m_optm(pf, sim)
    m_optm.run()
    m_optm.Results(pf, sim)

    logging.info("Optimization model results: %s", m_optm)
    return m_optm

def extract_data_fullsim(pf):
    logging.info("Extracting full simulation data")
    pf_time, freq_values = pf.extract_data('Term_10_4.ElmTerm', 'm:fehz', return_time=True)

    # Extrer datos de cada generador
    cont=0
    return_time = False
    list_gens_values = np.zeros((pf.ngen+pf.ngenstat+pf.ngenpv, len(pf_time)))
    #list_pt_values = np.zeros((pf.ngen, len(pf_time)))
    #list_pgt_values = np.zeros((pf.ngen, len(pf_time)))
    for gen in pf.all_gen:
        var_values = pf.extract_data(gen+'.ElmSym', 'm:Psum:bus1', return_time)
        list_gens_values[cont] = var_values
        #list_pt_values[cont] = pf.extract_data(gen+'.ElmSym', 's:pt', return_time)
        #list_pgt_values[cont] = pf.extract_data(gen+'.ElmSym', 's:pgt', return_time)
        cont+=1

    for gen in pf.all_genstat:
        var_values = pf.extract_data(gen+'.ElmGenstat', 'm:Psum:bus1', return_time)
        list_gens_values[cont] = var_values
        cont+=1

    for gen in pf.all_pv:
        var_values = pf.extract_data(gen+'.ElmPvsys', 'm:Psum:bus1', return_time)
        list_gens_values[cont] = var_values
        cont+=1


    # Extraer datos de todos las demandas
    cargas_name = list()
    list_load_values = np.zeros((len(pf.cargas), len(pf_time)))
    cont=0
    for carga in pf.cargas:
        var_values = pf.extract_data(carga.loc_name+'.ElmLod', 'm:Psum:bus1', return_time)
        cargas_name.append(carga.loc_name)
        list_load_values[cont] = var_values
        cont+=1

    # Extraer datos de todas las líneas
    lineas_name = list()
    list_line_values = np.zeros((len(pf.lineas), len(pf_time)))
    cont=0
    for linea in pf.lineas:
        var_values = pf.extract_data(linea.loc_name+'.ElmLne', 'c:Losses', return_time)
        lineas_name.append(linea.loc_name)
        list_line_values[cont] = np.array(var_values)/1000
        cont+=1

    trafos_name = list()
    list_trafos_values= np.zeros((len(pf.trafos), len(pf_time)))
    cont=0
    for trafo in pf.trafos:
        trafos_name.append(trafo.loc_name)
        var_hv = pf.extract_data(trafo.loc_name+'.ElmTr2', 'm:P:bushv', return_time)
        var_lv = pf.extract_data(trafo.loc_name+'.ElmTr2', 'm:P:buslv', return_time)
        list_trafos_values[cont] = abs(np.array(var_hv)+np.array(var_lv))
        cont+=1

    sum_gens = np.zeros(len(pf_time))
    for tiempo in range(len(pf_time)):
        sum_gens[tiempo] = list_gens_values[:,tiempo].sum()

    sum_loads = np.zeros(len(pf_time))
    for tiempo in range(len(pf_time)):
        sum_loads[tiempo] = list_line_values[:,tiempo].sum()
        sum_loads[tiempo] += list_trafos_values[:,tiempo].sum()
        sum_loads[tiempo] += list_load_values[:,tiempo].sum()


    data = pd.DataFrame(np.vstack((pf_time, freq_values, list_gens_values, list_load_values, list_line_values, list_trafos_values, sum_loads)).T, 
                        columns=['Time', 'Freq'] + pf.all_gen + pf.all_genstat + pf.all_pv + cargas_name + lineas_name + trafos_name + ['Sum Total Loads'])


    data.to_excel('generator_load_data.xlsx', index=False)
logging.info("Data extraction complete, saved to generator_load_data.xlsx")


def Variabilidad(pf, sim, m_optm, t_initial, tstop_cpf, t_final):
    logging.info("Starting variability analysis")
    # Tiempo entre cada simulación
    t_step = 4
    n_step = len(range(tstop_cpf+t_step, t_final, t_step))

    part_factors = np.zeros((pf.n_gen_agc,n_step+1,pf.Nt))


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

            ### NUEVO SF con TRANSMISSION SWITCHING
            new_SF = pfsim.new_SF(pf, m_optm, ti, scen)

            ### Eventos de desconexiones de líneas al inicio del CSF
            pfsim.CreateEvents_line(pf, m_optm, scen, ti, tstop_cpf)
            #pfsim.Set_param_agc(pf, tstop_cpf, m_optm.part_factors[:,0,0])


            events_folder = pf.IntEvt.GetContents()

            evt = events_folder[sim.name_events.index('Salida Gen')]
            evt.outserv = 0
            evt.time = t_initial
            evt.p_target = pf.generadores[list(pf.dict_gen).index(gen_out)]


            Pout_origin = sim.P_out[scen,ti]-pf.dda_barra[:,ti].sum() + sim.D_pfc[:,scen,ti].sum() - sim.PL_pre_line[ti]
            Pout_real = sim.P_out[0,0]-pf.dda_barra[:,0].sum()

            P_ens = np.zeros(n_step)
            
            Pb_genAGC = pf.Pgen_pre[m_optm.pos_gen_csf,ti].sum()
            Pb_genNOagc = -sim.P_out[0,0] + pf.Pgen_pre[:,ti].sum() - Pb_genAGC + pf.Pgenstat_pre[:,ti].sum() + pf.Ppv_pre[:,ti].sum()
            Pb_load = sim.D_pfc[:,scen,ti].sum()
            Pb_losses = sim.PL_pre_line[ti]

            Precp_original = np.zeros(n_step)
            Precp = np.zeros(n_step)
            P_change_list = np.zeros(n_step)
            
            Pn_genAGC_list= np.zeros(n_step)
            Pn_genNOagc_list= np.zeros(n_step)
            Pn_load_list= np.zeros(n_step)
            Pn_losses_list= np.zeros(n_step)
            freq_list = np.zeros(n_step)
            j = -1
            previous_part_factors = np.zeros(pf.n_gen_agc)
            sum_p_model_list = np.zeros(n_step)
            for t_int in range(tstop_cpf+t_step, t_final, t_step):
                j+=1

                print('Simulación - t_int: ' + str(t_int))
                ShortSim = pfsim.ShortSim(pf, m_optm, sim, new_SF, gen_out, t_int)

                #t_sim, freq_values = pf.extract_data('Term_10_4.ElmTerm', 'm:fehz', return_time=True, start_time=t_int)

                freq = pf.app.GetCalcRelevantObjects('Term_10_4.ElmTerm')[0].GetAttribute('m:fehz')

                freq_list[j] = freq

                Pn_genAGC = ShortSim.p_genAGC
                Pn_genNOagc = ShortSim.p_gen.sum() - Pn_genAGC + ShortSim.p_genstat.sum() + ShortSim.pv_gen.sum()
                Pn_load = ShortSim.D_t.sum()
                Pn_losses = ShortSim.PL.sum()

                Pn_genAGC_list[j] = Pn_genAGC
                Pn_genNOagc_list[j] = Pn_genNOagc
                Pn_load_list[j] = Pn_load
                Pn_losses_list[j] = Pn_losses

                Precp_original[j] = Pout_origin + Pn_losses
                Precp[j] = Pout_origin + Pn_losses - Pn_genAGC 


                P_change = -Pn_genNOagc + Pb_genNOagc
                P_change_list[j] = P_change

                P_change2 = Pout_real - P_change


                PartialModel = optm.PartialModel(pf, m_optm, sim, ShortSim, new_SF, P_change ,scen, ti, t_int)
                
                print('P_partialModel = ' + str(-PartialModel.pg_dec.x.sum()+PartialModel.pg_inc.x.sum()))


                part_factors[:,int((t_int-tstop_cpf)/t_step),ti] = PartialModel.part_factors
                #print('Part Factors: ' + str([factor for factor in PartialModel.part_factors if factor != 0]))

                sum_p_model = PartialModel.pg_inc.x.sum() - PartialModel.pg_dec.x.sum() + PartialModel.p_ens.x.sum()    
                #print('P_ens: ' + str(PartialModel.p_ens.x.sum()))
                if PartialModel.p_ens.x.sum() > 0:
                    print('Error: P_ens > 0')
                    P_ens[j] = PartialModel.p_ens.x.sum()

                pfsim.Set_param_agc(pf, t_int, PartialModel.part_factors, previous_part_factors)

                previous_part_factors = PartialModel.part_factors.copy()


            # Control Terciario de frecuencia

            # Realizar optimización

            # Simulación

    data = pd.DataFrame(np.vstack((range(tstop_cpf+t_step, t_final, t_step), freq_list, Pn_genAGC_list, Pn_genNOagc_list, Pn_load_list, Pn_losses_list, sum_p_model_list, P_change_list)).T, 
                        columns=['Time', 'Freq', 'SumAGC', 'SumNOagc', 'SumLoad', 'SumLoss', 'SumP', 'P_Change'])

    data.to_excel('data_short_sim.xlsx', index=False)
    logging.info("Variability analysis complete, saved to data_short_sim.xlsx")




    # Plotting the results in 4 figures
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))

    axs[0, 0].plot(Pn_genAGC_list, label='Pb_genAGC')
    axs[0, 0].axhline(y=Pb_genAGC, color='r', linestyle='--', label='Pb_genAGC Line')
    axs[0, 0].set_title('Pb_genAGC')
    axs[0, 0].set_xlabel('Time Step')
    axs[0, 0].set_ylabel('Power (MW)')
    axs[0, 0].grid(True)

    axs[0, 1].plot(Pn_genNOagc_list, label='Pb_genNOagc')
    axs[0, 1].axhline(y=Pb_genNOagc, color='r', linestyle='--', label='Pb_genAGC Line')
    axs[0, 1].set_title('Pb_genNOagc')
    axs[0, 1].set_xlabel('Time Step')
    axs[0, 1].set_ylabel('Power (MW)')
    axs[0, 1].grid(True)

    axs[1, 0].plot(Pn_load_list, label='Pb_load')
    axs[1, 0].axhline(y=Pb_load, color='r', linestyle='--', label='Pb_genAGC Line')
    axs[1, 0].set_title('Pb_load')
    axs[1, 0].set_xlabel('Time Step')
    axs[1, 0].set_ylabel('Power (MW)')
    axs[1, 0].grid(True)

    axs[1, 1].plot(Pn_losses_list, label='Pb_losses')
    axs[1, 1].axhline(y=Pb_losses, color='r', linestyle='--', label='Pb_genAGC Line')
    axs[1, 1].set_title('Pb_losses')
    axs[1, 1].set_xlabel('Time Step')
    axs[1, 1].set_ylabel('Power (MW)')
    axs[1, 1].grid(True)

    # Save the figure
    plt.tight_layout()
    plt.savefig('power_generation_and_load.png')
    plt.close()



def main():
    logging.info("Starting main function")
    estudio = 6
    nt = 1  # N° de escenarios
    params = setup_simulation(estudio)

    logging.info("Starting simulation")
    t0 = time.time()
    project_file = 'Articulo2 - TS'
    pf = pfsim.PowerFactorySim(project_file)

    lines = ['Line3_1']
    values = [500]
    change_max_line_values(pf, lines, values)

    pf.TS = params['TS']
    pf.Nt = nt
    export_simulation_data(pf, estudio)

    dict_lineas = pf.get_data()
    t_initial = 0
    t_final = 200 + t_initial
    tstop_cpf = 30 + t_initial

    sim = pfsim.Simulacion(pf, t_initial, tstop_cpf)
    t1 = time.time()

    m_optm = run_optimization_model(pf, sim, params)
    t2 = time.time()
    logging.info("Optimization model formulated in %.4f seconds", t2 - t1)

    t3 = time.time()
    logging.info("Optimization model solved in %.4f seconds", t3 - t2)

    t4 = time.time()
    logging.info("Results extracted in %.4f seconds", t4 - t3)

    print('=> Extract data time: %.4f (s)' % (t1 - t0))
    print('=> Formulation time: %.4f (s)' % (t2 - t1))
    print('=> Solution time: %.4f (s)' % (t3 - t2))

    pf.app.ResetCalculation()

    if False:
        pfsim.Set_param_agc(pf, tstop_cpf, m_optm.part_factors[:,0,0])

        # Do simulación
        Monitored_vars = {'*.ElmSym': ['m:Psum:bus1'],
                          '*.ElmGenstat': ['m:Psum:bus1'],
                          '*.ElmPvsys': ['m:Psum:bus1'],
                          '*.ElmLod': ['m:Psum:bus1'],
                          '*.ElmLne': ['c:Losses'],
                          '*.ElmTr2':['m:P:bushv','m:P:buslv']}

        pf.prepare_dynamic_sim(Monitored_vars, 'rms')
        pf.run_dynamic_sim(end_sim=t_final)

        # Extract data
        # Extraer frecuencia
    pf_time, freq_values = pf.extract_data('Term_10_4.ElmTerm', 'm:fehz', return_time=True)

    # Extrer datos de cada generador
    cont=0
    return_time = False
    list_gens_values = np.zeros((pf.ngen+pf.ngenstat+pf.ngenpv, len(pf_time)))
    #list_pt_values = np.zeros((pf.ngen, len(pf_time)))
    #list_pgt_values = np.zeros((pf.ngen, len(pf_time)))
    for gen in pf.all_gen:
        var_values = pf.extract_data(gen+'.ElmSym', 'm:Psum:bus1', return_time)
        list_gens_values[cont] = var_values
        #list_pt_values[cont] = pf.extract_data(gen+'.ElmSym', 's:pt', return_time)
        #list_pgt_values[cont] = pf.extract_data(gen+'.ElmSym', 's:pgt', return_time)
        cont+=1

    for gen in pf.all_genstat:
        var_values = pf.extract_data(gen+'.ElmGenstat', 'm:Psum:bus1', return_time)
        list_gens_values[cont] = var_values
        cont+=1

    for gen in pf.all_pv:
        var_values = pf.extract_data(gen+'.ElmPvsys', 'm:Psum:bus1', return_time)
        list_gens_values[cont] = var_values
        cont+=1


    # Extraer datos de todos las demandas
    cargas_name = list()
    list_load_values = np.zeros((len(pf.cargas), len(pf_time)))
    cont=0
    for carga in pf.cargas:
        var_values = pf.extract_data(carga.loc_name+'.ElmLod', 'm:Psum:bus1', return_time)
        cargas_name.append(carga.loc_name)
        list_load_values[cont] = var_values
        cont+=1

    # Extraer datos de todas las líneas
    lineas_name = list()
    list_line_values = np.zeros((len(pf.lineas), len(pf_time)))
    cont=0
    for linea in pf.lineas:
        var_values = pf.extract_data(linea.loc_name+'.ElmLne', 'c:Losses', return_time)
        lineas_name.append(linea.loc_name)
        list_line_values[cont] = np.array(var_values)/1000
        cont+=1
    
    trafos_name = list()
    list_trafos_values= np.zeros((len(pf.trafos), len(pf_time)))
    cont=0
    for trafo in pf.trafos:
        trafos_name.append(trafo.loc_name)
        trafo_name = trafo.loc_name
        var_hv = pf.extract_data(trafo_name +'.ElmTr2', 'm:P:bushv', return_time)
        var_lv = pf.extract_data(trafo_name +'.ElmTr2', 'm:P:buslv', return_time)
        list_trafos_values[cont] = abs(np.array(var_hv)+np.array(var_lv))
        cont+=1

    sum_gens = np.zeros(len(pf_time))
    for tiempo in range(len(pf_time)):
        sum_gens[tiempo] = list_gens_values[:,tiempo].sum()

    sum_loads = np.zeros(len(pf_time))
    for tiempo in range(len(pf_time)):
        sum_loads[tiempo] = list_line_values[:,tiempo].sum()
        sum_loads[tiempo] += list_trafos_values[:,tiempo].sum()
        sum_loads[tiempo] += list_load_values[:,tiempo].sum()


    data = pd.DataFrame(np.vstack((pf_time, freq_values, list_gens_values, list_load_values, list_line_values, list_trafos_values, sum_loads)).T, 
                        columns=['Time', 'Freq'] + pf.all_gen + pf.all_genstat + pf.all_pv + cargas_name + lineas_name + trafos_name + ['Sum Total Loads'])
    

    data.to_excel('generator_load_data.xlsx', index=False)



# Extrer datos de cada generador
#return_time = True
#for gen in pf.all_gen:
#    if return_time:
#        pf_time, var_values = pf.extract_data(gen.loc_name, 'm:Psum:bus1', return_time)
#    else:
#        var_values = pf.extract_data(gen.loc_name, 'm:Psum:bus1', return_time)
            
#return_time = True
#for gen in pf.all_genstat:
#    if return_time:
#        pf_time, var_values = pf.extract_data(gen.loc_name, 'm:Psum:bus1', return_time)
#    else:
#        var_values = pf.extract_data(gen.loc_name, 'm:Psum:bus1', return_time)#

#return_time = True
#for gen in pf.all_pv:
#    if return_time:
#        pf_time, var_values = pf.extract_data(gen.loc_name, 'm:Psum:bus1', return_time)
#    else:
#        var_values = pf.extract_data(gen.loc_name, 'm:Psum:bus1', return_time)
#

    if True:
        t4 = time.time()
        ### Chekeos
        #1
        all_gen = np.array(pf.all_gen)
        PgenPre_all = list()
        NamePPre_all =  ['P valle', 'P media', 'P peak']
        NamePPre_all = NamePPre_all[0:nt]
        for ti in range(nt):
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
        all_P_out = np.zeros(nt*pf.Ns)
        for ti in range(nt):
            all_P_out[ti*pf.Ns:ti*pf.Ns+pf.Ns] = sim.P_out[:,ti]
        dddf = np.vstack((range(1,pf.Ns*nt+1),np.tile(pf.Gen_Outages,nt),all_P_out,np.tile(sim.Barra_gen_out,nt)))
        Res_escenarios = pd.DataFrame(dddf.T,index=range(1,pf.Ns*nt+1), columns=['Escenario', 'Gen', 'Potencia ', 'Posición'])
        
        
        #df_genOnAGC = pd.DataFrame(np.vstack((Gen_outages, Gen_on_AGC[ti,:,:])).T, columns=np.insert(all_gen,0,'Gen', axis=0))
        
        # %%
        with pd.ExcelWriter('00Resultados'+ project_file+ '.xlsx') as writer:
            Res_escenarios.to_excel(writer, sheet_name='Escenarios - Gen OUT', index=False)
            datos_gen.to_excel(writer,sheet_name='DatosGenAGC', index = False)
        
            for ti in range(nt):
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
        
        


        
        # Simulación PowerFactory
        results_line = np.zeros((len(dict_lineas),pf.Ns,nt))
        result_gen = np.zeros((pf.ngen,pf.Ns,nt))
        if params['Sim_PF']:
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


            for ti in range(nt):
                pf.scens[ti].Activate()
                print('Escenario: ' + str(ti+1))
                for scen in range(pf.Ns):
                    for e in pf.events_folder:
                        e.outserv = 1   
                    print('Progreso: ' + str(scen+1+pf.Ns*ti) + '/' + str(pf.Ns*nt))
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
                    for obj in pf.indices_obj:
                        for line in pf.lineas:
                            if obj in pf.TS_lines and line.loc_name == obj:
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
            for ti in range(nt):    
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
            for ti in range(nt):
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
                    elif m_optm.flujos:
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
            


    logging.info("Simulation complete")
    print('finish!')
#  %% Comprobación de resultados

if __name__ == "__main__":
    main()
    
    