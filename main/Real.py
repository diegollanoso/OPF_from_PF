import sys
import pfsim
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import optm
import logging

# logging.basicConfig(filename='report.log',
#                     filemode='a',
#                     level=logging.INFO,
#                     format='%(asctime)s - %(levelname)s - %(message)s')


class PrintToLogger:
    """
    Custom class to redirect stdout to logger
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        #self.stdout = sys.stdout
        self.terminal = sys.__stdout__

    def write(self, buf):
        self.terminal.write(buf)
        # For each line in the buffer, strip whitespace and log if not empty
        for line in buf.rstrip().splitlines():
            if line:
                self.logger.log(self.log_level, line.rstrip())
    
    def flush(self):
        self.terminal.flush()


def setup_logging(log_file = 'report.log'):
    """
    Configure logging to output to both file and console with different formats
    """

    # Create formatters for different outputs
    file_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_formatter = logging.Formatter(
        '%(levelname)s - %(message)s'
    )

    # Configure file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    # Get root logger and configure it
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    sys.stdout = PrintToLogger(root_logger, logging.INFO)


def setup_simulation(estudio):
    logging.info("Setting up simulation with estudio: %d", estudio)
    params = {
        'VOLL': 300,
        'Pot_Down': 0,
        'Flujos': 1,
        'Perdidas': 0,
        'PerdidasPlus': 0,
        'Sim_PF': True,
        'Costo_ts': 0,
        'TS': False
    }

    if estudio == 1:  # Sin sistema de tranmisión
        params.update({'Flujos': 0})
    elif estudio == 2: # Sistema de transmisión sin perdidas y con efecto de voltage
        pass
    elif estudio == 3:  # Sistema de transmisión con perdidas
        params.update({'Perdidas': 1, 
                       'PerdidasPlus': 1})
    elif estudio == 4:  # No se utilizan las condiciones de complementaridad / exclusividad
        params.update({'Perdidas': 1})
    elif estudio == 5:  # Se puede bajar la potencia de los generadores
        params.update({'Pot_Down': 1, 
                       'Perdidas': 1, 
                       'PerdidasPlus': 1})
    elif estudio == 6:  # Transmission Switching con perdidas
        params.update({'VOLL': 300, 
                       'Pot_Down': 1, 
                       'Perdidas': 1,
                       'PerdidasPlus': 0,
                       'Costo_ts': 10, 
                       'TS': 1})

    logging.info("Simulation parameters: %s", params)
    return params






def run_optimization_model(pf, sim, params):
    logging.info("Running optimization model with params: %s", params)
    t1 = time.time()
    m_optm = optm.Modelo()
    m_optm.pot_down = params['Pot_Down']
    m_optm.flujos = params['Flujos']
    m_optm.losses = params['Perdidas']
    m_optm.losses_plus = params['PerdidasPlus']
    m_optm.Voll = params['VOLL']
    m_optm.TS = params['TS']
    m_optm.costo_ts = params['Costo_ts']

    m_optm(pf, sim)
    t2 = time.time()
    m_optm.run()
    t3 = time.time()
    m_optm.Results(pf, sim)


    print('=> Formulation time: %.4f (s)' % (t2-t1))
    print('=> Solution time: %.4f (s)' % (t3-t2))

    logging.info("Optimization model results: %s", m_optm)
    return m_optm

# Do a UC and set values for all generatos
def unit_commitment(pf):
    model_uc = optm.Model_UC(pf)

    for gen,value in zip(pf.generadores, 100*model_uc.pg.x):
        gen.pgini = value
    for gen,value in zip(pf.genstate, 100*model_uc.pg.x[pf.ngen:pf.ngen+pf.ngenstat]/pf.ngenstat_par):
        gen.pgini = value
    for gen,value in zip(pf.pvsys, 100*model_uc.pg.x[pf.ngen+pf.ngenstat:pf.ngen+pf.ngenstat+pf.ngenpv]/pf.ngenpv_par):
        gen.pgini = value

def Alternative_solutions(pf, m_optm):
    if m_optm.m.SolCount > 1:
        logging.info("Extracting %d alternatives solutions", m_optm.m.SolCount)
        for e in range(1,m_optm.m.SolCount):
            m_optm.m.setParam(m_optm.gp.GRB.Param.SolutionNumber,e)
            print('Another Solution: ' + str(m_optm.m.PoolObjVal))

            if m_optm.pot_down:
                Csfc = pf.Ccte_gen[pf.pos_gen_agc_list] @ (m_optm.vg_inc.xn + m_optm.vg_dec.xn)
            else:
                Csfc = pf.Ccte_gen[pf.pos_gen_agc_list] @ (m_optm.vg_inc.xn)

            Cts = 0
            Cop = 0
            Cue = 0
            Cpl = 0

            for ti in range(pf.Nt):
                Cop_s = 0
                Cue_s = 0
                Cpl_s = 0
                Cts_s = 0
                for s in range(pf.Ns):
                    if m_optm.pot_down:
                        inc = m_optm.pg_inc[:,s,ti].xn + m_optm.pg_dec[:,s,ti].xn
                    else:
                        inc = m_optm.pg_inc[:,s,ti].xn

                    Cop_s += pf.Cvar_gen[pf.pos_gen_agc_list]*pf.Sb @ inc
                    Cue_s += (m_optm.Voll*m_optm.p_ens[:,s,ti].xn).sum()*pf.Sb
                    if m_optm.losses and not m_optm.losses_plus:
                        Cpl_s += (m_optm.Voll*m_optm.ploss[:,s,ti].xn).sum()*pf.Sb
                    if m_optm.TS:
                        Cts_s += m_optm.costo_ts * (1-m_optm.s_ts.x[:,s,ti]).sum()

                Cop += Cop_s
                Cue += Cue_s
                Cpl += Cpl_s
                Cts += Cts_s

            print_cts = ''
            if m_optm.TS:
                print_cts = ' + Cts = %.2f' % (Cts/pf.Ns)
            if m_optm.losses and not m_optm.losses_plus:
                print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($) + Cpl = %.2f ($)' % (m_optm.m.PoolObjVal,Csfc,Cop/pf.Ns,Cue/pf.Ns,Cpl/pf.Ns) + print_cts)
            else:
                print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($)' % (m_optm.m.PoolObjVal,Csfc,Cop/pf.Ns,Cue/pf.Ns) + print_cts)


def Variabilidad(pf, t_initial, t_final, tstop_cpf, sim, m_optm):
    # Tiempo entre cada simulación
    t_step = 4
    n_step = len(range(tstop_cpf+t_step, t_final, t_step))

    part_factors = np.zeros((pf.n_gen_agc,n_step+1,pf.Nt))

    logging.info("Running simulation every %d seconds", t_step)
    # Variabilidad
    for ti in range(pf.Nt):
        pf.scens[ti].Activate()
        logging.info("Running simulation for scenario %d", ti)
        scen = -1
        for gen_out in pf.Gen_Outages:
            logging.info("Running simulation for generator %s", gen_out)
            scen += 1

            #Desactivar todos los eventos
            events_folder = pf.IntEvt.GetContents()
            for e in events_folder:
                    e.outserv = 1

            ### NUEVO SF con TRANSMISSION SWITCHING
            new_SF = pfsim.new_SF(pf, m_optm, ti, scen)

            ### Eventos de desconexiones de líneas al inicio del CSF
            pfsim.CreateEvents_line(pf, m_optm, scen, ti, tstop_cpf)
            
            ### Set factores de participación en AGC
            #pfsim.Set_param_agc(pf, tstop_cpf, m_optm.part_factors[:,scen,ti])

            ## Generadores que participan en el CSF
            pos_gen_csf = m_optm.pos_gen_csf
            gen_csf = np.array(pf.all_gen)[pos_gen_csf]
            if gen_out in gen_csf:
                pos = np.where(gen_csf == gen_out)
                gen_csf = np.delete(gen_csf, pos)
                pos_gen_csf = np.delete(pos_gen_csf, pos)
            logging.info("Generadores que participan en el CSF: %s", gen_csf)

            events_folder = pf.IntEvt.GetContents()
            evt = events_folder[sim.name_events.index('Salida Gen')]
            evt.outserv = 0
            evt.time = t_initial
            evt.p_target = pf.generadores[list(pf.dict_gen).index(gen_out)]

            # Potencia a recuperar
            # le faltan las perdidas actuales 
            Pout_origin = sim.P_out[scen,ti]-pf.dda_barra[:,ti].sum() + sim.D_pfc[:,scen,ti].sum() - sim.PL_pre_line[ti]

            ### POTENCIAS PREVIAS
            # Potencia previa de generadores participantes del CSF
            Pb_genAGC = pf.Pgen_pre[pos_gen_csf,ti].sum()
            # Potencia previa del resto de generadores
            Pb_genNOagc = -sim.P_out[scen,ti] + pf.Pgen_pre[:,ti].sum() - Pb_genAGC + pf.Pgenstat_pre[:,ti].sum() + pf.Ppv_pre[:,ti].sum()
            Pb_genVar = pf.Pgenstat_pre[:,ti].sum() + pf.Ppv_pre[:,ti].sum()
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
            previous_part_factors = np.zeros(pf.n_gen_agc)
            sum_pinc_model = np.zeros(n_step)
            sum_pdec_model = np.zeros(n_step)
            sum_ploss_model = np.zeros(n_step)
            P_faltante = np.zeros(n_step)
            P_ens = np.zeros(n_step)
            j = -1
            for t_int in range(tstop_cpf+t_step, t_final, t_step):
                j+=1

                print('Simulación - t_int: ' + str(t_int))
                ShortSim = pfsim.ShortSim(pf, gen_csf, sim, new_SF, gen_out, t_int)

                #t_sim, freq_values = pf.extract_data('Term_10_4.ElmTerm', 'm:fehz', return_time=True, start_time=t_int)

                freq_list[j] = pf.app.GetCalcRelevantObjects('Term_10_4.ElmTerm')[0].GetAttribute('m:fehz')

                #if freq_list[j] > 49.99:
                #    break


                ### POTENCIAS NUEVAS
                Pn_genAGC = ShortSim.p_genAGC
                Pn_genVar = ShortSim.p_genstat.sum() + ShortSim.pv_gen.sum()
                Pn_genNOagc = ShortSim.p_gen.sum() - Pn_genAGC + ShortSim.p_genstat.sum() + ShortSim.pv_gen.sum()
                Pn_load = ShortSim.D_t.sum()
                Pn_losses = ShortSim.PL.sum()

                Pn_genAGC_list[j] = Pn_genAGC
                Pn_genNOagc_list[j] = Pn_genNOagc
                Pn_load_list[j] = Pn_load
                Pn_losses_list[j] = Pn_losses

                # Potencia recuperada SOLO por el AGC
                Precp[j] = Pb_genAGC - Pn_genAGC
 

                # Variación de potencia de ERV y demanda variable 
                #P_change = -Pn_genNOagc + Pb_genNOagc
                P_change_ERV = -Pn_genVar + Pb_genVar
                P_change_load = -Pn_load + Pb_load
                P_change_list[j] = P_change_ERV


                #P_change2 = Pout_real - P_change


                PartialModel = optm.PartialModel(pf, m_optm, sim, ShortSim, new_SF, P_change_list[j] + Precp[j] ,scen, ti, t_int)

                # Potencia Faltante

                P_faltante[j] = P_change_list[j] + Precp[j] + sim.P_out[scen,ti] - pf.dda_barra[:,ti].sum() + ShortSim.D_t.sum() - sim.PL_pre_line[ti] + PartialModel.ploss.x.sum()

                logging.debug("Potencia AGC por recuperar: %s" % (PartialModel.pg_inc.x.sum() - PartialModel.pg_dec.x.sum() + PartialModel.p_ens.x.sum())), 

                part_factors[:,int((t_int-tstop_cpf)/t_step),ti] = PartialModel.part_factors
                #print('Part Factors: ' + str([factor for factor in PartialModel.part_factors if factor != 0]))

                sum_pinc_model[j] = PartialModel.pg_inc.x.sum()
                sum_pdec_model[j] = PartialModel.pg_dec.x.sum() 
                sum_ploss_model[j] = PartialModel.ploss.x.sum()
                #print('P_ens: ' + str(PartialModel.p_ens.x.sum()))
                if PartialModel.p_ens.x.sum() > 0:
                    logging.debug('Error: P_ens > 0')
                    P_ens[j] = PartialModel.p_ens.x.sum()

                pfsim.Set_param_agc(pf, t_int, PartialModel.part_factors, previous_part_factors)

                previous_part_factors = PartialModel.part_factors.copy()


            # Control Terciario de frecuencia

            # Realizar optimización

            # Simulación

            logging.info("Simulation finished for contigency: %s", gen_out)
            data = pd.DataFrame(np.vstack((range(tstop_cpf+t_step, t_final, t_step), freq_list, Pn_genAGC_list, Pn_genNOagc_list, Pn_load_list, Pn_losses_list,P_faltante, P_change_list + Precp, sum_pinc_model, sum_pdec_model, sum_ploss_model)).T, 
                                columns=['Time', 'Freq', 'SumAGC', 'SumNOagc', 'SumLoad', 'SumLoss', 'P_faltante', 'P_Change', 'SumPinc', 'SumPdec', 'SumLoss'])

            data.to_excel(f'data_short_sim_{gen_out}_{time.strftime("%Y%m%d-%H%M%S")}.xlsx', index=False)
            logging.info("Data saved to data_short_sim.xlsx")


            # Plotting the results in 4 figures
            fig, axs = plt.subplots(2, 2, figsize=(15, 10))

            axs[0, 0].plot(Pn_genAGC_list, label='Pb_genAGC')
            axs[0, 0].axhline(y=Pb_genAGC, color='r', linestyle='--', label='Pb_genAGC Line')
            axs[0, 0].axhline(y=Pn_genAGC, color='b', linestyle='--', label='Pn_genAGC Line')
            axs[0, 0].set_title('Pb_genAGC')
            axs[0, 0].set_xlabel('Time Step')
            axs[0, 0].set_ylabel('Power (MW)')
            axs[0, 0].grid(True)

            axs[0, 1].plot(Pn_genNOagc_list, label='Pb_genNOagc')
            axs[0, 1].axhline(y=Pb_genNOagc, color='r', linestyle='--', label='Pb_genNOagc Line')
            axs[0, 1].axhline(y=Pn_genNOagc, color='b', linestyle='--', label='Pn_genNOagc Line')
            axs[0, 1].set_title('Pb_genNOagc')
            axs[0, 1].set_xlabel('Time Step')
            axs[0, 1].set_ylabel('Power (MW)')
            axs[0, 1].grid(True)

            axs[1, 0].plot(Pn_load_list, label='Pb_load')
            axs[1, 0].axhline(y=Pb_load, color='r', linestyle='--', label='Pb_load Line')
            axs[1, 0].axhline(y=Pn_load, color='b', linestyle='--', label='Pn_load Line')
            axs[1, 0].set_title('Pb_load')
            axs[1, 0].set_xlabel('Time Step')
            axs[1, 0].set_ylabel('Power (MW)')
            axs[1, 0].grid(True)

            axs[1, 1].plot(Pn_losses_list, label='Pb_losses')
            axs[1, 1].axhline(y=Pb_losses, color='r', linestyle='--', label='Pb_losses Line')
            axs[1, 1].axhline(y=Pn_losses, color='b', linestyle='--', label='Pn_losses Line')
            axs[1, 1].set_title('Pb_losses')
            axs[1, 1].set_xlabel('Time Step')
            axs[1, 1].set_ylabel('Power (MW)')
            axs[1, 1].grid(True)

            # Save the figure
            plt.tight_layout()
            plt.savefig(f'power_generation_and_load_{gen_out}_{time.strftime("%Y%m%d-%H%M%S")}.png')
            plt.close()
            logging.info("Power generation and load plot saved to power_generation_and_load.png")



def Revision_big_optm(params, project_file, pf, dict_lineas, t_initial, tstop_cpf, t_final, sim, m_optm):
    logging.info("Revision big optm with PowerFactory")
    t4 = time.time()
    ### Chekeos
    #1
    all_gen = np.array(pf.all_gen)
    PgenPre_all = list()
    NamePPre_all =  ['P valle', 'P media', 'P peak']
    NamePPre_all = NamePPre_all[0:pf.Nt]
    for ti in range(pf.Nt):
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
    all_P_out = np.zeros(pf.Nt*pf.Ns)
    for ti in range(pf.Nt):
        all_P_out[ti*pf.Ns:ti*pf.Ns+pf.Ns] = sim.P_out[:,ti]
    dddf = np.vstack((range(1,pf.Ns*pf.Nt+1),np.tile(pf.Gen_Outages,pf.Nt),all_P_out,np.tile(sim.Barra_gen_out,pf.Nt)))
    Res_escenarios = pd.DataFrame(dddf.T,index=range(1,pf.Ns*pf.Nt+1), columns=['Escenario', 'Gen', 'Potencia ', 'Posición'])


    #df_genOnAGC = pd.DataFrame(np.vstack((Gen_outages, Gen_on_AGC[ti,:,:])).T, columns=np.insert(all_gen,0,'Gen', axis=0))

    # %%
    with pd.ExcelWriter(f'00Resultados{project_file}_{time.strftime("%Y%m%d-%H%M%S")}.xlsx') as writer:
        Res_escenarios.to_excel(writer, sheet_name='Escenarios - Gen OUT', index=False)
        datos_gen.to_excel(writer,sheet_name='DatosGenAGC', index = False)

        for ti in range(pf.Nt):
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
    logging.info("Results saved to Excel file %s", f'00Resultados{project_file}_{time.strftime("%Y%m%d-%H%M%S")}.xlsx')




    # Simulación PowerFactory
    results_line = np.zeros((len(dict_lineas),pf.Ns,pf.Nt))
    result_gen = np.zeros((pf.ngen,pf.Ns,pf.Nt))
    if params['Sim_PF']:            
        #signal_list = list()
        #for gen in pf.gen_agc:
        #    comp = gen.c_pmod
        #    for slot in comp.pblk:
        #        ## Slot pcu
        #        if slot.loc_name[:9] == 'pcu Slot' or slot.loc_name[:9] == 'Gov Slot':
        #            val = comp.pelm[comp.pblk.index(slot)].signal[3][-2:]
        #            if not val[0].isdigit():
        #                val = val[1]
        #            signal_list.append('Evento Gamma ' + val)


        for ti in range(pf.Nt):
            pf.scens[ti].Activate()
            print('Escenario: ' + str(ti+1))
            for scen in range(pf.Ns):
                events_folder = pf.app.GetFromStudyCase('IntEvt').GetContents()
                for e in events_folder:
                    e.outserv = 1   
                print('Progreso: ' + str(scen+1+pf.Ns*ti) + '/' + str(pf.Ns*pf.Nt))
                print('Generador Out: ' + pf.Gen_Outages[scen])
                evt = pf.events_folder[sim.name_events.index('Salida Gen')]
                evt.outserv = 0
                evt.time = t_initial
                evt.p_target = pf.generadores[list(pf.dict_gen).index(pf.Gen_Outages[scen])]
                
                pfsim.Set_param_agc(pf, tstop_cpf, m_optm.part_factors[:,scen,ti])
                #count=0
                #for gen in pf.Gen_AGC:
                #    evt2 = pf.events_folder[sim.name_events.index(signal_list[count])]
                #    evt2.outserv = 0
                #    evt2.value = str(m_optm.part_factors[:,scen,ti][pf.name_gen_agc_list.index(gen)])
                #    evt2.time = t_initial
                #    count += 1

                line_switch = pfsim.CreateEvents_line(pf, m_optm, scen, ti, tstop_cpf)
                #cont = 0
                #line_sim_out = list()
                #for obj in pf.indices_obj:
                #    for line in pf.lineas:
                #        if obj in pf.TS_lines and line.loc_name == obj:
                #            if m_optm.s_ts.x[cont,scen,ti] == 0:
                #                line.outserv = 1
                #                line_sim_out.append(line)
                #            cont += 1

                pf.prepare_dynamic_sim({}, 'rms')
                pf.run_dynamic_sim(end_sim=t_final)

                for gen in pf.generadores:
                    if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                        result_gen[list(all_gen).index(gen.loc_name),scen,ti] = gen.GetAttribute('m:Psum:bus1')

                for line in pf.lineas:
                    if line.loc_name not in line_switch:
                        results_line[pf.all_line.index(line.loc_name),scen,ti] = (line.GetAttribute('m:Psum:bus1') - line.GetAttribute('m:Psum:bus2'))/2
                #results_line[:,scen] = np.array(list(map(lambda x: (x.GetAttribute('m:Psum:bus1') + x.GetAttribute('m:Psum:bus2'))/2 , pf.lineas)))  


                #for line in line_sim_out:
                #    line.outserv = 0

    with pd.ExcelWriter(f'01Comparaciones{project_file}_{time.strftime("%Y%m%d-%H%M%S")}.xlsx') as writer:
        for ti in range(pf.Nt):    
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

    # Define column names for different scenarios
    escenarios = ['peak', 'media', 'valle']
    escenariosPy_col = [f'Py_{esc}' for esc in escenarios]
    escenariosPF_col = [f'PF_{esc}' for esc in escenarios]

    with pd.ExcelWriter(f'03Casos{project_file}_{time.strftime("%Y%m%d-%H%M%S")}.xlsx') as writer:
        for ti in range(pf.Nt):
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
            diff_gen_pf = np.round(diff_gen_pf,1)
            results_pf = pd.DataFrame(np.vstack((pf.Gen_Outages,diff_gen_pf)).T, columns = [escenariosPF_col[ti]] + pf.name_gen_agc_list)
            results_pf.to_excel(writer, index = False, sheet_name='Case', startrow= 12*(3+ti))

            diff_line = pd.DataFrame(np.vstack((m_optm.f.x[:,:,ti])).T, columns = list(pf.indices_obj))
            diff_line.to_excel(writer, index = False, sheet_name='Case', startrow= ti*12, startcol=len(columns_diffpower)+2)


    t5 = time.time()

    print('=> Check time: %.4f (s)' % (t5-t4))
    logging.info("Check completed")


#def extract_data_fullsim(pf):
#    logging.info("Extracting data from PowerFactory")


def main():
    setup_logging()
    
    logging.info("Starting main function")
    estudio = 6
    nt = 1 # N° de escenarios
    params = setup_simulation(estudio)

    logging.info("Start up PowerFactory")
    t0 = time.time()

    #project_file = '3Bus_TS'
    #project_file = 'Articulo2'
    project_file = 'Articulo2 - TS'
    pf = pfsim.PowerFactorySim(project_file)

    # Valores originales
    #lines = ['Line3_1', 'Line3_3', 'Line3_4']
    #val = [600, 725, 590.2]

    # Cambio para que ocurra el TS
    lines = ['Line3_1']
    val = [500]

    pf.ChangeMaxLine(lines, val)

    pf.TS = params['TS']
    pf.Nt = nt # N° de escenarios
    pf.export_csv('AC')
    
    #SF, indices_obj = pfsim.ShiftFactors()

    pf.use_sf = 1 # 1 = SF calculados con matrices;  0 = SF de PF

    dict_lineas = pf.get_data()    
    

    # Dynamic Simulation until 30 seg - Finishing CPF
    t_initial = 0 # Tiempo Perturbación
    t_final = 200 + t_initial # Tiempo total

    tstop_cpf = 30 + t_initial# Tiempo hasta el CPF


    sim = pfsim.Simulacion(pf, t_initial, tstop_cpf)

    t1 = time.time()
    
    


    m_optm = run_optimization_model(pf, sim, params)
    print('=> Extract data time: %.4f (s)'% (t1-t0))

    Alternative_solutions(pf, m_optm)

    Revision_big_optm(params, project_file, pf, dict_lineas, t_initial, tstop_cpf, t_final, sim, m_optm)


    pf.app.ResetCalculation()


    # Extract data in each step of simulation
    # only for one scenario
    # Takes too long
    #pf.extract_data_each_step(pf, t_final, tstop_cpf, m_optm.part_factors[:,0,0])


    # VARIABILIDAD
    Variabilidad(pf, t_initial, t_final, tstop_cpf, sim, m_optm)




    if False:
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


        logging.info('Completed!')



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



if __name__ == "__main__":
    main()
