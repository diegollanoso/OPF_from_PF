import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import *
from scipy.sparse import csr_matrix as sparse, identity as sparseI
import logging


# Importar DIgSILENT
import sys
sys.path.append(r'C:\Program Files\DIgSILENT\PowerFactory 2023 SP3A\Python\3.11')
import powerfactory as pf

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# La matriz de shift factors se limpia y ordena
# Return matriz de shif factors y indices con nombres de elementos
def ShiftFactors(ptdf_dataframe):
    ejes = ptdf_dataframe.axes
    ejes_corregidos = list()

    # Elementos del SEP (NO Buses!) se encuentran dentro de otros objetos
    # ej: Talca - Itahue 66 kV L1\Tap San Rafael - Panguilemo  66 kV_bhk1
    # La línea se encuentra dentro de Talca - Itahue 66 kV L1
    # Ejes corregidos es el nombre del elemento sin incluir la carpeta de origen

    for i in ejes[1]:
        if '\\' in i:
            ejes_corregidos.append(i.split('\\')[1])
        else:
            ejes_corregidos.append(i)
    del ejes_corregidos[0:3]


    df = pd.DataFrame()
    j=0
    for i in ejes_corregidos:
        df[ejes_corregidos[j]] = (ptdf_dataframe[ejes_corregidos[j]]-ptdf_dataframe[ejes_corregidos[j+1]])/2
        j += 2
        if len(ejes_corregidos) == j:
            break


    ptdf_array_T = np.asarray(df)
    #ptdf_array_T = np.delete(ptdf_array_T ,[0,1,2], axis = 1)


    # PROBABLEMENTE!! Bastaría con este código
    #ptdf_array_T = np.asarray(ptdf_dataframe)
    #ptdf_array_T = np.delete(ptdf_array_T ,[0,1,2], axis = 1)
    
    return np.transpose(ptdf_array_T).astype(float), df.axes[1]


class PowerFactorySim(object):
    def __init__(self, project_name='Project'):
        logging.info("Initializing PowerFactory simulation for project: %s", project_name)
        # Start PowerFactory
        self.app = pf.GetApplication()
        self.app.Show()
        # Activate project
        self.app.ActivateProject(project_name)
        folder_scens = self.app.GetProjectFolder('scen')
        self.scens = folder_scens.GetContents()
        self.scens.sort(key=lambda x: x.loc_name[0])
        if self.scens:
            self.scens[0].Activate()

        # Objetos relevantes
        self.lineas = self.app.GetCalcRelevantObjects('*.ElmLne')
        self.generadores = self.app.GetCalcRelevantObjects('*.ElmSym') 
        self.cargas = self.app.GetCalcRelevantObjects('*.ElmLod')
        self.barras = self.app.GetCalcRelevantObjects('*.ElmTerm')
        self.trafos = self.app.GetCalcRelevantObjects('*.ElmTr2')
        self.trafos3 = self.app.GetCalcRelevantObjects('*.ElmTr3')
        self.capacitores = self.app.GetCalcRelevantObjects('*.ElmScap')
        self.inductores = self.app.GetCalcRelevantObjects('*.ElmSind')
        self.genstate = self.app.GetCalcRelevantObjects('*.ElmGenstat')
        self.asincronicos = self.app.GetCalcRelevantObjects('*.ElmAsm')
        self.shunt = self.app.GetCalcRelevantObjects('*.ElmShnt')
        self.pvsys = self.app.GetCalcRelevantObjects('*.ElmPvsys')


        self.ldf = self.app.GetFromStudyCase('ComLdf')
        self.ldf.iopt_net = 0
        
        if 'Articulo2' in project_name:
            self.dsl_agc_bloques = self.app.GetCalcRelevantObjects('agc_bloques.ElmDsl')[0]
        elif '3Bus_TS' == project_name:
            self.dsl_agc_bloques = self.app.GetCalcRelevantObjects('AGC Control.ElmDsl')[0]

        self.IntEvt = self.app.GetFromStudyCase('IntEvt')
        self.res = self.app.GetFromStudyCase('Results.ElmRes')
        self.events_folder = self.IntEvt.GetContents()

        try:
            study_case = self.app.GetProjectFolder('study').GetContents('Base Case.IntCase')[0]
        except:
            try:
                study_case = self.app.GetProjectFolder('study').GetContents('Study Case.IntCase')[0]
            except:
                logging.error('Error in StudyCase name')


        self.gen_agc = study_case.GetContents('Gen AGC.SetSelect')[0].All()
        self.Gen_AGC = list(map(lambda x: x.loc_name, self.gen_agc))

        gen_out = study_case.GetContents('Gen OUT.SetSelect')[0].All()
        self.Gen_Outages = list(map(lambda x: x.loc_name, gen_out))
        #self.Gen_Outages = list(map(lambda x: x.loc_name, gen_out))[:1]
        self.Ns = len(self.Gen_Outages)

        self.use_sf = False # 1 = SF calculados con matrices;  0 = SF de PF

        self.TS = False
        self.Nt = 3

        self.potencia_ac = 'm:Psum:bus1'
        self.potencia_ac2 = 'm:Psum:bus2'
        
        try:
            line_ts = study_case.GetContents('Line TS.SetSelect')[0].All()
            self.TS_lines = list(map(lambda x: x.loc_name, line_ts))
        except:
            print('Not Found: Lineas candidatas a TS')


        # Lista con N° gamma de cada generador en agc
        self.signal_list = list()
        for gen in self.gen_agc:
            comp = gen.c_pmod
            for slot in comp.pblk:
                ## Slot pcu
                if slot.loc_name[:9] == 'pcu Slot' or slot.loc_name[:9] == 'Gov Slot':
                    val = comp.pelm[comp.pblk.index(slot)].signal[3][-2:]
                    if not val[0].isdigit():
                        val = val[1]
                    self.signal_list.append('gamma' + val)

        self.Sb = 100
        self.raiz3 = 1.73205080757
        self.kV_fm = 10            # Tensiones menores, tendrán Fm = 1e9
        self.L = 11

    # Return pandas con matriz de Shift-factors
    # Return lista ordenada con barras 
    def export_csv(self, flujo:str = 'AC'):
        logging.info("Exporting CSV with flujo: %s", flujo)
        sfactors = self.app.GetFromStudyCase('ComVstab')
        if flujo == 'AC':
            sfactors.iopt_method = 0 # 0 = AC; 2 = DC
        if flujo == 'DC':
            sfactors.iopt_method = 2 # 0 = AC; 2 = DC
            self.ldf.iopt_net = 2
        #Results = sfactors.pResult
        #for i in Results.GetContents():
        #    i.Delete()
        #Line_Results = Results.CreateObject('IntMon','Line')
        #Line_Results.classnm = 'ElmLne'
        #Line_Results.AddVar('m:dPdP:bus1')
        #Line_Results.AddVar('m:dPdP:bus2')
        sfactors.frmElmFilt4Res = 0
        sfactors.frmLimitsBrc = 1
        sfactors.cdpflim = 0 
        sfactors.frmLimitsBus = 0
        
        sfactors.Execute()

        Results = self.app.GetFromStudyCase(sfactors.pResult.loc_name+'.ElmRes') # Archivo Resultados
        #print(Results.GetContents())
        export = self.app.GetFromStudyCase('ComRes') # Método para exportar 
        export.pResult = Results
        export.iopt_vars = 0
        export.iopt_tsel = 0
        export.iopt_exp = 6 #tipo de archivo a exportar
        export.f_name = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\SF_datos.csv'
        export.iopt_csel = 0 
        export.locn = 1
        export.iopt_sep=0
        export.iopt_locn = 1 
        export.ciopt_head = 0
        export.dec_Sep = '.'
        export.col_Sep = ';'
        export.Execute()

        ptdf_df = pd.read_csv(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\SF_datos.csv',skiprows = 0,delimiter=';')
        self.indices_bus = sfactors.p_bus.All()

        ejes = ptdf_df.axes
        ejes_corregidos = list()


        for i in ejes[1]:
            if '\\' in i:
                ejes_corregidos.append(i.split('\\')[1])
            else:
                ejes_corregidos.append(i)
        del ejes_corregidos[0:3]


        df = pd.DataFrame()
        j=0
        for i in ejes_corregidos:
            df[ejes_corregidos[j]] = (ptdf_df[ejes_corregidos[j]]-ptdf_df[ejes_corregidos[j+1]])/2
            j += 2
            if len(ejes_corregidos) == j:
                break

        ptdf_array_T = np.asarray(df)
        self.SF = np.transpose(ptdf_array_T).astype(float)
        self.indices_obj = df.axes[1]
        logging.info("CSV export complete")
        
    # Return diccionarios de los elementos 
    # bus, carga, linea, trafo, gen, genstat
    def get_data(self):
        logging.info("Getting data from PowerFactory simulation")
        self.ldf.Execute()
        self.dict_barras = dict()
        cont = 0
        for bus in self.indices_bus:
            self.dict_barras[bus.loc_name] = cont
            cont += 1

        # CARGAS - en servicio y switch cerrado
        # name carga = (N° barra, NombreBarra, Potencia peak MW, Potencia mean, Potencia valle)
        dict_cargas = dict()                              
        for c in self.cargas:
            if c.outserv == 0 and c.bus1.cpCB.on_off == 1:                  # carga en servicio y switch conectado
                # name carga = (N° barra, Potencia MW, NombreBarra)
                potencia = c.GetAttribute(self.potencia_ac)
                dict_cargas[c.loc_name] = [self.dict_barras[c.bus1.cterm.loc_name], c.bus1.cterm.loc_name, potencia, 0, 0]


        # LINEAS - en servicio
        # # name linea = (R, X, fmax, N° paralelas, from, to)
        dict_lineas = dict()
        for line in self.lineas:
            if line.outserv == 0:  #Linea en servicio
                #Zb = round(line.bus1.cterm.uknom*line.bus2.cterm.uknom/self.Sb,6)     
                Zb = line.bus1.cterm.uknom*line.bus2.cterm.uknom/self.Sb
                # Escoger flujo máximo de la línea cuando es menor a kV_fm
                if line.bus1.cterm.uknom < self.kV_fm:
                    fmax_kV = 1e9
                else:
                    fmax_kV = line.Inom_a*line.bus1.cterm.uknom*self.raiz3
                # name linea = (R, X, fmax, N° paralelas, bus1, bus2)
                dict_lineas[line.loc_name]=(line.R1/Zb, line.X1/Zb, fmax_kV, line.nlnum, line.bus1.cterm.loc_name, line.bus2.cterm.loc_name)

        # TRAFOS 2 DEV - en servicio
        # name trafo = (Xd, R, N° paralelo, fmax, from HV, to LV)
        dict_trafos = dict()
        for t in self.trafos:
            if t.typ_id.r1pu == 0:      # Checkear R1 trafo 2 dev
                r = 1e-8
            else:
                r = t.typ_id.r1pu/1000  # Cooper losses en kW
            
            if t.outserv == 0:      # Trafo en servicio
                if t.bushv.cterm.uknom < self.kV_fm:
                    fmax_kV = 1e9
                else:
                    fmax_kV = t.Snom_a
                
                Zb = t.bushv.cterm.uknom * t.bushv.cterm.uknom/self.Sb
                Zbt = t.bushv.cterm.uknom * t.bushv.cterm.uknom/ (t.Snom_a/t.ntnum)

                # name trafo = (Xd, R, N° paralelo, fmax)
                dict_trafos[t.loc_name] = (3*t.typ_id.x1pu*(Zbt/Zb), 3*r*(Zbt/Zb), t.ntnum, fmax_kV, t.bushv.cterm.loc_name, t.buslv.cterm.loc_name)

        # GENERADORES
        #name generador = (N° Barra, N° Gen paralelo,
        #                   Outserv, Pmin, Pmax, Pref,
        #                   costos var, costo fijo, rampa)
        self.dict_gen = dict()
        self.Name_all_gen = dict()
        for gen in self.generadores:
            self.Name_all_gen[gen.loc_name] = gen
            if gen.ip_ctrl == 1:    # Buscar barra slack
                self.bus_slack = gen.bus1.cterm.loc_name
            if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                potencia = gen.GetAttribute(self.potencia_ac)
                # name generador = (N° Barra, N° Gen paralelo,
                #                   Outserv, Pmin, Pmax,
                #                   costos var, costo fijo, rampa)
                self.dict_gen[gen.loc_name] = [self.dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          0, gen.Pmin_uc, gen.Pmax_uc,
                                          gen.penaltyCosts, gen.fixedCosts, gen.limRampUp,
                                          potencia, 0, 0]
            else:
                self.dict_gen[gen.loc_name] = [self.dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          1, gen.Pmin_uc, gen.Pmax_uc,
                                          gen.penaltyCosts, gen.fixedCosts, gen.limRampUp,
                                          0, 0, 0]

        # GENERADORES ESTATICOS
        # name generador = (N° Barra, N° Gen paralelo,
        #                   Outserv, Pmin, Pmax, Pref,
        #                   costos var, costo fijo, rampa)
        self.dict_genstat = dict()
        self.Name_all_genstat = dict()
        for gen in self.genstate:
            self.Name_all_genstat[gen.loc_name] = gen
            #if gen.ip_ctrl == 1:    # Buscar barra slack
            #    self.bus_slack = gen.bus1.cterm.loc_name
            if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                potencia = gen.GetAttribute(self.potencia_ac)
                # name generador = (N° Barra, N° Gen paralelo,
                #                   Outserv, Pmin, Pmax,
                #                   costos var, costo fijo, rampa)
                self.dict_genstat[gen.loc_name] = [self.dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          0, gen.Pmin_uc, gen.Pmax_uc,
                                          gen.penaltyCosts, gen.fixedCosts, gen.limRampUp,
                                          potencia, 0, 0]
            else:
                self.dict_genstat[gen.loc_name] = [self.dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          1, gen.Pmin_uc, gen.Pmax_uc,
                                          gen.penaltyCosts, gen.fixedCosts, gen.limRampUp,
                                          0, 0, 0]


        # PV generadores
        self.dict_pv = dict()
        self.Name_all_pv = dict()
        for pv in self.pvsys:
            self.Name_all_pv[pv.loc_name] = pv
            if pv.outserv == 0 and  pv.bus1.cpCB.on_off == 1:
                potencia = pv.GetAttribute(self.potencia_ac)
                # name generador = (N° Barra, N° Gen paralelo,
                #                   Outserv, Pmin, Pmax,
                #                   costos var, costo fijo, rampa)
                self.dict_pv[pv.loc_name] = [self.dict_barras[pv.bus1.cterm.loc_name], pv.ngnum,
                                          0, pv.Pmin_uc, pv.Pmax_uc,
                                          pv.penaltyCosts, pv.fixedCosts, pv.limRampUp,
                                          potencia, 0, 0]
            else:
                self.dict_pv[pv.loc_name] = [self.dict_barras[pv.bus1.cterm.loc_name], pv.ngnum,
                                          1, pv.Pmin_uc, pv.Pmax_uc,
                                          pv.penaltyCosts, pv.fixedCosts, pv.limRampUp,
                                          0, 0, 0]
        #FALTAN
        # TRAFOS 3 DEVANADOS
        # CAPACITORES SERIE
        # INDUCTORES SERIE
        # ASINCRONICOS
        # GEN VIRTUALES


        #### Listas de elementos
        
        nom_trf = list()
        nom_lin = list()

        nom_trf3 = list()
        nom_cap = list()
        nom_ind = list()
        



        for t in self.trafos:
            nom_trf.append(t.loc_name)

        for l in self.lineas:
            nom_lin.append(l.loc_name)

        for t in self.trafos3:
            nom_trf3.append(t.loc_name)

        for c in self.capacitores:
            nom_cap.append(c.loc_name)

        for i in self.inductores:
            nom_ind.append(i.loc_name)



        for sc in range(1,self.Nt):
            self.scens[sc].Activate()
            self.ldf.Execute()

                #Agregar potencias de carga en otros escenarios.

            for carga in self.cargas:
                dict_cargas[carga.loc_name][sc+2] = carga.GetAttribute(self.potencia_ac)

            for gen in self.generadores:
                if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                    potencia = gen.GetAttribute(self.potencia_ac)
                else:
                    potencia = 0
                self.dict_gen[gen.loc_name][sc+8] = potencia

            for gen in self.genstate:
                if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                    potencia = gen.GetAttribute(self.potencia_ac)
                else:
                    potencia = 0
                self.dict_genstat[gen.loc_name][sc+8] = potencia

            for pv in self.pvsys:
                if pv.outserv == 0 and  pv.bus1.cpCB.on_off == 1:
                    potencia = pv.GetAttribute(self.potencia_ac)
                else:
                    potencia = 0
                self.dict_pv[pv.loc_name][sc+8] = potencia


        # Demanda por barra
        self.Nb = len(self.indices_bus)
        self.dda_barra = np.zeros((self.Nb, self.Nt))
        for load in dict_cargas:
            for ti in range(self.Nt):
                self.dda_barra[dict_cargas[load][0], ti] += dict_cargas[load][ti+2]/self.Sb


        # Generadores
        self.ngen=len(self.dict_gen)
        self.n_gen_agc = len(self.Gen_AGC)

        self.pos_gen = np.zeros(self.ngen).astype(int)
        self.ngen_par = np.zeros(self.ngen)
        self.Cvar_gen = np.zeros(self.ngen)
        self.Ccte_gen = np.zeros(self.ngen)
        self.Pmin_gen = np.zeros(self.ngen)
        self.Pmax_gen =np.zeros(self.ngen)
        self.Gen_on_AGC = np.zeros(self.ngen)
        self.pos_gen_agc_list = list()
        self.name_gen_agc_list = list()
        self.Pgen_pre = np.zeros((self.ngen,self.Nt))
        self.Ramp_gen = np.zeros(self.ngen)
        self.all_gen = list()
        cont=0
        cont2 = 0
        for gen in self.dict_gen:
            self.all_gen.append(gen)
            self.pos_gen[cont] = self.dict_gen[gen][0]
            self.ngen_par[cont] = self.dict_gen[gen][1]
            self.Pmin_gen[cont] = self.dict_gen[gen][3]/self.Sb
            self.Pmax_gen[cont] = self.dict_gen[gen][4]/self.Sb
            self.Cvar_gen[cont] = self.dict_gen[gen][5]
            self.Ccte_gen[cont] = self.dict_gen[gen][6]
            self.Ramp_gen[cont] = self.dict_gen[gen][7]/self.Sb
            for ti in range(self.Nt):
                self.Pgen_pre[cont,ti] = self.dict_gen[gen][8+ti]/self.Sb
            if gen in self.Gen_AGC:
                self.Gen_on_AGC[cont] = 1
                self.pos_gen_agc_list.append(cont)
                self.name_gen_agc_list.append(gen)

                cont2 += 1

            cont+=1


        # Generadores estaticos

        self.ngenstat = len(self.dict_genstat)
        self.pos_genstat = np.zeros(self.ngenstat).astype(int)
        self.ngenstat_par = np.zeros(self.ngenstat)
        self.Cvar_genstat = np.zeros(self.ngenstat)
        self.Ccte_genstat = np.zeros(self.ngenstat)
        self.Pmin_genstat = np.zeros(self.ngenstat)
        self.Pmax_genstat =np.zeros(self.ngenstat)
        self.Pgenstat_pre = np.zeros((self.ngenstat,self.Nt))
        self.Ramp_genstat = np.zeros(self.ngenstat)
        self.all_genstat = list()
        cont=0
        for gen in self.dict_genstat:
            self.all_genstat.append(gen)
            self.pos_genstat[cont] = self.dict_genstat[gen][0]
            self.ngenstat_par[cont] = self.dict_genstat[gen][1]
            self.Pmin_genstat[cont] = self.dict_genstat[gen][3]/self.Sb
            self.Pmax_genstat[cont] = self.dict_genstat[gen][4]/self.Sb
            self.Cvar_genstat[cont] = self.dict_genstat[gen][5]
            self.Ccte_genstat[cont] = self.dict_genstat[gen][6]
            self.Ramp_genstat[cont] = self.dict_genstat[gen][7]/self.Sb
            for ti in range(self.Nt):
                self.Pgenstat_pre[cont,ti] = self.dict_genstat[gen][8+ti]/self.Sb
            cont+=1

        # PV Generadores

        self.ngenpv = len(self.dict_pv)
        self.pos_pv = np.zeros(self.ngenpv).astype(int)
        self.ngenpv_par = np.zeros(self.ngenpv)
        self.Cvar_pv = np.zeros(self.ngenpv)
        self.Ccte_pv = np.zeros(self.ngenpv)
        self.Pmin_pv = np.zeros(self.ngenpv)
        self.Pmax_pv =np.zeros(self.ngenpv)
        self.Ppv_pre = np.zeros((self.ngenpv,self.Nt))
        self.Ramp_pv = np.zeros(self.ngenpv)
        self.all_pv = list()
        cont=0
        for pv in self.dict_pv:
            self.all_pv.append(pv)
            self.pos_pv[cont] = self.dict_pv[pv][0]
            self.ngenpv_par[cont] = self.dict_pv[pv][1]
            self.Pmin_pv[cont] = self.dict_pv[pv][3]/self.Sb
            self.Pmax_pv[cont] = self.dict_pv[pv][4]/self.Sb
            self.Cvar_pv[cont] = self.dict_pv[pv][5]
            self.Ccte_pv[cont] = self.dict_pv[pv][6]
            self.Ramp_pv[cont] = self.dict_pv[pv][7]/self.Sb
            for ti in range(self.Nt):
                self.Ppv_pre[cont,ti] = self.dict_pv[pv][8+ti]/self.Sb
            cont+=1

        self.n_elem = len(self.indices_obj)   # N° de líneas + 'trf2'
        # Obtención de fmax, R y X
        dict_full = dict()
        self.pos_line = list()
        self.all_line = list()
        self.FMax = np.zeros(self.n_elem)
        self.R = np.zeros(self.n_elem)
        self.X = np.zeros(self.n_elem)
        self.i_buses = np.zeros(self.n_elem).astype(int)
        self.j_buses = np.zeros(self.n_elem).astype(int)
        self.pos_ts = list()
        self.all_branch = list()
        cont = -1
        for i in self.indices_obj:
            cont += 1
            if i in dict_lineas:
                self.all_line.append(i)
                dict_full[i] = dict_lineas[i]
                Fmax_i = dict_lineas[i][2]
                R_i = dict_lineas[i][0]
                X_i = dict_lineas[i][1]
                i_bus = self.dict_barras[dict_lineas[i][4]]
                j_bus = self.dict_barras[dict_lineas[i][5]]
                self.pos_line.append(cont)
                if self.TS:
                    if i in self.TS_lines:
                        self.pos_ts.append(cont)

            elif i in dict_trafos:
                dict_full[i] = dict_trafos[i]
                Fmax_i = dict_trafos[i][3]
                R_i = dict_trafos[i][1]
                X_i = dict_trafos[i][0]
                i_bus = self.dict_barras[dict_trafos[i][4]]
                j_bus = self.dict_barras[dict_trafos[i][5]]
            else:
                print('Elemento no encontrado en diccs: ' + i)
            self.all_branch.append(i)
            self.FMax[cont] = Fmax_i/self.Sb
            self.R[cont] = R_i
            self.X[cont] = X_i
            self.i_buses[cont] = i_bus
            self.j_buses[cont] = j_bus
        
        if self.TS:
            self.pos_nots = [pos for pos in list(range(len(self.indices_obj))) if pos not in self.pos_ts]
            # Lineas Candidatass TS
            self.i_buses_ts = self.i_buses[self.pos_ts]
            self.j_buses_ts = self.j_buses[self.pos_ts]

        I = np.r_[range(self.n_elem), range(self.n_elem)]
        S = sparse((np.r_[np.ones(self.n_elem), -np.ones(self.n_elem)], (I, np.r_[self.i_buses, self.j_buses])), (self.n_elem, self.Nb))
        self.A = np.array(S.todense())


        val_min = 1e-24 #Valor mínimo para evitar valores nulos en matriz X
        yprim = np.zeros(self.n_elem).astype(complex)
        for i in range(self.n_elem):
            if self.X[i] == 0:
                self.X[i] = val_min
            yprim[i] = 1/(complex(self.R[i],self.X[i]))
        self.G = np.real(yprim)
        self.B = np.imag(yprim)

        ## CALCULO DE SHIFT FACTORS
        if self.use_sf:
            n_bus_slack = list(self.dict_barras.keys()).index(self.bus_slack)
            self.noslack = np.flatnonzero(np.arange(self.Nb) != n_bus_slack)

            BfR = sparse((np.r_[np.imag(yprim), -np.imag(yprim)], (I, np.r_[self.i_buses, self.j_buses])), (self.n_elem,self.Nb))
            BbusR = S.T * BfR
            SFR = np.zeros((self.n_elem,self.Nb))
            SFR[:,self.noslack] = BfR[:, self.noslack].todense()*np.linalg.inv(BbusR[np.ix_(self.noslack, self.noslack)].todense())    

            self.SF = SFR

            #SF_pf = pd.DataFrame(np.vstack((self.name_gen_agc_list,self.SF[:,self.pos_gen[self.pos_gen_agc_list]])).T, columns=np.insert(self.indices_obj, 0, 'x', axis=0))
            #SF_pf.to_excel(r"C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\SF_gens.xlsx", index = False)


        #return [nom_bus, nom_trf, nom_lin, nom_trf3, nom_cap, nom_ind]
        logging.info("Data extraction complete")
        return (dict_lineas)



    def prepare_dynamic_sim(self, monitored_variables,
                            sim_type='rms', start_time=0.0,
                            step_size=0.01):
        logging.info("Preparing dynamic simulation with sim_type: %s", sim_type)
        # get result file
        # select results variables to monitor
        for elm_name, var_names in monitored_variables.items():
            # get all network elements that match 'elm_name'
            elements = self.app.GetCalcRelevantObjects(elm_name)
            # select variables to monitor for each element
            for element in elements:
                self.res.AddVars(element, *var_names)
        # retrieve initial conditions and time domain sim. objects
        self.inc = self.app.GetFromStudyCase('ComInc')
        self.sim = self.app.GetFromStudyCase('ComSim')
        # set simulation type: 'rms' or 'ins' (for EMT)
        self.inc.iopt_sim = sim_type
        # set start time, step size and end time
        self.inc.tstart = start_time
        self.inc.dtgrd = step_size
        # set initial conditions
        self.inc.Execute()
        logging.info("Dynamic simulation prepared")
    
    def run_dynamic_sim(self, end_sim = 49.9):
        self.sim.tstop = end_sim
        logging.info("Running dynamic simulation until time: %.2f", end_sim)
        result = bool(self.sim.Execute())
        logging.info("Dynamic simulation complete with result: %s", result)
        return result


    def get_dynamic_results(self, elm_name, var_name):
        logging.info("Getting dynamic results for element: %s, variable: %s", elm_name, var_name)
        # get network element of interest
        element = self.app.GetCalcRelevantObjects(elm_name)[0]
        # load results from file
        self.app.ResLoadData(self.res)
        # find colum in results file that holds result of interest
        col_index = self.app.ResGetIndex(
            self.res, element, var_name)
        # get number of rows (points in time) in the result file
        n_rows = self.app.ResGetValueCount(self.res, 0)
        # read results and time and store them in lists
        time = []
        var_values = []
        for i in range(n_rows):
            time.append(self.app.ResGetData(self.res, i, -1)[1])
            var_values.append(
                self.app.ResGetData(self.res, i, col_index)[1])
        return time, var_values


    # Check results
    def check_results(self, gen_eff, genstat_eff, p_g, p_gstat, f_line):
        logging.info("Checking results")
        cont = 0
        for gen in self.generadores:
            if gen.loc_name in gen_eff:
                gen.pgini = p_g[cont].x*self.Sb
                print('%s => %.2f (MW)' % (gen.loc_name,gen.ngnum*p_g[cont].X*self.Sb))
                cont += 1

        cont = 0
        for gen in self.genstate:
            if gen.loc_name in genstat_eff:
                gen.pgini = p_gstat[cont].x*self.Sb
                print('%s => %.2f (MW)' % (gen.loc_name,gen.ngnum*p_gstat[cont].X*self.Sb))
                cont += 1

        
        self.ldf.Execute()
        cont=0
        for line in self.lineas:
            if line.loc_name in f_line[0]:
                print('%s ; Gurobi (DC) = %.3f [MW]; DIgSILENT (AC) = %.3f [MW]' % (line.loc_name, float(f_line[1][np.where(f_line[0] == line.loc_name)])*self.Sb, (-line.GetAttribute(self.potencia_ac2)+line.GetAttribute(self.potencia_ac))/2))
                cont += 1
        logging.info("Results check complete")


    def ChangeMaxLine(self, obj:list, value:list):
        logging.info("Changing max line values for objects: %s with values: %s", obj, value)
        for line in self.lineas:
            if line.loc_name in obj:
                line.typ_id.sline = value[obj.index(line.loc_name)]/1000
                

    def extract_data(self, elem: str, variable: str, return_time: bool = True, start_time: float = None):
        logging.info("Extracting data for element: %s, variable: %s", elem, variable)
        element = self.app.GetCalcRelevantObjects(elem)[0]
        self.app.ResLoadData(self.res)
        col_index = self.app.ResGetIndex(self.res, element, variable)
        n_rows = self.app.ResGetValueCount(self.res, 0)
        time = []
        var_values = []
        for i in range(n_rows):
            current_time = self.app.ResGetData(self.res, i, -1)[1]
            # skip data points before start_time
            if (start_time is not None) and (current_time < start_time):
                continue
            if return_time:
                time.append(current_time)
            var_values.append(self.app.ResGetData(self.res, i, col_index)[1])
        logging.info("Data extraction complete")
        #if variable == 's:pt':
        #    p_nominal = element.typ_id.sgn
        #    var_values = np.array(var_values)*p_nominal

        if return_time:
            return time, var_values
        else:
            return var_values


class Simulacion(object):
    def __init__(self, data, t_initial=0, tstop_cpf=30):
        logging.info("Initializing simulation with t_initial: %d, tstop_cpf: %d", t_initial, tstop_cpf)
        # Elementos a monitorear
        bus_freq = data.app.GetCalcRelevantObjects('Term_10_4.ElmTerm')
        #Inercia
        H=0.94
    
        # Switch event
        # Gen Out
        self.name_events = list()
        for e in data.events_folder:
            self.name_events.append(e.loc_name)
            e.outserv = 1   

        self.D_pfc = np.zeros((data.Nb, data.Ns, data.Nt)) # Demanda al finaliza el CPF
        self.PL_pre_line = np.zeros(data.Nt) # Perdidas antes de la contingencia
        #PL_pre_trafo2 = np.zeros((len(dict_trafos),Ns))
        self.Pgen_pfc = np.zeros((data.ngen, data.Ns, data.Nt))
        self.Pgenstat_pfc = np.zeros((data.ngenstat, data.Ns, data.Nt))
        self.PgenPV_pfc = np.zeros((data.ngenpv, data.Ns, data.Nt))
        self.Gen_on_AGC = np.tile(data.Gen_on_AGC,(data.Ns,1)).T
        self.Gen_on_AGC = np.stack([self.Gen_on_AGC]*3,axis=0)
        #self.Genstat_on_AGC = np.tile(data.Genstat_on_AGC,(Ns,1)).T
        self.P_out = np.zeros((data.Ns, data.Nt))
        self.P_out_f = np.zeros((data.Ns, data.Nt))
        self.Barra_gen_out = np.zeros((data.Ns))
        for ti in range(data.Nt):
            data.scens[ti].Activate()
            print('Escenario: ' + str(ti+1))
            cont = -1
            for gen_out in data.Gen_Outages:
                cont += 1
                sys.stdout.write('\r')
                sys.stdout.write('Progreso: ' + str(cont+1+data.Ns*ti) + '/' + str(data.Ns*data.Nt))
                sys.stdout.write(' - Generador Out: ' + gen_out)
                sys.stdout.flush()
                
                self.P_out[cont,ti] = data.dict_gen[gen_out][8+ti]/data.Sb
                self.Barra_gen_out[cont] = data.dict_gen[gen_out][0]
                
                evt = data.events_folder[self.name_events.index('Salida Gen')]
                evt.outserv = 0
                evt.time = t_initial
                evt.p_target = data.generadores[list(data.dict_gen).index(gen_out)]

                data.prepare_dynamic_sim({}, 'rms')
                
                if self.PL_pre_line[ti] == 0: # Obtener perdidas solo una vez para cada contingencia y cada operación demanda
                    self.PL_pre_line[ti] = np.array(list(map(lambda x: x.GetAttribute('c:Losses'), data.lineas))).sum()/1000/data.Sb
                
                data.run_dynamic_sim(end_sim=tstop_cpf)

                for load in list(map(lambda x: (x.GetAttribute(data.potencia_ac),data.dict_barras[x.bus1.cterm.loc_name]), data.cargas)):  # Demanda al final del CPF por carga
                    self.D_pfc[load[1], cont, ti] += load[0]/data.Sb

                #PL_pre_trafo2[:,cont] = np.array(list(map(lambda x: x.GetAttribute('c:Ploss'), data.trafos)))  


                for gen in data.generadores:
                    if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                        potencia = gen.GetAttribute(data.potencia_ac)/data.Sb
                    else:
                        potencia = 0
                    self.Pgen_pfc[list(data.all_gen).index(gen.loc_name),cont, ti] = potencia

                for gen in data.genstate:
                    if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                        potencia = gen.GetAttribute(data.potencia_ac)/data.Sb
                    else:
                        potencia = 0
                    self.Pgenstat_pfc[list(data.all_genstat).index(gen.loc_name),cont, ti] = potencia

                for pv in data.pvsys:
                    if pv.outserv == 0 and  pv.bus1.cpCB.on_off == 1:
                        potencia = pv.GetAttribute(data.potencia_ac)/data.Sb
                    else:
                        potencia = 0
                    self.PgenPV_pfc[list(data.all_pv).index(pv.loc_name),cont, ti] = potencia

                #Pgen_pfc[:,cont] = np.array(list(map(lambda x: x.GetAttribute(data.potencia_ac), data.generadores)))/data.Sb  # Potencia generadores al final del CPF
                #Pgenstat_pfc[:,cont,ti] = np.array(list(map(lambda x: x.GetAttribute(data.potencia_ac), data.genstate)))/data.Sb  # Potencia generadores al final del CPF

                #if (gen_out in data.Gen_AGC) and (gen_out in data.dict_gen):
                #    self.Gen_on_AGC[ti, list(data.all_gen).index(gen_out),cont] = 0

                #elif (gen_out in Gen_AGC) and (gen_out in dict_genstat_eff):
                #    self.Genstat_on_AGC[ti, list(dict_genstat_eff.keys()).index(gen_out),cont] = 0
                
                
                #data.app.ResLoadData(data.res)
                #row_time = data.app.ResGetValueCount(data.res,0)
                #col_freq = data.app.ResGetIndex(data.res, bus_freq[0], 'm:fehz')

                #time = []
                #freq_value = []
                #for i in range(row_time):
                #    time.append(data.app.ResGetData(data.res, i, -1)[1])
                #    freq_value.append(data.app.ResGetData(data.res, i, col_freq)[1])


                #step_size = data.app.GetFromStudyCase('ComInc').dtgrd

                #Rocof a las 500ms
                #delta_t = int(0.5/step_size)
                #rocof = (freq_value[t_initial+delta_t]-freq_value[t_initial])/0.5

                #self.P_out_f[cont,ti] = rocof

                #self.P_out_f[cont,ti] = -2*H*rocof
        logging.info("Simulation initialized")






class ShortSim(object):
    def __init__(self, data, big_optm, simm, new_SF, gen_out, ti):
        logging.info("Initializing short simulation for gen_out: %s, ti: %d", gen_out, ti)
        #count=0
        #for gen in pf.Gen_AGC:
        #    evt2 = simm.events_folder[simm.name_events.index(signal_list[count])]
        #    evt2.outserv = 0
        #    evt2.value = str(big_optm.part_factors[:,scen,ti][pf.name_gen_agc_list.index(gen)])
        #    count += 1

        data.prepare_dynamic_sim({}, 'rms')
        data.run_dynamic_sim(end_sim=ti)

        # Potencia en cada generador
        self.p_gen = np.zeros(data.ngen)
        self.p_genAGC = 0
        i=0
        for gen in data.generadores:
            self.p_gen[i] = gen.GetAttribute('m:Psum:bus1')/data.Sb
            i += 1
            if gen.loc_name in big_optm.gen_csf and gen.loc_name != gen_out:
                self.p_genAGC += gen.GetAttribute('m:Psum:bus1')/data.Sb

        self.p_genstat = np.zeros(data.ngenstat)
        i=0
        for gen in data.genstate:
            self.p_genstat[i] = gen.GetAttribute('m:Psum:bus1')/data.Sb
            i += 1
        
        self.pv_gen = np.zeros(data.ngenpv)
        i=0
        for pv in data.pvsys:
            self.pv_gen[i] = pv.GetAttribute('m:Psum:bus1')/data.Sb
            i += 1


        # Demanda en cada barra
        self.D_t = np.zeros(data.Nb)
        for load in list(map(lambda x: (x.GetAttribute('m:Psum:bus1'),data.dict_barras[x.bus1.cterm.loc_name]), data.cargas)):  # Demanda al final del CPF por carga
            self.D_t[load[1]] += load[0]/data.Sb

        # Perdidas

        self.PL = np.zeros(new_SF.new_n_elem)
        for line in data.app.GetCalcRelevantObjects('*.ElmLne'):
            if line.loc_name in new_SF.all_branch:
                perdida = line.GetAttribute('c:Losses')/1000/data.Sb
                self.PL[new_SF.all_branch.index(line.loc_name)] = perdida
        for trafo in data.app.GetCalcRelevantObjects('*.ElmTr2'):
            p_hv = trafo.GetAttribute('m:P:bushv')/data.Sb
            p_lv = trafo.GetAttribute('m:P:buslv')/data.Sb
            perdida = abs(p_hv + p_lv)
            self.PL[new_SF.all_branch.index(trafo.loc_name)] = perdida
        logging.info("Short simulation initialized")


class new_SF(object):
    def __init__(self, pf, m_optm, ti, scen):
        logging.info("Initializing new shift factors for ti: %d, scen: %d", ti, scen)
        N_ts = int(sum(1-m_optm.s_ts.x[:,scen,ti]))

        self.new_n_elem = pf.n_elem - N_ts # N° de líneas + trf2 - Líneas TS
        new_R = np.zeros(self.new_n_elem)
        new_X = np.zeros(self.new_n_elem)
        new_i_buses = np.zeros(self.new_n_elem)
        new_j_buses = np.zeros(self.new_n_elem)
        self.new_FMax = np.zeros(self.new_n_elem)

        self.all_branch = list()

        cont=0
        j=-1
        for i in range(pf.n_elem):
            if i in pf.pos_line:
                j += 1
                if m_optm.s_ts[j,scen,ti].x:
                    new_R[cont] = pf.R[i]
                    new_X[cont] = pf.X[i]
                    new_i_buses[cont] = pf.i_buses[i]
                    new_j_buses[cont] = pf.j_buses[i]
                    self.new_FMax[cont] = pf.FMax[i]
                    self.all_branch.append(pf.all_branch[i])
                    cont+=1
            else:
                new_R[cont] = pf.R[i]
                new_X[cont] = pf.X[i]
                new_i_buses[cont] = pf.i_buses[i]
                new_j_buses[cont] = pf.j_buses[i]
                self.new_FMax[cont] = pf.FMax[i]
                self.all_branch.append(pf.all_branch[i])
                cont+=1

        I = np.r_[range(self.new_n_elem), range(self.new_n_elem)]
        S = sparse((np.r_[np.ones(self.new_n_elem), -np.ones(self.new_n_elem)], (I, np.r_[new_i_buses, new_j_buses])), (self.new_n_elem, pf.Nb))
        self.A = np.array(S.todense())


        val_min = 1e-24 #Valor mínimo para evitar valores nulos en matriz X
        yprim = np.zeros(self.new_n_elem).astype(complex)
        for i in range(self.new_n_elem):
            if new_X[i] == 0:
                new_X[i] = val_min
            yprim[i] = 1/(complex(new_R[i],new_X[i]))
        self.G = np.real(yprim)
        self.B = np.imag(yprim)


        BfR = sparse((np.r_[np.imag(yprim), -np.imag(yprim)], (I, np.r_[new_i_buses, new_j_buses])), (self.new_n_elem,pf.Nb))
        BbusR = S.T * BfR
        SFR = np.zeros((self.new_n_elem,pf.Nb))
        SFR[:,pf.noslack] = BfR[:, pf.noslack].todense()*np.linalg.inv(BbusR[np.ix_(pf.noslack, pf.noslack)].todense())    

        self.SF = SFR
        logging.info("New shift factors initialized")

def CreateEvents_line(data, optm, scen, ti, tstop_cpf):
    logging.info("Creating events for lines with scen: %d, ti: %d", scen, ti)
    for line in data.lineas:
        if line.loc_name in data.TS_lines:
            index = data.all_line.index(line.loc_name)
            if optm.s_ts.x[index,scen,ti] == 0:
                event_folder = data.IntEvt.GetContents()
                lista_eventos = list(map(lambda x: x.loc_name , event_folder))
                name_evt = 'Evt '  + str(line.loc_name)
                if name_evt in lista_eventos:
                    position = lista_eventos.index(name_evt)
                    evento = data.events_folder[position]
                    evento.outserv = 0
                    evento.time = tstop_cpf
                else:
                    evento = data.IntEvt.CreateObject('EvtOutage', 'Evt ' + str(line.loc_name))   
                    evento.p_target = line
                    evento.i_what = 0
                    evento.time = tstop_cpf 
    logging.info("Events for lines created")


def CreateEvents_gen(data, optm, gen, value, ti):
    logging.info("Creating events for generator: %s, ti: %d", gen, ti)
    cont=0
    for gen in data.name_gen_agc_list:
        name_event =  'Evt Gamma ' + str(gen) + ' - ' + str(ti)
        if not name_event in data.IntEvt.GetContents():
            evento = data.IntEvt.CreateObject('EvtParam', 'Evt Gamma ' + str(gen) + ' - ' + str(ti))    
            evento.time = ti
            evento.mtime = ti//60
            print('Evento en minuto: ' + str(ti//60))
            if ti//60 != 0:
                evento.time = ti%60
                print('Evento en segundo: ' + str(ti%60))
            print('Evento en segundo: ' + str(ti%60))

            evento.p_target = data.Name_all_gen[gen]
        cont+=1
    logging.info("Events for generator created")

def Set_param_agc(pf, t_int, part_factors, previous_part_factors=None):
    logging.info("Setting AGC parameters for t_int: %d", t_int)
    if previous_part_factors is None:
        previous_part_factors = [0.0] * len(part_factors)

    cont=0
    for gen in pf.gen_agc:
        name_event =  'Evt Gamma ' + str(gen.loc_name) + ' - ' + str(t_int)
        posicion = pf.name_gen_agc_list.index(gen.loc_name)
        event_folder = pf.IntEvt.GetContents()
        lista_eventos = list(map(lambda x: x.loc_name , event_folder))
        if  (float(part_factors[posicion]) != 0.0 and
            float(part_factors[posicion]) != float(previous_part_factors[posicion])):

            if (name_event not in lista_eventos):
                evento = pf.IntEvt.CreateObject('EvtParam', 'Evt Gamma ' + str(gen.loc_name) + ' - ' + str(t_int))    
                evento.variable = pf.signal_list[cont]

            elif (name_event in lista_eventos):
                evento = event_folder[lista_eventos.index(name_event)]
                evento.outserv = 0
            
            evento.variable = pf.signal_list[cont]
            evento.time = t_int
            evento.mtime = t_int//60
            if t_int//60 != 0:
                evento.time = t_int%60
            evento.p_target = pf.dsl_agc_bloques
        cont+=1
    logging.info("AGC parameters set")

#class LpGurobi(object):
#    def __init__(self):
#        self.m = gp.Model('Modelo 1')
#        self.m.Params.MIPGap = 1e-5
#        self.m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

