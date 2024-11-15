import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import *
from scipy.sparse import csr_matrix as sparse, identity as sparseI


# Importar DIgSILENT
import sys
sys.path.append(r'C:\Program Files\DIgSILENT\PowerFactory 2023 SP3A\Python\3.11')
import powerfactory as pf


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
        # Start PowerFactory
        self.app = pf.GetApplication()
        #self.app.Show()
        # Activate project
        self.app.ActivateProject(project_name)
        folder_scens = self.app.GetProjectFolder('scen')
        self.scens = folder_scens.GetContents()
        self.scens.sort(key=lambda x: x.loc_name[0])
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
        
        self.ldf = self.app.GetFromStudyCase('ComLdf')
        self.ldf.iopt_net = 0

        self.dsl_agc_bloques = self.app.GetCalcRelevantObjects('agc_bloques.ElmDsl')[0]
        self.IntEvt = self.app.GetFromStudyCase('IntEvt')
        self.res = self.app.GetFromStudyCase('Results.ElmRes')
        self.events_folder = self.IntEvt.GetContents()

        try:
            study_case = self.app.GetProjectFolder('study').GetContents('Base Case.IntCase')[0]
        except:
            try:
                study_case = self.app.GetProjectFolder('study').GetContents('Study Case.IntCase')[0]
            except:
                print('Error en nombre StudyCase')


        self.gen_agc = study_case.GetContents('Gen AGC.SetSelect')[0].All()
        self.Gen_AGC = list(map(lambda x: x.loc_name, self.gen_agc))

        gen_out = study_case.GetContents('Gen OUT.SetSelect')[0].All()
        self.Gen_Outages = list(map(lambda x: x.loc_name, gen_out))
        #self.Gen_Outages = list(map(lambda x: x.loc_name, gen_out))[-1:]
        self.Ns = len(self.Gen_Outages)

        self.use_sf = False

        self.TS = False
        self.Nt = 3

        self.flujo_dc = 0
        self.potencia_dc = 'm:P:bus1'
        self.potencia_dc2 = 'm:P:bus2'
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
        self.kV_fm = 110            # Tensiones menores, tendrán Fm = 1e9


    # Return pandas con matriz de Shift-factors
    # Return lista ordenada con barras 
    def export_csv(self, flujo:str = 'AC'):
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
        
    # Return diccionarios de los elementos 
    # bus, carga, linea, trafo, gen, genstat
    def get_data(self):
        if self.flujo_dc:
            self.potencia_ac = self.potencia_dc
            self.potencia_ac2 = self.potencia_dc2
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
                r = t.typ_id.r1pu
            
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
        dict_genstat = dict()
        for gen in self.genstate:
            if gen.ip_ctrl == 1:    # Buscar barra slack
                self.bus_slack = gen.bus1.cterm.loc_name
            if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                # name generador = (N° Barra, N° Gen paralelo,
                #                   Outserv, Pmin, Pmax, Pref,
                #                   costos var, costo fijo, rampa)
                dict_genstat[gen.loc_name] = (self.dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          0, gen.Pmin_uc, gen.Pmax_uc, gen.pgini,
                                          gen.penaltyCosts, gen.fixedCosts, gen.limRampUp)
            else:
                dict_genstat[gen.loc_name] = (self.dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          1, gen.Pmin_uc, gen.Pmax_uc, 0,
                                          gen.penaltyCosts, gen.fixedCosts, gen.limRampUp)

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
        self.Cvar_gen = np.zeros(self.n_gen_agc)
        self.Ccte_gen = np.zeros(self.n_gen_agc)
        self.Pmin_gen = np.zeros(self.n_gen_agc)
        self.Pmax_gen =np.zeros(self.n_gen_agc)
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
            self.Ramp_gen[cont] = self.dict_gen[gen][7]/self.Sb
            for ti in range(self.Nt):
                self.Pgen_pre[cont,ti] = self.dict_gen[gen][8+ti]/self.Sb
            if gen in self.Gen_AGC:
                self.Gen_on_AGC[cont] = 1
                self.pos_gen_agc_list.append(cont)
                self.name_gen_agc_list.append(gen)
                self.Pmin_gen[cont2] = self.dict_gen[gen][3]/self.Sb
                self.Pmax_gen[cont2] = self.dict_gen[gen][4]/self.Sb
                self.Cvar_gen[cont2] = self.dict_gen[gen][5]
                self.Ccte_gen[cont2] = self.dict_gen[gen][6]
                cont2 += 1

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
            if self.flujo_dc:
                self.R[cont] = 0
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
        return (dict_lineas)



    def prepare_dynamic_sim(self, monitored_variables,
                            sim_type='rms', start_time=0.0,
                            step_size=0.01):
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
    
    def run_dynamic_sim(self, end_sim = 49.9):
        self.sim.tstop = end_sim
        return bool(self.sim.Execute())


    def get_dynamic_results(self, elm_name, var_name):
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


    def ChangeMaxLine(self, obj:list, value:list):
        for line in self.lineas:
            if line.loc_name in obj:
                line.typ_id.sline = value[obj.index(line.loc_name)]/1000
                

class Simulacion(object):
    def __init__(self, data, t_initial=0, tstop_cpf=30):
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
        self.PL_pre_line = np.zeros(data.Nt) # Perdidas al finalizar el CPF
        #PL_pre_trafo2 = np.zeros((len(dict_trafos),Ns))
        self.Pgen_pfc = np.zeros((data.ngen, data.Ns, data.Nt))
        #self.Pgenstat_pfc = np.zeros((ngenstat_eff, Ns, Nt))
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
                
                if data.flujo_dc:
                    continue

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
                    self.Pgenstat_pfc[list(data.all_gen).index(gen.loc_name),cont, ti] = gen.GetAttribute(data.potencia_ac)/data.Sb

                #Pgen_pfc[:,cont] = np.array(list(map(lambda x: x.GetAttribute(data.potencia_ac), data.generadores)))/data.Sb  # Potencia generadores al final del CPF
                #Pgenstat_pfc[:,cont,ti] = np.array(list(map(lambda x: x.GetAttribute(data.potencia_ac), data.genstate)))/data.Sb  # Potencia generadores al final del CPF

                if (gen_out in data.Gen_AGC) and (gen_out in data.dict_gen):
                    self.Gen_on_AGC[ti, list(data.all_gen).index(gen_out),cont] = 0

                #elif (gen_out in Gen_AGC) and (gen_out in dict_genstat_eff):
                #    self.Genstat_on_AGC[ti, list(dict_genstat_eff.keys()).index(gen_out),cont] = 0
                data.app.ResLoadData(data.res)
                row_time = data.app.ResGetValueCount(data.res,0)
                col_freq = data.app.ResGetIndex(data.res, bus_freq[0], 'm:fehz')

                time = []
                freq_value = []
                for i in range(row_time):
                    time.append(data.app.ResGetData(data.res, i, -1)[1])
                    freq_value.append(data.app.ResGetData(data.res, i, col_freq)[1])


                step_size = data.app.GetFromStudyCase('ComInc').dtgrd

                #Rocof a las 500ms
                delta_t = int(0.5/step_size)
                rocof = (freq_value[t_initial+delta_t]-freq_value[t_initial])/0.5

                self.P_out_f[cont,ti] = rocof

                #self.P_out_f[cont,ti] = -2*H*rocof

                print('x')

class ShortSim(object):
    def __init__(self, data, big_optm, simm, gen_out, ti):
        #count=0
        #for gen in pf.Gen_AGC:
        #    evt2 = simm.events_folder[simm.name_events.index(signal_list[count])]
        #    evt2.outserv = 0
        #    evt2.value = str(big_optm.part_factors[:,scen,ti][pf.name_gen_agc_list.index(gen)])
        #    count += 1

        data.prepare_dynamic_sim({}, 'rms')
        data.run_dynamic_sim(end_sim=ti)

        # Potencia en cada generador
        self.p_gen = np.zeros(len(data.generadores))
        i=0
        for gen in data.generadores:
            self.p_gen[i] = gen.GetAttribute('m:Psum:bus1')/data.Sb
            i += 1

        # Demanda en cada barra
        self.D_t = np.zeros(data.Nb)
        for load in list(map(lambda x: (x.GetAttribute('m:Psum:bus1'),data.dict_barras[x.bus1.cterm.loc_name]), data.cargas)):  # Demanda al final del CPF por carga
            self.D_t[load[1]] += load[0]/data.Sb




class new_SF(object):
    def __init__(self, pf, m_optm, ti, scen):
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

def CreateEvents_line(data, optm, scen,ti):
    for line in data.lineas:
        if line.loc_name in data.TS_lines:
            index = data.all_line.index(line.loc_name)
            if optm.s_ts.x[index,scen,ti] == 0:
                event_folder = data.IntEvt.GetContents()
                lista_eventos = list(map(lambda x: x.loc_name , event_folder))
                name_evt = 'Evt '  + str(line.loc_name)
                if name_evt in lista_eventos:
                    position = lista_eventos.index(name_evt)
                    data.events_folder[position].outserv = 0
                else:
                    evento = data.IntEvt.CreateObject('EvtOutage', 'Evt ' + str(line.loc_name))   
                    evento.p_target = line
                    evento.i_what = 0
                    evento.time = 30 


def CreateEvents_gen(data, optm, gen, value, ti):
    cont=0
    for gen in data.name_gen_agc_list:
        name_event =  'Evt Gamma ' + str(gen) + ' - ' + str(ti)
        if not name_event in data.IntEvt.GetContents():
            evento = data.IntEvt.CreateObject('EvtParam', 'Evt Gamma ' + str(gen) + ' - ' + str(ti))    
            evento.time = ti
            evento.p_target = data.Name_all_gen[gen]
        cont+=1


#class LpGurobi(object):
#    def __init__(self):
#        self.m = gp.Model('Modelo 1')
#        self.m.Params.MIPGap = 1e-5
#        self.m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

