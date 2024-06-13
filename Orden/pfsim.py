import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import *

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

    if False:
        eliminar_col = list()
        corregir = list()
        [nom_bus, nom_trf, nom_lin, nom_trf3, nom_cap, nom_ind] = lista_nom


        #PROBABLEMENTE!!! NO ES NECESARIO SI SE MODIFICA LAS VARIABLES DE SENSIBILIDAD EN EL POWERFACTORY

        # Buscar en los nombres obtenidos del excel y comparar con los obtenidos directamente desde PF.
        # Se guardan en 'eliminar_col', el numero de la columna que se desea eliminar
        # Esto se hace para eliminar los Factores de sensibilidad que no se necesitan
        # Para las líneas solo se utiliza dP/dP:bus1
        # Para los trf se utilizan dP/dP:bushv  y  dP/dP:buslv
        # Para los trf de 3dev se utilizan dP/dP:bushv   dP/dP:busmv   y dP/dP:buslv
        # Para los inductores solo se utiliza dP/dP:bus1
        # Para los capacitores solo se utiliza dP/dP:bus1
        for i in ejes_corregidos:  
            if i in nom_lin:
                continue
            elif i in nom_trf:
                corregir.append(ejes_corregidos.index(i)+1)
                continue
            elif i in nom_trf3:
                corregir.append(ejes_corregidos.index(i)+1)
                corregir.append(ejes_corregidos.index(i)+2)
                continue
            elif i in nom_ind:
                continue
            elif i in nom_cap:
                continue
            else:
                eliminar_col.append(ejes_corregidos.index(i))

# 
        for i in corregir:
            print(i)
            eliminar_col.remove(i)


    ptdf_array_T = np.asarray(df)
    #ptdf_array_T = np.delete(ptdf_array_T ,[0,1,2], axis = 1)


    # PROBABLEMENTE!! Bastaría con este código
    #ptdf_array_T = np.asarray(ptdf_dataframe)
    #ptdf_array_T = np.delete(ptdf_array_T ,[0,1,2], axis = 1)
    
    return np.transpose(ptdf_array_T).astype(float), df.axes[1]


class PowerFactorySim(object):
    def __init__(self, project_name='Project'):
        # start PowerFactory
        self.app = pf.GetApplication()
        self.app.Show()
        #activate project
        self.app.ActivateProject(project_name)
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
        self.IntEvt = self.app.GetFromStudyCase('IntEvt')
        self.res = self.app.GetFromStudyCase('*.ElmRes')

        self.Sb = 100
        self.raiz3 = 1.73205080757
        self.kV_fm = 110            # Tensiones menores, tendrán Fm = 1e9

    # Return pandas con matriz de Shift-factors
    # Return lista ordenada con barras 
    def export_csv(self, flujo):
        sfactors = self.app.GetFromStudyCase('ComVstab')
        if flujo == 'AC':
            sfactors.iopt_method = 0 # 0 = AC; 2 = DC
        if flujo == 'DC':
            sfactors.iopt_method = 2 # 0 = AC; 2 = DC

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

        return (pd.read_csv(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\SF_datos.csv',skiprows = 0,delimiter=';'), sfactors.p_bus.All())
        
    # Return diccionarios de los elementos 
    # bus, carga, linea, trafo, gen, genstat
    def get_data(self, buses, state='post'):
        dict_barras = dict()
        cont = 0
        for bus in buses:
            dict_barras[bus.loc_name] = cont
            cont += 1

        # CARGAS - en servicio y switch cerrado
        # name carga = (N° barra, Potencia MW, NombreBarra)
        dict_cargas = dict()                              
        for c in self.cargas:
            if c.outserv == 0 and c.bus1.cpCB.on_off == 1:                  # carga en servicio y switch conectado
                # name carga = (N° barra, Potencia MW, NombreBarra)
                if state == 'post':
                    potencia = c.GetAttribute('m:Psum:bus1')
                else:
                    potencia = c.plini
                dict_cargas[c.loc_name] = (dict_barras[c.bus1.cterm.loc_name],potencia ,c.bus1.cterm.loc_name)


        # LINEAS - en servicio
        # # name linea = (R, X, fmax, N° paralelas)
        dict_lineas = dict()
        for line in self.lineas:
            if line.outserv == 0:  #Linea en servicio
                Zb = round(line.bus1.cterm.uknom*line.bus2.cterm.uknom/self.Sb,6)     
                # Escoger flujo máximo de la línea cuando es menor a kV_fm
                if line.bus1.cterm.uknom < self.kV_fm:
                    fmax_kV = 1e9
                else:
                    fmax_kV = round(line.Inom_a*line.bus1.cterm.uknom*self.raiz3,4)
                # name linea = (R, X, fmax, N° paralelas)
                dict_lineas[line.loc_name]=(round(line.R1/Zb,6), round(line.X1/Zb,6), fmax_kV, line.nlnum)

        # TRAFOS 2 DEV - en servicio
        # name trafo = (Xd, R, N° paralelo, fmax)
        dict_trafos = dict()
        for t in self.trafos:
            if t.typ_id.r1pu == 0:      # Checkear R1 trafo 2 dev
                r = 1e-8
            else:
                r = t.typ_id.r1pu
            
            if t.outserv == 0:      # Trafo en servicio
                if t.bushv.cterm.uknom < self.kV_fm and t.buslv.cterm.uknom:
                    fmax_kV = 1e9
                else:
                    fmax_kV = t.Snom_a
                
                Zb = t.bushv.cterm.uknom * t.bushv.cterm.uknom/self.Sb
                Zbt = t.bushv.cterm.uknom * t.bushv.cterm.uknom/ (t.Snom_a/t.ntnum)

                # name trafo = (Xd, R, N° paralelo, fmax)
                dict_trafos[t.loc_name] = (3*t.typ_id.x1pu*(Zbt/Zb), 3*r*(Zbt/Zb), t.ntnum, fmax_kV)

        # GENERADORES
        #name generador = (N° Barra, N° Gen paralelo,
        #                   Outserv, Pmin, Pmax, Pref,
        #                   costos var, costo fijo)
        dict_gen = dict()
        for gen in self.generadores:
            if gen.ip_ctrl == 1:    # Buscar barra slack
                bus_slack = gen.bus1.cterm.loc_name
            if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                if state == 'post':
                    potencia = gen.GetAttribute('m:Psum:bus1')
                else:
                    potencia = gen.pgini
                # name generador = (N° Barra, N° Gen paralelo,
                #                   Outserv, Pmin, Pmax, Pref,
                #                   costos var, costo fijo)
                dict_gen[gen.loc_name] = (dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          0, gen.Pmin_uc, gen.Pmax_uc, potencia,
                                          gen.penaltyCosts, gen.fixedCosts)
            else:
                dict_gen[gen.loc_name] = (dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          1, gen.Pmin_uc, gen.Pmax_uc, 0,
                                          gen.penaltyCosts, gen.fixedCosts)

        # GENERADORES ESTATICOS
        # name generador = (N° Barra, N° Gen paralelo,
        #                   Outserv, Pmin, Pmax, Pref,
        #                   costos var, costo fijo)
        dict_genstat = dict()
        for gen in self.genstate:
            if gen.ip_ctrl == 1:    # Buscar barra slack
                bus_slack = gen.bus1.cterm.loc_name
            if gen.outserv == 0 and  gen.bus1.cpCB.on_off == 1:
                # name generador = (N° Barra, N° Gen paralelo,
                #                   Outserv, Pmin, Pmax, Pref,
                #                   costos var, costo fijo)
                dict_genstat[gen.loc_name] = (dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          0, gen.Pmin_uc, gen.Pmax_uc, gen.pgini,
                                          gen.penaltyCosts, gen.fixedCosts)
            else:
                dict_genstat[gen.loc_name] = (dict_barras[gen.bus1.cterm.loc_name], gen.ngnum,
                                          1, gen.Pmin_uc, gen.Pmax_uc, 0,
                                          gen.penaltyCosts, gen.fixedCosts)

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

        #return [nom_bus, nom_trf, nom_lin, nom_trf3, nom_cap, nom_ind]
        return (dict_barras, dict_cargas, dict_lineas, dict_trafos, dict_gen, dict_genstat)
        ## Lista de cargas
        #cargas = self.app.GetCalcRelevantObjects('*.ElmLod')
        #nc = len(cargas) # N° de cargas
        #nl = len(lineas) # N° de lineas
        #nb = len(barras) # N° de barras
        #nt = len(trafos) # N° de transformadores
        #nt3 = len(trafos3) # N° de transformadores 3 devanados
        #ncap = len(capacitores) #N° de capacitores series
        #nind = len(inductores) #N° de inductores series
        # # Elementos activos
        #ng = len(generadores)
        #ngs = len(genstate)
        #nga = len(asincronicos)

    def prepare_dynamic_sim(self, monitored_variables,
                            sim_type='rms', start_time=0.0,
                            step_size=0.01, end_time=10.0):
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
        self.sim.tstop = end_time
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
                print('%s ; Gurobi (DC) = %.3f [MW]; DIgSILENT (AC) = %.3f [MW]' % (line.loc_name, float(f_line[1][np.where(f_line[0] == line.loc_name)])*self.Sb, (-line.GetAttribute('m:Psum:bus2')+line.GetAttribute('m:Psum:bus1'))/2))
                cont += 1


#class LpGurobi(object):
#    def __init__(self):
#        self.m = gp.Model('Modelo 1')
#        self.m.Params.MIPGap = 1e-5
#        self.m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

