import pandas as pd
import numpy as np

# Importar DIgSILENT
import sys
sys.path.append(r'C:\Program Files\DIgSILENT\PowerFactory 2023 SP3A\Python\3.11')
import powerfactory as pf


# La matriz de shift factors se limpia y ordena
# Return matriz de shif factors
def ShiftFactors(ptdf_dataframe, lista_nom):
    ejes = ptdf_dataframe.axes
    ejes_corregidos = list()

    #Líneas se encuentran dentro de otros objetos
    # ej: Talca - Itahue 66 kV L1\Tap San Rafael - Panguilemo  66 kV_bhk1
    # La línea se encuentra dentro de Talca - Itahue 66 kV L1
    # Ejes corregidos es el nombre del elemento sin incluir la carpeta de origen

    for i in ejes[1]:
        if '\\' in i:
            ejes_corregidos.append(i.split('\\')[1])
        else:
            ejes_corregidos.append(i)

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


    ptdf_array_T = np.asarray(ptdf_dataframe)
    ptdf_array_T = np.delete(ptdf_array_T ,eliminar_col, axis = 1)

    # PROBABLEMENTE!! Bastaría con este código
    #ptdf_array_T = np.asarray(ptdf_dataframe)
    #ptdf_array_T = np.delete(ptdf_array_T ,[0,1,2], axis = 1)
    
    return np.transpose(ptdf_array_T).astype(float)


class PowerFactorySim(object):
    def __init__(self, project_name='Project'):
        # start PowerFactory
        self.app = pf.GetApplication()
        #self.app.Show()
        #activate project
        self.app.ActivateProject(project_name)

    # Return pandas con matriz de Shift-factors
    # Return lista ordenada con barras 
    def export_csv(self):
        sfactors = self.app.GetFromStudyCase('ComVstab')
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
        export.f_name = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.csv'
        export.iopt_csel = 0 
        export.locn = 1
        export.iopt_sep=0
        export.iopt_locn = 1 
        export.ciopt_head = 0
        export.dec_Sep = '.'
        export.col_Sep = ';'
        export.Execute()

        return (pd.read_csv(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.csv',skiprows = 0,delimiter=';'), sfactors.pbus.All())
        
    # Return 6 Listas con los nombres de cada uno de los elementos 
    # bus, trafo, linea, trfo 3dev, cap, ind
    def get_data(self):
        generadores = self.app.GetCalcRelevantObjects('*.ElmSym') 
        cargas = self.app.GetCalcRelevantObjects('*.ElmLod')
        lineas = self.app.GetCalcRelevantObjects('*.ElmLne')
        barras = self.app.GetCalcRelevantObjects('*.ElmTerm')
        trafos = self.app.GetCalcRelevantObjects('*.ElmTr2')
        trafos3 = self.app.GetCalcRelevantObjects('*.ElmTr3')
        capacitores = self.app.GetCalcRelevantObjects('*.ElmScap')
        inductores = self.app.GetCalcRelevantObjects('*.ElmSind')
        genstate = self.app.GetCalcRelevantObjects('*.ElmGenstat')
        asincronicos = self.app.GetCalcRelevantObjects('*.ElmAsm')
        shunt = self.app.GetCalcRelevantObjects('*.ElmShnt')
        
        #### Listas de elementos
        nom_bus = list()
        nom_trf = list()
        nom_lin = list()
        lineas_error = list()

        nom_trf3 = list()
        nom_trf3_eq = ['Rapel equivalente 1_acr1','Rapel equivalente 2_acs0']
        nom_cap = list()
        nom_ind = list()
        nom_cargas = list()
        

        for b in barras:
            nom_bus.append(b.loc_name)

        for t in trafos:
            nom_trf.append(t.loc_name)

        for l in lineas:
            nom_lin.append(l.loc_name)

        for t in trafos3:
            nom_trf3.append(t.loc_name)

        for c in capacitores:
            nom_cap.append(c.loc_name)

        for i in inductores:
            nom_ind.append(i.loc_name)

        return [nom_bus, nom_trf, nom_lin, nom_trf3, nom_cap, nom_ind]
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