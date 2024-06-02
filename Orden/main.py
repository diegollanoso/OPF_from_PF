import pfsim
import pandas as pd
import numpy as np

sim = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#sim = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')




ptdf_dataframe, indices_bus = sim.export_csv()
#ptdf_dataframe = pd.read_csv(r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Tesis\Code\SF_datos.csv',skiprows = 0,delimiter=';')

shift_factors, indices_obj = pfsim.ShiftFactors(ptdf_dataframe)

(dict_barras, dict_cargas, dict_lineas, dict_trafos, dict_gen, dict_genstat) = sim.get_data(indices_bus)

Cv = np.zeros(len(dict_gen))



#indices = pd.read_csv(r'C:\Users\lldie\Desktop\SF_Indices.csv',skiprows = 0,delimiter=';', usecols=[2], encoding= 'ISO-8859-1')

#array_indices = indices.to_numpy()[:,0]
#array_indices_corregidos = list()
#for i in array_indices:
#    if '\\' in i:
#        array_indices_corregidos.append(i.split('\\')[1])
#    else:
#        array_indices_corregidos.append(i)
#