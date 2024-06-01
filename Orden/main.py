import pfsim
import pandas as pd

sim = pfsim.PowerFactorySim('IEEE 14-bus con PE y PV')
#sim = pfsim.PowerFactorySim('Taller_AGClisto2207-BD-OP-COORD-DMAP')

ptdf_dataframe = sim.export_csv()

shift_factors, indices = pfsim.ShiftFactors(ptdf_dataframe)


lista_nom = sim.get_data()

#indices = pd.read_csv(r'C:\Users\lldie\Desktop\SF_Indices.csv',skiprows = 0,delimiter=';', usecols=[2], encoding= 'ISO-8859-1')

#array_indices = indices.to_numpy()[:,0]
#array_indices_corregidos = list()
#for i in array_indices:
#    if '\\' in i:
#        array_indices_corregidos.append(i.split('\\')[1])
#    else:
#        array_indices_corregidos.append(i)
#