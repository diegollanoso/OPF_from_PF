import sys
import pfsimV2 as pfsim
import pandas as pd
import numpy as np
import time
import gurobipy as gp
from gurobipy import *


t0  = time.time()
project_file = 'IEEE 14-bus AGC(2)'
pf = pfsim.PowerFactorySim(project_file, False)
pf.Nt = 1

pf.export_csv('AC')

dict_lineas = pf.get_data()

Monitored_vars = {'*.ElmSym': ['m:Psum:bus1'],
                  '*.ElmLod': ['m:Psum:bus1'],
                  '*.ElmLne': ['c:Losses'],
                  '*.ElmTr2':['m:P:bushv','m:P:buslv']}

pf.prepare_dynamic_sim(Monitored_vars, 'rms')
pf.run_dynamic_sim(end_sim=150)

pf_time, freq_values = pf.extract_data('Bus 2.ElmTerm', 'm:fehz', return_time=True)

# Extrer datos de cada generador
cont = 0
return_time = False
list_gens_values = np.zeros((pf.ngen + pf.ngenstat + pf.ngenpv, len(pf_time)))
for gen in pf.all_gen:
    var_values = pf.extract_data(gen + '.ElmSym', 'm:Psum:bus1', return_time)
    list_gens_values[cont] = var_values
    cont += 1

for gen in pf.all_genstat:
    var_values = pf.extract_data(gen + '.ElmGenstat', 'm:Psum:bus1', return_time)
    list_gens_values[cont] = var_values
    cont += 1

for gen in pf.all_pv:
    var_values = pf.extract_data(gen + '.ElmPvsys', 'm:Psum:bus1', return_time)
    list_gens_values[cont] = var_values
    cont += 1

# Extraer datos de todos las demandas
cargas_name = list()
list_load_values = np.zeros((len(pf.cargas), len(pf_time)))
cont = 0
for carga in pf.cargas:
    var_values = pf.extract_data(carga.loc_name + '.ElmLod', 'm:Psum:bus1', return_time)
    cargas_name.append(carga.loc_name)
    list_load_values[cont] = var_values
    cont += 1

# Extraer datos de todas las líneas
lineas_name = list()
list_line_values = np.zeros((len(pf.lineas), len(pf_time)))
cont = 0
for linea in pf.lineas:
    var_values = pf.extract_data(linea.loc_name + '.ElmLne', 'c:Losses', return_time)
    lineas_name.append(linea.loc_name)
    list_line_values[cont] = np.array(var_values) / 1000
    cont += 1

trafos_name = list()
list_trafos_values = np.zeros((len(pf.trafos), len(pf_time)))
cont = 0
for trafo in pf.trafos:
    trafos_name.append(trafo.loc_name)
    var_hv = pf.extract_data(trafo.loc_name + '.ElmTr2', 'm:P:bushv', return_time)
    var_lv = pf.extract_data(trafo.loc_name + '.ElmTr2', 'm:P:buslv', return_time)
    list_trafos_values[cont] = abs(np.array(var_hv) + np.array(var_lv))
    cont += 1

sum_gens = np.zeros(len(pf_time))
for tiempo in range(len(pf_time)):
    sum_gens[tiempo] = list_gens_values[:, tiempo].sum()

sum_loads = np.zeros(len(pf_time))
for tiempo in range(len(pf_time)):
    sum_loads[tiempo] = list_line_values[:, tiempo].sum()
    sum_loads[tiempo] += list_trafos_values[:, tiempo].sum()
    sum_loads[tiempo] += list_load_values[:, tiempo].sum()

data = pd.DataFrame(
    np.vstack(
        (
            pf_time,
            freq_values,
            list_gens_values,
            list_load_values,
            list_line_values,
            list_trafos_values,
            sum_gens,
            sum_loads,
        )
    ).T,
    columns=['Time','Freq']
    + pf.all_gen
    + pf.all_genstat
    + pf.all_pv
    + cargas_name
    + lineas_name
    + trafos_name
    + ['Sum Gens', 'Sum Total Loads'],
    )

data.to_excel('data_IEEE14Bus_emt.xlsx', index=False)

print('Finish')