# Code

### CASOS

- A:Sin sistema de transmisión
    VOLL = 150000
    No hay penalización por perdidas

- B:Sistema de transmisión
    Restricciones de complementaridad y adyacencia
    VOLL = 150000
    No hay penalización por perdidas

### TBD

- Añadir contigencias de carga para variabilidad

- Probar en el otro archivo con gen renovable
- Añadir variabiliadd ERV y demanda
- Hacer ajuste de filtro de AGC
- AGC se activa solo para ciertas frecuencias

- Añadir gráfico que muestre la frecuencia y la variación de potencia en los Gens, en particular, los del gen del csf

- Se añadio restricción ENS = 0

- Al añadir contingencia de carga, revisar la lista optm.gen_csf, crear dos listas. Una para vg_inc y otra para vg_dec

- Fijar los factores de participación utilizando el modelo original en t=30s

- Para las potencias de los generadores, extraer solo los que participan en el AGC. 

- No convergencia dynamic model

- P_agc = P_out + P_variacones

- P_out debe disminuir

- Revisar el valor de P_out de optm.py
(Variaciones)

- Variabilidad en demanda (alta, media y baja variación), definir por áres o mezclar

- Datos de generación renovable

- Añadir control terciario de frecuencia utilizando UC

- Arreglar muestreo al AGC

- Revisar restriccion de suma de dpk = flujos

- Realizar simulaciones para confirmar valores de flujo DC

- Añadir factor de seguridad a FMax

- Revisar restricciones de complementaridad y exclusividad

- Falta añadir la probabilidad de ocurrencia de las contigencnias en la función objetivo

- Que tomar como perdidas? Solo líneas o tambien los trafos?? Se consideran ambos! No se pueden añadir perdidas en simulaciones 

- Adecuar restricciones de rampa cuando existe desconexión de carga.

- ¿las perdidas PL_PRE son obtenidas antes de realizar la contingencia?

- Idea: Solo extraer la información de los generadores que participan en el AGC.

- Eliminar carpetas y archivos innecesarios.

### Considerar

- Cuántos generadores estaticos? 50%
aumentar la penetración con 2 + 1 + 1 ERV en cada área

- Se tuvo que ajustar el parámetro IntFeasTol del modelo para que las restricciones Big-M 
https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/numerics_guide.html#dealing-with-big-m-constraints

- Se asume que todos las contigencias afectan al sistema, es decir, antes de la salida de cargas o generadores se encuentran consumiendo o entregando potencia.

- Añadir escenarios de demanda distintos al modelo (valle, promedio, peak)
    - Posibles soluciones: 
        1.- Usar listas para las cargas y generadores.
        2.- Exportar datos de potencia de cargas y generadores e importarlos como tuplas.

- Diferencias en las cargas, al final del CPF y se obtiene la demanda de las cargas durante la simulación, no se colocan condiciones. Mientras, que en el diccionario original se colocan 2 resitriccion, Out of Service / Circuit Breaker.

### questions

- Que hacer con la parte reactiva de las cargas. El txt solo tiene un dato, que considere como potencia activa

### Progreso

- Función eventos gamma generadores

- Revisar Flujos de perdidas PartialModel

- 9/9: Código de OPF diapositiva clases

- 8/9: Código de AGC realizado

- 3/9: Añadir flujos DC

- 2/9: Código más ordenado con clases para powerfactory, simulación y gurobi. 

- 24/8: Código más ordenado con clases para simulación y modelo gurobi

- 16/8: Se añaden la función objetivo y restricciones para los distintos escenarios de demanda.

- 13/8: Se solucionan los SF, problemas en el control de frecuencia - potencia del generador. Se añaden los escenarios valle y mean.

- 12/8: Se obtiene SF de forma algebraica, existen diferencias con el obtenido por PowerFactory

- 11/8: Se soluciono el problema de rampas
- 11/8: Se añadieron Excels con los resultados.

- 06/8: Añadir variación de demanda y perdidas por la variación de tensión. 
    - No olvidar la restricción de la potencia ENS <= DemandaPFC

- 05/8: Se añadio la variación de demanda por la diferencia de tensión en las cargas.

- 04/8: Se termino de añadir el efecto de disminución de potencia en unidades participantes del AGC de la formulación.

- 30/7: Se comenzo a añadir el efecto de disminución de potencia en unidades participantes del AGC de la formulación.

## Documentation

### Real.py
This script sets up and runs a power system simulation using PowerFactory and an optimization model. It includes functions for:
- Setting up the simulation parameters (`setup_simulation`)
- Changing line values (`change_max_line_values`)
- Exporting data (`export_simulation_data`)
- Running the optimization model (`run_optimization_model`)

The `main` function orchestrates the entire process, including initializing the simulation, running the optimization model, and handling the results.

### optm.py
This file contains the optimization model classes used in the simulation:
- `Modelo`: Main optimization model for the power system simulation.
- `PartialModel`: A partial optimization model used for specific scenarios.
- `Model_UC`: Unit commitment model.
- `Model_CTF`: Control Terciario de Frecuencia model.

Each class includes methods for setting up the model, defining the objective function, adding constraints, running the optimization, and extracting results.

### pfsim.py
This file contains the PowerFactory simulation classes and functions:
- `ShiftFactors`: Function to clean and order the shift factors matrix.
- `PowerFactorySim`: Class to handle the PowerFactory simulation, including methods for exporting data, getting simulation data, preparing and running dynamic simulations, and extracting results.
- `Simulacion`: Class to handle the simulation process, including monitoring variables and handling events.
- `ShortSim`: Class for short-term simulations.
- `new_SF`: Class to handle new shift factors with transmission switching.
- `CreateEvents_line` and `CreateEvents_gen`: Functions to create events for line and generator changes.
- `Set_param_agc`: Function to set AGC parameters.

The classes and functions in this file are used to manage the PowerFactory simulation, including setting up the simulation environment, running simulations, and extracting and processing results.