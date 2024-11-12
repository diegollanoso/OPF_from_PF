import gurobipy as gp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters for the synthetic ACE data
time = np.arange(0, 300, 1)  # Time from 0 to 300 seconds
noise_level = 0.2  # Level of noise to introduce
trend = 50  # General trend for the ACE data

# Create synthetic ACE data
np.random.seed(42)  # For reproducibility
noise = noise_level * np.random.randn(len(time))  # Random noise
#fluctuations = 2 * np.sin(0.05 * time)  # Low-frequency fluctuations
#ace_data = trend + fluctuations + noise  # Combine trend, fluctuations, and noise
ace_data = trend + noise  # Combine trend, fluctuations, and noise

# Create a DataFrame for easier manipulation and potential further analysis
#df = pd.DataFrame({"Time (s)": time, "ACE Crudo (MW)": ace_data})

# Plotting the synthetic ACE data
plt.figure(figsize=(10, 6))
plt.plot(time, ace_data, label="ACE crudo (sintético)", color="blue", alpha=0.7)
plt.xlabel("Tiempo (s)")
plt.ylabel("ACE (MW)")
plt.legend()
plt.grid()
plt.title("Datos sintéticos para ACE crudo")
plt.show()


# Parámetros iniciales
delt = 4  # Intervalo de tiempo entre ejecuciones del AGC
const = 0.1  # Valor de suavizado que deseas mantener constante
max_n = 50
max_T = 100

ACE_array = ace_data  # Array de ACE crudo
p = len(ACE_array)  # Número de puntos en el array de ACE crudo

# Create a new Gurobi model
model = gp.Model("optimization_model")

# Add variables to the model
G_avg = model.addVar(vtype=gp.GRB.CONTINUOUS,  lb=-gp.GRB.INFINITY, ub=gp.GRB.INFINITY, name="G_ave")
G_exp = model.addVar(vtype=gp.GRB.CONTINUOUS,  lb=-gp.GRB.INFINITY, ub=gp.GRB.INFINITY, name="G_exp")
n = model.addVar(vtype=gp.GRB.INTEGER, lb=1, ub=max_n, name="n")
T = model.addVar(vtype=gp.GRB.CONTINUOUS, lb=1, ub=max_T, name="T")

# Variables auxiliares
ACE_avg = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="ACE_avg")
ACE_exp = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="ACE_exp")
FACE = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="FACE")

# Variables auxiliares para los inversos
inv_n = model.addVar(vtype=gp.GRB.CONTINUOUS, name="inv_n", lb=1/max_n, ub=1)
inv_T = model.addVar(vtype=gp.GRB.CONTINUOUS, name="inv_T", lb=1/max_T, ub=1)

# Restricciones para los inversos
model.addConstr(inv_n * n == 1, "Inv_n_Constraint")
model.addConstr(inv_T * T == 1, "Inv_T_Constraint")

# Variables auxiliares
ACE_avg = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="ACE_avg")
ACE_exp = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="ACE_exp")
FACE = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="FACE")
abs_diff = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="abs_diff")
abs_smooth = model.addVars(p, vtype=gp.GRB.CONTINUOUS, name="abs_smooth")

# Restricciones para ACE_avg (Promedio Móvil)
for t in range(p):
    model.addConstr(
        ACE_avg[t] == gp.quicksum(ACE_array[max(0, t - i)] for i in range(max_n)) * inv_n,
        name=f"ACE_avg_{t}")

# Restricciones para ACE_exp (Filtro Exponencial)
model.addConstr(ACE_exp[0] == ACE_array[0], "Condición_inicial")
for t in range(1, p):
    model.addConstr(
        ACE_exp[t] == ACE_exp[t-1] + delt * inv_T * (ACE_array[t] - ACE_exp[t-1]),
        name=f"ACE_exp_{t}"
    )

# Restricciones para FACE
for t in range(p):
    model.addConstr(
        FACE[t] == G_avg * ACE_avg[t] + G_exp * ACE_exp[t],
        name=f"FACE_{t}"
    )

# Función de costo a minimizar usando variables auxiliares para el valor absoluto
for t in range(p):
    model.addConstr(abs_diff[t] >= ACE_array[t] - FACE[t], name=f"AbsDiff_Pos_{t}")
    model.addConstr(abs_diff[t] >= FACE[t] - ACE_array[t], name=f"AbsDiff_Neg_{t}")

# Función objetivo
cost = gp.quicksum(abs_diff[t] for t in range(p)) / p
model.setObjective(cost, gp.GRB.MINIMIZE)

# Restricción de suavizado usando variables auxiliares para el valor absoluto
for t in range(2, p):
    smooth_diff = FACE[t] - 2 * FACE[t-1] + FACE[t-2]
    model.addConstr(abs_smooth[t] >= smooth_diff, name=f"Smooth_Pos_{t}")
    model.addConstr(abs_smooth[t] >= -smooth_diff, name=f"Smooth_Neg_{t}")

smo = gp.quicksum(abs_smooth[t] for t in range(2, p)) / p
model.addConstr(smo == const, "Suavizado")

# Optimización
model.optimize()

# Resultados
if model.status == gp.GRB.OPTIMAL:
    print(f"Solución Óptima:")
    print(f"n = {n.X}")
    print(f"T = {T.X}")
    print(f"G_avg = {G_avg.X}")
    print(f"G_exp = {G_exp.X}")
    print(f"Costo = {model.objVal}")
else:
    model.computeIIS() 
    model.write("GTCEP_OPTM-agc.ilp")
    print("No se encontró solución óptima.")

print('Optimal objective: %g' % model.objVal)