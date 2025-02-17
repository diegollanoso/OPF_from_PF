import matplotlib.pyplot as plt
import control
import numpy as np
import os

#System Dynamic
H = 10
D = 0.1

ps = control.TransferFunction([1], [H, D])

#Turbine and generator 
K = 1.3
F = 0.01
R = 0.05
T = 1.7

Tg = control.TransferFunction([K*F*T, K], [T*R, R])

tf = control.feedback(ps,Tg, sign=-1)

#print(tf)



t0 = 0
t1 = 60
dt = 0.001
nt = int((t1-t0)/dt) + 1
t = np.linspace(t0, t1, nt)
u = 10*np.ones(nt)

(t, y) = control.forced_response(tf, t, -u)


plt.plot(t,y+50)
output_dir = r'c:/Users/lldie/OneDrive - Universidad Técnica Federico Santa María/Universidad/Memoria/Code/Graphs'

output_path = output_dir + '/response.svg'
plt.grid()
#plt.show()
plt.savefig(output_path)
#os.startfile(output_path)