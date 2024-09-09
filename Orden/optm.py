import gurobipy as gp
import numpy as np

class Modelo():
    def __init__(self):
        self.gp = gp
        self.m = gp.Model('Modelo AGC')
        #m.setParam('DualReductions', 0)
        self.m.Params.MIPGap = 1e-7
        self.m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi

        self.T_Sfc = 15
        self.pot_down = 0
        self.flujos = 0
        self.losses = 0
        self.losses_plus = 0
        self.L = 11
        self.Voll = 500

    def __call__(self,data,simm):
        self.pg_inc = self.m.addMVar((data.n_gen_agc, data.Ns, data.Nt), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_inc')
        self.vg_inc = self.m.addMVar(data.n_gen_agc, vtype=self.gp.GRB.BINARY, name='vg_inc')

        self.p_ens = self.m.addMVar((data.Nb, data.Ns, data.Nt), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='P_ens')

        if self.pot_down:
            self.vg_dec = self.m.addMVar(data.n_gen_agc, vtype=self.gp.GRB.BINARY, name='vg_dec')
            self.pg_dec = self.m.addMVar((data.n_gen_agc, data.Ns, data.Nt), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_dec')
        
        self.f = self.m.addMVar((data.n_elem, data.Ns, data.Nt),vtype=self.gp.GRB.CONTINUOUS, lb=-self.gp.GRB.INFINITY, ub=self.gp.GRB.INFINITY, name='flujo') # Flujo por las líneas
        
        if self.losses:
            fp = self.m.addMVar((data.n_elem, data.Ns, data.Nt),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fp') # Flujo-p por las líneas
            fn = self.m.addMVar((data.n_elem, data.Ns, data.Nt),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fn') # Flujo-n por las líneas

            self.ploss = self.m.addMVar((data.n_elem, data.Ns, data.Nt),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='losses') # Flujo por las líneas
            dpk = self.m.addMVar((data.n_elem, self.L, data.Ns, data.Nt), vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo

            n_l = self.m.addMVar((data.n_elem, data.Ns, data.Nt), vtype=self.gp.GRB.BINARY, name='n_l')                                          # variable binaria complentaridad
            n_a = self.m.addMVar((data.n_elem, data.Ns, data.Nt, self.L), vtype=self.gp.GRB.BINARY, name='n_a')                                          # variable binaria adyacencia


        self.ObjFunction(data)

        for ti  in range(data.Nt):
            for s in range(data.Ns):
                self.Balance(data, simm, s, ti)
                self.R_Generadores(data, simm, ti, s)

                if self.losses:
                    self.R_FlujoLosses(data, simm, ti, s, fp, fn, dpk)
                    if self.losses_plus:
                        self.R_FlujoLossesPlus(data, fp, fn, n_l, n_a , dpk, ti, s)

                elif self.flujos:
                    self.R_FlujoLosses(data, simm, ti, s, fp=0,fn=0,dpk=0)
                


        self.m.write('OPF.lp')


    def ObjFunction(self, data):
        if self.pot_down:
            self.Csfc = data.Ccte_gen[data.pos_gen_agc_list] @ (self.vg_inc + self.vg_dec)            

        else:
            self.Csfc = data.Ccte_gen[data.pos_gen_agc_list] @ (self.vg_inc)

        f_obj = 0

        #Csfc += Ccte_genstat[pos_genstat_agc_list] @ (vstatg_inc + vstatg_dec)
        
        self.Cop = 0
        self.Cue = 0
        self.Cpl = 0
        self.Cts = 0
        for ti in range(data.Nt):
            Cop_s = 0
            Cue_s = 0
            Cpl_s = 0
            Cts_s = 0
            for s in range(data.Ns):
                if self.pot_down:
                    inc = self.pg_inc[:,s,ti] + self.pg_dec[:,s,ti]
                else:
                    inc = self.pg_inc[:,s,ti]

                #Cop_s += ss.Gen_on_AGC[ti,:,s][ss.pos_gen_agc_list] * ss.Cvar_gen[ss.pos_gen_agc_list]*Sb @ (inc)
                Cop_s += data.Cvar_gen[data.pos_gen_agc_list]*data.Sb @ (inc)
                #if np.size(Cvar_genstat) != 0:
                #    Cop_s += Genstat_on_AGC[:,s] * Cvar_genstat*sim.Sb @ pstatg_inc[:,s]

                Cue_s += (self.Voll * self.p_ens[:,s,ti]).sum() * data.Sb
                if self.losses and not self.losses_plus:
                    Cpl_s += (self.Voll * self.ploss[:,s,ti]).sum() * data.Sb
                #Cts_s += Costo_ts*(1-s_ts[:,s,ti]).sum()

            self.Cop += Cop_s
            self.Cue += Cue_s
            self.Cpl += Cpl_s
            self.Cts += Cts_s

        f_obj = self.Csfc + self.Cop + self.Cue + self.Cpl + self.Cts

        self.m.setObjective(f_obj,self.gp.GRB.MINIMIZE)

    def Balance(self, data, simm, s, ti):    
        # Con perdidas y gen pueden bajar potencia
        if self.flujos and self.losses and self.pot_down:
            self.m.addConstr(self.pg_inc[:,s,ti].sum() - self.pg_dec[:,s,ti].sum() + self.p_ens[:,s].sum() == simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum() - simm.PL_pre_line[ti] + self.ploss[:,s,ti].sum(), name ='Balance')
        
        # Con perdidas 
        elif self.flujos and self.losses:
            self.m.addConstr(self.pg_inc[:,s,ti].sum() + self.p_ens[:,s,ti].sum() == simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum()- simm.PL_pre_line[ti] + self.ploss[:,s,ti].sum(), name ='Balance')
        
        # Sin Perdidas y con sistema de transmisión
        elif self.flujos and not data.flujo_dc:
            self.m.addConstr(self.pg_inc[:,s,ti].sum() + self.p_ens[:,s,ti].sum() == simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum(), name ='Balance')

        # Sin Perdidas / Inc potencia
        else:
            self.m.addConstr(self.pg_inc[:,s,ti].sum() + self.p_ens[:,s,ti].sum() == simm.P_out[s,ti], name ='Balance')


    def R_FlujoLosses(self, data, simm, ti, s, fp, fn, dpk):

        if self.pot_down:
            inc = self.pg_inc[:,s,ti] - self.pg_dec[:,s,ti]

        else:
            inc = self.pg_inc[:,s,ti]


        f_gen =  data.SF[:,data.pos_gen] @ data.Pgen_pre[:,ti]
        f_genstat = 0
        #f_genstat =  SF[:,simm.pos_genstat] @ simm.Pgenstat_pre
        f_ens =  data.SF @ self.p_ens[:,s,ti]
        f_gen_agc = data.SF[:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        #f_gen_agc = SF[:,simm.pos_g2] @ (inc)
        f_genstat_agc = 0
        #if np.size(Cvar_genstat) != 0:
        #    f_genstat_agc = Genstat_on_AGC[:,s] * SF[:,pos_genstat] @ (pstatg_inc[:,s] - pstatg_dec[:,s])
        f_gen_out = data.SF[:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        #Flujo_dda = data.SF @ simm.D_pfc[:,s,ti]
        Flujo_dda = data.SF @ simm.D_pfc[:,s,ti]
        if self.losses:
            f_loss = 0.5 * data.SF @ abs(data.A.T) @ self.ploss[:,s,ti]
        else:
            f_loss = 0
        
        self.m.addConstr(self.f[:,s,ti] == f_gen + f_genstat + f_ens + f_gen_agc + f_genstat_agc - Flujo_dda - f_gen_out - f_loss, name='flujos') 

        if self.losses:
            self.m.addConstr(self.f[:,s,ti] == fp[:,s,ti] - fn[:,s,ti], name = 'f')
            self.m.addConstr(fp[:,s,ti] + fn[:,s,ti] == dpk[:,:,s,ti].sum(1), name='SumaDpk')
        
            kl = np.zeros((data.n_elem, self.L))
            for l in range(self.L):
                kl[:,l] = (2*(l+1)-1)*(data.FMax)/self.L
            self.m.addConstr(self.ploss[:,s,ti] == data.G/(data.B**2)*(self.gp.quicksum(kl[:,i]*dpk[:,i,s,ti] for i in range(self.L))), name = 'Ploss')  
        
            self.m.addConstr(-fp[:,s,ti] - fn[:,s,ti] - 0.5*self.ploss[:,s,ti] >= -data.FMax, name = 'fp s='+str(s)+' c='+str(ti))
            self.m.addConstr(-fp[:,s,ti] >= -data.FMax, name = 'fp+  s='+str(s)+' c='+str(ti))
            self.m.addConstr(-fn[:,s,ti] >= -data.FMax, name = 'fn+  s='+str(s)+' c='+str(ti))
            for l in range(self.L):
                self.m.addConstr(-dpk[:,l,s,ti] >= -data.FMax/self.L, name = 'LimiteDpk')

        else:
            self.m.addConstr(-self.f[:,s,ti] >= -data.FMax, name = 'fMax+ s='+str(s)+' c='+str(ti))
            self.m.addConstr(self.f[:,s,ti] >= -data.FMax, name = 'fMax- s='+str(s)+' c='+str(ti))


    def R_FlujoLossesPlus(self, data, fp, fn, n_l, n_a , dpk, ti, s):

        self.m.addConstr(-fp[:,s,ti] >= -n_l[:,s,ti]*data.FMax, name='fp_n  s='+str(s)+' c='+str(ti))   #flujo positvo-restriccion de complementaridad
        self.m.addConstr(-fn[:,s,ti] >= (1-n_l[:,s,ti])*(-data.FMax), name='fn_n  s='+str(s)+' c='+str(ti)) #flujo nefativo-restriccion de complementaridad
        
        for l in range(self.L): 
            if l==0:
                self.m.addConstr(-dpk[:,l,s,ti]>= -data.FMax/self.L, name='d_f_Res_max_A_l')
                self.m.addConstr(dpk[:,l,s,ti] >= n_a[:,s,ti,l]*(data.FMax/self.L), name='d_f_Res_min_A_l')
            elif l==self.L-1:
                self.m.addConstr(-dpk[:,l,s,ti] >=-n_a[:,s,ti,l-1]*data.FMax/self.L, name='d_f_Res_max_A_L')
                self.m.addConstr(dpk[:,l,s,ti] >=0, name='d_f_Res_min_A_L')
            else:
                self.m.addConstr(-dpk[:,l,s,ti] >=-n_a[:,s,ti,l-1]*data.FMax/self.L, name='d_f_Res_max_A_L-1')
                self.m.addConstr(dpk[:,l,s,ti] >=n_a[:,s,ti,l]*data.FMax/self.L, name='d_f_Res_min_A_L-1')



    def R_Generadores(self, data, simm, ti, s):

        self.m.addConstr(-self.pg_inc[:,s,ti] >= -self.vg_inc * (data.Pmax_gen[data.pos_gen_agc_list] * data.ngen_par[data.pos_gen_agc_list] - data.Pgen_pre[data.pos_gen_agc_list,ti]), name = 'PMax  s='+str(s)+' c='+str(ti))
        if self.pot_down:
            self.m.addConstr(-self.pg_dec[:,s,ti] >= -self.vg_dec * (data.Pgen_pre[data.pos_gen_agc_list,ti] - data.Pmin_gen[data.pos_gen_agc_list] * data.ngen_par[data.pos_gen_agc_list]), name = 'PMin  s='+str(s)+' c='+str(ti))
        #Nueva Lista que solo tiene unidades participantes en el AGC / se quita de unidades participantes la unidad que sale de servicio
        pos_part_gen_agc_list = data.pos_gen_agc_list[:]
        x = list(range(data.n_gen_agc))
        if data.Gen_Outages[s] in data.name_gen_agc_list:
            pos_part_gen_agc_list.remove(data.all_gen.index(data.Gen_Outages[s]))
            x.pop(data.name_gen_agc_list.index(data.Gen_Outages[s]))
            self.m.addConstr(self.pg_inc[data.name_gen_agc_list.index(data.Gen_Outages[s]), s, ti] == 0, name = 'GenOut+  s='+str(s)+' c='+str(ti))
            if self.pot_down:
                self.m.addConstr(self.pg_dec[data.name_gen_agc_list.index(data.Gen_Outages[s]), s, ti] == 0, name = 'GenOut-  s='+str(s)+' c='+str(ti))
        if data.flujo_dc:
            self.m.addConstr(-self.pg_inc[x, s, ti] >= -self.T_Sfc * data.Ramp_gen[pos_part_gen_agc_list] - (data.Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + data.Gen_Outages[s] +'_Ramp+  s='+str(s)+' c='+str(ti))
        else:
            self.m.addConstr(-self.pg_inc[x, s, ti] >= -self.T_Sfc * data.Ramp_gen[pos_part_gen_agc_list] - (simm.Pgen_pfc[pos_part_gen_agc_list,s,ti] - data.Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + data.Gen_Outages[s] +'_Ramp+  s='+str(s)+' c='+str(ti))
        if self.pot_down:
            self.m.addConstr(-self.pg_dec[x, s, ti] >= -self.T_Sfc * data.Ramp_gen[pos_part_gen_agc_list] + (simm.Pgen_pfc[pos_part_gen_agc_list,s,ti] - data.Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + data.Gen_Outages[s] +'_Ramp-  s='+str(s)+' c='+str(ti))
        #if np.size(Cvar_genstat) != 0:
            #self.m.addConstr(-pg_inc[pos_genstat_agc_list, s, ti] >= -vstatg_inc * (Pmax_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list] - Pgen_pre[pos_genstat_agc_list,ti]))
            #self.m.addConstr(-pg_dec[pos_genstat_agc_list, s, ti] >= -vstatg_dec * (Pgen_pre[pos_genstat_agc_list,ti] - Pmin_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list]))
            #self.m.addConstr(-pg_inc[pos_genstat_agc_list,s] >= -T_Sfc * Ramp_genstat[pos_genstat_agc_list] - (Pgenstat_pfc[pos_genstat_agc_list, s] - Pgenstat_pre[pos_genstat_agc_list]))
            #self.m.addConstr(-pg_dec[pos_genstat_agc_list,s] >= -T_Sfc * Ramp_genstat[pos_genstat_agc_list] + (Pgenstat_pfc[pos_genstat_agc_list, s] - Pgenstat_pre[pos_genstat_agc_list]))


        self.m.addConstr(-self.p_ens[:,s,ti] >= -simm.D_pfc[:,s,ti], 'LimENS  s='+str(s)+' c='+str(ti))

    def run(self):
        self.m.optimize()

    def Results(self, data, simm):
        status = self.m.Status
        if status == gp.GRB.Status.OPTIMAL:
            print('-----------------------------')
            #print('La demanda total del sistema es: %.2f (MW)' % (pf.dda_barra[:,ti].sum()*pf.Sb))
            
            if self.losses and not self.losses_plus:
                print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($) + Cpl = %.2f ($)' % (self.m.objVal,self.Csfc.getValue(),self.Cop.getValue(),self.Cue.getValue(),self.Cpl.getValue()))
            
            else:
                print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($)' % (self.m.objVal,self.Csfc.getValue(),self.Cop.getValue(),self.Cue.getValue()))


            print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (self.m.NumVars,self.m.NumConstrs,self.m.DNumNZs)) #print('num_Vars =  %d / num_Const =  %d' % (len(m.getVars()), len(m.getConstrs())))      
            #print ('Total P_loss = %.2f [MW]'%(pf.Sb*ploss.sum().getValue()))
            print('=> Solver time: %.4f (s)' % (self.m.Runtime))
            print ('Costo => %.2f ($/h)' % self.m.objVal) 
            #print ('Las perdidas son %.2f (MW)' % sum(pk_loss[i].X for i in range(len(pk_loss)))) 
        
            Post_gen_out = np.zeros((data.Ns,data.Nt))
    
            self.part_factors = np.zeros((data.n_gen_agc,data.Ns,data.Nt))
            for gen, s, ti in np.ndindex(self.pg_inc.x.shape):
                
                if self.losses:
                    Post_gen_out[s, ti] = simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum() - simm.PL_pre_line[ti] + self.ploss[:,s,ti].x.sum()

                elif self.flujos:
                    Post_gen_out[s, ti] = simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum()
                
                else:
                    Post_gen_out[s, ti] = simm.P_out[s,ti]

                
                if self.pot_down:
                    if self.pg_inc.x[gen,s,ti] != 0 or self.pg_dec.x[gen,s,ti] != 0:
                        self.part_factors[gen,s,ti] = (self.pg_inc.x[gen,s,ti]-self.pg_dec.x[gen,s,ti])/Post_gen_out[s,ti]
                
                else:
                    if self.pg_inc.x[gen,s,ti] != 0:
                        self.part_factors[gen,s,ti] = (self.pg_inc.x[gen,s,ti])/Post_gen_out[s,ti]
                





    
        elif status == gp.GRB.Status.INF_OR_UNBD or \
            status == gp.GRB.Status.INFEASIBLE  or \
            status == gp.GRB.Status.UNBOUNDED:
            print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
            self.m.computeIIS() 
            self.m.write("GTCEP.ilp")

