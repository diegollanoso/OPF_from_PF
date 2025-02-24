import gurobipy as gp
import numpy as np
import logging


class GurobiFilter(logging.Filter):
    def __init__(self, name="GurobiFilter"):
        super().__init__(name)

    def filter(self, record):
        return False

grbfilter = GurobiFilter()

grblogger = logging.getLogger('gurobipy')
if grblogger is not None:
    grblogger.addFilter(grbfilter)
    grblogger = grblogger.getChild('gurobipy')
    if grblogger is not None:
        grblogger.addFilter(grbfilter)



class Modelo():
    def __init__(self):
        logging.info('Initializing Model')
        self.gp = gp
        self.m = gp.Model('Modelo AGC')
        #m.setParam('DualReductions', 0)
        self.m.Params.MIPGap = 1e-2
        #self.m.Params.PoolSolutions = 3
        #self.m.Params.OutputFlag = 0 # eliminar mensajes adicioneales Gurobi
        #self.m.Params.IntFeasTol = 1e-6

        self.T_Sfc = 15     # 15 minutos
        self.pot_down = 0
        self.flujos = 0
        self.losses = 0
        self.losses_plus = 0
        self.Voll = 500
        self.TS = False
        self.costo_ts = 0
        logging.info("Optimization model initialized")

    def __call__(self,data,simm):
        logging.info('Formulating Model')
        self.L = 11
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

        if self.losses_plus:
            n_l = self.m.addMVar((data.n_elem, data.Ns, data.Nt), vtype=self.gp.GRB.BINARY, name='n_l')                                          # variable binaria complentaridad
            n_a = self.m.addMVar((data.n_elem, data.Ns, data.Nt, self.L), vtype=self.gp.GRB.BINARY, name='n_a')                                          # variable binaria adyacencia

        if self.TS:
            self.s_ts = self.m.addMVar((len(data.pos_ts),data.Ns, data.Nt), vtype=self.gp.GRB.BINARY, name='s_ts')                                    # variable binaria de TS
            self.f_ts = self.m.addMVar((len(data.pos_ts),data.Ns, data.Nt),vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=-self.gp.GRB.INFINITY, name='f_ts')


        self.ObjFunction(data)
        logging.info('Formulating Objective Function')
        logging.info('Adding Constraints')
        for ti  in range(data.Nt):
            for s in range(data.Ns):
                self.Balance(data, simm, s, ti)
                self.R_Generadores(data, simm, ti, s)
                logging.info('Adding Balance and Generators Constraints for Scenario %d and Time %d' % (s,ti))
                if self.TS:
                    #self.Flujos_TS(data,simm, ti, s, f_ts, dpk, fp, fn)
                    self.Flujos_TS(data,simm, ti, s, dpk, fp, fn)
                    logging.info('Adding TS Constraints for Scenario %d and Time %d' % (s,ti))
                elif self.losses:
                    self.R_FlujoLosses(data, simm, ti, s, fp, fn, dpk)
                    logging.info('Adding Losses Constraints for Scenario %d and Time %d' % (s,ti))
                elif self.flujos:
                    self.R_FlujoLosses(data, simm, ti, s, fp=0,fn=0,dpk=0)
                    logging.info('Adding Flows Constraints for Scenario %d and Time %d' % (s,ti))

                if self.losses_plus and (self.TS or self.losses):
                    self.R_FlujoLossesPlus(data, fp, fn, n_l, n_a , dpk, ti, s)
                    logging.info('Adding Losses-Plus Constraints for Scenario %d and Time %d' % (s,ti))

        self.m.write('OPF_DownPower.lp')


    def ObjFunction(self, data):
        if self.pot_down:
            self.Csfc = data.Ccte_gen[data.pos_gen_agc_list] @ (self.vg_inc + self.vg_dec)            

        else:
            self.Csfc = data.Ccte_gen[data.pos_gen_agc_list] @ (self.vg_inc)

        f_obj = 0

        self.Cts = 0 
        self.Cop = 0
        self.Cue = 0
        self.Cpl = 0
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

                Cop_s += data.Cvar_gen[data.pos_gen_agc_list]*data.Sb @ (inc)
                #if np.size(Cvar_genstat) != 0:
                #    Cop_s += Genstat_on_AGC[:,s] * Cvar_genstat*sim.Sb @ pstatg_inc[:,s]

                Cue_s += (self.Voll * self.p_ens[:,s,ti]).sum() * data.Sb
                if self.losses and not self.losses_plus:
                    #Cpl_s += (self.Voll * self.ploss[:,s,ti]).sum() * data.Sb
                    Cpl_s += (self.Voll * self.ploss[:,s,ti]).sum() * data.Sb

                if self.TS:
                    Cts_s += self.costo_ts * (1 - self.s_ts[:,s,ti]).sum()

            self.Cts += Cts_s
            self.Cop += Cop_s
            self.Cue += Cue_s
            self.Cpl += Cpl_s

        f_obj = self.Csfc + (self.Cop + self.Cue + self.Cpl + self.Cts)/data.Ns

        self.m.setObjective(f_obj,self.gp.GRB.MINIMIZE)
        logging.info("Objective function set")

    def Balance(self, data, simm, s, ti):    
        logging.info("Adding balance constraint for s: %d, ti: %d", s, ti)
        # Con perdidas y gen pueden bajar potencia
        if self.pot_down:
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
        if  data.ngenstat != 0:
            f_genstat = data.SF[:,data.pos_genstat] @ data.Pgenstat_pre[:,ti]
        f_pv = 0
        if data.ngenpv != 0:
            f_pv = data.SF[:,data.pos_pv] @ data.Ppv_pre[:,ti]
        f_ens =  data.SF @ self.p_ens[:,s,ti]
        f_gen_agc = data.SF[:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        f_gen_out = data.SF[:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        Flujo_dda = data.SF @ simm.D_pfc[:,s,ti]

                
        #print('f_gen = ' + str(f_gen[0]))
        #print('Flujo_dda = ' + str(Flujo_dda[0]))
        #print('f_gen_out = ' + str(f_gen_out[0]))


        if self.losses:
            f_loss = 0.5 * data.SF @ abs(data.A.T) @ self.ploss[:,s,ti]
        else:
            f_loss = 0
        
        self.m.addConstr(self.f[:,s,ti] == f_gen + f_genstat + f_pv + f_ens + f_gen_agc - Flujo_dda - f_gen_out - f_loss, name='flujos') 

        if self.losses:
            self.m.addConstr(self.f[:,s,ti] == fp[:,s,ti] - fn[:,s,ti], name = 'f')
            self.m.addConstr(fp[:,s,ti] + fn[:,s,ti] == dpk[:,:,s,ti].sum(1), name='SumaDpk')
        
            kl = np.zeros((data.n_elem, self.L))
            for l in range(self.L):
                kl[:,l] = (2*(l+1)-1)*(data.FMax)/self.L
            self.m.addConstr(self.ploss[:,s,ti] == data.G/(data.B**2)*(self.gp.quicksum(kl[:,i]*dpk[:,i,s,ti] for i in range(self.L))), name = 'Ploss')  
        
            self.m.addConstr(-fp[:,s,ti] - fn[:,s,ti] - 0.5*self.ploss[:,s,ti] >= -1.1*data.FMax, name = 'fp s='+str(s)+' c='+str(ti))
            self.m.addConstr(-fp[:,s,ti] >= -data.FMax, name = 'fp+  s='+str(s)+' c='+str(ti))
            self.m.addConstr(-fn[:,s,ti] >= -data.FMax, name = 'fn+  s='+str(s)+' c='+str(ti))
            for l in range(self.L):
                self.m.addConstr(-dpk[:,l,s,ti] >= -data.FMax/self.L, name = 'LimiteDpk')

        else:
            self.m.addConstr(-self.f[:,s,ti] >= -data.FMax, name = 'fMax+ s='+str(s)+' c='+str(ti))
            self.m.addConstr(self.f[:,s,ti] >= -data.FMax, name = 'fMax- s='+str(s)+' c='+str(ti))

    #def Flujos_TS(self ,data ,simm, ti, s, f_ts, dpk, fp, fn, M = 1e6):
    def Flujos_TS(self ,data ,simm, ti, s, dpk, fp, fn, M = 100):

        if self.pot_down:
            inc = self.pg_inc[:,s,ti] - self.pg_dec[:,s,ti]

        else:
            inc = self.pg_inc[:,s,ti]


        f_loss_c = 0
        f_loss_nc = 0

        if self.losses:
            kl = np.zeros((data.n_elem, self.L))
            for l in range(self.L):
                kl[:,l] = (2*(l+1)-1)*(data.FMax)/self.L
            
            self.m.addConstr(self.ploss[:,s,ti] == data.G/(data.B**2)*(self.gp.quicksum(kl[:,i]*dpk[:,i,s,ti] for i in range(self.L))), name = 'TS_Ploss')  
        
        
            f_loss_nc = 0.5 * data.SF[data.pos_nots,:] @ abs(data.A.T) @ self.ploss[:,s,ti]
            f_loss_c = 0.5 * data.SF[data.pos_ts,:] @ abs(data.A.T) @ self.ploss[:,s,ti]

            for l in range(self.L):
                self.m.addConstr(-dpk[:,l,s,ti] >= -data.FMax/self.L, name = 'LimiteDpk')



        #Líneas NO candidatas

        fnots_gen = data.SF[data.pos_nots,:][:,data.pos_gen] @ data.Pgen_pre[:,ti]
        fnots_genstat = 0
        if  data.ngenstat != 0:
            fnots_genstat = data.SF[data.pos_nots,:][:,data.pos_genstat] @ data.Pgenstat_pre[:,ti]
        fnots_pv = 0
        if  data.ngenpv != 0:
            fnots_pv = data.SF[data.pos_nots,:][:,data.pos_pv] @ data.Ppv_pre[:,ti]
        fnots_ens = data.SF[data.pos_nots,:] @ self.p_ens[:,s,ti]
        fnots_gen_agc = 0
        if np.size(data.SF[data.pos_nots,:][:,data.pos_gen[data.pos_gen_agc_list]]) !=0:
            fnots_gen_agc = data.SF[data.pos_nots,:][:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        fnots_gen_out = 0
        if np.size(data.SF[data.pos_nots,:][:,int(simm.Barra_gen_out[s])]) != 0:
            fnots_gen_out = data.SF[data.pos_nots,:][:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        fnots_dda = data.SF[data.pos_nots,:] @ simm.D_pfc[:,s,ti]

        #print('fnots_gen = ' + str(fnots_gen[44]))
        #print('fnots_genstat = ' + str(fnots_genstat[44]))
        #print('fnots_pv = ' + str(fnots_pv[44]))
        #print('fnots_dda = ' + str(fnots_dda[44]))
        #print('fnots_gen_out = ' + str(fnots_gen_out[44]))



        fe = fnots_gen + fnots_genstat + fnots_pv + fnots_ens + fnots_gen_agc - fnots_dda - fnots_gen_out

        fv = (data.SF[data.pos_nots,:] @ data.A[data.pos_ts,:].T) @ self.f_ts[:,s,ti]

        #print('fe = ' + str(fe[44]))
        #print('fv = ' + str(fv[44]))
        #self.m.addConstr(-(fe-f_loss_nc+fv) >= -data.FMax[data.pos_nots], name = 'fe_p')
        #self.m.addConstr(fe-f_loss_nc+fv >= -data.FMax[data.pos_nots], name = 'fe_n')

        #Líneas candidatas
        fts_gen = data.SF[data.pos_ts,:][:,data.pos_gen] @ data.Pgen_pre[:,ti]
        fts_genstat = 0
        if data.ngenstat != 0:
            fts_genstat = data.SF[data.pos_ts,:][:,data.pos_genstat] @ data.Pgenstat_pre[:,ti]
        fts_pv = 0
        if data.ngenpv != 0:
            fts_pv = data.SF[data.pos_ts,:][:,data.pos_pv] @ data.Ppv_pre[:,ti]
        fts_ens = data.SF[data.pos_ts,:] @ self.p_ens[:,s,ti]        
        fts_gen_agc = data.SF[data.pos_ts,:][:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        fts_gen_out = data.SF[data.pos_ts,:][:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        fts_dda = data.SF[data.pos_ts,:] @ simm.D_pfc[:,s,ti]


        f1 = fts_gen + fts_genstat + fts_pv + fts_ens + fts_gen_agc - fts_dda - fts_gen_out
        f2 = self.f_ts[:,s,ti] - (data.SF[data.pos_ts,:] @ data.A[data.pos_ts,:].T) @ self.f_ts[:,s,ti]

        self.m.addConstr(self.f_ts[:,s,ti] <= M*(1 - self.s_ts[:,s,ti]), name = 'fs2_p') # 2
        self.m.addConstr(self.f_ts[:,s,ti] >= -M*(1 - self.s_ts[:,s,ti]), name = 'fs2_n') # 2

        # NO LOSSES
        #self.m.addConstr(f1-f2 <= 1.1*(data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1p_lc') # 1
        #self.m.addConstr(f1-f2 >= -1.1*(data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1n_lc')

        #self.m.addConstr(fe+fv <= 1.1*data.FMax[data.pos_nots], name = 'fs1p_lnc') # 1
        #self.m.addConstr(fe+fv >= -1.1*data.FMax[data.pos_nots], name = 'fs1n_lnc') # 1

        # LOSSES
        self.m.addConstr(f1-f_loss_c-f2 <= 1.1*(data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1_p') # 1
        self.m.addConstr(f1-f_loss_c-f2 >= -1.1*(data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1_n')

        self.m.addConstr(-fp[data.pos_ts,s,ti] - fn[data.pos_ts,s,ti] - 0.5*self.ploss[data.pos_ts,s,ti] >= -1.1*(data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1')
        self.m.addConstr(-fp[data.pos_nots,s,ti] - fn[data.pos_nots,s,ti] - 0.5*self.ploss[data.pos_nots,s,ti] >= -1.1*data.FMax[data.pos_nots], name = 'suma_f')

        self.m.addConstr(dpk[:,:,s,ti].sum(1) == fp[:,s,ti] + fn[:,s,ti], 'dpk')

        self.m.addConstr(self.f[:,s,ti] == fp[:,s,ti] - fn[:,s,ti], name = 'flujo')

        self.m.addConstr(self.f[data.pos_ts,s,ti] == f1-f_loss_c-f2)
        self.m.addConstr(self.f[data.pos_nots,s,ti] == fe-f_loss_nc+fv)


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

        self.m.addConstr(self.p_ens.sum()== 0)

        self.m.addConstr(-self.pg_inc[:,s,ti] >= -self.vg_inc * (data.Pmax_gen[data.pos_gen_agc_list] * data.ngen_par[data.pos_gen_agc_list] - data.Pgen_pre[data.pos_gen_agc_list,ti]), name = 'PMax_s='+str(s)+'_c='+str(ti))
        if self.pot_down:
            self.m.addConstr(-self.pg_dec[:,s,ti] >= -self.vg_dec * (data.Pgen_pre[data.pos_gen_agc_list,ti] - data.Pmin_gen[data.pos_gen_agc_list]  * data.ngen_par[data.pos_gen_agc_list]), name = 'PMin_s='+str(s)+'_c='+str(ti))
        #Nueva Lista que solo tiene unidades participantes en el AGC / se quita de unidades participantes la unidad que sale de servicio
        pos_part_gen_agc_list = data.pos_gen_agc_list[:]
        x = list(range(data.n_gen_agc))
        if data.Gen_Outages[s] in data.name_gen_agc_list:
            pos_part_gen_agc_list.remove(data.all_gen.index(data.Gen_Outages[s]))
            x.pop(data.name_gen_agc_list.index(data.Gen_Outages[s]))
            self.m.addConstr(self.pg_inc[data.name_gen_agc_list.index(data.Gen_Outages[s]), s, ti] == 0, name = 'GenOut+_s='+str(s)+'_c='+str(ti))
            if self.pot_down:
                self.m.addConstr(self.pg_dec[data.name_gen_agc_list.index(data.Gen_Outages[s]), s, ti] == 0, name = 'GenOut-_s='+str(s)+'_c='+str(ti))
        
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
            
            print_cts = ''
            if self.TS:
                print_cts = ' + Cts = %.2f' % (self.Cts.getValue()/data.Nt)


            if self.losses and not self.losses_plus:
                print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($) + Cpl = %.2f ($)' % (self.m.objVal,self.Csfc.getValue()/data.Nt,self.Cop.getValue()/data.Nt,self.Cue.getValue()/data.Nt,self.Cpl.getValue()/data.Nt) + print_cts)
            
            else:
                print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($)' % (self.m.objVal,self.Csfc.getValue()/data.Nt,self.Cop.getValue()/data.Nt,self.Cue.getValue()/data.Nt) + print_cts)


            print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (self.m.NumVars,self.m.NumConstrs,self.m.DNumNZs)) #print('num_Vars =  %d / num_Const =  %d' % (len(m.getVars()), len(m.getConstrs())))      
            #print ('Total P_loss = %.2f [MW]'%(pf.Sb*ploss.sum().getValue()))
            print('=> Solver time: %.4f (s)' % (self.m.Runtime))
            #print ('Costo => %.2f ($/h)' % self.m.objVal) 
            #print ('Las perdidas son %.2f (MW)' % sum(pk_loss[i].X for i in range(len(pk_loss)))) 

            self.Post_gen_out = np.zeros((data.Ns,data.Nt))
    
            self.part_factors = np.zeros((data.n_gen_agc,data.Ns,data.Nt))
            for gen, s, ti in np.ndindex(self.pg_inc.x.shape):
                
                if self.losses:
                    self.Post_gen_out[s, ti] = simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum() - simm.PL_pre_line[ti] + self.ploss[:,s,ti].x.sum()

                elif self.flujos:
                    self.Post_gen_out[s, ti] = simm.P_out[s,ti] - data.dda_barra[:,ti].sum() + simm.D_pfc[:,s,ti].sum()
                
                else:
                    self.Post_gen_out[s, ti] = simm.P_out[s,ti]

                
                if self.pot_down:
                    if self.pg_inc.x[gen,s,ti] != 0 or self.pg_dec.x[gen,s,ti] != 0:
                        self.part_factors[gen,s,ti] = (self.pg_inc.x[gen,s,ti]-self.pg_dec.x[gen,s,ti])/self.Post_gen_out[s,ti]
                
                else:
                    if self.pg_inc.x[gen,s,ti] != 0:
                        self.part_factors[gen,s,ti] = (self.pg_inc.x[gen,s,ti])/self.Post_gen_out[s,ti]
                   
            
            self.gen_csf = list()
            self.pos_gen_csf = list()
            cont=-1
            for gen in data.name_gen_agc_list:
                cont+=1
                if float(self.vg_inc[cont].x + self.vg_dec[cont].x) != 0:
                    self.gen_csf.append(gen)
                    self.pos_gen_csf.append(data.pos_gen_agc_list[cont])
        
        
        elif status == gp.GRB.Status.INF_OR_UNBD or \
            status == gp.GRB.Status.INFEASIBLE  or \
            status == gp.GRB.Status.UNBOUNDED:
            logging.error("Model is infeasible or unbounded")
            print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
            self.m.computeIIS() 
            self.m.write("GTCEP.ilp")
        logging.info("Results extraction complete")


class PartialModel(object):
    def __init__(self, data, big_optm, simm, short_sim, data_SF, P_change, scen, ti, t_int):
        self.gp = gp
        self.m = gp.Model('Modelo Variable AGC')
        self.m.Params.MIPGap = 1e-7

        self.vg_inc = big_optm.vg_inc.x
        #self.vg_inc = np.array([0., 0., 1., 1., 0., 0., 0., 1., 1., 0., 0., 0., 0.])
        #self.vg_inc = np.array([1, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        self.vg_dec = big_optm.vg_dec.x

        self.pg_inc = self.m.addMVar((data.n_gen_agc), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_inc')

        self.p_ens = self.m.addMVar((data.Nb), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='P_ens')

        if big_optm.pot_down:
            self.pg_dec = self.m.addMVar((data.n_gen_agc), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_dec')
        
        self.f = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=-self.gp.GRB.INFINITY, ub=self.gp.GRB.INFINITY, name='flujo') # Flujo por las líneas
        
        if big_optm.losses:
            fp = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fp') # Flujo-p por las líneas
            fn = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fn') # Flujo-n por las líneas

            self.ploss = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='losses') # Flujo por las líneas
            dpk = self.m.addMVar((data_SF.new_n_elem, big_optm.L), vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo

        if big_optm.losses_plus:
            n_l = self.m.addMVar((data_SF.new_n_elem), vtype=self.gp.GRB.BINARY, name='n_l')                                          # variable binaria complentaridad
            n_a = self.m.addMVar((data_SF.new_n_elem, big_optm.L), vtype=self.gp.GRB.BINARY, name='n_a')                                          # variable binaria adyacencia

        ## RESTRICCIONES

        self.ObjFunction(data, big_optm)

        self.Balance(data, big_optm, simm, short_sim, P_change, scen, ti)

        self.R_FlujoLosses(data, data_SF, big_optm, simm, short_sim, fp, fn, dpk, scen, ti)
        
        if big_optm.losses_plus:
            self.R_FlujoLossesPlus(data_SF, big_optm, fp, fn, n_l, n_a, dpk)

        self.R_Generadores(data, big_optm, simm, scen, ti)

        #####
        self.m.write('OPF_partial_OPT_' + str(t_int) + '.lp')

        ## OPTMIZAR
        self.m.optimize()


        ## Resulatdos 
        self.Results(data, big_optm, simm)


    def ObjFunction(self, data, big_optm):
        f_obj = 0

        self.Cpl = 0

        if big_optm.pot_down:
            inc = self.pg_inc + self.pg_dec
        else:
            inc = self.pg_inc

        self.Cop = data.Cvar_gen[data.pos_gen_agc_list]*data.Sb @ (inc)
        #if np.size(Cvar_genstat) != 0:
        #    Cop_s += Genstat_on_AGC[:,s] * Cvar_genstat*sim.Sb @ pstatg_inc[:,s]

        self.Cue = (big_optm.Voll * self.p_ens).sum() * data.Sb
        if big_optm.losses and not big_optm.losses_plus:
            self.Cpl = (big_optm.Voll * self.ploss).sum() * data.Sb

        f_obj = self.Cop + self.Cue + self.Cpl

        self.m.setObjective(f_obj,self.gp.GRB.MINIMIZE)

    def Balance(self, data, big_optm, simm, short_sim, P_change, scen, ti):    
        # Con perdidas y gen pueden bajar potenciati
        if big_optm.pot_down:
            #self.P_agc = self.pg_inc.sum() - self.pg_dec.sum() + self.p_ens.sum()
            self.m.addConstr(self.pg_inc.sum() - self.pg_dec.sum() + self.p_ens.sum() == P_change + simm.P_out[scen,ti] - data.dda_barra[:,ti].sum() + short_sim.D_t.sum() - simm.PL_pre_line[ti] + self.ploss.sum(), name ='Balance')
            #print('P_out = ' + str(simm.P_out[scen,ti]))
            #print('dda_barra = ' + str(data.dda_barra[:,ti].sum()))
            #print('D_pfc = ' + str(simm.D_pfc[:,scen,ti].sum()))        
            #print('PL_pre_line = ' + str(simm.PL_pre_line[ti]))        
        
        
        # Con perdidas 
        elif big_optm.flujos and big_optm.losses:
            self.m.addConstr(self.pg_inc.sum() + self.p_ens.sum() == simm.P_out[scen] - data.dda_barra.sum() + simm.D_pfc[:,scen,ti].sum()- simm.PL_pre_line[ti] + self.ploss.sum(), name ='Balance')
        
        # Sin Perdidas y con sistema de transmisión
        elif big_optm.flujos:
            self.m.addConstr(self.pg_inc.sum() + self.p_ens.sum() == simm.P_out[scen] - data.dda_barra.sum() + simm.D_pfc[:,scen,ti].sum(), name ='Balance')

        # Sin Perdidas / Inc potencia
        else:
            self.m.addConstr(self.pg_inc.sum() + self.p_ens.sum() == simm.P_out[scen], name ='Balance')

    def R_FlujoLosses(self, data, data_SF, big_optm, simm, short_sim, fp, fn, dpk, s, ti):

        if big_optm.pot_down:
            inc = self.pg_inc - self.pg_dec

        else:
            inc = self.pg_inc


        #f_gen =  data_SF.SF[:,data.pos_gen] @ data.Pgen_pre[:,ti]
        f_gen = data_SF.SF[:,data.pos_gen] @ short_sim.p_gen
        f_genstat = 0
        if data.ngenstat != 0:
            f_genstat =  data_SF.SF[:,data.pos_genstat] @ short_sim.p_genstat
        f_pv = 0
        if data.ngenpv != 0:
            f_pv = data_SF.SF[:,data.pos_pv] @ short_sim.pv_gen
        #f_genstat =  SF[:,simm.pos_genstat] @ simm.Pgenstat_pre
        f_ens =  data_SF.SF @ self.p_ens
        f_gen_agc = data_SF.SF[:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        f_genstat_agc = 0
        #f_gen_out = data_SF.SF[:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        #Flujo_dda = data.SF @ simm.D_pfc[:,s,ti]
        Flujo_dda = data_SF.SF @ short_sim.D_t
        if big_optm.losses:
            f_loss = 0.5 * data_SF.SF @ abs(data_SF.A.T) @ self.ploss
        else:
            f_loss = 0
        
        #print('f_gen = ' + str(f_gen[0]))
        #print('f_genstat = ' + str(f_genstat[0]))
        #print('f_pv = ' + str(f_pv[0]))
        #print('Flujo_dda = ' + str(Flujo_dda[0]))


        #print('f_gen = ' + str(f_gen[79]))
        #print('f_genstat = ' + str(f_genstat[79]))
        #print('f_pv = ' + str(f_pv[79]))
        #print('Flujo_dda = ' + str(Flujo_dda[79]))



        #self.m.addConstr(self.f == f_gen + f_genstat + f_pv + f_ens + f_gen_agc + f_genstat_agc - Flujo_dda - f_gen_out - f_loss, name='flujos') 
        self.m.addConstr(self.f == f_gen + f_genstat + f_pv + f_ens + f_gen_agc + f_genstat_agc - Flujo_dda - f_loss, name='flujos') 

        #if False:
        if big_optm.losses:
            self.m.addConstr(self.f == fp - fn, name = 'f')
            self.m.addConstr(fp + fn == dpk.sum(1), name='SumaDpk')
            #self.m.addConstr(self.p_ens.sum() == 0)

            kl = np.zeros((data_SF.new_n_elem, big_optm.L))
            for l in range(big_optm.L):
                kl[:,l] = (2*(l+1)-1)*(data_SF.new_FMax)/big_optm.L
            self.m.addConstr(self.ploss == data_SF.G/(data_SF.B**2)*(self.gp.quicksum(kl[:,i]*dpk[:,i] for i in range(big_optm.L))), name = 'Ploss')  
        
            self.m.addConstr(-fp - fn- 0.5*self.ploss >= -data_SF.new_FMax, name = 'fp+fn')
            self.m.addConstr(-fp >= -data_SF.new_FMax, name = 'fp+')
            self.m.addConstr(-fn >= -data_SF.new_FMax, name = 'fn+')
            for l in range(big_optm.L):
                self.m.addConstr(-dpk[:,l] >= -data_SF.new_FMax/big_optm.L, name = 'LimiteDpk')

        else:
            self.m.addConstr(-self.f >= -data_SF.new_FMax, name = 'fMax+')
            self.m.addConstr(self.f >= -data_SF.new_FMax, name = 'fMax-')

    def R_FlujoLossesPlus(self, data_SF, big_optm, fp, fn, n_l, n_a , dpk):

        self.m.addConstr(-fp >= -n_l*data_SF.new_FMax, name='fp_n')   #flujo positvo-restriccion de complementaridad
        self.m.addConstr(-fn >= (1-n_l)*(-data_SF.new_FMax), name='fn_n') #flujo nefativo-restriccion de complementaridad
        
        for l in range(big_optm.L): 
            if l==0:
                self.m.addConstr(-dpk[:,l]>= -data_SF.new_FMax/big_optm.L, name='d_f_Res_max_A_l')
                self.m.addConstr(dpk[:,l]>= n_a[:,l]*(data_SF.new_FMax/big_optm.L), name='d_f_Res_min_A_l')
            elif l==big_optm.L-1:
                self.m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*data_SF.new_FMax/big_optm.L, name='d_f_Res_max_A_L')
                self.m.addConstr(dpk[:,l] >=0, name='d_f_Res_min_A_L')
            else:
                self.m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*data_SF.new_FMax/big_optm.L, name='d_f_Res_max_A_L-1')
                self.m.addConstr(dpk[:,l] >=n_a[:,l]*data_SF.new_FMax/big_optm.L, name='d_f_Res_min_A_L-1')

    def R_Generadores(self, data, big_optm, simm, s, ti):

        self.m.addConstr(-self.pg_inc >= -self.vg_inc * (data.Pmax_gen [data.pos_gen_agc_list]* data.ngen_par[data.pos_gen_agc_list] - data.Pgen_pre[data.pos_gen_agc_list,ti]), name = 'PMax')
        if big_optm.pot_down:
            self.m.addConstr(-self.pg_dec >= -self.vg_dec * (data.Pgen_pre[data.pos_gen_agc_list,ti] - data.Pmin_gen [data.pos_gen_agc_list]* data.ngen_par[data.pos_gen_agc_list]), name = 'PMin')
        #Nueva Lista que solo tiene unidades participantes en el AGC / se quita de unidades participantes la unidad que sale de servicio
        pos_part_gen_agc_list = data.pos_gen_agc_list[:]
        x = list(range(data.n_gen_agc))
        if data.Gen_Outages[s] in data.name_gen_agc_list:
            pos_part_gen_agc_list.remove(data.all_gen.index(data.Gen_Outages[s]))
            x.pop(data.name_gen_agc_list.index(data.Gen_Outages[s]))
            self.m.addConstr(self.pg_inc[data.name_gen_agc_list.index(data.Gen_Outages[s])]== 0, name = 'GenOut+')
            if big_optm.pot_down:
                self.m.addConstr(self.pg_dec[data.name_gen_agc_list.index(data.Gen_Outages[s])] == 0, name = 'GenOut-')
        self.m.addConstr(-self.pg_inc[x] >= -big_optm.T_Sfc * data.Ramp_gen[pos_part_gen_agc_list] - (simm.Pgen_pfc[pos_part_gen_agc_list,s,ti] - data.Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + data.Gen_Outages[s] +'_Ramp+')
        if big_optm.pot_down:
            self.m.addConstr(-self.pg_dec[x] >= -big_optm.T_Sfc * data.Ramp_gen[pos_part_gen_agc_list] + (simm.Pgen_pfc[pos_part_gen_agc_list,s,ti] - data.Pgen_pre[pos_part_gen_agc_list,ti]), name = 'E' + str(ti+1) + '-' + data.Gen_Outages[s] +'_Ramp-')
        #if np.size(Cvar_genstat) != 0:
            #self.m.addConstr(-pg_inc[pos_genstat_agc_list, s, ti] >= -vstatg_inc * (Pmax_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list] - Pgen_pre[pos_genstat_agc_list,ti]))
            #self.m.addConstr(-pg_dec[pos_genstat_agc_list, s, ti] >= -vstatg_dec * (Pgen_pre[pos_genstat_agc_list,ti] - Pmin_gen[pos_genstat_agc_list] * ngen_par[pos_genstat_agc_list]))
            #self.m.addConstr(-pg_inc[pos_genstat_agc_list,s] >= -T_Sfc * Ramp_genstat[pos_genstat_agc_list] - (Pgenstat_pfc[pos_genstat_agc_list, s] - Pgenstat_pre[pos_genstat_agc_list]))
            #self.m.addConstr(-pg_dec[pos_genstat_agc_list,s] >= -T_Sfc * Ramp_genstat[pos_genstat_agc_list] + (Pgenstat_pfc[pos_genstat_agc_list, s] - Pgenstat_pre[pos_genstat_agc_list]))


        self.m.addConstr(-self.p_ens >= -simm.D_pfc[:,s,ti], 'LimENS')


    def Results(self, data, big_optm, simm):
        status = self.m.Status
        if status == gp.GRB.Status.OPTIMAL:
            print('-----------------------------')

            print('Costo = %.2f' % (self.m.ObjVal))

            self.part_factors = np.zeros(data.n_gen_agc)
            for i in range(data.n_gen_agc):
                if big_optm.pot_down:
                    inc = self.pg_inc[i].x - self.pg_dec[i].x
                    self.P_agc = self.pg_inc.x.sum() - self.pg_dec.x.sum() + self.p_ens.x.sum()
                    if inc != 0:
                        self.part_factors[i] = inc/self.P_agc
                else:
                    if self.pg_inc[i].x != 0:
                        self.P_agc = self.pg_inc.x.sum() + self.p_ens.x.sum()
                        self.part_factors[i] = self.pg_inc[i].x/self.P_agc

        elif status == gp.GRB.Status.INF_OR_UNBD or \
            status == gp.GRB.Status.INFEASIBLE  or \
            status == gp.GRB.Status.UNBOUNDED:
            print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
            self.m.computeIIS() 
            self.m.write("GTCEP_partial.ilp")


  
class Model_UC(object):
    def __init__(self, data):
        self.gp = gp
        self.m = gp.Model('Modelo UC')
        self.m.Params.MIPGap = 1e-7

        self.pg = self.m.addMVar((data.ngen + data.ngenstat + data.ngenpv), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_inc')
        self.f = self.m.addMVar((data.n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=-self.gp.GRB.INFINITY, ub=self.gp.GRB.INFINITY, name='flujo') # Flujo por las líneas

        fp = self.m.addMVar((data.n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fp') # Flujo-p por las líneas
        fn = self.m.addMVar((data.n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fn') # Flujo-n por las líneas

        dpk = self.m.addMVar((data.n_elem, data.L), vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='dpk')   # Perdidas por cada línea por cada tramo
        self.ploss = self.m.addMVar((data.n_elem), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Ploss')

        n_l = self.m.addMVar((data.n_elem), vtype=self.gp.GRB.BINARY, name='n_l')   # variable binaria complentaridad
        n_a = self.m.addMVar((data.n_elem, data.L), vtype=self.gp.GRB.BINARY, name='n_a')  # variable binaria adyacencia

        f_obj = 0
        Cop = 0
        Cop += data.Cvar_gen*data.Sb @ self.pg[:data.ngen]
        Cop += data.Cvar_genstat*data.Sb @ self.pg[data.ngen:data.ngen+data.ngenstat]
        Cop += data.Cvar_pv*data.Sb @ self.pg[data.ngen+data.ngenstat:data.ngen+data.ngenstat+data.ngenpv]
        f_obj += Cop

        self.m.setObjective(f_obj,self.gp.GRB.MINIMIZE)
        
        #Balance
        self.m.addConstr(self.pg.sum() == data.dda_barra[:,0].sum() + self.ploss.sum(), name ='Balance')

        #Flujos
        f_gen = data.SF[:,data.pos_gen] @ self.pg[:data.ngen]
        f_genstat = data.SF[:,data.pos_genstat] @ self.pg[data.ngen:data.ngen+data.ngenstat]
        f_pv = data.SF[:,data.pos_pv] @ self.pg[data.ngen+data.ngenstat:data.ngen+data.ngenstat+data.ngenpv]
        Flujo_dda = data.SF @ data.dda_barra[:,0]
        f_loss = 0.5 * data.SF @ abs(data.A.T) @ self.ploss

        self.m.addConstr(self.f == f_gen + f_genstat + f_pv - Flujo_dda - f_loss, name='flujos')

        self.m.addConstr(self.f == fp - fn, name = 'f')
        self.m.addConstr(fp + fn == dpk.sum(1), name='SumaDpk')
        self.m.addConstr(self.ploss == data.G/(data.B**2)*(self.gp.quicksum((2*(l+1)-1)*(data.FMax)/data.L*dpk[:,l] for l in range(data.L))), name = 'Ploss')

        self.m.addConstr(-fp - fn- 0.5*self.ploss >= -0.9*data.FMax, name = 'fp+fn')
        self.m.addConstr(-fp >= -data.FMax, name = 'fp+')
        self.m.addConstr(-fn >= -data.FMax, name = 'fn+')

        for l in range(data.L):
            self.m.addConstr(-dpk[:,l] >= -data.FMax/data.L, name = 'LimiteDpk')

        self.m.addConstr(-fp >= -n_l*data.FMax, name='fp_n')   #flujo positvo-restriccion de complementaridad
        self.m.addConstr(-fn >= (1-n_l)*(-data.FMax), name='fn_n') #flujo nefativo-restriccion de complementaridad

        for l in range(data.L):
            if l==0:
                self.m.addConstr(-dpk[:,l]>= -data.FMax/data.L, name='d_f_Res_max_A_l')
                self.m.addConstr(dpk[:,l]>= n_a[:,l]*(data.FMax/data.L), name='d_f_Res_min_A_l')
            elif l==data.L-1:
                self.m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*data.FMax/data.L, name='d_f_Res_max_A_L')
                self.m.addConstr(dpk[:,l] >=0, name='d_f_Res_min_A_L')
            else:
                self.m.addConstr(-dpk[:,l] >=-n_a[:,l-1]*data.FMax/data.L, name='d_f_Res_max_A_L-1')
                self.m.addConstr(dpk[:,l] >=n_a[:,l]*data.FMax/data.L, name='d_f_Res_min_A_L-1')

        self.m.addConstr(self.pg[:data.ngen] >= data.Pmin_gen*data.ngen_par, name = 'PMin')
        self.m.addConstr(self.pg[:data.ngen] <= data.Pmax_gen*data.ngen_par, name = 'PMax')

        self.m.addConstr(self.pg[data.ngen:data.ngen+data.ngenstat] >= data.Pmin_genstat*data.ngenstat_par, name = 'PMinStat')
        self.m.addConstr(self.pg[data.ngen:data.ngen+data.ngenstat] <= data.Pmax_genstat*data.ngenstat_par, name = 'PMaxStat')

        self.m.addConstr(self.pg[data.ngen+data.ngenstat:data.ngen+data.ngenstat+data.ngenpv] >= data.Pmin_pv*data.ngenpv_par, name = 'PMinPV')
        self.m.addConstr(self.pg[data.ngen+data.ngenstat:data.ngen+data.ngenstat+data.ngenpv] <= data.Pmax_pv*data.ngenpv_par, name = 'PMaxPV')

        self.m.optimize()

        status = self.m.Status
        if status == gp.GRB.Status.OPTIMAL:
            print('-----------------------------')
            print('Costo = %.2f' % (self.m.ObjVal))
            print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (self.m.NumVars,self.m.NumConstrs,self.m.DNumNZs))
            print('=> Solver time: %.4f (s)' % (self.m.Runtime))
        elif status == gp.GRB.Status.INF_OR_UNBD or \
            status == gp.GRB.Status.INFEASIBLE  or \
            status == gp.GRB.Status.UNBOUNDED:
            print('The model cannot be solved because it is infeasible or unbounded => status "%d"' % status)
            self.m.computeIIS() 
            self.m.write("GTCEP_UC.ilp")


class Model_CTF(object):
    def __init__(self, data, simm, big_optm, short_sim, data_SF, scen, ti):
        self.gp = gp
        self.m = gp.Model('Modelo UC CTF')
        self.m.Params.MIPGap = 1e-7

        self.vg_inc = big_optm.vg_inc.x
        self.vg_dec = big_optm.vg_dec.x

        self.pg = self.m.addMVar((data.n_gen_agc), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_inc')

        self.p_ens = self.m.addMVar((data.Nb), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='P_ens')

        if big_optm.pot_down:
            self.pg_dec = self.m.addMVar((data.n_gen_agc), vtype=self.gp.GRB.CONTINUOUS, ub=self.gp.GRB.INFINITY, lb=0, name='Pg_dec')
        
        self.f = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=-self.gp.GRB.INFINITY, ub=self.gp.GRB.INFINITY, name='flujo') # Flujo por las líneas
        
        if big_optm.losses:
            fp = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fp') # Flujo-p por las líneas
            fn = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='fn') # Flujo-n por las líneas

            self.ploss = self.m.addMVar((data_SF.new_n_elem),vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self.gp.GRB.INFINITY, name='losses') # Flujo por las líneas
            dpk = self.m.addMVar((data_SF.new_n_elem, big_optm.L), vtype=self.gp.GRB.CONTINUOUS, lb=0, ub=self)





if False:
    def MoreSolutions(self,data,simm):
        print('-----------------------------')
        #print('La demanda total del sistema es: %.2f (MW)' % (pf.dda_barra[:,ti].sum()*pf.Sb))
            
        print_cts = ''
        if self.TS:
            print_cts = ' + Cts = %.2f' % (self.Cts.getValue())


        if self.losses and not self.losses_plus:
            print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($) + Cpl = %.2f ($)' % (self.m.objVal,self.Csfc.getValue(),self.Cop.getValue(),self.Cue.getValue(),self.Cpl.getValue()) + print_cts)
            
        else:
            print ('Cost = %.2f ($) => Csfc = %.2f ($) + Cop = %.2f ($) + Cue = %.2f ($)' % (self.m.objVal,self.Csfc.getValue(),self.Cop.getValue(),self.Cue.getValue()) + print_cts)


        print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (self.m.NumVars,self.m.NumConstrs,self.m.DNumNZs)) #print('num_Vars =  %d / num_Const =  %d' % (len(m.getVars()), len(m.getConstrs())))      
        #print ('Total P_loss = %.2f [MW]'%(pf.Sb*ploss.sum().getValue()))
        print('=> Solver time: %.4f (s)' % (self.m.Runtime))
        print ('Costo => %.2f ($/h)' % self.m.objVal) 
        #print ('Las perdidas son %.2f (MW)' % sum(pk_loss[i].X for i in range(len(pk_loss)))




if False:
     def Flujos_TS(self,data,simm, ti, s, f_ts, dpk, fp, fn, M = 1e6):

        if self.pot_down:
            inc = self.pg_inc[:,s,ti] - self.pg_dec[:,s,ti]

        else:
            inc = self.pg_inc[:,s,ti]


        f_loss_c = 0
        f_loss_nc = 0

        if self.losses:
            kl = np.zeros((data.n_elem, self.L))
            for l in range(self.L):
                kl[:,l] = (2*(l+1)-1)*(data.FMax)/self.L
            
            self.m.addConstr(self.ploss[:,s,ti] == data.G/(data.B**2)*(gp.quicksum(kl[:,i]*dpk[:,i,s,ti] for i in range(self.L))), name = 'Ploss')  
        
        
            f_loss_nc = 0.5 * data.SF[data.pos_nots,:] @ abs(data.A[data.pos_nots,:].T) @ self.ploss[data.pos_nots,s,ti]
            f_loss_c = 0.5 * data.SF[data.pos_ts,:] @ abs(data.A[data.pos_ts,:].T) @ self.ploss[data.pos_ts,s,ti]



        #Líneas NO candidatas

        fnots_gen = data.SF[data.pos_nots,:][:,data.pos_gen] @ data.Pgen_pre[:,ti]
        fnots_ens = data.SF[data.pos_nots,:] @ self.p_ens[:,s,ti]
        if np.size(data.SF[data.pos_nots,:][:,data.pos_gen[data.pos_gen_agc_list]]) !=0:
            fnots_gen_agc = data.SF[data.pos_nots,:][:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        else:
            fnots_gen_agc= 0
        if np.size(data.SF[data.pos_nots,:][:,int(simm.Barra_gen_out[s])]) != 0:
            fnots_gen_out = data.SF[data.pos_nots,:][:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        else:
            fnots_gen_out =0
        fnots_dda = data.SF[data.pos_nots,:] @ simm.D_pfc[:,s,ti]

        fe = fnots_gen + fnots_ens + fnots_gen_agc - fnots_dda - fnots_gen_out

        fv = (data.SF[data.pos_nots,:] @ data.A[data.pos_ts,:].T) @ f_ts[:,s,ti]

        self.m.addConstr(-(fe-f_loss_nc+fv) >= -data.FMax[data.pos_nots], name = 'fe_p')
        self.m.addConstr(fe-f_loss_nc+fv >= -data.FMax[data.pos_nots], name = 'fe_n')

        #Líneas candidatas
        fts_gen = data.SF[data.pos_ts,:][:,data.pos_gen] @ data.Pgen_pre[:,ti]
        fts_ens = data.SF[data.pos_ts,:] @ self.p_ens[:,s,ti]        
        fts_gen_agc = data.SF[data.pos_ts,:][:,data.pos_gen[data.pos_gen_agc_list]] @ (inc)
        fts_gen_out = data.SF[data.pos_ts,:][:,int(simm.Barra_gen_out[s])]*simm.P_out[s,ti]
        fts_dda = data.SF[data.pos_ts,:] @ simm.D_pfc[:,s,ti]


        f1 = fts_gen + fts_ens + fts_gen_agc - fts_dda - fts_gen_out

        f2 = f_ts[:,s,ti] - (data.SF[data.pos_ts,:] @ data.A[data.pos_ts,:].T) @ f_ts[:,s,ti]

        self.m.addConstr(f1-f_loss_c-f2 <= (data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1_p') # 1
        self.m.addConstr(f1-f_loss_c-f2 >= -(data.FMax[data.pos_ts]) * self.s_ts[:,s,ti], name = 'fs1_n')

        self.m.addConstr(f_ts[:,s,ti] <= M*(1 - self.s_ts[:,s,ti]), name = 'fs2_p') # 2
        self.m.addConstr(f_ts[:,s,ti] >= -M*(1 - self.s_ts[:,s,ti]), name = 'fs2_n') # 2

        self.m.addConstr(dpk[:,:,s,ti].sum(1) == fp[:,s,ti] + fn[:,s,ti], 'dpk')

        self.m.addConstr(self.f[:,s,ti] == fp[:,s,ti] - fn[:,s,ti], name = 'flujo')

        self.m.addConstr(self.f[data.pos_ts,s,ti] == f1-f_loss_c-f2)
        self.m.addConstr(self.f[data.pos_nots,s,ti] == fe-f_loss_nc+fv)

        self.m.addConstr(-fp[:,s,ti] - fn[:,s,ti] - 0.5*self.ploss[:,s,ti] >= -data.FMax, name = 'suma_f')
