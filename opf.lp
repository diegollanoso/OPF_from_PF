\ Model Modelo 1
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  80 Pg[2] + 100 Pg[3] + 100 Pg[4] + 150 Pg[5] + 150 Pg[6] + 200 Constant
Subject To
 Pmin: Pg[0] >= 0
 Pmax: - Pg[0] >= 0
 Pmin: Pg[1] >= 0
 Pmax: - Pg[1] >= 0
 Pmin: Pg[2] >= 100
 Pmax: - Pg[2] >= -492
 Pmin: Pg[3] >= 100
 Pmax: - Pg[3] >= -500
 Pmin: Pg[4] >= 100
 Pmax: - Pg[4] >= -500
 Pmin: Pg[5] >= 100
 Pmax: - Pg[5] >= -492
 Pmin: Pg[6] >= 100
 Pmax: - Pg[6] >= -492
 Pmin_pv: Pg_pv[0] >= 0
 Pmax_pv: - Pg_pv[0] >= -110
 Pmin_pe: Pg_pe[0] >= 0
 Pmax_pe: - Pg_pe[0] >= -2
 Balance: Pg[0] + Pg[1] + Pg[2] + Pg[3] + Pg[4] + Pg[5] + Pg[6] + Pg_pv[0]
   + 25 Pg_pe[0] = 760.8899912834167
Bounds
 Constant = 1
End
