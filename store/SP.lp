\ Model sub problem
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  22 x_0_0 + 33 x_0_1 + 24 x_0_2 + 33 x_1_0 + 23 x_1_1 + 30 x_1_2 + 20 x_2_0
   + 25 x_2_1 + 27 x_2_2
Subject To
 sub_capacity_1_z_0: x_0_0 + x_0_1 + x_0_2 <= 772
 sub_capacity_1_z_1: x_1_0 + x_1_1 + x_1_2 <= 0
 sub_capacity_1_z_2: x_2_0 + x_2_1 + x_2_2 <= 0
 R3: - d_0 + x_0_0 + x_1_0 + x_2_0 >= 0
 R4: - d_1 + x_0_1 + x_1_1 + x_2_1 >= 0
 R5: - d_2 + x_0_2 + x_1_2 + x_2_2 >= 0
 R6: pi_0 - theta_0 <= 22
 R7: pi_0 - theta_1 <= 33
 R8: pi_0 - theta_2 <= 24
 R9: pi_1 - theta_0 <= 33
 R10: pi_1 - theta_1 <= 23
 R11: pi_1 - theta_2 <= 30
 R12: pi_2 - theta_0 <= 20
 R13: pi_2 - theta_1 <= 25
 R14: pi_2 - theta_2 <= 27
 R15: - 40 g_0 + d_0 = 206
 R16: - 40 g_1 + d_1 = 274
 R17: - 40 g_2 + d_2 = 220
 R18: g_0 + g_1 + g_2 <= 1.8
 R19: g_0 + g_1 <= 1.2
 R20: - pi_0 - 100000 v_0 <= 0
 sub_capacity_2_z_0: 100000 v_0 - x_0_0 - x_0_1 - x_0_2 <= 99228
 R22: - pi_1 - 100000 v_1 <= 0
 sub_capacity_2_z_1: 100000 v_1 - x_1_0 - x_1_1 - x_1_2 <= 100000
 R24: - pi_2 - 100000 v_2 <= 0
 sub_capacity_2_z_2: 100000 v_2 - x_2_0 - x_2_1 - x_2_2 <= 100000
 R26: - 100000 w_0 - theta_0 <= 0
 R27: 100000 w_0 - d_0 + x_0_0 + x_1_0 + x_2_0 <= 100000
 R28: - 100000 w_1 - theta_1 <= 0
 R29: 100000 w_1 - d_1 + x_0_1 + x_1_1 + x_2_1 <= 100000
 R30: - 100000 w_2 - theta_2 <= 0
 R31: 100000 w_2 - d_2 + x_0_2 + x_1_2 + x_2_2 <= 100000
 R32: - 100000 h_0_0 + x_0_0 <= 0
 R33: - pi_0 + theta_0 + 100000 h_0_0 <= 99978
 R34: - 100000 h_1_0 + x_1_0 <= 0
 R35: - pi_1 + theta_0 + 100000 h_1_0 <= 99967
 R36: - 100000 h_2_0 + x_2_0 <= 0
 R37: - pi_2 + theta_0 + 100000 h_2_0 <= 99980
 R38: - 100000 h_0_1 + x_0_1 <= 0
 R39: - pi_0 + theta_1 + 100000 h_0_1 <= 99967
 R40: - 100000 h_1_1 + x_1_1 <= 0
 R41: - pi_1 + theta_1 + 100000 h_1_1 <= 99977
 R42: - 100000 h_2_1 + x_2_1 <= 0
 R43: - pi_2 + theta_1 + 100000 h_2_1 <= 99975
 R44: - 100000 h_0_2 + x_0_2 <= 0
 R45: - pi_0 + theta_2 + 100000 h_0_2 <= 99976
 R46: - 100000 h_1_2 + x_1_2 <= 0
 R47: - pi_1 + theta_2 + 100000 h_1_2 <= 99970
 R48: - 100000 h_2_2 + x_2_2 <= 0
 R49: - pi_2 + theta_2 + 100000 h_2_2 <= 99973
Bounds
 -infinity <= pi_0 <= 0
 -infinity <= pi_1 <= 0
 -infinity <= pi_2 <= 0
 g_0 <= 1
 -infinity <= theta_0 <= 0
 g_1 <= 1
 -infinity <= theta_1 <= 0
 g_2 <= 1
 -infinity <= theta_2 <= 0
Binaries
 v_0 v_1 v_2 w_0 w_1 w_2 h_0_0 h_0_1 h_0_2 h_1_0 h_1_1 h_1_2 h_2_0 h_2_1
 h_2_2
End
