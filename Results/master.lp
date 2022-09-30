\ Model master_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  _Theta
Subject To
 R0: _f[9,Do_nothing] + _f[9,R1] + _f[9,R2] + _f[9,R3] = 1
 R1: _f[6,Do_nothing] + _f[6,R1] + _f[6,R2] + _f[6,R3] = 1
 R2: _f[2,Do_nothing] + _f[2,R1] + _f[2,R2] + _f[2,R3] = 1
 R3: _f[10,Do_nothing] + _f[10,R1] + _f[10,R2] + _f[10,R3] = 1
 R4: _f[5,Do_nothing] + _f[5,R1] + _f[5,R2] + _f[5,R3] = 1
 R5: _f[1,Do_nothing] + _f[1,R1] + _f[1,R2] + _f[1,R3] = 1
 R6: _f[4,Do_nothing] + _f[4,R1] + _f[4,R2] + _f[4,R3] = 1
 R7: _f[8,Do_nothing] + _f[8,R1] + _f[8,R2] + _f[8,R3] = 1
 R8: _f[3,Do_nothing] + _f[3,R1] + _f[3,R2] + _f[3,R3] = 1
 R9: _f[7,Do_nothing] + _f[7,R1] + _f[7,R2] + _f[7,R3] = 1
 R10: _Theta >= 0
 R11: 1.330991e+06 _f[9,R1] + 1.672772e+06 _f[9,R2] + 3.002586e+06 _f[9,R3]
   + 745011 _f[6,R1] + 945943 _f[6,R2] + 1.697949e+06 _f[6,R3]
   + 1.470866e+06 _f[2,R1] + 1.84857e+06 _f[2,R2] + 3.318145e+06 _f[2,R3]
   + 1.841752e+06 _f[10,R1] + 2.314703e+06 _f[10,R2]
   + 4.154833e+06 _f[10,R3] + 1.427401e+06 _f[5,R1] + 1.808916e+06 _f[5,R2]
   + 3.246972e+06 _f[5,R3] + 754217 _f[1,R1] + 988682 _f[1,R2]
   + 1.774661e+06 _f[1,R3] + 1.317138e+06 _f[4,R1] + 1.655368e+06 _f[4,R2]
   + 2.971353e+06 _f[4,R3] + 2.894596e+06 _f[8,R1] + 3.637895e+06 _f[8,R2]
   + 6.529949e+06 _f[8,R3] + 1.800986e+06 _f[3,R1] + 2.263469e+06 _f[3,R2]
   + 4.062874e+06 _f[3,R3] + 527129 _f[7,R1] + 662486 _f[7,R2]
   + 1.189149e+06 _f[7,R3] + 1.67294641e+08 _r0[9,Do_nothing,Recover]
   + 1.6724473e+08 _r0[9,R1,Recover] + 1.19991452e+08 _r0[9,R2,Recover]
   + 7.3054158e+07 _r0[9,R3,Recover]
   + 1.12287189e+08 _r0[6,Do_nothing,Recover]
   + 1.12253953e+08 _r0[6,R1,Recover] + 8.0694437e+07 _r0[6,R2,Recover]
   + 4.93764e+07 _r0[6,R3,Recover]
   + 2.34957414e+08 _r0[2,Do_nothing,Recover]
   + 2.34887316e+08 _r0[2,R1,Recover] + 1.68522322e+08 _r0[2,R2,Recover]
   + 1.02601111e+08 _r0[2,R3,Recover]
   + 2.68081439e+08 _r0[10,Do_nothing,Recover]
   + 2.68001459e+08 _r0[10,R1,Recover] + 1.92280404e+08 _r0[10,R2,Recover]
   + 1.17065697e+08 _r0[10,R3,Recover]
   + 3.10887702e+08 _r0[5,Do_nothing,Recover]
   + 3.10795373e+08 _r0[5,R1,Recover] + 2.23233935e+08 _r0[5,R2,Recover]
   + 1.36306708e+08 _r0[5,R3,Recover]
   + 1.77010253e+08 _r0[1,Do_nothing,Recover]
   + 1.76958131e+08 _r0[1,R1,Recover] + 1.27368662e+08 _r0[1,R2,Recover]
   + 7.819011e+07 _r0[1,R3,Recover]
   + 1.76857603e+08 _r0[4,Do_nothing,Recover]
   + 1.76804839e+08 _r0[4,R1,Recover] + 1.26850451e+08 _r0[4,R2,Recover]
   + 7.7230108e+07 _r0[4,R3,Recover]
   + 3.64132206e+08 _r0[8,Do_nothing,Recover]
   + 3.6402357e+08 _r0[8,R1,Recover] + 2.61172456e+08 _r0[8,R2,Recover]
   + 1.59009108e+08 _r0[8,R3,Recover]
   + 3.96428872e+08 _r0[3,Do_nothing,Recover]
   + 3.963106e+08 _r0[3,R1,Recover] + 2.84337119e+08 _r0[3,R2,Recover]
   + 1.73112404e+08 _r0[3,R3,Recover]
   + 1.46873379e+08 _r0[7,Do_nothing,Recover]
   + 1.46829561e+08 _r0[7,R1,Recover] + 1.05344379e+08 _r0[7,R2,Recover]
   + 6.4136609e+07 _r0[7,R3,Recover] <= 3e+07
 R12: - _f[9,Do_nothing] + _r0[9,Do_nothing,Do_nothing]
   + _r0[9,Do_nothing,Recover] = 0
 R13: - _f[9,R1] + _r0[9,R1,Do_nothing] + _r0[9,R1,Recover] = 0
 R14: - _f[9,R2] + _r0[9,R2,Do_nothing] + _r0[9,R2,Recover] = 0
 R15: - _f[9,R3] + _r0[9,R3,Do_nothing] + _r0[9,R3,Recover] = 0
 R16: - _f[6,Do_nothing] + _r0[6,Do_nothing,Do_nothing]
   + _r0[6,Do_nothing,Recover] = 0
 R17: - _f[6,R1] + _r0[6,R1,Do_nothing] + _r0[6,R1,Recover] = 0
 R18: - _f[6,R2] + _r0[6,R2,Do_nothing] + _r0[6,R2,Recover] = 0
 R19: - _f[6,R3] + _r0[6,R3,Do_nothing] + _r0[6,R3,Recover] = 0
 R20: - _f[2,Do_nothing] + _r0[2,Do_nothing,Do_nothing]
   + _r0[2,Do_nothing,Recover] = 0
 R21: - _f[2,R1] + _r0[2,R1,Do_nothing] + _r0[2,R1,Recover] = 0
 R22: - _f[2,R2] + _r0[2,R2,Do_nothing] + _r0[2,R2,Recover] = 0
 R23: - _f[2,R3] + _r0[2,R3,Do_nothing] + _r0[2,R3,Recover] = 0
 R24: - _f[10,Do_nothing] + _r0[10,Do_nothing,Do_nothing]
   + _r0[10,Do_nothing,Recover] = 0
 R25: - _f[10,R1] + _r0[10,R1,Do_nothing] + _r0[10,R1,Recover] = 0
 R26: - _f[10,R2] + _r0[10,R2,Do_nothing] + _r0[10,R2,Recover] = 0
 R27: - _f[10,R3] + _r0[10,R3,Do_nothing] + _r0[10,R3,Recover] = 0
 R28: - _f[5,Do_nothing] + _r0[5,Do_nothing,Do_nothing]
   + _r0[5,Do_nothing,Recover] = 0
 R29: - _f[5,R1] + _r0[5,R1,Do_nothing] + _r0[5,R1,Recover] = 0
 R30: - _f[5,R2] + _r0[5,R2,Do_nothing] + _r0[5,R2,Recover] = 0
 R31: - _f[5,R3] + _r0[5,R3,Do_nothing] + _r0[5,R3,Recover] = 0
 R32: - _f[1,Do_nothing] + _r0[1,Do_nothing,Do_nothing]
   + _r0[1,Do_nothing,Recover] = 0
 R33: - _f[1,R1] + _r0[1,R1,Do_nothing] + _r0[1,R1,Recover] = 0
 R34: - _f[1,R2] + _r0[1,R2,Do_nothing] + _r0[1,R2,Recover] = 0
 R35: - _f[1,R3] + _r0[1,R3,Do_nothing] + _r0[1,R3,Recover] = 0
 R36: - _f[4,Do_nothing] + _r0[4,Do_nothing,Do_nothing]
   + _r0[4,Do_nothing,Recover] = 0
 R37: - _f[4,R1] + _r0[4,R1,Do_nothing] + _r0[4,R1,Recover] = 0
 R38: - _f[4,R2] + _r0[4,R2,Do_nothing] + _r0[4,R2,Recover] = 0
 R39: - _f[4,R3] + _r0[4,R3,Do_nothing] + _r0[4,R3,Recover] = 0
 R40: - _f[8,Do_nothing] + _r0[8,Do_nothing,Do_nothing]
   + _r0[8,Do_nothing,Recover] = 0
 R41: - _f[8,R1] + _r0[8,R1,Do_nothing] + _r0[8,R1,Recover] = 0
 R42: - _f[8,R2] + _r0[8,R2,Do_nothing] + _r0[8,R2,Recover] = 0
 R43: - _f[8,R3] + _r0[8,R3,Do_nothing] + _r0[8,R3,Recover] = 0
 R44: - _f[3,Do_nothing] + _r0[3,Do_nothing,Do_nothing]
   + _r0[3,Do_nothing,Recover] = 0
 R45: - _f[3,R1] + _r0[3,R1,Do_nothing] + _r0[3,R1,Recover] = 0
 R46: - _f[3,R2] + _r0[3,R2,Do_nothing] + _r0[3,R2,Recover] = 0
 R47: - _f[3,R3] + _r0[3,R3,Do_nothing] + _r0[3,R3,Recover] = 0
 R48: - _f[7,Do_nothing] + _r0[7,Do_nothing,Do_nothing]
   + _r0[7,Do_nothing,Recover] = 0
 R49: - _f[7,R1] + _r0[7,R1,Do_nothing] + _r0[7,R1,Recover] = 0
 R50: - _f[7,R2] + _r0[7,R2,Do_nothing] + _r0[7,R2,Recover] = 0
 R51: - _f[7,R3] + _r0[7,R3,Do_nothing] + _r0[7,R3,Recover] = 0
 R52: _Theta - 5069 _r1[9,Do_nothing,Do_nothing]
   - 5023 _r1[9,Do_nothing,Recover] - 5069 _r1[9,R1,Do_nothing]
   - 5022 _r1[9,R1,Recover] - 4475 _r1[9,R2,Do_nothing]
   - 3834 _r1[9,R2,Recover] - 4161 _r1[9,R3,Do_nothing]
   - 3208 _r1[9,R3,Recover] - 7616 _r1[2,Do_nothing,Do_nothing]
   - 7547 _r1[2,Do_nothing,Recover] - 7615 _r1[2,R1,Do_nothing]
   - 7546 _r1[2,R1,Recover] - 6723 _r1[2,R2,Do_nothing]
   - 5761 _r1[2,R2,Recover] - 6252 _r1[2,R3,Do_nothing]
   - 4819 _r1[2,R3,Recover] - 9778 _r1[8,Do_nothing,Do_nothing]
   - 9689 _r1[8,Do_nothing,Recover] - 9778 _r1[8,R1,Do_nothing]
   - 9688 _r1[8,R1,Recover] - 8632 _r1[8,R2,Do_nothing]
   - 7396 _r1[8,R2,Recover] - 8027 _r1[8,R3,Do_nothing]
   - 6188 _r1[8,R3,Recover] - 5364 _r1[3,Do_nothing,Do_nothing]
   - 5316 _r1[3,Do_nothing,Recover] - 5364 _r1[3,R1,Do_nothing]
   - 5315 _r1[3,R1,Recover] - 4735 _r1[3,R2,Do_nothing]
   - 4058 _r1[3,R2,Recover] - 4404 _r1[3,R3,Do_nothing]
   - 3395 _r1[3,R3,Recover] >= 0
 R53: 1.330991e+06 _f[9,R1] + 1.672772e+06 _f[9,R2] + 3.002586e+06 _f[9,R3]
   + 745011 _f[6,R1] + 945943 _f[6,R2] + 1.697949e+06 _f[6,R3]
   + 1.470866e+06 _f[2,R1] + 1.84857e+06 _f[2,R2] + 3.318145e+06 _f[2,R3]
   + 1.841752e+06 _f[10,R1] + 2.314703e+06 _f[10,R2]
   + 4.154833e+06 _f[10,R3] + 1.427401e+06 _f[5,R1] + 1.808916e+06 _f[5,R2]
   + 3.246972e+06 _f[5,R3] + 754217 _f[1,R1] + 988682 _f[1,R2]
   + 1.774661e+06 _f[1,R3] + 1.317138e+06 _f[4,R1] + 1.655368e+06 _f[4,R2]
   + 2.971353e+06 _f[4,R3] + 2.894596e+06 _f[8,R1] + 3.637895e+06 _f[8,R2]
   + 6.529949e+06 _f[8,R3] + 1.800986e+06 _f[3,R1] + 2.263469e+06 _f[3,R2]
   + 4.062874e+06 _f[3,R3] + 527129 _f[7,R1] + 662486 _f[7,R2]
   + 1.189149e+06 _f[7,R3] + 1.67294641e+08 _r1[9,Do_nothing,Recover]
   + 1.6724473e+08 _r1[9,R1,Recover] + 1.19991452e+08 _r1[9,R2,Recover]
   + 7.3054158e+07 _r1[9,R3,Recover]
   + 1.12287189e+08 _r1[6,Do_nothing,Recover]
   + 1.12253953e+08 _r1[6,R1,Recover] + 8.0694437e+07 _r1[6,R2,Recover]
   + 4.93764e+07 _r1[6,R3,Recover]
   + 2.34957414e+08 _r1[2,Do_nothing,Recover]
   + 2.34887316e+08 _r1[2,R1,Recover] + 1.68522322e+08 _r1[2,R2,Recover]
   + 1.02601111e+08 _r1[2,R3,Recover]
   + 2.68081439e+08 _r1[10,Do_nothing,Recover]
   + 2.68001459e+08 _r1[10,R1,Recover] + 1.92280404e+08 _r1[10,R2,Recover]
   + 1.17065697e+08 _r1[10,R3,Recover]
   + 3.10887702e+08 _r1[5,Do_nothing,Recover]
   + 3.10795373e+08 _r1[5,R1,Recover] + 2.23233935e+08 _r1[5,R2,Recover]
   + 1.36306708e+08 _r1[5,R3,Recover]
   + 1.77010253e+08 _r1[1,Do_nothing,Recover]
   + 1.76958131e+08 _r1[1,R1,Recover] + 1.27368662e+08 _r1[1,R2,Recover]
   + 7.819011e+07 _r1[1,R3,Recover]
   + 1.76857603e+08 _r1[4,Do_nothing,Recover]
   + 1.76804839e+08 _r1[4,R1,Recover] + 1.26850451e+08 _r1[4,R2,Recover]
   + 7.7230108e+07 _r1[4,R3,Recover]
   + 3.64132206e+08 _r1[8,Do_nothing,Recover]
   + 3.6402357e+08 _r1[8,R1,Recover] + 2.61172456e+08 _r1[8,R2,Recover]
   + 1.59009108e+08 _r1[8,R3,Recover]
   + 3.96428872e+08 _r1[3,Do_nothing,Recover]
   + 3.963106e+08 _r1[3,R1,Recover] + 2.84337119e+08 _r1[3,R2,Recover]
   + 1.73112404e+08 _r1[3,R3,Recover]
   + 1.46873379e+08 _r1[7,Do_nothing,Recover]
   + 1.46829561e+08 _r1[7,R1,Recover] + 1.05344379e+08 _r1[7,R2,Recover]
   + 6.4136609e+07 _r1[7,R3,Recover] <= 3e+07
 R54: - _f[9,Do_nothing] + _r1[9,Do_nothing,Do_nothing]
   + _r1[9,Do_nothing,Recover] = 0
 R55: - _f[9,R1] + _r1[9,R1,Do_nothing] + _r1[9,R1,Recover] = 0
 R56: - _f[9,R2] + _r1[9,R2,Do_nothing] + _r1[9,R2,Recover] = 0
 R57: - _f[9,R3] + _r1[9,R3,Do_nothing] + _r1[9,R3,Recover] = 0
 R58: - _f[6,Do_nothing] + _r1[6,Do_nothing,Do_nothing]
   + _r1[6,Do_nothing,Recover] = 0
 R59: - _f[6,R1] + _r1[6,R1,Do_nothing] + _r1[6,R1,Recover] = 0
 R60: - _f[6,R2] + _r1[6,R2,Do_nothing] + _r1[6,R2,Recover] = 0
 R61: - _f[6,R3] + _r1[6,R3,Do_nothing] + _r1[6,R3,Recover] = 0
 R62: - _f[2,Do_nothing] + _r1[2,Do_nothing,Do_nothing]
   + _r1[2,Do_nothing,Recover] = 0
 R63: - _f[2,R1] + _r1[2,R1,Do_nothing] + _r1[2,R1,Recover] = 0
 R64: - _f[2,R2] + _r1[2,R2,Do_nothing] + _r1[2,R2,Recover] = 0
 R65: - _f[2,R3] + _r1[2,R3,Do_nothing] + _r1[2,R3,Recover] = 0
 R66: - _f[10,Do_nothing] + _r1[10,Do_nothing,Do_nothing]
   + _r1[10,Do_nothing,Recover] = 0
 R67: - _f[10,R1] + _r1[10,R1,Do_nothing] + _r1[10,R1,Recover] = 0
 R68: - _f[10,R2] + _r1[10,R2,Do_nothing] + _r1[10,R2,Recover] = 0
 R69: - _f[10,R3] + _r1[10,R3,Do_nothing] + _r1[10,R3,Recover] = 0
 R70: - _f[5,Do_nothing] + _r1[5,Do_nothing,Do_nothing]
   + _r1[5,Do_nothing,Recover] = 0
 R71: - _f[5,R1] + _r1[5,R1,Do_nothing] + _r1[5,R1,Recover] = 0
 R72: - _f[5,R2] + _r1[5,R2,Do_nothing] + _r1[5,R2,Recover] = 0
 R73: - _f[5,R3] + _r1[5,R3,Do_nothing] + _r1[5,R3,Recover] = 0
 R74: - _f[1,Do_nothing] + _r1[1,Do_nothing,Do_nothing]
   + _r1[1,Do_nothing,Recover] = 0
 R75: - _f[1,R1] + _r1[1,R1,Do_nothing] + _r1[1,R1,Recover] = 0
 R76: - _f[1,R2] + _r1[1,R2,Do_nothing] + _r1[1,R2,Recover] = 0
 R77: - _f[1,R3] + _r1[1,R3,Do_nothing] + _r1[1,R3,Recover] = 0
 R78: - _f[4,Do_nothing] + _r1[4,Do_nothing,Do_nothing]
   + _r1[4,Do_nothing,Recover] = 0
 R79: - _f[4,R1] + _r1[4,R1,Do_nothing] + _r1[4,R1,Recover] = 0
 R80: - _f[4,R2] + _r1[4,R2,Do_nothing] + _r1[4,R2,Recover] = 0
 R81: - _f[4,R3] + _r1[4,R3,Do_nothing] + _r1[4,R3,Recover] = 0
 R82: - _f[8,Do_nothing] + _r1[8,Do_nothing,Do_nothing]
   + _r1[8,Do_nothing,Recover] = 0
 R83: - _f[8,R1] + _r1[8,R1,Do_nothing] + _r1[8,R1,Recover] = 0
 R84: - _f[8,R2] + _r1[8,R2,Do_nothing] + _r1[8,R2,Recover] = 0
 R85: - _f[8,R3] + _r1[8,R3,Do_nothing] + _r1[8,R3,Recover] = 0
 R86: - _f[3,Do_nothing] + _r1[3,Do_nothing,Do_nothing]
   + _r1[3,Do_nothing,Recover] = 0
 R87: - _f[3,R1] + _r1[3,R1,Do_nothing] + _r1[3,R1,Recover] = 0
 R88: - _f[3,R2] + _r1[3,R2,Do_nothing] + _r1[3,R2,Recover] = 0
 R89: - _f[3,R3] + _r1[3,R3,Do_nothing] + _r1[3,R3,Recover] = 0
 R90: - _f[7,Do_nothing] + _r1[7,Do_nothing,Do_nothing]
   + _r1[7,Do_nothing,Recover] = 0
 R91: - _f[7,R1] + _r1[7,R1,Do_nothing] + _r1[7,R1,Recover] = 0
 R92: - _f[7,R2] + _r1[7,R2,Do_nothing] + _r1[7,R2,Recover] = 0
 R93: - _f[7,R3] + _r1[7,R3,Do_nothing] + _r1[7,R3,Recover] = 0
Bounds
 _Theta free
Binaries
 _f[9,Do_nothing] _f[9,R1] _f[9,R2] _f[9,R3] _f[6,Do_nothing] _f[6,R1]
 _f[6,R2] _f[6,R3] _f[2,Do_nothing] _f[2,R1] _f[2,R2] _f[2,R3]
 _f[10,Do_nothing] _f[10,R1] _f[10,R2] _f[10,R3] _f[5,Do_nothing] _f[5,R1]
 _f[5,R2] _f[5,R3] _f[1,Do_nothing] _f[1,R1] _f[1,R2] _f[1,R3]
 _f[4,Do_nothing] _f[4,R1] _f[4,R2] _f[4,R3] _f[8,Do_nothing] _f[8,R1]
 _f[8,R2] _f[8,R3] _f[3,Do_nothing] _f[3,R1] _f[3,R2] _f[3,R3]
 _f[7,Do_nothing] _f[7,R1] _f[7,R2] _f[7,R3] _r0[9,Do_nothing,Do_nothing]
 _r0[9,Do_nothing,Recover] _r0[9,R1,Do_nothing] _r0[9,R1,Recover]
 _r0[9,R2,Do_nothing] _r0[9,R2,Recover] _r0[9,R3,Do_nothing]
 _r0[9,R3,Recover] _r0[6,Do_nothing,Do_nothing] _r0[6,Do_nothing,Recover]
 _r0[6,R1,Do_nothing] _r0[6,R1,Recover] _r0[6,R2,Do_nothing]
 _r0[6,R2,Recover] _r0[6,R3,Do_nothing] _r0[6,R3,Recover]
 _r0[2,Do_nothing,Do_nothing] _r0[2,Do_nothing,Recover]
 _r0[2,R1,Do_nothing] _r0[2,R1,Recover] _r0[2,R2,Do_nothing]
 _r0[2,R2,Recover] _r0[2,R3,Do_nothing] _r0[2,R3,Recover]
 _r0[10,Do_nothing,Do_nothing] _r0[10,Do_nothing,Recover]
 _r0[10,R1,Do_nothing] _r0[10,R1,Recover] _r0[10,R2,Do_nothing]
 _r0[10,R2,Recover] _r0[10,R3,Do_nothing] _r0[10,R3,Recover]
 _r0[5,Do_nothing,Do_nothing] _r0[5,Do_nothing,Recover]
 _r0[5,R1,Do_nothing] _r0[5,R1,Recover] _r0[5,R2,Do_nothing]
 _r0[5,R2,Recover] _r0[5,R3,Do_nothing] _r0[5,R3,Recover]
 _r0[1,Do_nothing,Do_nothing] _r0[1,Do_nothing,Recover]
 _r0[1,R1,Do_nothing] _r0[1,R1,Recover] _r0[1,R2,Do_nothing]
 _r0[1,R2,Recover] _r0[1,R3,Do_nothing] _r0[1,R3,Recover]
 _r0[4,Do_nothing,Do_nothing] _r0[4,Do_nothing,Recover]
 _r0[4,R1,Do_nothing] _r0[4,R1,Recover] _r0[4,R2,Do_nothing]
 _r0[4,R2,Recover] _r0[4,R3,Do_nothing] _r0[4,R3,Recover]
 _r0[8,Do_nothing,Do_nothing] _r0[8,Do_nothing,Recover]
 _r0[8,R1,Do_nothing] _r0[8,R1,Recover] _r0[8,R2,Do_nothing]
 _r0[8,R2,Recover] _r0[8,R3,Do_nothing] _r0[8,R3,Recover]
 _r0[3,Do_nothing,Do_nothing] _r0[3,Do_nothing,Recover]
 _r0[3,R1,Do_nothing] _r0[3,R1,Recover] _r0[3,R2,Do_nothing]
 _r0[3,R2,Recover] _r0[3,R3,Do_nothing] _r0[3,R3,Recover]
 _r0[7,Do_nothing,Do_nothing] _r0[7,Do_nothing,Recover]
 _r0[7,R1,Do_nothing] _r0[7,R1,Recover] _r0[7,R2,Do_nothing]
 _r0[7,R2,Recover] _r0[7,R3,Do_nothing] _r0[7,R3,Recover]
 _r1[9,Do_nothing,Do_nothing] _r1[9,Do_nothing,Recover]
 _r1[9,R1,Do_nothing] _r1[9,R1,Recover] _r1[9,R2,Do_nothing]
 _r1[9,R2,Recover] _r1[9,R3,Do_nothing] _r1[9,R3,Recover]
 _r1[6,Do_nothing,Do_nothing] _r1[6,Do_nothing,Recover]
 _r1[6,R1,Do_nothing] _r1[6,R1,Recover] _r1[6,R2,Do_nothing]
 _r1[6,R2,Recover] _r1[6,R3,Do_nothing] _r1[6,R3,Recover]
 _r1[2,Do_nothing,Do_nothing] _r1[2,Do_nothing,Recover]
 _r1[2,R1,Do_nothing] _r1[2,R1,Recover] _r1[2,R2,Do_nothing]
 _r1[2,R2,Recover] _r1[2,R3,Do_nothing] _r1[2,R3,Recover]
 _r1[10,Do_nothing,Do_nothing] _r1[10,Do_nothing,Recover]
 _r1[10,R1,Do_nothing] _r1[10,R1,Recover] _r1[10,R2,Do_nothing]
 _r1[10,R2,Recover] _r1[10,R3,Do_nothing] _r1[10,R3,Recover]
 _r1[5,Do_nothing,Do_nothing] _r1[5,Do_nothing,Recover]
 _r1[5,R1,Do_nothing] _r1[5,R1,Recover] _r1[5,R2,Do_nothing]
 _r1[5,R2,Recover] _r1[5,R3,Do_nothing] _r1[5,R3,Recover]
 _r1[1,Do_nothing,Do_nothing] _r1[1,Do_nothing,Recover]
 _r1[1,R1,Do_nothing] _r1[1,R1,Recover] _r1[1,R2,Do_nothing]
 _r1[1,R2,Recover] _r1[1,R3,Do_nothing] _r1[1,R3,Recover]
 _r1[4,Do_nothing,Do_nothing] _r1[4,Do_nothing,Recover]
 _r1[4,R1,Do_nothing] _r1[4,R1,Recover] _r1[4,R2,Do_nothing]
 _r1[4,R2,Recover] _r1[4,R3,Do_nothing] _r1[4,R3,Recover]
 _r1[8,Do_nothing,Do_nothing] _r1[8,Do_nothing,Recover]
 _r1[8,R1,Do_nothing] _r1[8,R1,Recover] _r1[8,R2,Do_nothing]
 _r1[8,R2,Recover] _r1[8,R3,Do_nothing] _r1[8,R3,Recover]
 _r1[3,Do_nothing,Do_nothing] _r1[3,Do_nothing,Recover]
 _r1[3,R1,Do_nothing] _r1[3,R1,Recover] _r1[3,R2,Do_nothing]
 _r1[3,R2,Recover] _r1[3,R3,Do_nothing] _r1[3,R3,Recover]
 _r1[7,Do_nothing,Do_nothing] _r1[7,Do_nothing,Recover]
 _r1[7,R1,Do_nothing] _r1[7,R1,Recover] _r1[7,R2,Do_nothing]
 _r1[7,R2,Recover] _r1[7,R3,Do_nothing] _r1[7,R3,Recover]
End
