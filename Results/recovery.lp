\ Model recoveryProblem_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5069 _r[9,Do_nothing,Do_nothing] + 5023 _r[9,Do_nothing,Recover]
   + 5069 _r[9,R1,Do_nothing] + 5022 _r[9,R1,Recover]
   + 4475 _r[9,R2,Do_nothing] + 3834 _r[9,R2,Recover]
   + 4161 _r[9,R3,Do_nothing] + 3208 _r[9,R3,Recover]
   + 7616 _r[2,Do_nothing,Do_nothing] + 7547 _r[2,Do_nothing,Recover]
   + 7615 _r[2,R1,Do_nothing] + 7546 _r[2,R1,Recover]
   + 6723 _r[2,R2,Do_nothing] + 5761 _r[2,R2,Recover]
   + 6252 _r[2,R3,Do_nothing] + 4819 _r[2,R3,Recover]
   + 9778 _r[8,Do_nothing,Do_nothing] + 9689 _r[8,Do_nothing,Recover]
   + 9778 _r[8,R1,Do_nothing] + 9688 _r[8,R1,Recover]
   + 8632 _r[8,R2,Do_nothing] + 7396 _r[8,R2,Recover]
   + 8027 _r[8,R3,Do_nothing] + 6188 _r[8,R3,Recover]
   + 5364 _r[3,Do_nothing,Do_nothing] + 5316 _r[3,Do_nothing,Recover]
   + 5364 _r[3,R1,Do_nothing] + 5315 _r[3,R1,Recover]
   + 4735 _r[3,R2,Do_nothing] + 4058 _r[3,R2,Recover]
   + 4404 _r[3,R3,Do_nothing] + 3395 _r[3,R3,Recover]
Subject To
 R0: 1.67294641e+08 _r[9,Do_nothing,Recover]
   + 1.6724473e+08 _r[9,R1,Recover] + 1.19991452e+08 _r[9,R2,Recover]
   + 7.3054158e+07 _r[9,R3,Recover]
   + 1.12287189e+08 _r[6,Do_nothing,Recover]
   + 1.12253953e+08 _r[6,R1,Recover] + 8.0694437e+07 _r[6,R2,Recover]
   + 4.93764e+07 _r[6,R3,Recover] + 2.34957414e+08 _r[2,Do_nothing,Recover]
   + 2.34887316e+08 _r[2,R1,Recover] + 1.68522322e+08 _r[2,R2,Recover]
   + 1.02601111e+08 _r[2,R3,Recover]
   + 2.68081439e+08 _r[10,Do_nothing,Recover]
   + 2.68001459e+08 _r[10,R1,Recover] + 1.92280404e+08 _r[10,R2,Recover]
   + 1.17065697e+08 _r[10,R3,Recover]
   + 3.10887702e+08 _r[5,Do_nothing,Recover]
   + 3.10795373e+08 _r[5,R1,Recover] + 2.23233935e+08 _r[5,R2,Recover]
   + 1.36306708e+08 _r[5,R3,Recover]
   + 1.77010253e+08 _r[1,Do_nothing,Recover]
   + 1.76958131e+08 _r[1,R1,Recover] + 1.27368662e+08 _r[1,R2,Recover]
   + 7.819011e+07 _r[1,R3,Recover]
   + 1.76857603e+08 _r[4,Do_nothing,Recover]
   + 1.76804839e+08 _r[4,R1,Recover] + 1.26850451e+08 _r[4,R2,Recover]
   + 7.7230108e+07 _r[4,R3,Recover]
   + 3.64132206e+08 _r[8,Do_nothing,Recover]
   + 3.6402357e+08 _r[8,R1,Recover] + 2.61172456e+08 _r[8,R2,Recover]
   + 1.59009108e+08 _r[8,R3,Recover]
   + 3.96428872e+08 _r[3,Do_nothing,Recover]
   + 3.963106e+08 _r[3,R1,Recover] + 2.84337119e+08 _r[3,R2,Recover]
   + 1.73112404e+08 _r[3,R3,Recover]
   + 1.46873379e+08 _r[7,Do_nothing,Recover]
   + 1.46829561e+08 _r[7,R1,Recover] + 1.05344379e+08 _r[7,R2,Recover]
   + 6.4136609e+07 _r[7,R3,Recover] <= 1.3086446e+07
 R1: _r[9,Do_nothing,Do_nothing] + _r[9,Do_nothing,Recover] = 0
 R2: _r[9,R1,Do_nothing] + _r[9,R1,Recover] = 0
 R3: _r[9,R2,Do_nothing] + _r[9,R2,Recover] = 0
 R4: _r[9,R3,Do_nothing] + _r[9,R3,Recover] = 1
 R5: _r[6,Do_nothing,Do_nothing] + _r[6,Do_nothing,Recover] = 1
 R6: _r[6,R1,Do_nothing] + _r[6,R1,Recover] = 0
 R7: _r[6,R2,Do_nothing] + _r[6,R2,Recover] = 0
 R8: _r[6,R3,Do_nothing] + _r[6,R3,Recover] = 0
 R9: _r[2,Do_nothing,Do_nothing] + _r[2,Do_nothing,Recover] = 0
 R10: _r[2,R1,Do_nothing] + _r[2,R1,Recover] = 0
 R11: _r[2,R2,Do_nothing] + _r[2,R2,Recover] = 0
 R12: _r[2,R3,Do_nothing] + _r[2,R3,Recover] = 1
 R13: _r[10,Do_nothing,Do_nothing] + _r[10,Do_nothing,Recover] = 1
 R14: _r[10,R1,Do_nothing] + _r[10,R1,Recover] = 0
 R15: _r[10,R2,Do_nothing] + _r[10,R2,Recover] = 0
 R16: _r[10,R3,Do_nothing] + _r[10,R3,Recover] = 0
 R17: _r[5,Do_nothing,Do_nothing] + _r[5,Do_nothing,Recover] = 1
 R18: _r[5,R1,Do_nothing] + _r[5,R1,Recover] = 0
 R19: _r[5,R2,Do_nothing] + _r[5,R2,Recover] = 0
 R20: _r[5,R3,Do_nothing] + _r[5,R3,Recover] = 0
 R21: _r[1,Do_nothing,Do_nothing] + _r[1,Do_nothing,Recover] = 1
 R22: _r[1,R1,Do_nothing] + _r[1,R1,Recover] = 0
 R23: _r[1,R2,Do_nothing] + _r[1,R2,Recover] = 0
 R24: _r[1,R3,Do_nothing] + _r[1,R3,Recover] = 0
 R25: _r[4,Do_nothing,Do_nothing] + _r[4,Do_nothing,Recover] = 1
 R26: _r[4,R1,Do_nothing] + _r[4,R1,Recover] = 0
 R27: _r[4,R2,Do_nothing] + _r[4,R2,Recover] = 0
 R28: _r[4,R3,Do_nothing] + _r[4,R3,Recover] = 0
 R29: _r[8,Do_nothing,Do_nothing] + _r[8,Do_nothing,Recover] = 0
 R30: _r[8,R1,Do_nothing] + _r[8,R1,Recover] = 0
 R31: _r[8,R2,Do_nothing] + _r[8,R2,Recover] = 0
 R32: _r[8,R3,Do_nothing] + _r[8,R3,Recover] = 1
 R33: _r[3,Do_nothing,Do_nothing] + _r[3,Do_nothing,Recover] = 0
 R34: _r[3,R1,Do_nothing] + _r[3,R1,Recover] = 0
 R35: _r[3,R2,Do_nothing] + _r[3,R2,Recover] = 0
 R36: _r[3,R3,Do_nothing] + _r[3,R3,Recover] = 1
 R37: _r[7,Do_nothing,Do_nothing] + _r[7,Do_nothing,Recover] = 1
 R38: _r[7,R1,Do_nothing] + _r[7,R1,Recover] = 0
 R39: _r[7,R2,Do_nothing] + _r[7,R2,Recover] = 0
 R40: _r[7,R3,Do_nothing] + _r[7,R3,Recover] = 0
Bounds
Binaries
 _r[9,Do_nothing,Do_nothing] _r[9,Do_nothing,Recover] _r[9,R1,Do_nothing]
 _r[9,R1,Recover] _r[9,R2,Do_nothing] _r[9,R2,Recover] _r[9,R3,Do_nothing]
 _r[9,R3,Recover] _r[6,Do_nothing,Do_nothing] _r[6,Do_nothing,Recover]
 _r[6,R1,Do_nothing] _r[6,R1,Recover] _r[6,R2,Do_nothing] _r[6,R2,Recover]
 _r[6,R3,Do_nothing] _r[6,R3,Recover] _r[2,Do_nothing,Do_nothing]
 _r[2,Do_nothing,Recover] _r[2,R1,Do_nothing] _r[2,R1,Recover]
 _r[2,R2,Do_nothing] _r[2,R2,Recover] _r[2,R3,Do_nothing] _r[2,R3,Recover]
 _r[10,Do_nothing,Do_nothing] _r[10,Do_nothing,Recover]
 _r[10,R1,Do_nothing] _r[10,R1,Recover] _r[10,R2,Do_nothing]
 _r[10,R2,Recover] _r[10,R3,Do_nothing] _r[10,R3,Recover]
 _r[5,Do_nothing,Do_nothing] _r[5,Do_nothing,Recover] _r[5,R1,Do_nothing]
 _r[5,R1,Recover] _r[5,R2,Do_nothing] _r[5,R2,Recover] _r[5,R3,Do_nothing]
 _r[5,R3,Recover] _r[1,Do_nothing,Do_nothing] _r[1,Do_nothing,Recover]
 _r[1,R1,Do_nothing] _r[1,R1,Recover] _r[1,R2,Do_nothing] _r[1,R2,Recover]
 _r[1,R3,Do_nothing] _r[1,R3,Recover] _r[4,Do_nothing,Do_nothing]
 _r[4,Do_nothing,Recover] _r[4,R1,Do_nothing] _r[4,R1,Recover]
 _r[4,R2,Do_nothing] _r[4,R2,Recover] _r[4,R3,Do_nothing] _r[4,R3,Recover]
 _r[8,Do_nothing,Do_nothing] _r[8,Do_nothing,Recover] _r[8,R1,Do_nothing]
 _r[8,R1,Recover] _r[8,R2,Do_nothing] _r[8,R2,Recover] _r[8,R3,Do_nothing]
 _r[8,R3,Recover] _r[3,Do_nothing,Do_nothing] _r[3,Do_nothing,Recover]
 _r[3,R1,Do_nothing] _r[3,R1,Recover] _r[3,R2,Do_nothing] _r[3,R2,Recover]
 _r[3,R3,Do_nothing] _r[3,R3,Recover] _r[7,Do_nothing,Do_nothing]
 _r[7,Do_nothing,Recover] _r[7,R1,Do_nothing] _r[7,R1,Recover]
 _r[7,R2,Do_nothing] _r[7,R2,Recover] _r[7,R3,Do_nothing] _r[7,R3,Recover]
End
