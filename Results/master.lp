\ Model master_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  _Theta
Subject To
 R0: _f[32,Do_nothing] + _f[32,R1] + _f[32,R2] + _f[32,R3] = 1
 R1: _f[70,Do_nothing] + _f[70,R1] + _f[70,R2] + _f[70,R3] = 1
 R2: _f[5,Do_nothing] + _f[5,R1] + _f[5,R2] + _f[5,R3] = 1
 R3: _f[85,Do_nothing] + _f[85,R1] + _f[85,R2] + _f[85,R3] = 1
 R4: _f[71,Do_nothing] + _f[71,R1] + _f[71,R2] + _f[71,R3] = 1
 R5: _f[92,Do_nothing] + _f[92,R1] + _f[92,R2] + _f[92,R3] = 1
 R6: _f[38,Do_nothing] + _f[38,R1] + _f[38,R2] + _f[38,R3] = 1
 R7: _f[21,Do_nothing] + _f[21,R1] + _f[21,R2] + _f[21,R3] = 1
 R8: _f[41,Do_nothing] + _f[41,R1] + _f[41,R2] + _f[41,R3] = 1
 R9: _f[56,Do_nothing] + _f[56,R1] + _f[56,R2] + _f[56,R3] = 1
 R10: _f[80,Do_nothing] + _f[80,R1] + _f[80,R2] + _f[80,R3] = 1
 R11: _f[62,Do_nothing] + _f[62,R1] + _f[62,R2] + _f[62,R3] = 1
 R12: _f[20,Do_nothing] + _f[20,R1] + _f[20,R2] + _f[20,R3] = 1
 R13: _f[10,Do_nothing] + _f[10,R1] + _f[10,R2] + _f[10,R3] = 1
 R14: _f[55,Do_nothing] + _f[55,R1] + _f[55,R2] + _f[55,R3] = 1
 R15: _f[36,Do_nothing] + _f[36,R1] + _f[36,R2] + _f[36,R3] = 1
 R16: _f[44,Do_nothing] + _f[44,R1] + _f[44,R2] + _f[44,R3] = 1
 R17: _f[48,Do_nothing] + _f[48,R1] + _f[48,R2] + _f[48,R3] = 1
 R18: _f[8,Do_nothing] + _f[8,R1] + _f[8,R2] + _f[8,R3] = 1
 R19: _f[15,Do_nothing] + _f[15,R1] + _f[15,R2] + _f[15,R3] = 1
 R20: _f[61,Do_nothing] + _f[61,R1] + _f[61,R2] + _f[61,R3] = 1
 R21: _f[78,Do_nothing] + _f[78,R1] + _f[78,R2] + _f[78,R3] = 1
 R22: _f[34,Do_nothing] + _f[34,R1] + _f[34,R2] + _f[34,R3] = 1
 R23: _f[100,Do_nothing] + _f[100,R1] + _f[100,R2] + _f[100,R3] = 1
 R24: _f[90,Do_nothing] + _f[90,R1] + _f[90,R2] + _f[90,R3] = 1
 R25: _f[64,Do_nothing] + _f[64,R1] + _f[64,R2] + _f[64,R3] = 1
 R26: _f[79,Do_nothing] + _f[79,R1] + _f[79,R2] + _f[79,R3] = 1
 R27: _f[1,Do_nothing] + _f[1,R1] + _f[1,R2] + _f[1,R3] = 1
 R28: _f[22,Do_nothing] + _f[22,R1] + _f[22,R2] + _f[22,R3] = 1
 R29: _f[31,Do_nothing] + _f[31,R1] + _f[31,R2] + _f[31,R3] = 1
 R30: _f[54,Do_nothing] + _f[54,R1] + _f[54,R2] + _f[54,R3] = 1
 R31: _f[82,Do_nothing] + _f[82,R1] + _f[82,R2] + _f[82,R3] = 1
 R32: _f[11,Do_nothing] + _f[11,R1] + _f[11,R2] + _f[11,R3] = 1
 R33: _f[75,Do_nothing] + _f[75,R1] + _f[75,R2] + _f[75,R3] = 1
 R34: _f[50,Do_nothing] + _f[50,R1] + _f[50,R2] + _f[50,R3] = 1
 R35: _f[19,Do_nothing] + _f[19,R1] + _f[19,R2] + _f[19,R3] = 1
 R36: _f[99,Do_nothing] + _f[99,R1] + _f[99,R2] + _f[99,R3] = 1
 R37: _f[72,Do_nothing] + _f[72,R1] + _f[72,R2] + _f[72,R3] = 1
 R38: _f[47,Do_nothing] + _f[47,R1] + _f[47,R2] + _f[47,R3] = 1
 R39: _f[49,Do_nothing] + _f[49,R1] + _f[49,R2] + _f[49,R3] = 1
 R40: _f[77,Do_nothing] + _f[77,R1] + _f[77,R2] + _f[77,R3] = 1
 R41: _f[95,Do_nothing] + _f[95,R1] + _f[95,R2] + _f[95,R3] = 1
 R42: _f[98,Do_nothing] + _f[98,R1] + _f[98,R2] + _f[98,R3] = 1
 R43: _f[53,Do_nothing] + _f[53,R1] + _f[53,R2] + _f[53,R3] = 1
 R44: _f[65,Do_nothing] + _f[65,R1] + _f[65,R2] + _f[65,R3] = 1
 R45: _f[63,Do_nothing] + _f[63,R1] + _f[63,R2] + _f[63,R3] = 1
 R46: _f[14,Do_nothing] + _f[14,R1] + _f[14,R2] + _f[14,R3] = 1
 R47: _f[23,Do_nothing] + _f[23,R1] + _f[23,R2] + _f[23,R3] = 1
 R48: _f[93,Do_nothing] + _f[93,R1] + _f[93,R2] + _f[93,R3] = 1
 R49: _f[28,Do_nothing] + _f[28,R1] + _f[28,R2] + _f[28,R3] = 1
 R50: _f[2,Do_nothing] + _f[2,R1] + _f[2,R2] + _f[2,R3] = 1
 R51: _f[81,Do_nothing] + _f[81,R1] + _f[81,R2] + _f[81,R3] = 1
 R52: _f[45,Do_nothing] + _f[45,R1] + _f[45,R2] + _f[45,R3] = 1
 R53: _f[88,Do_nothing] + _f[88,R1] + _f[88,R2] + _f[88,R3] = 1
 R54: _f[89,Do_nothing] + _f[89,R1] + _f[89,R2] + _f[89,R3] = 1
 R55: _f[25,Do_nothing] + _f[25,R1] + _f[25,R2] + _f[25,R3] = 1
 R56: _f[91,Do_nothing] + _f[91,R1] + _f[91,R2] + _f[91,R3] = 1
 R57: _f[43,Do_nothing] + _f[43,R1] + _f[43,R2] + _f[43,R3] = 1
 R58: _f[76,Do_nothing] + _f[76,R1] + _f[76,R2] + _f[76,R3] = 1
 R59: _f[42,Do_nothing] + _f[42,R1] + _f[42,R2] + _f[42,R3] = 1
 R60: _f[7,Do_nothing] + _f[7,R1] + _f[7,R2] + _f[7,R3] = 1
 R61: _f[17,Do_nothing] + _f[17,R1] + _f[17,R2] + _f[17,R3] = 1
 R62: _f[37,Do_nothing] + _f[37,R1] + _f[37,R2] + _f[37,R3] = 1
 R63: _f[26,Do_nothing] + _f[26,R1] + _f[26,R2] + _f[26,R3] = 1
 R64: _f[9,Do_nothing] + _f[9,R1] + _f[9,R2] + _f[9,R3] = 1
 R65: _f[60,Do_nothing] + _f[60,R1] + _f[60,R2] + _f[60,R3] = 1
 R66: _f[16,Do_nothing] + _f[16,R1] + _f[16,R2] + _f[16,R3] = 1
 R67: _f[83,Do_nothing] + _f[83,R1] + _f[83,R2] + _f[83,R3] = 1
 R68: _f[13,Do_nothing] + _f[13,R1] + _f[13,R2] + _f[13,R3] = 1
 R69: _f[67,Do_nothing] + _f[67,R1] + _f[67,R2] + _f[67,R3] = 1
 R70: _f[59,Do_nothing] + _f[59,R1] + _f[59,R2] + _f[59,R3] = 1
 R71: _f[58,Do_nothing] + _f[58,R1] + _f[58,R2] + _f[58,R3] = 1
 R72: _f[4,Do_nothing] + _f[4,R1] + _f[4,R2] + _f[4,R3] = 1
 R73: _f[68,Do_nothing] + _f[68,R1] + _f[68,R2] + _f[68,R3] = 1
 R74: _f[87,Do_nothing] + _f[87,R1] + _f[87,R2] + _f[87,R3] = 1
 R75: _f[40,Do_nothing] + _f[40,R1] + _f[40,R2] + _f[40,R3] = 1
 R76: _f[30,Do_nothing] + _f[30,R1] + _f[30,R2] + _f[30,R3] = 1
 R77: _f[27,Do_nothing] + _f[27,R1] + _f[27,R2] + _f[27,R3] = 1
 R78: _f[73,Do_nothing] + _f[73,R1] + _f[73,R2] + _f[73,R3] = 1
 R79: _f[84,Do_nothing] + _f[84,R1] + _f[84,R2] + _f[84,R3] = 1
 R80: _f[35,Do_nothing] + _f[35,R1] + _f[35,R2] + _f[35,R3] = 1
 R81: _f[46,Do_nothing] + _f[46,R1] + _f[46,R2] + _f[46,R3] = 1
 R82: _f[39,Do_nothing] + _f[39,R1] + _f[39,R2] + _f[39,R3] = 1
 R83: _f[51,Do_nothing] + _f[51,R1] + _f[51,R2] + _f[51,R3] = 1
 R84: _f[52,Do_nothing] + _f[52,R1] + _f[52,R2] + _f[52,R3] = 1
 R85: _f[6,Do_nothing] + _f[6,R1] + _f[6,R2] + _f[6,R3] = 1
 R86: _f[96,Do_nothing] + _f[96,R1] + _f[96,R2] + _f[96,R3] = 1
 R87: _f[74,Do_nothing] + _f[74,R1] + _f[74,R2] + _f[74,R3] = 1
 R88: _f[86,Do_nothing] + _f[86,R1] + _f[86,R2] + _f[86,R3] = 1
 R89: _f[29,Do_nothing] + _f[29,R1] + _f[29,R2] + _f[29,R3] = 1
 R90: _f[18,Do_nothing] + _f[18,R1] + _f[18,R2] + _f[18,R3] = 1
 R91: _f[24,Do_nothing] + _f[24,R1] + _f[24,R2] + _f[24,R3] = 1
 R92: _f[94,Do_nothing] + _f[94,R1] + _f[94,R2] + _f[94,R3] = 1
 R93: _f[33,Do_nothing] + _f[33,R1] + _f[33,R2] + _f[33,R3] = 1
 R94: _f[12,Do_nothing] + _f[12,R1] + _f[12,R2] + _f[12,R3] = 1
 R95: _f[57,Do_nothing] + _f[57,R1] + _f[57,R2] + _f[57,R3] = 1
 R96: _f[97,Do_nothing] + _f[97,R1] + _f[97,R2] + _f[97,R3] = 1
 R97: _f[69,Do_nothing] + _f[69,R1] + _f[69,R2] + _f[69,R3] = 1
 R98: _f[3,Do_nothing] + _f[3,R1] + _f[3,R2] + _f[3,R3] = 1
 R99: _f[66,Do_nothing] + _f[66,R1] + _f[66,R2] + _f[66,R3] = 1
Bounds
 _Theta free
Binaries
 _f[32,Do_nothing] _f[32,R1] _f[32,R2] _f[32,R3] _f[70,Do_nothing]
 _f[70,R1] _f[70,R2] _f[70,R3] _f[5,Do_nothing] _f[5,R1] _f[5,R2] _f[5,R3]
 _f[85,Do_nothing] _f[85,R1] _f[85,R2] _f[85,R3] _f[71,Do_nothing]
 _f[71,R1] _f[71,R2] _f[71,R3] _f[92,Do_nothing] _f[92,R1] _f[92,R2]
 _f[92,R3] _f[38,Do_nothing] _f[38,R1] _f[38,R2] _f[38,R3]
 _f[21,Do_nothing] _f[21,R1] _f[21,R2] _f[21,R3] _f[41,Do_nothing]
 _f[41,R1] _f[41,R2] _f[41,R3] _f[56,Do_nothing] _f[56,R1] _f[56,R2]
 _f[56,R3] _f[80,Do_nothing] _f[80,R1] _f[80,R2] _f[80,R3]
 _f[62,Do_nothing] _f[62,R1] _f[62,R2] _f[62,R3] _f[20,Do_nothing]
 _f[20,R1] _f[20,R2] _f[20,R3] _f[10,Do_nothing] _f[10,R1] _f[10,R2]
 _f[10,R3] _f[55,Do_nothing] _f[55,R1] _f[55,R2] _f[55,R3]
 _f[36,Do_nothing] _f[36,R1] _f[36,R2] _f[36,R3] _f[44,Do_nothing]
 _f[44,R1] _f[44,R2] _f[44,R3] _f[48,Do_nothing] _f[48,R1] _f[48,R2]
 _f[48,R3] _f[8,Do_nothing] _f[8,R1] _f[8,R2] _f[8,R3] _f[15,Do_nothing]
 _f[15,R1] _f[15,R2] _f[15,R3] _f[61,Do_nothing] _f[61,R1] _f[61,R2]
 _f[61,R3] _f[78,Do_nothing] _f[78,R1] _f[78,R2] _f[78,R3]
 _f[34,Do_nothing] _f[34,R1] _f[34,R2] _f[34,R3] _f[100,Do_nothing]
 _f[100,R1] _f[100,R2] _f[100,R3] _f[90,Do_nothing] _f[90,R1] _f[90,R2]
 _f[90,R3] _f[64,Do_nothing] _f[64,R1] _f[64,R2] _f[64,R3]
 _f[79,Do_nothing] _f[79,R1] _f[79,R2] _f[79,R3] _f[1,Do_nothing] _f[1,R1]
 _f[1,R2] _f[1,R3] _f[22,Do_nothing] _f[22,R1] _f[22,R2] _f[22,R3]
 _f[31,Do_nothing] _f[31,R1] _f[31,R2] _f[31,R3] _f[54,Do_nothing]
 _f[54,R1] _f[54,R2] _f[54,R3] _f[82,Do_nothing] _f[82,R1] _f[82,R2]
 _f[82,R3] _f[11,Do_nothing] _f[11,R1] _f[11,R2] _f[11,R3]
 _f[75,Do_nothing] _f[75,R1] _f[75,R2] _f[75,R3] _f[50,Do_nothing]
 _f[50,R1] _f[50,R2] _f[50,R3] _f[19,Do_nothing] _f[19,R1] _f[19,R2]
 _f[19,R3] _f[99,Do_nothing] _f[99,R1] _f[99,R2] _f[99,R3]
 _f[72,Do_nothing] _f[72,R1] _f[72,R2] _f[72,R3] _f[47,Do_nothing]
 _f[47,R1] _f[47,R2] _f[47,R3] _f[49,Do_nothing] _f[49,R1] _f[49,R2]
 _f[49,R3] _f[77,Do_nothing] _f[77,R1] _f[77,R2] _f[77,R3]
 _f[95,Do_nothing] _f[95,R1] _f[95,R2] _f[95,R3] _f[98,Do_nothing]
 _f[98,R1] _f[98,R2] _f[98,R3] _f[53,Do_nothing] _f[53,R1] _f[53,R2]
 _f[53,R3] _f[65,Do_nothing] _f[65,R1] _f[65,R2] _f[65,R3]
 _f[63,Do_nothing] _f[63,R1] _f[63,R2] _f[63,R3] _f[14,Do_nothing]
 _f[14,R1] _f[14,R2] _f[14,R3] _f[23,Do_nothing] _f[23,R1] _f[23,R2]
 _f[23,R3] _f[93,Do_nothing] _f[93,R1] _f[93,R2] _f[93,R3]
 _f[28,Do_nothing] _f[28,R1] _f[28,R2] _f[28,R3] _f[2,Do_nothing] _f[2,R1]
 _f[2,R2] _f[2,R3] _f[81,Do_nothing] _f[81,R1] _f[81,R2] _f[81,R3]
 _f[45,Do_nothing] _f[45,R1] _f[45,R2] _f[45,R3] _f[88,Do_nothing]
 _f[88,R1] _f[88,R2] _f[88,R3] _f[89,Do_nothing] _f[89,R1] _f[89,R2]
 _f[89,R3] _f[25,Do_nothing] _f[25,R1] _f[25,R2] _f[25,R3]
 _f[91,Do_nothing] _f[91,R1] _f[91,R2] _f[91,R3] _f[43,Do_nothing]
 _f[43,R1] _f[43,R2] _f[43,R3] _f[76,Do_nothing] _f[76,R1] _f[76,R2]
 _f[76,R3] _f[42,Do_nothing] _f[42,R1] _f[42,R2] _f[42,R3] _f[7,Do_nothing]
 _f[7,R1] _f[7,R2] _f[7,R3] _f[17,Do_nothing] _f[17,R1] _f[17,R2] _f[17,R3]
 _f[37,Do_nothing] _f[37,R1] _f[37,R2] _f[37,R3] _f[26,Do_nothing]
 _f[26,R1] _f[26,R2] _f[26,R3] _f[9,Do_nothing] _f[9,R1] _f[9,R2] _f[9,R3]
 _f[60,Do_nothing] _f[60,R1] _f[60,R2] _f[60,R3] _f[16,Do_nothing]
 _f[16,R1] _f[16,R2] _f[16,R3] _f[83,Do_nothing] _f[83,R1] _f[83,R2]
 _f[83,R3] _f[13,Do_nothing] _f[13,R1] _f[13,R2] _f[13,R3]
 _f[67,Do_nothing] _f[67,R1] _f[67,R2] _f[67,R3] _f[59,Do_nothing]
 _f[59,R1] _f[59,R2] _f[59,R3] _f[58,Do_nothing] _f[58,R1] _f[58,R2]
 _f[58,R3] _f[4,Do_nothing] _f[4,R1] _f[4,R2] _f[4,R3] _f[68,Do_nothing]
 _f[68,R1] _f[68,R2] _f[68,R3] _f[87,Do_nothing] _f[87,R1] _f[87,R2]
 _f[87,R3] _f[40,Do_nothing] _f[40,R1] _f[40,R2] _f[40,R3]
 _f[30,Do_nothing] _f[30,R1] _f[30,R2] _f[30,R3] _f[27,Do_nothing]
 _f[27,R1] _f[27,R2] _f[27,R3] _f[73,Do_nothing] _f[73,R1] _f[73,R2]
 _f[73,R3] _f[84,Do_nothing] _f[84,R1] _f[84,R2] _f[84,R3]
 _f[35,Do_nothing] _f[35,R1] _f[35,R2] _f[35,R3] _f[46,Do_nothing]
 _f[46,R1] _f[46,R2] _f[46,R3] _f[39,Do_nothing] _f[39,R1] _f[39,R2]
 _f[39,R3] _f[51,Do_nothing] _f[51,R1] _f[51,R2] _f[51,R3]
 _f[52,Do_nothing] _f[52,R1] _f[52,R2] _f[52,R3] _f[6,Do_nothing] _f[6,R1]
 _f[6,R2] _f[6,R3] _f[96,Do_nothing] _f[96,R1] _f[96,R2] _f[96,R3]
 _f[74,Do_nothing] _f[74,R1] _f[74,R2] _f[74,R3] _f[86,Do_nothing]
 _f[86,R1] _f[86,R2] _f[86,R3] _f[29,Do_nothing] _f[29,R1] _f[29,R2]
 _f[29,R3] _f[18,Do_nothing] _f[18,R1] _f[18,R2] _f[18,R3]
 _f[24,Do_nothing] _f[24,R1] _f[24,R2] _f[24,R3] _f[94,Do_nothing]
 _f[94,R1] _f[94,R2] _f[94,R3] _f[33,Do_nothing] _f[33,R1] _f[33,R2]
 _f[33,R3] _f[12,Do_nothing] _f[12,R1] _f[12,R2] _f[12,R3]
 _f[57,Do_nothing] _f[57,R1] _f[57,R2] _f[57,R3] _f[97,Do_nothing]
 _f[97,R1] _f[97,R2] _f[97,R3] _f[69,Do_nothing] _f[69,R1] _f[69,R2]
 _f[69,R3] _f[3,Do_nothing] _f[3,R1] _f[3,R2] _f[3,R3] _f[66,Do_nothing]
 _f[66,R1] _f[66,R2] _f[66,R3]
End
