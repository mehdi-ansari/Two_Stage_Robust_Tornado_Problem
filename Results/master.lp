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
 R100: _Theta >= 0
 R101: 184945 _f[32,R1] + 232438 _f[32,R2] + 417219 _f[32,R3]
   + 77012 _f[70,R1] + 96789 _f[70,R2] + 173733 _f[70,R3] + 42150 _f[5,R1]
   + 52975 _f[5,R2] + 95089 _f[5,R3] + 83838 _f[85,R1] + 105366 _f[85,R2]
   + 189129 _f[85,R3] + 19605 _f[71,R1] + 24639 _f[71,R2] + 44226 _f[71,R3]
   + 112508 _f[92,R1] + 141402 _f[92,R2] + 253810 _f[92,R3]
   + 158010 _f[38,R1] + 198586 _f[38,R2] + 356456 _f[38,R3]
   + 227840 _f[21,R1] + 286345 _f[21,R2] + 513985 _f[21,R3]
   + 134812 _f[41,R1] + 169431 _f[41,R2] + 304131 _f[41,R3]
   + 125910 _f[56,R1] + 158240 _f[56,R2] + 284041 _f[56,R3]
   + 111072 _f[80,R1] + 139593 _f[80,R2] + 250567 _f[80,R3]
   + 88700 _f[62,R1] + 121094 _f[62,R2] + 217364 _f[62,R3]
   + 274315 _f[20,R1] + 344756 _f[20,R2] + 618825 _f[20,R3]
   + 96147 _f[10,R1] + 120838 _f[10,R2] + 216902 _f[10,R3]
   + 141199 _f[55,R1] + 177452 _f[55,R2] + 318522 _f[55,R3]
   + 90202 _f[36,R1] + 113364 _f[36,R2] + 203487 _f[36,R3]
   + 136344 _f[44,R1] + 171354 _f[44,R2] + 307578 _f[44,R3]
   + 69011 _f[48,R1] + 86734 _f[48,R2] + 155683 _f[48,R3] + 216395 _f[8,R1]
   + 271965 _f[8,R2] + 488172 _f[8,R3] + 108443 _f[15,R1]
   + 136293 _f[15,R2] + 244641 _f[15,R3] + 28666 _f[61,R1]
   + 36027 _f[61,R2] + 64668 _f[61,R3] + 139534 _f[78,R1]
   + 175362 _f[78,R2] + 314775 _f[78,R3] + 163384 _f[34,R1]
   + 205344 _f[34,R2] + 368587 _f[34,R3] + 76581 _f[100,R1]
   + 96246 _f[100,R2] + 172762 _f[100,R3] + 182173 _f[90,R1]
   + 228954 _f[90,R2] + 410968 _f[90,R3] + 253339 _f[64,R1]
   + 318395 _f[64,R2] + 571513 _f[64,R3] + 84551 _f[79,R1]
   + 106262 _f[79,R2] + 190738 _f[79,R3] + 104250 _f[1,R1]
   + 131021 _f[1,R2] + 235179 _f[1,R3] + 86007 _f[22,R1] + 108093 _f[22,R2]
   + 194022 _f[22,R3] + 149381 _f[31,R1] + 187743 _f[31,R2]
   + 336989 _f[31,R3] + 232480 _f[54,R1] + 292181 _f[54,R2]
   + 524455 _f[54,R3] + 76086 _f[82,R1] + 95624 _f[82,R2]
   + 171642 _f[82,R3] + 104028 _f[11,R1] + 130741 _f[11,R2]
   + 234681 _f[11,R3] + 139646 _f[75,R1] + 175502 _f[75,R2]
   + 315022 _f[75,R3] + 588581 _f[50,R1] + 739723 _f[50,R2]
   + 1.327789e+06 _f[50,R3] + 106626 _f[19,R1] + 134005 _f[19,R2]
   + 240534 _f[19,R3] + 113557 _f[99,R1] + 142720 _f[99,R2]
   + 256178 _f[99,R3] + 129156 _f[72,R1] + 162324 _f[72,R2]
   + 291373 _f[72,R3] + 88597 _f[47,R1] + 111344 _f[47,R2]
   + 199866 _f[47,R3] + 123902 _f[49,R1] + 155714 _f[49,R2]
   + 279506 _f[49,R3] + 21608 _f[77,R1] + 27157 _f[77,R2] + 48746 _f[77,R3]
   + 124290 _f[95,R1] + 156209 _f[95,R2] + 280390 _f[95,R3]
   + 112314 _f[98,R1] + 141154 _f[98,R2] + 253370 _f[98,R3]
   + 109154 _f[53,R1] + 137182 _f[53,R2] + 246237 _f[53,R3]
   + 89855 _f[65,R1] + 112926 _f[65,R2] + 202701 _f[65,R3]
   + 15197 _f[63,R1] + 19099 _f[63,R2] + 34283 _f[63,R3] + 82896 _f[14,R1]
   + 104182 _f[14,R2] + 187009 _f[14,R3] + 171814 _f[23,R1]
   + 215936 _f[23,R2] + 387595 _f[23,R3] + 82531 _f[93,R1]
   + 103723 _f[93,R2] + 186181 _f[93,R3] + 139133 _f[28,R1]
   + 174860 _f[28,R2] + 313866 _f[28,R3] + 266795 _f[2,R1]
   + 350290 _f[2,R2] + 628762 _f[2,R3] + 65502 _f[81,R1] + 82323 _f[81,R2]
   + 147768 _f[81,R3] + 95888 _f[45,R1] + 120511 _f[45,R2]
   + 216314 _f[45,R3] + 141611 _f[88,R1] + 177974 _f[88,R2]
   + 319460 _f[88,R3] + 48010 _f[89,R1] + 60340 _f[89,R2]
   + 108310 _f[89,R3] + 36212 _f[25,R1] + 45512 _f[25,R2] + 81689 _f[25,R3]
   + 130813 _f[91,R1] + 164403 _f[91,R2] + 295099 _f[91,R3]
   + 459976 _f[43,R1] + 578092 _f[43,R2] + 1.037663e+06 _f[43,R3]
   + 172013 _f[76,R1] + 216187 _f[76,R2] + 388050 _f[76,R3]
   + 170159 _f[42,R1] + 213855 _f[42,R2] + 383865 _f[42,R3]
   + 67468 _f[7,R1] + 84797 _f[7,R2] + 152208 _f[7,R3] + 124758 _f[17,R1]
   + 156793 _f[17,R2] + 281442 _f[17,R3] + 55522 _f[37,R1]
   + 69779 _f[37,R2] + 125252 _f[37,R3] + 169821 _f[26,R1]
   + 213431 _f[26,R2] + 383103 _f[26,R3] + 393339 _f[9,R1]
   + 494348 _f[9,R2] + 887341 _f[9,R3] + 120732 _f[60,R1]
   + 151737 _f[60,R2] + 272364 _f[60,R3] + 152263 _f[16,R1]
   + 191364 _f[16,R2] + 343496 _f[16,R3] + 99675 _f[83,R1]
   + 125272 _f[83,R2] + 224860 _f[83,R3] + 193418 _f[13,R1]
   + 243086 _f[13,R2] + 436338 _f[13,R3] + 127707 _f[67,R1]
   + 201287 _f[67,R2] + 361306 _f[67,R3] + 176600 _f[59,R1]
   + 221950 _f[59,R2] + 398395 _f[59,R3] + 100371 _f[58,R1]
   + 126143 _f[58,R2] + 226422 _f[58,R3] + 206083 _f[4,R1]
   + 259003 _f[4,R2] + 464904 _f[4,R3] + 50733 _f[68,R1] + 63759 _f[68,R2]
   + 114447 _f[68,R3] + 45515 _f[87,R1] + 57203 _f[87,R2]
   + 102678 _f[87,R3] + 205522 _f[40,R1] + 258296 _f[40,R2]
   + 463642 _f[40,R3] + 92954 _f[30,R1] + 116822 _f[30,R2]
   + 209693 _f[30,R3] + 67691 _f[27,R1] + 85070 _f[27,R2]
   + 152699 _f[27,R3] + 878240 _f[73,R1] + 1.103759e+06 _f[73,R2]
   + 1.981229e+06 _f[73,R3] + 200216 _f[84,R1] + 251629 _f[84,R2]
   + 451670 _f[84,R3] + 52951 _f[35,R1] + 66548 _f[35,R2]
   + 119449 _f[35,R3] + 191369 _f[46,R1] + 240512 _f[46,R2]
   + 431707 _f[46,R3] + 165478 _f[39,R1] + 207972 _f[39,R2]
   + 373305 _f[39,R3] + 890435 _f[51,R1] + 1.119091e+06 _f[51,R2]
   + 2.008749e+06 _f[51,R3] + 210223 _f[52,R1] + 264209 _f[52,R2]
   + 474249 _f[52,R3] + 100228 _f[6,R1] + 125965 _f[6,R2] + 226105 _f[6,R3]
   + 156136 _f[96,R1] + 196233 _f[96,R2] + 352232 _f[96,R3]
   + 26744 _f[74,R1] + 33612 _f[74,R2] + 60332 _f[74,R3] + 22382 _f[86,R1]
   + 28129 _f[86,R2] + 50491 _f[86,R3] + 59968 _f[29,R1] + 75367 _f[29,R2]
   + 135283 _f[29,R3] + 78936 _f[18,R1] + 99206 _f[18,R2]
   + 178071 _f[18,R3] + 118132 _f[24,R1] + 148466 _f[24,R2]
   + 266493 _f[24,R3] + 46689 _f[94,R1] + 58681 _f[94,R2]
   + 105329 _f[94,R3] + 63134 _f[33,R1] + 79346 _f[33,R2]
   + 142423 _f[33,R3] + 33082 _f[12,R1] + 41577 _f[12,R2] + 74630 _f[12,R3]
   + 221052 _f[57,R1] + 277816 _f[57,R2] + 498675 _f[57,R3]
   + 49387 _f[97,R1] + 62067 _f[97,R2] + 111410 _f[97,R3] + 41699 _f[69,R1]
   + 52406 _f[69,R2] + 94068 _f[69,R3] + 11993 _f[3,R1] + 15072 _f[3,R2]
   + 27055 _f[3,R3] + 36687 _f[66,R1] + 46107 _f[66,R2] + 82761 _f[66,R3]
   + 1.6250306e+07 _r0[32,Do_nothing,Recover]
   + 1.6245457e+07 _r0[32,R1,Recover] + 1.1655471e+07 _r0[32,R2,Recover]
   + 7.096177e+06 _r0[32,R3,Recover]
   + 8.395899e+06 _r0[70,Do_nothing,Recover]
   + 8.393394e+06 _r0[70,R1,Recover] + 6.021927e+06 _r0[70,R2,Recover]
   + 3.666318e+06 _r0[70,R3,Recover]
   + 3.681994e+06 _r0[5,Do_nothing,Recover]
   + 3.680895e+06 _r0[5,R1,Recover] + 2.640896e+06 _r0[5,R2,Recover]
   + 1.607852e+06 _r0[5,R3,Recover]
   + 2.0876574e+07 _r0[85,Do_nothing,Recover]
   + 2.0870346e+07 _r0[85,R1,Recover] + 1.4973644e+07 _r0[85,R2,Recover]
   + 9.116374e+06 _r0[85,R3,Recover]
   + 1.926495e+06 _r0[71,Do_nothing,Recover]
   + 1.92592e+06 _r0[71,R1,Recover] + 1.381771e+06 _r0[71,R2,Recover]
   + 841261 _r0[71,R3,Recover] + 1.3779239e+07 _r0[92,Do_nothing,Recover]
   + 1.3775128e+07 _r0[92,R1,Recover] + 9.883107e+06 _r0[92,R2,Recover]
   + 6.017113e+06 _r0[92,R3,Recover]
   + 1.6740245e+07 _r0[38,Do_nothing,Recover]
   + 1.6735251e+07 _r0[38,R1,Recover] + 1.2006878e+07 _r0[38,R2,Recover]
   + 7.310124e+06 _r0[38,R3,Recover]
   + 1.697441e+07 _r0[21,Do_nothing,Recover]
   + 1.6969346e+07 _r0[21,R1,Recover] + 1.2174832e+07 _r0[21,R2,Recover]
   + 7.412379e+06 _r0[21,R3,Recover]
   + 1.3483233e+07 _r0[41,Do_nothing,Recover]
   + 1.347921e+07 _r0[41,R1,Recover] + 9.670798e+06 _r0[41,R2,Recover]
   + 5.887853e+06 _r0[41,R3,Recover]
   + 1.5514733e+07 _r0[56,Do_nothing,Recover]
   + 1.5510104e+07 _r0[56,R1,Recover] + 1.1127884e+07 _r0[56,R2,Recover]
   + 6.774967e+06 _r0[56,R3,Recover]
   + 3.3460116e+07 _r0[80,Do_nothing,Recover]
   + 3.3450133e+07 _r0[80,R1,Recover] + 2.3999142e+07 _r0[80,R2,Recover]
   + 1.461135e+07 _r0[80,R3,Recover]
   + 1.3793197e+07 _r0[62,Do_nothing,Recover]
   + 1.3789358e+07 _r0[62,R1,Recover] + 1.0056913e+07 _r0[62,R2,Recover]
   + 6.381213e+06 _r0[62,R3,Recover]
   + 2.6317165e+07 _r0[20,Do_nothing,Recover]
   + 2.6309313e+07 _r0[20,R1,Recover] + 1.8875887e+07 _r0[20,R2,Recover]
   + 1.1492169e+07 _r0[20,R3,Recover]
   + 2.3539324e+07 _r0[10,Do_nothing,Recover]
   + 2.3532301e+07 _r0[10,R1,Recover] + 1.6883492e+07 _r0[10,R2,Recover]
   + 1.0279143e+07 _r0[10,R3,Recover]
   + 2.0671468e+07 _r0[55,Do_nothing,Recover]
   + 2.06653e+07 _r0[55,R1,Recover] + 1.4826533e+07 _r0[55,R2,Recover]
   + 9.026808e+06 _r0[55,R3,Recover]
   + 1.026872e+07 _r0[36,Do_nothing,Recover]
   + 1.0265657e+07 _r0[36,R1,Recover] + 7.365201e+06 _r0[36,R2,Recover]
   + 4.484141e+06 _r0[36,R3,Recover]
   + 2.556896e+06 _r0[44,Do_nothing,Recover]
   + 2.556133e+06 _r0[44,R1,Recover] + 1.833924e+06 _r0[44,R2,Recover]
   + 1.116544e+06 _r0[44,R3,Recover]
   + 1.677052e+07 _r0[48,Do_nothing,Recover]
   + 1.6765517e+07 _r0[48,R1,Recover] + 1.2028593e+07 _r0[48,R2,Recover]
   + 7.323344e+06 _r0[48,R3,Recover]
   + 1.6412617e+07 _r0[8,Do_nothing,Recover]
   + 1.640772e+07 _r0[8,R1,Recover] + 1.1771888e+07 _r0[8,R2,Recover]
   + 7.167055e+06 _r0[8,R3,Recover]
   + 2.4310998e+07 _r0[15,Do_nothing,Recover]
   + 2.4303745e+07 _r0[15,R1,Recover] + 1.7436972e+07 _r0[15,R2,Recover]
   + 1.0616117e+07 _r0[15,R3,Recover]
   + 1.664175e+06 _r0[61,Do_nothing,Recover]
   + 1.663678e+06 _r0[61,R1,Recover] + 1.193623e+06 _r0[61,R2,Recover]
   + 726711 _r0[61,R3,Recover] + 2.1085368e+07 _r0[78,Do_nothing,Recover]
   + 2.1079077e+07 _r0[78,R1,Recover] + 1.5123401e+07 _r0[78,R2,Recover]
   + 9.20755e+06 _r0[78,R3,Recover]
   + 5.4495756e+07 _r0[34,Do_nothing,Recover]
   + 5.4479497e+07 _r0[34,R1,Recover] + 3.9086876e+07 _r0[34,R2,Recover]
   + 2.3797185e+07 _r0[34,R3,Recover]
   + 1.3726875e+07 _r0[100,Do_nothing,Recover]
   + 1.3722779e+07 _r0[100,R1,Recover] + 9.845549e+06 _r0[100,R2,Recover]
   + 5.994246e+06 _r0[100,R3,Recover]
   + 3.5521762e+07 _r0[90,Do_nothing,Recover]
   + 3.5511165e+07 _r0[90,R1,Recover] + 2.547785e+07 _r0[90,R2,Recover]
   + 1.5511629e+07 _r0[90,R3,Recover]
   + 1.471966e+07 _r0[64,Do_nothing,Recover]
   + 1.4715269e+07 _r0[64,R1,Recover] + 1.0557621e+07 _r0[64,R2,Recover]
   + 6.427775e+06 _r0[64,R3,Recover]
   + 3.5959429e+07 _r0[79,Do_nothing,Recover]
   + 3.5948701e+07 _r0[79,R1,Recover] + 2.5791766e+07 _r0[79,R2,Recover]
   + 1.5702749e+07 _r0[79,R3,Recover]
   + 3.911387e+07 _r0[1,Do_nothing,Recover]
   + 3.9102201e+07 _r0[1,R1,Recover] + 2.8054276e+07 _r0[1,R2,Recover]
   + 1.7080229e+07 _r0[1,R3,Recover]
   + 1.2599232e+07 _r0[22,Do_nothing,Recover]
   + 1.2595474e+07 _r0[22,R1,Recover] + 9.036752e+06 _r0[22,R2,Recover]
   + 5.501828e+06 _r0[22,R3,Recover]
   + 2.6726068e+07 _r0[31,Do_nothing,Recover]
   + 2.6718095e+07 _r0[31,R1,Recover] + 1.9169172e+07 _r0[31,R2,Recover]
   + 1.1670729e+07 _r0[31,R3,Recover]
   + 1.9445107e+07 _r0[54,Do_nothing,Recover]
   + 1.9439306e+07 _r0[54,R1,Recover] + 1.394693e+07 _r0[54,R2,Recover]
   + 8.491282e+06 _r0[54,R3,Recover]
   + 1.3534085e+07 _r0[82,Do_nothing,Recover]
   + 1.3530047e+07 _r0[82,R1,Recover] + 9.707272e+06 _r0[82,R2,Recover]
   + 5.910059e+06 _r0[82,R3,Recover]
   + 1.1946244e+07 _r0[11,Do_nothing,Recover]
   + 1.194268e+07 _r0[11,R1,Recover] + 8.568399e+06 _r0[11,R2,Recover]
   + 5.216681e+06 _r0[11,R3,Recover]
   + 1.9233548e+07 _r0[75,Do_nothing,Recover]
   + 1.922781e+07 _r0[75,R1,Recover] + 1.379519e+07 _r0[75,R2,Recover]
   + 8.398898e+06 _r0[75,R3,Recover]
   + 7.9890053e+07 _r0[50,Do_nothing,Recover]
   + 7.9866218e+07 _r0[50,R1,Recover] + 5.730084e+07 _r0[50,R2,Recover]
   + 3.4886357e+07 _r0[50,R3,Recover]
   + 3.7772492e+07 _r0[19,Do_nothing,Recover]
   + 3.7761223e+07 _r0[19,R1,Recover] + 2.7092178e+07 _r0[19,R2,Recover]
   + 1.6494477e+07 _r0[19,R3,Recover]
   + 1.401132e+07 _r0[99,Do_nothing,Recover]
   + 1.400714e+07 _r0[99,R1,Recover] + 1.0049567e+07 _r0[99,R2,Recover]
   + 6.118458e+06 _r0[99,R3,Recover]
   + 2.3511427e+07 _r0[72,Do_nothing,Recover]
   + 2.3504412e+07 _r0[72,R1,Recover] + 1.6863483e+07 _r0[72,R2,Recover]
   + 1.0266961e+07 _r0[72,R3,Recover]
   + 1.0302942e+07 _r0[47,Do_nothing,Recover]
   + 1.0299868e+07 _r0[47,R1,Recover] + 7.389746e+06 _r0[47,R2,Recover]
   + 4.499084e+06 _r0[47,R3,Recover]
   + 3.5335076e+07 _r0[49,Do_nothing,Recover]
   + 3.5324534e+07 _r0[49,R1,Recover] + 2.534395e+07 _r0[49,R2,Recover]
   + 1.5430107e+07 _r0[49,R3,Recover]
   + 2.275954e+06 _r0[77,Do_nothing,Recover]
   + 2.275275e+06 _r0[77,R1,Recover] + 1.632419e+06 _r0[77,R2,Recover]
   + 993863 _r0[77,R3,Recover] + 3.0225358e+07 _r0[95,Do_nothing,Recover]
   + 3.0216341e+07 _r0[95,R1,Recover] + 2.1679025e+07 _r0[95,R2,Recover]
   + 1.3198797e+07 _r0[95,R3,Recover]
   + 8.36295e+06 _r0[98,Do_nothing,Recover]
   + 8.360455e+06 _r0[98,R1,Recover] + 5.998294e+06 _r0[98,R2,Recover]
   + 3.65193e+06 _r0[98,R3,Recover]
   + 2.399908e+07 _r0[53,Do_nothing,Recover]
   + 2.399192e+07 _r0[53,R1,Recover] + 1.721325e+07 _r0[53,R2,Recover]
   + 1.0479909e+07 _r0[53,R3,Recover]
   + 1.3037471e+07 _r0[65,Do_nothing,Recover]
   + 1.3033581e+07 _r0[65,R1,Recover] + 9.351077e+06 _r0[65,R2,Recover]
   + 5.693198e+06 _r0[65,R3,Recover]
   + 5.549173e+06 _r0[63,Do_nothing,Recover]
   + 5.547518e+06 _r0[63,R1,Recover] + 3.980124e+06 _r0[63,R2,Recover]
   + 2.423211e+06 _r0[63,R3,Recover]
   + 2.8131869e+07 _r0[14,Do_nothing,Recover]
   + 2.8123476e+07 _r0[14,R1,Recover] + 2.0177477e+07 _r0[14,R2,Recover]
   + 1.2284613e+07 _r0[14,R3,Recover]
   + 1.4888056e+07 _r0[23,Do_nothing,Recover]
   + 1.4883614e+07 _r0[23,R1,Recover] + 1.0678402e+07 _r0[23,R2,Recover]
   + 6.50131e+06 _r0[23,R3,Recover]
   + 1.4395332e+07 _r0[93,Do_nothing,Recover]
   + 1.4391038e+07 _r0[93,R1,Recover] + 1.0324998e+07 _r0[93,R2,Recover]
   + 6.286148e+06 _r0[93,R3,Recover]
   + 1.3562886e+07 _r0[28,Do_nothing,Recover]
   + 1.3558839e+07 _r0[28,R1,Recover] + 9.727929e+06 _r0[28,R2,Recover]
   + 5.922636e+06 _r0[28,R3,Recover]
   + 2.598607e+07 _r0[2,Do_nothing,Recover]
   + 2.5978609e+07 _r0[2,R1,Recover] + 1.8811219e+07 _r0[2,R2,Recover]
   + 1.1725292e+07 _r0[2,R3,Recover]
   + 1.1245096e+07 _r0[81,Do_nothing,Recover]
   + 1.1241741e+07 _r0[81,R1,Recover] + 8.065503e+06 _r0[81,R2,Recover]
   + 4.910504e+06 _r0[81,R3,Recover]
   + 2.0483269e+07 _r0[45,Do_nothing,Recover]
   + 2.0477158e+07 _r0[45,R1,Recover] + 1.4691548e+07 _r0[45,R2,Recover]
   + 8.944626e+06 _r0[45,R3,Recover]
   + 1.7600312e+07 _r0[88,Do_nothing,Recover]
   + 1.7595061e+07 _r0[88,R1,Recover] + 1.2623758e+07 _r0[88,R2,Recover]
   + 7.685697e+06 _r0[88,R3,Recover]
   + 7.828206e+06 _r0[89,Do_nothing,Recover]
   + 7.82587e+06 _r0[89,R1,Recover] + 5.614751e+06 _r0[89,R2,Recover]
   + 3.418418e+06 _r0[89,R3,Recover]
   + 5.800854e+06 _r0[25,Do_nothing,Recover]
   + 5.799124e+06 _r0[25,R1,Recover] + 4.160641e+06 _r0[25,R2,Recover]
   + 2.533115e+06 _r0[25,R3,Recover]
   + 2.0526875e+07 _r0[91,Do_nothing,Recover]
   + 2.0520751e+07 _r0[91,R1,Recover] + 1.4722824e+07 _r0[91,R2,Recover]
   + 8.963668e+06 _r0[91,R3,Recover]
   + 4.9197858e+07 _r0[43,Do_nothing,Recover]
   + 4.918318e+07 _r0[43,R1,Recover] + 3.5286979e+07 _r0[43,R2,Recover]
   + 2.1483701e+07 _r0[43,R3,Recover]
   + 1.09461187e+08 _r0[76,Do_nothing,Recover]
   + 1.0942853e+08 _r0[76,R1,Recover] + 7.8510625e+07 _r0[76,R2,Recover]
   + 4.7799468e+07 _r0[76,R3,Recover]
   + 3.3561304e+07 _r0[42,Do_nothing,Recover]
   + 3.3551292e+07 _r0[42,R1,Recover] + 2.4071719e+07 _r0[42,R2,Recover]
   + 1.4655537e+07 _r0[42,R3,Recover]
   + 1.7648288e+07 _r0[7,Do_nothing,Recover]
   + 1.7643023e+07 _r0[7,R1,Recover] + 1.2658169e+07 _r0[7,R2,Recover]
   + 7.706648e+06 _r0[7,R3,Recover]
   + 1.9181497e+07 _r0[17,Do_nothing,Recover]
   + 1.9175774e+07 _r0[17,R1,Recover] + 1.3757856e+07 _r0[17,R2,Recover]
   + 8.376168e+06 _r0[17,R3,Recover]
   + 1.0390541e+07 _r0[37,Do_nothing,Recover]
   + 1.0387441e+07 _r0[37,R1,Recover] + 7.452577e+06 _r0[37,R2,Recover]
   + 4.537337e+06 _r0[37,R3,Recover]
   + 2.0918355e+07 _r0[26,Do_nothing,Recover]
   + 2.0912114e+07 _r0[26,R1,Recover] + 1.5003611e+07 _r0[26,R2,Recover]
   + 9.134619e+06 _r0[26,R3,Recover]
   + 5.8801993e+07 _r0[9,Do_nothing,Recover]
   + 5.878445e+07 _r0[9,R1,Recover] + 4.2175508e+07 _r0[9,R2,Recover]
   + 2.5677631e+07 _r0[9,R3,Recover]
   + 1.5265061e+07 _r0[60,Do_nothing,Recover]
   + 1.5260507e+07 _r0[60,R1,Recover] + 1.0948808e+07 _r0[60,R2,Recover]
   + 6.665941e+06 _r0[60,R3,Recover]
   + 1.5422118e+07 _r0[16,Do_nothing,Recover]
   + 1.5417517e+07 _r0[16,R1,Recover] + 1.1061456e+07 _r0[16,R2,Recover]
   + 6.734524e+06 _r0[16,R3,Recover]
   + 1.3153079e+07 _r0[83,Do_nothing,Recover]
   + 1.3149154e+07 _r0[83,R1,Recover] + 9.433996e+06 _r0[83,R2,Recover]
   + 5.743681e+06 _r0[83,R3,Recover]
   + 7.2047115e+07 _r0[13,Do_nothing,Recover]
   + 7.202562e+07 _r0[13,R1,Recover] + 5.1675522e+07 _r0[13,R2,Recover]
   + 3.1461506e+07 _r0[13,R3,Recover]
   + 5.5945218e+07 _r0[67,Do_nothing,Recover]
   + 5.5931138e+07 _r0[67,R1,Recover] + 4.1676653e+07 _r0[67,R2,Recover]
   + 2.7818295e+07 _r0[67,R3,Recover]
   + 2.1147689e+07 _r0[59,Do_nothing,Recover]
   + 2.114138e+07 _r0[59,R1,Recover] + 1.51681e+07 _r0[59,R2,Recover]
   + 9.234765e+06 _r0[59,R3,Recover]
   + 1.4856103e+07 _r0[58,Do_nothing,Recover]
   + 1.485167e+07 _r0[58,R1,Recover] + 1.0655484e+07 _r0[58,R2,Recover]
   + 6.487357e+06 _r0[58,R3,Recover]
   + 2.2818833e+07 _r0[4,Do_nothing,Recover]
   + 2.2812026e+07 _r0[4,R1,Recover] + 1.6366723e+07 _r0[4,R2,Recover]
   + 9.964519e+06 _r0[4,R3,Recover]
   + 1.8591493e+07 _r0[68,Do_nothing,Recover]
   + 1.8585947e+07 _r0[68,R1,Recover] + 1.3334679e+07 _r0[68,R2,Recover]
   + 8.118526e+06 _r0[68,R3,Recover]
   + 4.873506e+06 _r0[87,Do_nothing,Recover]
   + 4.872052e+06 _r0[87,R1,Recover] + 3.495504e+06 _r0[87,R2,Recover]
   + 2.128161e+06 _r0[87,R3,Recover]
   + 2.1259203e+07 _r0[40,Do_nothing,Recover]
   + 2.1252861e+07 _r0[40,R1,Recover] + 1.5248084e+07 _r0[40,R2,Recover]
   + 9.28346e+06 _r0[40,R3,Recover]
   + 9.148287e+06 _r0[30,Do_nothing,Recover]
   + 9.145558e+06 _r0[30,R1,Recover] + 6.561575e+06 _r0[30,R2,Recover]
   + 3.994871e+06 _r0[30,R3,Recover]
   + 1.2238711e+07 _r0[27,Do_nothing,Recover]
   + 1.2235059e+07 _r0[27,R1,Recover] + 8.778169e+06 _r0[27,R2,Recover]
   + 5.344395e+06 _r0[27,R3,Recover]
   + 3.1009091e+07 _r0[73,Do_nothing,Recover]
   + 3.099984e+07 _r0[73,R1,Recover] + 2.2241154e+07 _r0[73,R2,Recover]
   + 1.3541038e+07 _r0[73,R3,Recover]
   + 2.2537319e+07 _r0[84,Do_nothing,Recover]
   + 2.2530595e+07 _r0[84,R1,Recover] + 1.6164808e+07 _r0[84,R2,Recover]
   + 9.841588e+06 _r0[84,R3,Recover]
   + 4.660656e+06 _r0[35,Do_nothing,Recover]
   + 4.659265e+06 _r0[35,R1,Recover] + 3.342838e+06 _r0[35,R2,Recover]
   + 2.035213e+06 _r0[35,R3,Recover]
   + 8.4516505e+07 _r0[46,Do_nothing,Recover]
   + 8.4491291e+07 _r0[46,R1,Recover] + 6.0619146e+07 _r0[46,R2,Recover]
   + 3.6906634e+07 _r0[46,R3,Recover]
   + 2.0923204e+07 _r0[39,Do_nothing,Recover]
   + 2.0916962e+07 _r0[39,R1,Recover] + 1.500709e+07 _r0[39,R2,Recover]
   + 9.136737e+06 _r0[39,R3,Recover]
   + 1.9955349e+08 _r0[51,Do_nothing,Recover]
   + 1.99493955e+08 _r0[51,R1,Recover] + 1.43128991e+08 _r0[51,R2,Recover]
   + 8.714094e+07 _r0[51,R3,Recover]
   + 7.9890532e+07 _r0[52,Do_nothing,Recover]
   + 7.9866697e+07 _r0[52,R1,Recover] + 5.7301184e+07 _r0[52,R2,Recover]
   + 3.4886566e+07 _r0[52,R3,Recover]
   + 1.6577583e+07 _r0[6,Do_nothing,Recover]
   + 1.6572637e+07 _r0[6,R1,Recover] + 1.1890209e+07 _r0[6,R2,Recover]
   + 7.239093e+06 _r0[6,R3,Recover]
   + 1.2827184e+07 _r0[96,Do_nothing,Recover]
   + 1.2823357e+07 _r0[96,R1,Recover] + 9.200249e+06 _r0[96,R2,Recover]
   + 5.60137e+06 _r0[96,R3,Recover]
   + 1.0275635e+07 _r0[74,Do_nothing,Recover]
   + 1.0272569e+07 _r0[74,R1,Recover] + 7.37016e+06 _r0[74,R2,Recover]
   + 4.48716e+06 _r0[74,R3,Recover] + 880257 _r0[86,Do_nothing,Recover]
   + 879995 _r0[86,R1,Recover] + 631361 _r0[86,R2,Recover]
   + 384390 _r0[86,R3,Recover] + 5.748969e+06 _r0[29,Do_nothing,Recover]
   + 5.747254e+06 _r0[29,R1,Recover] + 4.123427e+06 _r0[29,R2,Recover]
   + 2.510458e+06 _r0[29,R3,Recover]
   + 5.511117e+06 _r0[18,Do_nothing,Recover]
   + 5.509472e+06 _r0[18,R1,Recover] + 3.952828e+06 _r0[18,R2,Recover]
   + 2.406592e+06 _r0[18,R3,Recover]
   + 2.1123609e+07 _r0[24,Do_nothing,Recover]
   + 2.1117307e+07 _r0[24,R1,Recover] + 1.5150829e+07 _r0[24,R2,Recover]
   + 9.224249e+06 _r0[24,R3,Recover]
   + 4.315511e+06 _r0[94,Do_nothing,Recover]
   + 4.314224e+06 _r0[94,R1,Recover] + 3.095284e+06 _r0[94,R2,Recover]
   + 1.884496e+06 _r0[94,R3,Recover]
   + 5.1681149e+07 _r0[33,Do_nothing,Recover]
   + 5.166573e+07 _r0[33,R1,Recover] + 3.706811e+07 _r0[33,R2,Recover]
   + 2.2568104e+07 _r0[33,R3,Recover]
   + 1.378918e+06 _r0[12,Do_nothing,Recover]
   + 1.378506e+06 _r0[12,R1,Recover] + 989024 _r0[12,R2,Recover]
   + 602145 _r0[12,R3,Recover] + 4.2913776e+07 _r0[57,Do_nothing,Recover]
   + 4.2900973e+07 _r0[57,R1,Recover] + 3.0779745e+07 _r0[57,R2,Recover]
   + 1.8739571e+07 _r0[57,R3,Recover]
   + 2.092752e+06 _r0[97,Do_nothing,Recover]
   + 2.092128e+06 _r0[97,R1,Recover] + 1.501018e+06 _r0[97,R2,Recover]
   + 913862 _r0[97,R3,Recover] + 1.461207e+06 _r0[69,Do_nothing,Recover]
   + 1.460771e+06 _r0[69,R1,Recover] + 1.048045e+06 _r0[69,R2,Recover]
   + 638079 _r0[69,R3,Recover] + 1.8390093e+07 _r0[3,Do_nothing,Recover]
   + 1.8384607e+07 _r0[3,R1,Recover] + 1.3190225e+07 _r0[3,R2,Recover]
   + 8.030579e+06 _r0[3,R3,Recover] + 396329 _r0[66,Do_nothing,Recover]
   + 396210 _r0[66,R1,Recover] + 284265 _r0[66,R2,Recover]
   + 173069 _r0[66,R3,Recover] <= 1.5e+07
 R102: - _f[32,Do_nothing] + _r0[32,Do_nothing,Do_nothing]
   + _r0[32,Do_nothing,Recover] = 0
 R103: - _f[32,R1] + _r0[32,R1,Do_nothing] + _r0[32,R1,Recover] = 0
 R104: - _f[32,R2] + _r0[32,R2,Do_nothing] + _r0[32,R2,Recover] = 0
 R105: - _f[32,R3] + _r0[32,R3,Do_nothing] + _r0[32,R3,Recover] = 0
 R106: - _f[70,Do_nothing] + _r0[70,Do_nothing,Do_nothing]
   + _r0[70,Do_nothing,Recover] = 0
 R107: - _f[70,R1] + _r0[70,R1,Do_nothing] + _r0[70,R1,Recover] = 0
 R108: - _f[70,R2] + _r0[70,R2,Do_nothing] + _r0[70,R2,Recover] = 0
 R109: - _f[70,R3] + _r0[70,R3,Do_nothing] + _r0[70,R3,Recover] = 0
 R110: - _f[5,Do_nothing] + _r0[5,Do_nothing,Do_nothing]
   + _r0[5,Do_nothing,Recover] = 0
 R111: - _f[5,R1] + _r0[5,R1,Do_nothing] + _r0[5,R1,Recover] = 0
 R112: - _f[5,R2] + _r0[5,R2,Do_nothing] + _r0[5,R2,Recover] = 0
 R113: - _f[5,R3] + _r0[5,R3,Do_nothing] + _r0[5,R3,Recover] = 0
 R114: - _f[85,Do_nothing] + _r0[85,Do_nothing,Do_nothing]
   + _r0[85,Do_nothing,Recover] = 0
 R115: - _f[85,R1] + _r0[85,R1,Do_nothing] + _r0[85,R1,Recover] = 0
 R116: - _f[85,R2] + _r0[85,R2,Do_nothing] + _r0[85,R2,Recover] = 0
 R117: - _f[85,R3] + _r0[85,R3,Do_nothing] + _r0[85,R3,Recover] = 0
 R118: - _f[71,Do_nothing] + _r0[71,Do_nothing,Do_nothing]
   + _r0[71,Do_nothing,Recover] = 0
 R119: - _f[71,R1] + _r0[71,R1,Do_nothing] + _r0[71,R1,Recover] = 0
 R120: - _f[71,R2] + _r0[71,R2,Do_nothing] + _r0[71,R2,Recover] = 0
 R121: - _f[71,R3] + _r0[71,R3,Do_nothing] + _r0[71,R3,Recover] = 0
 R122: - _f[92,Do_nothing] + _r0[92,Do_nothing,Do_nothing]
   + _r0[92,Do_nothing,Recover] = 0
 R123: - _f[92,R1] + _r0[92,R1,Do_nothing] + _r0[92,R1,Recover] = 0
 R124: - _f[92,R2] + _r0[92,R2,Do_nothing] + _r0[92,R2,Recover] = 0
 R125: - _f[92,R3] + _r0[92,R3,Do_nothing] + _r0[92,R3,Recover] = 0
 R126: - _f[38,Do_nothing] + _r0[38,Do_nothing,Do_nothing]
   + _r0[38,Do_nothing,Recover] = 0
 R127: - _f[38,R1] + _r0[38,R1,Do_nothing] + _r0[38,R1,Recover] = 0
 R128: - _f[38,R2] + _r0[38,R2,Do_nothing] + _r0[38,R2,Recover] = 0
 R129: - _f[38,R3] + _r0[38,R3,Do_nothing] + _r0[38,R3,Recover] = 0
 R130: - _f[21,Do_nothing] + _r0[21,Do_nothing,Do_nothing]
   + _r0[21,Do_nothing,Recover] = 0
 R131: - _f[21,R1] + _r0[21,R1,Do_nothing] + _r0[21,R1,Recover] = 0
 R132: - _f[21,R2] + _r0[21,R2,Do_nothing] + _r0[21,R2,Recover] = 0
 R133: - _f[21,R3] + _r0[21,R3,Do_nothing] + _r0[21,R3,Recover] = 0
 R134: - _f[41,Do_nothing] + _r0[41,Do_nothing,Do_nothing]
   + _r0[41,Do_nothing,Recover] = 0
 R135: - _f[41,R1] + _r0[41,R1,Do_nothing] + _r0[41,R1,Recover] = 0
 R136: - _f[41,R2] + _r0[41,R2,Do_nothing] + _r0[41,R2,Recover] = 0
 R137: - _f[41,R3] + _r0[41,R3,Do_nothing] + _r0[41,R3,Recover] = 0
 R138: - _f[56,Do_nothing] + _r0[56,Do_nothing,Do_nothing]
   + _r0[56,Do_nothing,Recover] = 0
 R139: - _f[56,R1] + _r0[56,R1,Do_nothing] + _r0[56,R1,Recover] = 0
 R140: - _f[56,R2] + _r0[56,R2,Do_nothing] + _r0[56,R2,Recover] = 0
 R141: - _f[56,R3] + _r0[56,R3,Do_nothing] + _r0[56,R3,Recover] = 0
 R142: - _f[80,Do_nothing] + _r0[80,Do_nothing,Do_nothing]
   + _r0[80,Do_nothing,Recover] = 0
 R143: - _f[80,R1] + _r0[80,R1,Do_nothing] + _r0[80,R1,Recover] = 0
 R144: - _f[80,R2] + _r0[80,R2,Do_nothing] + _r0[80,R2,Recover] = 0
 R145: - _f[80,R3] + _r0[80,R3,Do_nothing] + _r0[80,R3,Recover] = 0
 R146: - _f[62,Do_nothing] + _r0[62,Do_nothing,Do_nothing]
   + _r0[62,Do_nothing,Recover] = 0
 R147: - _f[62,R1] + _r0[62,R1,Do_nothing] + _r0[62,R1,Recover] = 0
 R148: - _f[62,R2] + _r0[62,R2,Do_nothing] + _r0[62,R2,Recover] = 0
 R149: - _f[62,R3] + _r0[62,R3,Do_nothing] + _r0[62,R3,Recover] = 0
 R150: - _f[20,Do_nothing] + _r0[20,Do_nothing,Do_nothing]
   + _r0[20,Do_nothing,Recover] = 0
 R151: - _f[20,R1] + _r0[20,R1,Do_nothing] + _r0[20,R1,Recover] = 0
 R152: - _f[20,R2] + _r0[20,R2,Do_nothing] + _r0[20,R2,Recover] = 0
 R153: - _f[20,R3] + _r0[20,R3,Do_nothing] + _r0[20,R3,Recover] = 0
 R154: - _f[10,Do_nothing] + _r0[10,Do_nothing,Do_nothing]
   + _r0[10,Do_nothing,Recover] = 0
 R155: - _f[10,R1] + _r0[10,R1,Do_nothing] + _r0[10,R1,Recover] = 0
 R156: - _f[10,R2] + _r0[10,R2,Do_nothing] + _r0[10,R2,Recover] = 0
 R157: - _f[10,R3] + _r0[10,R3,Do_nothing] + _r0[10,R3,Recover] = 0
 R158: - _f[55,Do_nothing] + _r0[55,Do_nothing,Do_nothing]
   + _r0[55,Do_nothing,Recover] = 0
 R159: - _f[55,R1] + _r0[55,R1,Do_nothing] + _r0[55,R1,Recover] = 0
 R160: - _f[55,R2] + _r0[55,R2,Do_nothing] + _r0[55,R2,Recover] = 0
 R161: - _f[55,R3] + _r0[55,R3,Do_nothing] + _r0[55,R3,Recover] = 0
 R162: - _f[36,Do_nothing] + _r0[36,Do_nothing,Do_nothing]
   + _r0[36,Do_nothing,Recover] = 0
 R163: - _f[36,R1] + _r0[36,R1,Do_nothing] + _r0[36,R1,Recover] = 0
 R164: - _f[36,R2] + _r0[36,R2,Do_nothing] + _r0[36,R2,Recover] = 0
 R165: - _f[36,R3] + _r0[36,R3,Do_nothing] + _r0[36,R3,Recover] = 0
 R166: - _f[44,Do_nothing] + _r0[44,Do_nothing,Do_nothing]
   + _r0[44,Do_nothing,Recover] = 0
 R167: - _f[44,R1] + _r0[44,R1,Do_nothing] + _r0[44,R1,Recover] = 0
 R168: - _f[44,R2] + _r0[44,R2,Do_nothing] + _r0[44,R2,Recover] = 0
 R169: - _f[44,R3] + _r0[44,R3,Do_nothing] + _r0[44,R3,Recover] = 0
 R170: - _f[48,Do_nothing] + _r0[48,Do_nothing,Do_nothing]
   + _r0[48,Do_nothing,Recover] = 0
 R171: - _f[48,R1] + _r0[48,R1,Do_nothing] + _r0[48,R1,Recover] = 0
 R172: - _f[48,R2] + _r0[48,R2,Do_nothing] + _r0[48,R2,Recover] = 0
 R173: - _f[48,R3] + _r0[48,R3,Do_nothing] + _r0[48,R3,Recover] = 0
 R174: - _f[8,Do_nothing] + _r0[8,Do_nothing,Do_nothing]
   + _r0[8,Do_nothing,Recover] = 0
 R175: - _f[8,R1] + _r0[8,R1,Do_nothing] + _r0[8,R1,Recover] = 0
 R176: - _f[8,R2] + _r0[8,R2,Do_nothing] + _r0[8,R2,Recover] = 0
 R177: - _f[8,R3] + _r0[8,R3,Do_nothing] + _r0[8,R3,Recover] = 0
 R178: - _f[15,Do_nothing] + _r0[15,Do_nothing,Do_nothing]
   + _r0[15,Do_nothing,Recover] = 0
 R179: - _f[15,R1] + _r0[15,R1,Do_nothing] + _r0[15,R1,Recover] = 0
 R180: - _f[15,R2] + _r0[15,R2,Do_nothing] + _r0[15,R2,Recover] = 0
 R181: - _f[15,R3] + _r0[15,R3,Do_nothing] + _r0[15,R3,Recover] = 0
 R182: - _f[61,Do_nothing] + _r0[61,Do_nothing,Do_nothing]
   + _r0[61,Do_nothing,Recover] = 0
 R183: - _f[61,R1] + _r0[61,R1,Do_nothing] + _r0[61,R1,Recover] = 0
 R184: - _f[61,R2] + _r0[61,R2,Do_nothing] + _r0[61,R2,Recover] = 0
 R185: - _f[61,R3] + _r0[61,R3,Do_nothing] + _r0[61,R3,Recover] = 0
 R186: - _f[78,Do_nothing] + _r0[78,Do_nothing,Do_nothing]
   + _r0[78,Do_nothing,Recover] = 0
 R187: - _f[78,R1] + _r0[78,R1,Do_nothing] + _r0[78,R1,Recover] = 0
 R188: - _f[78,R2] + _r0[78,R2,Do_nothing] + _r0[78,R2,Recover] = 0
 R189: - _f[78,R3] + _r0[78,R3,Do_nothing] + _r0[78,R3,Recover] = 0
 R190: - _f[34,Do_nothing] + _r0[34,Do_nothing,Do_nothing]
   + _r0[34,Do_nothing,Recover] = 0
 R191: - _f[34,R1] + _r0[34,R1,Do_nothing] + _r0[34,R1,Recover] = 0
 R192: - _f[34,R2] + _r0[34,R2,Do_nothing] + _r0[34,R2,Recover] = 0
 R193: - _f[34,R3] + _r0[34,R3,Do_nothing] + _r0[34,R3,Recover] = 0
 R194: - _f[100,Do_nothing] + _r0[100,Do_nothing,Do_nothing]
   + _r0[100,Do_nothing,Recover] = 0
 R195: - _f[100,R1] + _r0[100,R1,Do_nothing] + _r0[100,R1,Recover] = 0
 R196: - _f[100,R2] + _r0[100,R2,Do_nothing] + _r0[100,R2,Recover] = 0
 R197: - _f[100,R3] + _r0[100,R3,Do_nothing] + _r0[100,R3,Recover] = 0
 R198: - _f[90,Do_nothing] + _r0[90,Do_nothing,Do_nothing]
   + _r0[90,Do_nothing,Recover] = 0
 R199: - _f[90,R1] + _r0[90,R1,Do_nothing] + _r0[90,R1,Recover] = 0
 R200: - _f[90,R2] + _r0[90,R2,Do_nothing] + _r0[90,R2,Recover] = 0
 R201: - _f[90,R3] + _r0[90,R3,Do_nothing] + _r0[90,R3,Recover] = 0
 R202: - _f[64,Do_nothing] + _r0[64,Do_nothing,Do_nothing]
   + _r0[64,Do_nothing,Recover] = 0
 R203: - _f[64,R1] + _r0[64,R1,Do_nothing] + _r0[64,R1,Recover] = 0
 R204: - _f[64,R2] + _r0[64,R2,Do_nothing] + _r0[64,R2,Recover] = 0
 R205: - _f[64,R3] + _r0[64,R3,Do_nothing] + _r0[64,R3,Recover] = 0
 R206: - _f[79,Do_nothing] + _r0[79,Do_nothing,Do_nothing]
   + _r0[79,Do_nothing,Recover] = 0
 R207: - _f[79,R1] + _r0[79,R1,Do_nothing] + _r0[79,R1,Recover] = 0
 R208: - _f[79,R2] + _r0[79,R2,Do_nothing] + _r0[79,R2,Recover] = 0
 R209: - _f[79,R3] + _r0[79,R3,Do_nothing] + _r0[79,R3,Recover] = 0
 R210: - _f[1,Do_nothing] + _r0[1,Do_nothing,Do_nothing]
   + _r0[1,Do_nothing,Recover] = 0
 R211: - _f[1,R1] + _r0[1,R1,Do_nothing] + _r0[1,R1,Recover] = 0
 R212: - _f[1,R2] + _r0[1,R2,Do_nothing] + _r0[1,R2,Recover] = 0
 R213: - _f[1,R3] + _r0[1,R3,Do_nothing] + _r0[1,R3,Recover] = 0
 R214: - _f[22,Do_nothing] + _r0[22,Do_nothing,Do_nothing]
   + _r0[22,Do_nothing,Recover] = 0
 R215: - _f[22,R1] + _r0[22,R1,Do_nothing] + _r0[22,R1,Recover] = 0
 R216: - _f[22,R2] + _r0[22,R2,Do_nothing] + _r0[22,R2,Recover] = 0
 R217: - _f[22,R3] + _r0[22,R3,Do_nothing] + _r0[22,R3,Recover] = 0
 R218: - _f[31,Do_nothing] + _r0[31,Do_nothing,Do_nothing]
   + _r0[31,Do_nothing,Recover] = 0
 R219: - _f[31,R1] + _r0[31,R1,Do_nothing] + _r0[31,R1,Recover] = 0
 R220: - _f[31,R2] + _r0[31,R2,Do_nothing] + _r0[31,R2,Recover] = 0
 R221: - _f[31,R3] + _r0[31,R3,Do_nothing] + _r0[31,R3,Recover] = 0
 R222: - _f[54,Do_nothing] + _r0[54,Do_nothing,Do_nothing]
   + _r0[54,Do_nothing,Recover] = 0
 R223: - _f[54,R1] + _r0[54,R1,Do_nothing] + _r0[54,R1,Recover] = 0
 R224: - _f[54,R2] + _r0[54,R2,Do_nothing] + _r0[54,R2,Recover] = 0
 R225: - _f[54,R3] + _r0[54,R3,Do_nothing] + _r0[54,R3,Recover] = 0
 R226: - _f[82,Do_nothing] + _r0[82,Do_nothing,Do_nothing]
   + _r0[82,Do_nothing,Recover] = 0
 R227: - _f[82,R1] + _r0[82,R1,Do_nothing] + _r0[82,R1,Recover] = 0
 R228: - _f[82,R2] + _r0[82,R2,Do_nothing] + _r0[82,R2,Recover] = 0
 R229: - _f[82,R3] + _r0[82,R3,Do_nothing] + _r0[82,R3,Recover] = 0
 R230: - _f[11,Do_nothing] + _r0[11,Do_nothing,Do_nothing]
   + _r0[11,Do_nothing,Recover] = 0
 R231: - _f[11,R1] + _r0[11,R1,Do_nothing] + _r0[11,R1,Recover] = 0
 R232: - _f[11,R2] + _r0[11,R2,Do_nothing] + _r0[11,R2,Recover] = 0
 R233: - _f[11,R3] + _r0[11,R3,Do_nothing] + _r0[11,R3,Recover] = 0
 R234: - _f[75,Do_nothing] + _r0[75,Do_nothing,Do_nothing]
   + _r0[75,Do_nothing,Recover] = 0
 R235: - _f[75,R1] + _r0[75,R1,Do_nothing] + _r0[75,R1,Recover] = 0
 R236: - _f[75,R2] + _r0[75,R2,Do_nothing] + _r0[75,R2,Recover] = 0
 R237: - _f[75,R3] + _r0[75,R3,Do_nothing] + _r0[75,R3,Recover] = 0
 R238: - _f[50,Do_nothing] + _r0[50,Do_nothing,Do_nothing]
   + _r0[50,Do_nothing,Recover] = 0
 R239: - _f[50,R1] + _r0[50,R1,Do_nothing] + _r0[50,R1,Recover] = 0
 R240: - _f[50,R2] + _r0[50,R2,Do_nothing] + _r0[50,R2,Recover] = 0
 R241: - _f[50,R3] + _r0[50,R3,Do_nothing] + _r0[50,R3,Recover] = 0
 R242: - _f[19,Do_nothing] + _r0[19,Do_nothing,Do_nothing]
   + _r0[19,Do_nothing,Recover] = 0
 R243: - _f[19,R1] + _r0[19,R1,Do_nothing] + _r0[19,R1,Recover] = 0
 R244: - _f[19,R2] + _r0[19,R2,Do_nothing] + _r0[19,R2,Recover] = 0
 R245: - _f[19,R3] + _r0[19,R3,Do_nothing] + _r0[19,R3,Recover] = 0
 R246: - _f[99,Do_nothing] + _r0[99,Do_nothing,Do_nothing]
   + _r0[99,Do_nothing,Recover] = 0
 R247: - _f[99,R1] + _r0[99,R1,Do_nothing] + _r0[99,R1,Recover] = 0
 R248: - _f[99,R2] + _r0[99,R2,Do_nothing] + _r0[99,R2,Recover] = 0
 R249: - _f[99,R3] + _r0[99,R3,Do_nothing] + _r0[99,R3,Recover] = 0
 R250: - _f[72,Do_nothing] + _r0[72,Do_nothing,Do_nothing]
   + _r0[72,Do_nothing,Recover] = 0
 R251: - _f[72,R1] + _r0[72,R1,Do_nothing] + _r0[72,R1,Recover] = 0
 R252: - _f[72,R2] + _r0[72,R2,Do_nothing] + _r0[72,R2,Recover] = 0
 R253: - _f[72,R3] + _r0[72,R3,Do_nothing] + _r0[72,R3,Recover] = 0
 R254: - _f[47,Do_nothing] + _r0[47,Do_nothing,Do_nothing]
   + _r0[47,Do_nothing,Recover] = 0
 R255: - _f[47,R1] + _r0[47,R1,Do_nothing] + _r0[47,R1,Recover] = 0
 R256: - _f[47,R2] + _r0[47,R2,Do_nothing] + _r0[47,R2,Recover] = 0
 R257: - _f[47,R3] + _r0[47,R3,Do_nothing] + _r0[47,R3,Recover] = 0
 R258: - _f[49,Do_nothing] + _r0[49,Do_nothing,Do_nothing]
   + _r0[49,Do_nothing,Recover] = 0
 R259: - _f[49,R1] + _r0[49,R1,Do_nothing] + _r0[49,R1,Recover] = 0
 R260: - _f[49,R2] + _r0[49,R2,Do_nothing] + _r0[49,R2,Recover] = 0
 R261: - _f[49,R3] + _r0[49,R3,Do_nothing] + _r0[49,R3,Recover] = 0
 R262: - _f[77,Do_nothing] + _r0[77,Do_nothing,Do_nothing]
   + _r0[77,Do_nothing,Recover] = 0
 R263: - _f[77,R1] + _r0[77,R1,Do_nothing] + _r0[77,R1,Recover] = 0
 R264: - _f[77,R2] + _r0[77,R2,Do_nothing] + _r0[77,R2,Recover] = 0
 R265: - _f[77,R3] + _r0[77,R3,Do_nothing] + _r0[77,R3,Recover] = 0
 R266: - _f[95,Do_nothing] + _r0[95,Do_nothing,Do_nothing]
   + _r0[95,Do_nothing,Recover] = 0
 R267: - _f[95,R1] + _r0[95,R1,Do_nothing] + _r0[95,R1,Recover] = 0
 R268: - _f[95,R2] + _r0[95,R2,Do_nothing] + _r0[95,R2,Recover] = 0
 R269: - _f[95,R3] + _r0[95,R3,Do_nothing] + _r0[95,R3,Recover] = 0
 R270: - _f[98,Do_nothing] + _r0[98,Do_nothing,Do_nothing]
   + _r0[98,Do_nothing,Recover] = 0
 R271: - _f[98,R1] + _r0[98,R1,Do_nothing] + _r0[98,R1,Recover] = 0
 R272: - _f[98,R2] + _r0[98,R2,Do_nothing] + _r0[98,R2,Recover] = 0
 R273: - _f[98,R3] + _r0[98,R3,Do_nothing] + _r0[98,R3,Recover] = 0
 R274: - _f[53,Do_nothing] + _r0[53,Do_nothing,Do_nothing]
   + _r0[53,Do_nothing,Recover] = 0
 R275: - _f[53,R1] + _r0[53,R1,Do_nothing] + _r0[53,R1,Recover] = 0
 R276: - _f[53,R2] + _r0[53,R2,Do_nothing] + _r0[53,R2,Recover] = 0
 R277: - _f[53,R3] + _r0[53,R3,Do_nothing] + _r0[53,R3,Recover] = 0
 R278: - _f[65,Do_nothing] + _r0[65,Do_nothing,Do_nothing]
   + _r0[65,Do_nothing,Recover] = 0
 R279: - _f[65,R1] + _r0[65,R1,Do_nothing] + _r0[65,R1,Recover] = 0
 R280: - _f[65,R2] + _r0[65,R2,Do_nothing] + _r0[65,R2,Recover] = 0
 R281: - _f[65,R3] + _r0[65,R3,Do_nothing] + _r0[65,R3,Recover] = 0
 R282: - _f[63,Do_nothing] + _r0[63,Do_nothing,Do_nothing]
   + _r0[63,Do_nothing,Recover] = 0
 R283: - _f[63,R1] + _r0[63,R1,Do_nothing] + _r0[63,R1,Recover] = 0
 R284: - _f[63,R2] + _r0[63,R2,Do_nothing] + _r0[63,R2,Recover] = 0
 R285: - _f[63,R3] + _r0[63,R3,Do_nothing] + _r0[63,R3,Recover] = 0
 R286: - _f[14,Do_nothing] + _r0[14,Do_nothing,Do_nothing]
   + _r0[14,Do_nothing,Recover] = 0
 R287: - _f[14,R1] + _r0[14,R1,Do_nothing] + _r0[14,R1,Recover] = 0
 R288: - _f[14,R2] + _r0[14,R2,Do_nothing] + _r0[14,R2,Recover] = 0
 R289: - _f[14,R3] + _r0[14,R3,Do_nothing] + _r0[14,R3,Recover] = 0
 R290: - _f[23,Do_nothing] + _r0[23,Do_nothing,Do_nothing]
   + _r0[23,Do_nothing,Recover] = 0
 R291: - _f[23,R1] + _r0[23,R1,Do_nothing] + _r0[23,R1,Recover] = 0
 R292: - _f[23,R2] + _r0[23,R2,Do_nothing] + _r0[23,R2,Recover] = 0
 R293: - _f[23,R3] + _r0[23,R3,Do_nothing] + _r0[23,R3,Recover] = 0
 R294: - _f[93,Do_nothing] + _r0[93,Do_nothing,Do_nothing]
   + _r0[93,Do_nothing,Recover] = 0
 R295: - _f[93,R1] + _r0[93,R1,Do_nothing] + _r0[93,R1,Recover] = 0
 R296: - _f[93,R2] + _r0[93,R2,Do_nothing] + _r0[93,R2,Recover] = 0
 R297: - _f[93,R3] + _r0[93,R3,Do_nothing] + _r0[93,R3,Recover] = 0
 R298: - _f[28,Do_nothing] + _r0[28,Do_nothing,Do_nothing]
   + _r0[28,Do_nothing,Recover] = 0
 R299: - _f[28,R1] + _r0[28,R1,Do_nothing] + _r0[28,R1,Recover] = 0
 R300: - _f[28,R2] + _r0[28,R2,Do_nothing] + _r0[28,R2,Recover] = 0
 R301: - _f[28,R3] + _r0[28,R3,Do_nothing] + _r0[28,R3,Recover] = 0
 R302: - _f[2,Do_nothing] + _r0[2,Do_nothing,Do_nothing]
   + _r0[2,Do_nothing,Recover] = 0
 R303: - _f[2,R1] + _r0[2,R1,Do_nothing] + _r0[2,R1,Recover] = 0
 R304: - _f[2,R2] + _r0[2,R2,Do_nothing] + _r0[2,R2,Recover] = 0
 R305: - _f[2,R3] + _r0[2,R3,Do_nothing] + _r0[2,R3,Recover] = 0
 R306: - _f[81,Do_nothing] + _r0[81,Do_nothing,Do_nothing]
   + _r0[81,Do_nothing,Recover] = 0
 R307: - _f[81,R1] + _r0[81,R1,Do_nothing] + _r0[81,R1,Recover] = 0
 R308: - _f[81,R2] + _r0[81,R2,Do_nothing] + _r0[81,R2,Recover] = 0
 R309: - _f[81,R3] + _r0[81,R3,Do_nothing] + _r0[81,R3,Recover] = 0
 R310: - _f[45,Do_nothing] + _r0[45,Do_nothing,Do_nothing]
   + _r0[45,Do_nothing,Recover] = 0
 R311: - _f[45,R1] + _r0[45,R1,Do_nothing] + _r0[45,R1,Recover] = 0
 R312: - _f[45,R2] + _r0[45,R2,Do_nothing] + _r0[45,R2,Recover] = 0
 R313: - _f[45,R3] + _r0[45,R3,Do_nothing] + _r0[45,R3,Recover] = 0
 R314: - _f[88,Do_nothing] + _r0[88,Do_nothing,Do_nothing]
   + _r0[88,Do_nothing,Recover] = 0
 R315: - _f[88,R1] + _r0[88,R1,Do_nothing] + _r0[88,R1,Recover] = 0
 R316: - _f[88,R2] + _r0[88,R2,Do_nothing] + _r0[88,R2,Recover] = 0
 R317: - _f[88,R3] + _r0[88,R3,Do_nothing] + _r0[88,R3,Recover] = 0
 R318: - _f[89,Do_nothing] + _r0[89,Do_nothing,Do_nothing]
   + _r0[89,Do_nothing,Recover] = 0
 R319: - _f[89,R1] + _r0[89,R1,Do_nothing] + _r0[89,R1,Recover] = 0
 R320: - _f[89,R2] + _r0[89,R2,Do_nothing] + _r0[89,R2,Recover] = 0
 R321: - _f[89,R3] + _r0[89,R3,Do_nothing] + _r0[89,R3,Recover] = 0
 R322: - _f[25,Do_nothing] + _r0[25,Do_nothing,Do_nothing]
   + _r0[25,Do_nothing,Recover] = 0
 R323: - _f[25,R1] + _r0[25,R1,Do_nothing] + _r0[25,R1,Recover] = 0
 R324: - _f[25,R2] + _r0[25,R2,Do_nothing] + _r0[25,R2,Recover] = 0
 R325: - _f[25,R3] + _r0[25,R3,Do_nothing] + _r0[25,R3,Recover] = 0
 R326: - _f[91,Do_nothing] + _r0[91,Do_nothing,Do_nothing]
   + _r0[91,Do_nothing,Recover] = 0
 R327: - _f[91,R1] + _r0[91,R1,Do_nothing] + _r0[91,R1,Recover] = 0
 R328: - _f[91,R2] + _r0[91,R2,Do_nothing] + _r0[91,R2,Recover] = 0
 R329: - _f[91,R3] + _r0[91,R3,Do_nothing] + _r0[91,R3,Recover] = 0
 R330: - _f[43,Do_nothing] + _r0[43,Do_nothing,Do_nothing]
   + _r0[43,Do_nothing,Recover] = 0
 R331: - _f[43,R1] + _r0[43,R1,Do_nothing] + _r0[43,R1,Recover] = 0
 R332: - _f[43,R2] + _r0[43,R2,Do_nothing] + _r0[43,R2,Recover] = 0
 R333: - _f[43,R3] + _r0[43,R3,Do_nothing] + _r0[43,R3,Recover] = 0
 R334: - _f[76,Do_nothing] + _r0[76,Do_nothing,Do_nothing]
   + _r0[76,Do_nothing,Recover] = 0
 R335: - _f[76,R1] + _r0[76,R1,Do_nothing] + _r0[76,R1,Recover] = 0
 R336: - _f[76,R2] + _r0[76,R2,Do_nothing] + _r0[76,R2,Recover] = 0
 R337: - _f[76,R3] + _r0[76,R3,Do_nothing] + _r0[76,R3,Recover] = 0
 R338: - _f[42,Do_nothing] + _r0[42,Do_nothing,Do_nothing]
   + _r0[42,Do_nothing,Recover] = 0
 R339: - _f[42,R1] + _r0[42,R1,Do_nothing] + _r0[42,R1,Recover] = 0
 R340: - _f[42,R2] + _r0[42,R2,Do_nothing] + _r0[42,R2,Recover] = 0
 R341: - _f[42,R3] + _r0[42,R3,Do_nothing] + _r0[42,R3,Recover] = 0
 R342: - _f[7,Do_nothing] + _r0[7,Do_nothing,Do_nothing]
   + _r0[7,Do_nothing,Recover] = 0
 R343: - _f[7,R1] + _r0[7,R1,Do_nothing] + _r0[7,R1,Recover] = 0
 R344: - _f[7,R2] + _r0[7,R2,Do_nothing] + _r0[7,R2,Recover] = 0
 R345: - _f[7,R3] + _r0[7,R3,Do_nothing] + _r0[7,R3,Recover] = 0
 R346: - _f[17,Do_nothing] + _r0[17,Do_nothing,Do_nothing]
   + _r0[17,Do_nothing,Recover] = 0
 R347: - _f[17,R1] + _r0[17,R1,Do_nothing] + _r0[17,R1,Recover] = 0
 R348: - _f[17,R2] + _r0[17,R2,Do_nothing] + _r0[17,R2,Recover] = 0
 R349: - _f[17,R3] + _r0[17,R3,Do_nothing] + _r0[17,R3,Recover] = 0
 R350: - _f[37,Do_nothing] + _r0[37,Do_nothing,Do_nothing]
   + _r0[37,Do_nothing,Recover] = 0
 R351: - _f[37,R1] + _r0[37,R1,Do_nothing] + _r0[37,R1,Recover] = 0
 R352: - _f[37,R2] + _r0[37,R2,Do_nothing] + _r0[37,R2,Recover] = 0
 R353: - _f[37,R3] + _r0[37,R3,Do_nothing] + _r0[37,R3,Recover] = 0
 R354: - _f[26,Do_nothing] + _r0[26,Do_nothing,Do_nothing]
   + _r0[26,Do_nothing,Recover] = 0
 R355: - _f[26,R1] + _r0[26,R1,Do_nothing] + _r0[26,R1,Recover] = 0
 R356: - _f[26,R2] + _r0[26,R2,Do_nothing] + _r0[26,R2,Recover] = 0
 R357: - _f[26,R3] + _r0[26,R3,Do_nothing] + _r0[26,R3,Recover] = 0
 R358: - _f[9,Do_nothing] + _r0[9,Do_nothing,Do_nothing]
   + _r0[9,Do_nothing,Recover] = 0
 R359: - _f[9,R1] + _r0[9,R1,Do_nothing] + _r0[9,R1,Recover] = 0
 R360: - _f[9,R2] + _r0[9,R2,Do_nothing] + _r0[9,R2,Recover] = 0
 R361: - _f[9,R3] + _r0[9,R3,Do_nothing] + _r0[9,R3,Recover] = 0
 R362: - _f[60,Do_nothing] + _r0[60,Do_nothing,Do_nothing]
   + _r0[60,Do_nothing,Recover] = 0
 R363: - _f[60,R1] + _r0[60,R1,Do_nothing] + _r0[60,R1,Recover] = 0
 R364: - _f[60,R2] + _r0[60,R2,Do_nothing] + _r0[60,R2,Recover] = 0
 R365: - _f[60,R3] + _r0[60,R3,Do_nothing] + _r0[60,R3,Recover] = 0
 R366: - _f[16,Do_nothing] + _r0[16,Do_nothing,Do_nothing]
   + _r0[16,Do_nothing,Recover] = 0
 R367: - _f[16,R1] + _r0[16,R1,Do_nothing] + _r0[16,R1,Recover] = 0
 R368: - _f[16,R2] + _r0[16,R2,Do_nothing] + _r0[16,R2,Recover] = 0
 R369: - _f[16,R3] + _r0[16,R3,Do_nothing] + _r0[16,R3,Recover] = 0
 R370: - _f[83,Do_nothing] + _r0[83,Do_nothing,Do_nothing]
   + _r0[83,Do_nothing,Recover] = 0
 R371: - _f[83,R1] + _r0[83,R1,Do_nothing] + _r0[83,R1,Recover] = 0
 R372: - _f[83,R2] + _r0[83,R2,Do_nothing] + _r0[83,R2,Recover] = 0
 R373: - _f[83,R3] + _r0[83,R3,Do_nothing] + _r0[83,R3,Recover] = 0
 R374: - _f[13,Do_nothing] + _r0[13,Do_nothing,Do_nothing]
   + _r0[13,Do_nothing,Recover] = 0
 R375: - _f[13,R1] + _r0[13,R1,Do_nothing] + _r0[13,R1,Recover] = 0
 R376: - _f[13,R2] + _r0[13,R2,Do_nothing] + _r0[13,R2,Recover] = 0
 R377: - _f[13,R3] + _r0[13,R3,Do_nothing] + _r0[13,R3,Recover] = 0
 R378: - _f[67,Do_nothing] + _r0[67,Do_nothing,Do_nothing]
   + _r0[67,Do_nothing,Recover] = 0
 R379: - _f[67,R1] + _r0[67,R1,Do_nothing] + _r0[67,R1,Recover] = 0
 R380: - _f[67,R2] + _r0[67,R2,Do_nothing] + _r0[67,R2,Recover] = 0
 R381: - _f[67,R3] + _r0[67,R3,Do_nothing] + _r0[67,R3,Recover] = 0
 R382: - _f[59,Do_nothing] + _r0[59,Do_nothing,Do_nothing]
   + _r0[59,Do_nothing,Recover] = 0
 R383: - _f[59,R1] + _r0[59,R1,Do_nothing] + _r0[59,R1,Recover] = 0
 R384: - _f[59,R2] + _r0[59,R2,Do_nothing] + _r0[59,R2,Recover] = 0
 R385: - _f[59,R3] + _r0[59,R3,Do_nothing] + _r0[59,R3,Recover] = 0
 R386: - _f[58,Do_nothing] + _r0[58,Do_nothing,Do_nothing]
   + _r0[58,Do_nothing,Recover] = 0
 R387: - _f[58,R1] + _r0[58,R1,Do_nothing] + _r0[58,R1,Recover] = 0
 R388: - _f[58,R2] + _r0[58,R2,Do_nothing] + _r0[58,R2,Recover] = 0
 R389: - _f[58,R3] + _r0[58,R3,Do_nothing] + _r0[58,R3,Recover] = 0
 R390: - _f[4,Do_nothing] + _r0[4,Do_nothing,Do_nothing]
   + _r0[4,Do_nothing,Recover] = 0
 R391: - _f[4,R1] + _r0[4,R1,Do_nothing] + _r0[4,R1,Recover] = 0
 R392: - _f[4,R2] + _r0[4,R2,Do_nothing] + _r0[4,R2,Recover] = 0
 R393: - _f[4,R3] + _r0[4,R3,Do_nothing] + _r0[4,R3,Recover] = 0
 R394: - _f[68,Do_nothing] + _r0[68,Do_nothing,Do_nothing]
   + _r0[68,Do_nothing,Recover] = 0
 R395: - _f[68,R1] + _r0[68,R1,Do_nothing] + _r0[68,R1,Recover] = 0
 R396: - _f[68,R2] + _r0[68,R2,Do_nothing] + _r0[68,R2,Recover] = 0
 R397: - _f[68,R3] + _r0[68,R3,Do_nothing] + _r0[68,R3,Recover] = 0
 R398: - _f[87,Do_nothing] + _r0[87,Do_nothing,Do_nothing]
   + _r0[87,Do_nothing,Recover] = 0
 R399: - _f[87,R1] + _r0[87,R1,Do_nothing] + _r0[87,R1,Recover] = 0
 R400: - _f[87,R2] + _r0[87,R2,Do_nothing] + _r0[87,R2,Recover] = 0
 R401: - _f[87,R3] + _r0[87,R3,Do_nothing] + _r0[87,R3,Recover] = 0
 R402: - _f[40,Do_nothing] + _r0[40,Do_nothing,Do_nothing]
   + _r0[40,Do_nothing,Recover] = 0
 R403: - _f[40,R1] + _r0[40,R1,Do_nothing] + _r0[40,R1,Recover] = 0
 R404: - _f[40,R2] + _r0[40,R2,Do_nothing] + _r0[40,R2,Recover] = 0
 R405: - _f[40,R3] + _r0[40,R3,Do_nothing] + _r0[40,R3,Recover] = 0
 R406: - _f[30,Do_nothing] + _r0[30,Do_nothing,Do_nothing]
   + _r0[30,Do_nothing,Recover] = 0
 R407: - _f[30,R1] + _r0[30,R1,Do_nothing] + _r0[30,R1,Recover] = 0
 R408: - _f[30,R2] + _r0[30,R2,Do_nothing] + _r0[30,R2,Recover] = 0
 R409: - _f[30,R3] + _r0[30,R3,Do_nothing] + _r0[30,R3,Recover] = 0
 R410: - _f[27,Do_nothing] + _r0[27,Do_nothing,Do_nothing]
   + _r0[27,Do_nothing,Recover] = 0
 R411: - _f[27,R1] + _r0[27,R1,Do_nothing] + _r0[27,R1,Recover] = 0
 R412: - _f[27,R2] + _r0[27,R2,Do_nothing] + _r0[27,R2,Recover] = 0
 R413: - _f[27,R3] + _r0[27,R3,Do_nothing] + _r0[27,R3,Recover] = 0
 R414: - _f[73,Do_nothing] + _r0[73,Do_nothing,Do_nothing]
   + _r0[73,Do_nothing,Recover] = 0
 R415: - _f[73,R1] + _r0[73,R1,Do_nothing] + _r0[73,R1,Recover] = 0
 R416: - _f[73,R2] + _r0[73,R2,Do_nothing] + _r0[73,R2,Recover] = 0
 R417: - _f[73,R3] + _r0[73,R3,Do_nothing] + _r0[73,R3,Recover] = 0
 R418: - _f[84,Do_nothing] + _r0[84,Do_nothing,Do_nothing]
   + _r0[84,Do_nothing,Recover] = 0
 R419: - _f[84,R1] + _r0[84,R1,Do_nothing] + _r0[84,R1,Recover] = 0
 R420: - _f[84,R2] + _r0[84,R2,Do_nothing] + _r0[84,R2,Recover] = 0
 R421: - _f[84,R3] + _r0[84,R3,Do_nothing] + _r0[84,R3,Recover] = 0
 R422: - _f[35,Do_nothing] + _r0[35,Do_nothing,Do_nothing]
   + _r0[35,Do_nothing,Recover] = 0
 R423: - _f[35,R1] + _r0[35,R1,Do_nothing] + _r0[35,R1,Recover] = 0
 R424: - _f[35,R2] + _r0[35,R2,Do_nothing] + _r0[35,R2,Recover] = 0
 R425: - _f[35,R3] + _r0[35,R3,Do_nothing] + _r0[35,R3,Recover] = 0
 R426: - _f[46,Do_nothing] + _r0[46,Do_nothing,Do_nothing]
   + _r0[46,Do_nothing,Recover] = 0
 R427: - _f[46,R1] + _r0[46,R1,Do_nothing] + _r0[46,R1,Recover] = 0
 R428: - _f[46,R2] + _r0[46,R2,Do_nothing] + _r0[46,R2,Recover] = 0
 R429: - _f[46,R3] + _r0[46,R3,Do_nothing] + _r0[46,R3,Recover] = 0
 R430: - _f[39,Do_nothing] + _r0[39,Do_nothing,Do_nothing]
   + _r0[39,Do_nothing,Recover] = 0
 R431: - _f[39,R1] + _r0[39,R1,Do_nothing] + _r0[39,R1,Recover] = 0
 R432: - _f[39,R2] + _r0[39,R2,Do_nothing] + _r0[39,R2,Recover] = 0
 R433: - _f[39,R3] + _r0[39,R3,Do_nothing] + _r0[39,R3,Recover] = 0
 R434: - _f[51,Do_nothing] + _r0[51,Do_nothing,Do_nothing]
   + _r0[51,Do_nothing,Recover] = 0
 R435: - _f[51,R1] + _r0[51,R1,Do_nothing] + _r0[51,R1,Recover] = 0
 R436: - _f[51,R2] + _r0[51,R2,Do_nothing] + _r0[51,R2,Recover] = 0
 R437: - _f[51,R3] + _r0[51,R3,Do_nothing] + _r0[51,R3,Recover] = 0
 R438: - _f[52,Do_nothing] + _r0[52,Do_nothing,Do_nothing]
   + _r0[52,Do_nothing,Recover] = 0
 R439: - _f[52,R1] + _r0[52,R1,Do_nothing] + _r0[52,R1,Recover] = 0
 R440: - _f[52,R2] + _r0[52,R2,Do_nothing] + _r0[52,R2,Recover] = 0
 R441: - _f[52,R3] + _r0[52,R3,Do_nothing] + _r0[52,R3,Recover] = 0
 R442: - _f[6,Do_nothing] + _r0[6,Do_nothing,Do_nothing]
   + _r0[6,Do_nothing,Recover] = 0
 R443: - _f[6,R1] + _r0[6,R1,Do_nothing] + _r0[6,R1,Recover] = 0
 R444: - _f[6,R2] + _r0[6,R2,Do_nothing] + _r0[6,R2,Recover] = 0
 R445: - _f[6,R3] + _r0[6,R3,Do_nothing] + _r0[6,R3,Recover] = 0
 R446: - _f[96,Do_nothing] + _r0[96,Do_nothing,Do_nothing]
   + _r0[96,Do_nothing,Recover] = 0
 R447: - _f[96,R1] + _r0[96,R1,Do_nothing] + _r0[96,R1,Recover] = 0
 R448: - _f[96,R2] + _r0[96,R2,Do_nothing] + _r0[96,R2,Recover] = 0
 R449: - _f[96,R3] + _r0[96,R3,Do_nothing] + _r0[96,R3,Recover] = 0
 R450: - _f[74,Do_nothing] + _r0[74,Do_nothing,Do_nothing]
   + _r0[74,Do_nothing,Recover] = 0
 R451: - _f[74,R1] + _r0[74,R1,Do_nothing] + _r0[74,R1,Recover] = 0
 R452: - _f[74,R2] + _r0[74,R2,Do_nothing] + _r0[74,R2,Recover] = 0
 R453: - _f[74,R3] + _r0[74,R3,Do_nothing] + _r0[74,R3,Recover] = 0
 R454: - _f[86,Do_nothing] + _r0[86,Do_nothing,Do_nothing]
   + _r0[86,Do_nothing,Recover] = 0
 R455: - _f[86,R1] + _r0[86,R1,Do_nothing] + _r0[86,R1,Recover] = 0
 R456: - _f[86,R2] + _r0[86,R2,Do_nothing] + _r0[86,R2,Recover] = 0
 R457: - _f[86,R3] + _r0[86,R3,Do_nothing] + _r0[86,R3,Recover] = 0
 R458: - _f[29,Do_nothing] + _r0[29,Do_nothing,Do_nothing]
   + _r0[29,Do_nothing,Recover] = 0
 R459: - _f[29,R1] + _r0[29,R1,Do_nothing] + _r0[29,R1,Recover] = 0
 R460: - _f[29,R2] + _r0[29,R2,Do_nothing] + _r0[29,R2,Recover] = 0
 R461: - _f[29,R3] + _r0[29,R3,Do_nothing] + _r0[29,R3,Recover] = 0
 R462: - _f[18,Do_nothing] + _r0[18,Do_nothing,Do_nothing]
   + _r0[18,Do_nothing,Recover] = 0
 R463: - _f[18,R1] + _r0[18,R1,Do_nothing] + _r0[18,R1,Recover] = 0
 R464: - _f[18,R2] + _r0[18,R2,Do_nothing] + _r0[18,R2,Recover] = 0
 R465: - _f[18,R3] + _r0[18,R3,Do_nothing] + _r0[18,R3,Recover] = 0
 R466: - _f[24,Do_nothing] + _r0[24,Do_nothing,Do_nothing]
   + _r0[24,Do_nothing,Recover] = 0
 R467: - _f[24,R1] + _r0[24,R1,Do_nothing] + _r0[24,R1,Recover] = 0
 R468: - _f[24,R2] + _r0[24,R2,Do_nothing] + _r0[24,R2,Recover] = 0
 R469: - _f[24,R3] + _r0[24,R3,Do_nothing] + _r0[24,R3,Recover] = 0
 R470: - _f[94,Do_nothing] + _r0[94,Do_nothing,Do_nothing]
   + _r0[94,Do_nothing,Recover] = 0
 R471: - _f[94,R1] + _r0[94,R1,Do_nothing] + _r0[94,R1,Recover] = 0
 R472: - _f[94,R2] + _r0[94,R2,Do_nothing] + _r0[94,R2,Recover] = 0
 R473: - _f[94,R3] + _r0[94,R3,Do_nothing] + _r0[94,R3,Recover] = 0
 R474: - _f[33,Do_nothing] + _r0[33,Do_nothing,Do_nothing]
   + _r0[33,Do_nothing,Recover] = 0
 R475: - _f[33,R1] + _r0[33,R1,Do_nothing] + _r0[33,R1,Recover] = 0
 R476: - _f[33,R2] + _r0[33,R2,Do_nothing] + _r0[33,R2,Recover] = 0
 R477: - _f[33,R3] + _r0[33,R3,Do_nothing] + _r0[33,R3,Recover] = 0
 R478: - _f[12,Do_nothing] + _r0[12,Do_nothing,Do_nothing]
   + _r0[12,Do_nothing,Recover] = 0
 R479: - _f[12,R1] + _r0[12,R1,Do_nothing] + _r0[12,R1,Recover] = 0
 R480: - _f[12,R2] + _r0[12,R2,Do_nothing] + _r0[12,R2,Recover] = 0
 R481: - _f[12,R3] + _r0[12,R3,Do_nothing] + _r0[12,R3,Recover] = 0
 R482: - _f[57,Do_nothing] + _r0[57,Do_nothing,Do_nothing]
   + _r0[57,Do_nothing,Recover] = 0
 R483: - _f[57,R1] + _r0[57,R1,Do_nothing] + _r0[57,R1,Recover] = 0
 R484: - _f[57,R2] + _r0[57,R2,Do_nothing] + _r0[57,R2,Recover] = 0
 R485: - _f[57,R3] + _r0[57,R3,Do_nothing] + _r0[57,R3,Recover] = 0
 R486: - _f[97,Do_nothing] + _r0[97,Do_nothing,Do_nothing]
   + _r0[97,Do_nothing,Recover] = 0
 R487: - _f[97,R1] + _r0[97,R1,Do_nothing] + _r0[97,R1,Recover] = 0
 R488: - _f[97,R2] + _r0[97,R2,Do_nothing] + _r0[97,R2,Recover] = 0
 R489: - _f[97,R3] + _r0[97,R3,Do_nothing] + _r0[97,R3,Recover] = 0
 R490: - _f[69,Do_nothing] + _r0[69,Do_nothing,Do_nothing]
   + _r0[69,Do_nothing,Recover] = 0
 R491: - _f[69,R1] + _r0[69,R1,Do_nothing] + _r0[69,R1,Recover] = 0
 R492: - _f[69,R2] + _r0[69,R2,Do_nothing] + _r0[69,R2,Recover] = 0
 R493: - _f[69,R3] + _r0[69,R3,Do_nothing] + _r0[69,R3,Recover] = 0
 R494: - _f[3,Do_nothing] + _r0[3,Do_nothing,Do_nothing]
   + _r0[3,Do_nothing,Recover] = 0
 R495: - _f[3,R1] + _r0[3,R1,Do_nothing] + _r0[3,R1,Recover] = 0
 R496: - _f[3,R2] + _r0[3,R2,Do_nothing] + _r0[3,R2,Recover] = 0
 R497: - _f[3,R3] + _r0[3,R3,Do_nothing] + _r0[3,R3,Recover] = 0
 R498: - _f[66,Do_nothing] + _r0[66,Do_nothing,Do_nothing]
   + _r0[66,Do_nothing,Recover] = 0
 R499: - _f[66,R1] + _r0[66,R1,Do_nothing] + _r0[66,R1,Recover] = 0
 R500: - _f[66,R2] + _r0[66,R2,Do_nothing] + _r0[66,R2,Recover] = 0
 R501: - _f[66,R3] + _r0[66,R3,Do_nothing] + _r0[66,R3,Recover] = 0
 R502: _Theta - 803 _r1[32,Do_nothing,Do_nothing]
   - 795 _r1[32,Do_nothing,Recover] - 803 _r1[32,R1,Do_nothing]
   - 795 _r1[32,R1,Recover] - 709 _r1[32,R2,Do_nothing]
   - 607 _r1[32,R2,Recover] - 659 _r1[32,R3,Do_nothing]
   - 508 _r1[32,R3,Recover] - 135 _r1[5,Do_nothing,Do_nothing]
   - 134 _r1[5,Do_nothing,Recover] - 135 _r1[5,R1,Do_nothing]
   - 134 _r1[5,R1,Recover] - 119 _r1[5,R2,Do_nothing]
   - 102 _r1[5,R2,Recover] - 111 _r1[5,R3,Do_nothing]
   - 85 _r1[5,R3,Recover] - 526 _r1[92,Do_nothing,Do_nothing]
   - 521 _r1[92,Do_nothing,Recover] - 526 _r1[92,R1,Do_nothing]
   - 521 _r1[92,R1,Recover] - 465 _r1[92,R2,Do_nothing]
   - 398 _r1[92,R2,Recover] - 432 _r1[92,R3,Do_nothing]
   - 333 _r1[92,R3,Recover] - 432 _r1[21,Do_nothing,Do_nothing]
   - 428 _r1[21,Do_nothing,Recover] - 432 _r1[21,R1,Do_nothing]
   - 428 _r1[21,R1,Recover] - 381 _r1[21,R2,Do_nothing]
   - 327 _r1[21,R2,Recover] - 355 _r1[21,R3,Do_nothing]
   - 273 _r1[21,R3,Recover] - 474 _r1[56,Do_nothing,Do_nothing]
   - 469 _r1[56,Do_nothing,Recover] - 474 _r1[56,R1,Do_nothing]
   - 469 _r1[56,R1,Recover] - 418 _r1[56,R2,Do_nothing]
   - 358 _r1[56,R2,Recover] - 389 _r1[56,R3,Do_nothing]
   - 300 _r1[56,R3,Recover] - 955 _r1[80,Do_nothing,Do_nothing]
   - 947 _r1[80,Do_nothing,Recover] - 955 _r1[80,R1,Do_nothing]
   - 947 _r1[80,R1,Recover] - 843 _r1[80,R2,Do_nothing]
   - 723 _r1[80,R2,Recover] - 784 _r1[80,R3,Do_nothing]
   - 605 _r1[80,R3,Recover] - 691 _r1[55,Do_nothing,Do_nothing]
   - 684 _r1[55,Do_nothing,Recover] - 691 _r1[55,R1,Do_nothing]
   - 684 _r1[55,R1,Recover] - 610 _r1[55,R2,Do_nothing]
   - 522 _r1[55,R2,Recover] - 567 _r1[55,R3,Do_nothing]
   - 437 _r1[55,R3,Recover] - 657 _r1[8,Do_nothing,Do_nothing]
   - 651 _r1[8,Do_nothing,Recover] - 657 _r1[8,R1,Do_nothing]
   - 651 _r1[8,R1,Recover] - 580 _r1[8,R2,Do_nothing]
   - 497 _r1[8,R2,Recover] - 539 _r1[8,R3,Do_nothing]
   - 416 _r1[8,R3,Recover] - 1005 _r1[78,Do_nothing,Do_nothing]
   - 996 _r1[78,Do_nothing,Recover] - 1005 _r1[78,R1,Do_nothing]
   - 996 _r1[78,R1,Recover] - 887 _r1[78,R2,Do_nothing]
   - 760 _r1[78,R2,Recover] - 825 _r1[78,R3,Do_nothing]
   - 636 _r1[78,R3,Recover] - 1029 _r1[34,Do_nothing,Do_nothing]
   - 1019 _r1[34,Do_nothing,Recover] - 1029 _r1[34,R1,Do_nothing]
   - 1019 _r1[34,R1,Recover] - 908 _r1[34,R2,Do_nothing]
   - 778 _r1[34,R2,Recover] - 844 _r1[34,R3,Do_nothing]
   - 651 _r1[34,R3,Recover] - 466 _r1[22,Do_nothing,Do_nothing]
   - 462 _r1[22,Do_nothing,Recover] - 466 _r1[22,R1,Do_nothing]
   - 461 _r1[22,R1,Recover] - 411 _r1[22,R2,Do_nothing]
   - 352 _r1[22,R2,Recover] - 382 _r1[22,R3,Do_nothing]
   - 295 _r1[22,R3,Recover] - 794 _r1[54,Do_nothing,Do_nothing]
   - 787 _r1[54,Do_nothing,Recover] - 794 _r1[54,R1,Do_nothing]
   - 786 _r1[54,R1,Recover] - 701 _r1[54,R2,Do_nothing]
   - 600 _r1[54,R2,Recover] - 652 _r1[54,R3,Do_nothing]
   - 502 _r1[54,R3,Recover] - 639 _r1[75,Do_nothing,Do_nothing]
   - 633 _r1[75,Do_nothing,Recover] - 639 _r1[75,R1,Do_nothing]
   - 633 _r1[75,R1,Recover] - 564 _r1[75,R2,Do_nothing]
   - 483 _r1[75,R2,Recover] - 525 _r1[75,R3,Do_nothing]
   - 404 _r1[75,R3,Recover] - 831 _r1[72,Do_nothing,Do_nothing]
   - 824 _r1[72,Do_nothing,Recover] - 831 _r1[72,R1,Do_nothing]
   - 824 _r1[72,R1,Recover] - 734 _r1[72,R2,Do_nothing]
   - 629 _r1[72,R2,Recover] - 683 _r1[72,R3,Do_nothing]
   - 526 _r1[72,R3,Recover] - 699 _r1[53,Do_nothing,Do_nothing]
   - 692 _r1[53,Do_nothing,Recover] - 699 _r1[53,R1,Do_nothing]
   - 692 _r1[53,R1,Recover] - 617 _r1[53,R2,Do_nothing]
   - 528 _r1[53,R2,Recover] - 574 _r1[53,R3,Do_nothing]
   - 442 _r1[53,R3,Recover] - 875 _r1[2,Do_nothing,Do_nothing]
   - 867 _r1[2,Do_nothing,Recover] - 875 _r1[2,R1,Do_nothing]
   - 867 _r1[2,R1,Recover] - 775 _r1[2,R2,Do_nothing]
   - 667 _r1[2,R2,Recover] - 723 _r1[2,R3,Do_nothing]
   - 562 _r1[2,R3,Recover] - 237 _r1[89,Do_nothing,Do_nothing]
   - 235 _r1[89,Do_nothing,Recover] - 237 _r1[89,R1,Do_nothing]
   - 235 _r1[89,R1,Recover] - 209 _r1[89,R2,Do_nothing]
   - 179 _r1[89,R2,Recover] - 194 _r1[89,R3,Do_nothing]
   - 150 _r1[89,R3,Recover] - 1213 _r1[43,Do_nothing,Do_nothing]
   - 1202 _r1[43,Do_nothing,Recover] - 1213 _r1[43,R1,Do_nothing]
   - 1202 _r1[43,R1,Recover] - 1071 _r1[43,R2,Do_nothing]
   - 917 _r1[43,R2,Recover] - 996 _r1[43,R3,Do_nothing]
   - 768 _r1[43,R3,Recover] - 1405 _r1[9,Do_nothing,Do_nothing]
   - 1392 _r1[9,Do_nothing,Recover] - 1405 _r1[9,R1,Do_nothing]
   - 1392 _r1[9,R1,Recover] - 1240 _r1[9,R2,Do_nothing]
   - 1063 _r1[9,R2,Recover] - 1154 _r1[9,R3,Do_nothing]
   - 889 _r1[9,R3,Recover] - 923 _r1[40,Do_nothing,Do_nothing]
   - 914 _r1[40,Do_nothing,Recover] - 923 _r1[40,R1,Do_nothing]
   - 914 _r1[40,R1,Recover] - 814 _r1[40,R2,Do_nothing]
   - 698 _r1[40,R2,Recover] - 757 _r1[40,R3,Do_nothing]
   - 584 _r1[40,R3,Recover] - 348 _r1[27,Do_nothing,Do_nothing]
   - 345 _r1[27,Do_nothing,Recover] - 348 _r1[27,R1,Do_nothing]
   - 345 _r1[27,R1,Recover] - 307 _r1[27,R2,Do_nothing]
   - 263 _r1[27,R2,Recover] - 286 _r1[27,R3,Do_nothing]
   - 220 _r1[27,R3,Recover] - 1181 _r1[46,Do_nothing,Do_nothing]
   - 1171 _r1[46,Do_nothing,Recover] - 1181 _r1[46,R1,Do_nothing]
   - 1170 _r1[46,R1,Recover] - 1043 _r1[46,R2,Do_nothing]
   - 894 _r1[46,R2,Recover] - 970 _r1[46,R3,Do_nothing]
   - 748 _r1[46,R3,Recover] - 441 _r1[6,Do_nothing,Do_nothing]
   - 437 _r1[6,Do_nothing,Recover] - 441 _r1[6,R1,Do_nothing]
   - 437 _r1[6,R1,Recover] - 389 _r1[6,R2,Do_nothing]
   - 334 _r1[6,R2,Recover] - 362 _r1[6,R3,Do_nothing]
   - 279 _r1[6,R3,Recover] - 339 _r1[96,Do_nothing,Do_nothing]
   - 336 _r1[96,Do_nothing,Recover] - 339 _r1[96,R1,Do_nothing]
   - 336 _r1[96,R1,Recover] - 299 _r1[96,R2,Do_nothing]
   - 256 _r1[96,R2,Recover] - 278 _r1[96,R3,Do_nothing]
   - 214 _r1[96,R3,Recover] - 195 _r1[29,Do_nothing,Do_nothing]
   - 193 _r1[29,Do_nothing,Recover] - 195 _r1[29,R1,Do_nothing]
   - 193 _r1[29,R1,Recover] - 172 _r1[29,R2,Do_nothing]
   - 148 _r1[29,R2,Recover] - 160 _r1[29,R3,Do_nothing]
   - 124 _r1[29,R3,Recover] >= 0
 R503: 184945 _f[32,R1] + 232438 _f[32,R2] + 417219 _f[32,R3]
   + 77012 _f[70,R1] + 96789 _f[70,R2] + 173733 _f[70,R3] + 42150 _f[5,R1]
   + 52975 _f[5,R2] + 95089 _f[5,R3] + 83838 _f[85,R1] + 105366 _f[85,R2]
   + 189129 _f[85,R3] + 19605 _f[71,R1] + 24639 _f[71,R2] + 44226 _f[71,R3]
   + 112508 _f[92,R1] + 141402 _f[92,R2] + 253810 _f[92,R3]
   + 158010 _f[38,R1] + 198586 _f[38,R2] + 356456 _f[38,R3]
   + 227840 _f[21,R1] + 286345 _f[21,R2] + 513985 _f[21,R3]
   + 134812 _f[41,R1] + 169431 _f[41,R2] + 304131 _f[41,R3]
   + 125910 _f[56,R1] + 158240 _f[56,R2] + 284041 _f[56,R3]
   + 111072 _f[80,R1] + 139593 _f[80,R2] + 250567 _f[80,R3]
   + 88700 _f[62,R1] + 121094 _f[62,R2] + 217364 _f[62,R3]
   + 274315 _f[20,R1] + 344756 _f[20,R2] + 618825 _f[20,R3]
   + 96147 _f[10,R1] + 120838 _f[10,R2] + 216902 _f[10,R3]
   + 141199 _f[55,R1] + 177452 _f[55,R2] + 318522 _f[55,R3]
   + 90202 _f[36,R1] + 113364 _f[36,R2] + 203487 _f[36,R3]
   + 136344 _f[44,R1] + 171354 _f[44,R2] + 307578 _f[44,R3]
   + 69011 _f[48,R1] + 86734 _f[48,R2] + 155683 _f[48,R3] + 216395 _f[8,R1]
   + 271965 _f[8,R2] + 488172 _f[8,R3] + 108443 _f[15,R1]
   + 136293 _f[15,R2] + 244641 _f[15,R3] + 28666 _f[61,R1]
   + 36027 _f[61,R2] + 64668 _f[61,R3] + 139534 _f[78,R1]
   + 175362 _f[78,R2] + 314775 _f[78,R3] + 163384 _f[34,R1]
   + 205344 _f[34,R2] + 368587 _f[34,R3] + 76581 _f[100,R1]
   + 96246 _f[100,R2] + 172762 _f[100,R3] + 182173 _f[90,R1]
   + 228954 _f[90,R2] + 410968 _f[90,R3] + 253339 _f[64,R1]
   + 318395 _f[64,R2] + 571513 _f[64,R3] + 84551 _f[79,R1]
   + 106262 _f[79,R2] + 190738 _f[79,R3] + 104250 _f[1,R1]
   + 131021 _f[1,R2] + 235179 _f[1,R3] + 86007 _f[22,R1] + 108093 _f[22,R2]
   + 194022 _f[22,R3] + 149381 _f[31,R1] + 187743 _f[31,R2]
   + 336989 _f[31,R3] + 232480 _f[54,R1] + 292181 _f[54,R2]
   + 524455 _f[54,R3] + 76086 _f[82,R1] + 95624 _f[82,R2]
   + 171642 _f[82,R3] + 104028 _f[11,R1] + 130741 _f[11,R2]
   + 234681 _f[11,R3] + 139646 _f[75,R1] + 175502 _f[75,R2]
   + 315022 _f[75,R3] + 588581 _f[50,R1] + 739723 _f[50,R2]
   + 1.327789e+06 _f[50,R3] + 106626 _f[19,R1] + 134005 _f[19,R2]
   + 240534 _f[19,R3] + 113557 _f[99,R1] + 142720 _f[99,R2]
   + 256178 _f[99,R3] + 129156 _f[72,R1] + 162324 _f[72,R2]
   + 291373 _f[72,R3] + 88597 _f[47,R1] + 111344 _f[47,R2]
   + 199866 _f[47,R3] + 123902 _f[49,R1] + 155714 _f[49,R2]
   + 279506 _f[49,R3] + 21608 _f[77,R1] + 27157 _f[77,R2] + 48746 _f[77,R3]
   + 124290 _f[95,R1] + 156209 _f[95,R2] + 280390 _f[95,R3]
   + 112314 _f[98,R1] + 141154 _f[98,R2] + 253370 _f[98,R3]
   + 109154 _f[53,R1] + 137182 _f[53,R2] + 246237 _f[53,R3]
   + 89855 _f[65,R1] + 112926 _f[65,R2] + 202701 _f[65,R3]
   + 15197 _f[63,R1] + 19099 _f[63,R2] + 34283 _f[63,R3] + 82896 _f[14,R1]
   + 104182 _f[14,R2] + 187009 _f[14,R3] + 171814 _f[23,R1]
   + 215936 _f[23,R2] + 387595 _f[23,R3] + 82531 _f[93,R1]
   + 103723 _f[93,R2] + 186181 _f[93,R3] + 139133 _f[28,R1]
   + 174860 _f[28,R2] + 313866 _f[28,R3] + 266795 _f[2,R1]
   + 350290 _f[2,R2] + 628762 _f[2,R3] + 65502 _f[81,R1] + 82323 _f[81,R2]
   + 147768 _f[81,R3] + 95888 _f[45,R1] + 120511 _f[45,R2]
   + 216314 _f[45,R3] + 141611 _f[88,R1] + 177974 _f[88,R2]
   + 319460 _f[88,R3] + 48010 _f[89,R1] + 60340 _f[89,R2]
   + 108310 _f[89,R3] + 36212 _f[25,R1] + 45512 _f[25,R2] + 81689 _f[25,R3]
   + 130813 _f[91,R1] + 164403 _f[91,R2] + 295099 _f[91,R3]
   + 459976 _f[43,R1] + 578092 _f[43,R2] + 1.037663e+06 _f[43,R3]
   + 172013 _f[76,R1] + 216187 _f[76,R2] + 388050 _f[76,R3]
   + 170159 _f[42,R1] + 213855 _f[42,R2] + 383865 _f[42,R3]
   + 67468 _f[7,R1] + 84797 _f[7,R2] + 152208 _f[7,R3] + 124758 _f[17,R1]
   + 156793 _f[17,R2] + 281442 _f[17,R3] + 55522 _f[37,R1]
   + 69779 _f[37,R2] + 125252 _f[37,R3] + 169821 _f[26,R1]
   + 213431 _f[26,R2] + 383103 _f[26,R3] + 393339 _f[9,R1]
   + 494348 _f[9,R2] + 887341 _f[9,R3] + 120732 _f[60,R1]
   + 151737 _f[60,R2] + 272364 _f[60,R3] + 152263 _f[16,R1]
   + 191364 _f[16,R2] + 343496 _f[16,R3] + 99675 _f[83,R1]
   + 125272 _f[83,R2] + 224860 _f[83,R3] + 193418 _f[13,R1]
   + 243086 _f[13,R2] + 436338 _f[13,R3] + 127707 _f[67,R1]
   + 201287 _f[67,R2] + 361306 _f[67,R3] + 176600 _f[59,R1]
   + 221950 _f[59,R2] + 398395 _f[59,R3] + 100371 _f[58,R1]
   + 126143 _f[58,R2] + 226422 _f[58,R3] + 206083 _f[4,R1]
   + 259003 _f[4,R2] + 464904 _f[4,R3] + 50733 _f[68,R1] + 63759 _f[68,R2]
   + 114447 _f[68,R3] + 45515 _f[87,R1] + 57203 _f[87,R2]
   + 102678 _f[87,R3] + 205522 _f[40,R1] + 258296 _f[40,R2]
   + 463642 _f[40,R3] + 92954 _f[30,R1] + 116822 _f[30,R2]
   + 209693 _f[30,R3] + 67691 _f[27,R1] + 85070 _f[27,R2]
   + 152699 _f[27,R3] + 878240 _f[73,R1] + 1.103759e+06 _f[73,R2]
   + 1.981229e+06 _f[73,R3] + 200216 _f[84,R1] + 251629 _f[84,R2]
   + 451670 _f[84,R3] + 52951 _f[35,R1] + 66548 _f[35,R2]
   + 119449 _f[35,R3] + 191369 _f[46,R1] + 240512 _f[46,R2]
   + 431707 _f[46,R3] + 165478 _f[39,R1] + 207972 _f[39,R2]
   + 373305 _f[39,R3] + 890435 _f[51,R1] + 1.119091e+06 _f[51,R2]
   + 2.008749e+06 _f[51,R3] + 210223 _f[52,R1] + 264209 _f[52,R2]
   + 474249 _f[52,R3] + 100228 _f[6,R1] + 125965 _f[6,R2] + 226105 _f[6,R3]
   + 156136 _f[96,R1] + 196233 _f[96,R2] + 352232 _f[96,R3]
   + 26744 _f[74,R1] + 33612 _f[74,R2] + 60332 _f[74,R3] + 22382 _f[86,R1]
   + 28129 _f[86,R2] + 50491 _f[86,R3] + 59968 _f[29,R1] + 75367 _f[29,R2]
   + 135283 _f[29,R3] + 78936 _f[18,R1] + 99206 _f[18,R2]
   + 178071 _f[18,R3] + 118132 _f[24,R1] + 148466 _f[24,R2]
   + 266493 _f[24,R3] + 46689 _f[94,R1] + 58681 _f[94,R2]
   + 105329 _f[94,R3] + 63134 _f[33,R1] + 79346 _f[33,R2]
   + 142423 _f[33,R3] + 33082 _f[12,R1] + 41577 _f[12,R2] + 74630 _f[12,R3]
   + 221052 _f[57,R1] + 277816 _f[57,R2] + 498675 _f[57,R3]
   + 49387 _f[97,R1] + 62067 _f[97,R2] + 111410 _f[97,R3] + 41699 _f[69,R1]
   + 52406 _f[69,R2] + 94068 _f[69,R3] + 11993 _f[3,R1] + 15072 _f[3,R2]
   + 27055 _f[3,R3] + 36687 _f[66,R1] + 46107 _f[66,R2] + 82761 _f[66,R3]
   + 1.6250306e+07 _r1[32,Do_nothing,Recover]
   + 1.6245457e+07 _r1[32,R1,Recover] + 1.1655471e+07 _r1[32,R2,Recover]
   + 7.096177e+06 _r1[32,R3,Recover]
   + 8.395899e+06 _r1[70,Do_nothing,Recover]
   + 8.393394e+06 _r1[70,R1,Recover] + 6.021927e+06 _r1[70,R2,Recover]
   + 3.666318e+06 _r1[70,R3,Recover]
   + 3.681994e+06 _r1[5,Do_nothing,Recover]
   + 3.680895e+06 _r1[5,R1,Recover] + 2.640896e+06 _r1[5,R2,Recover]
   + 1.607852e+06 _r1[5,R3,Recover]
   + 2.0876574e+07 _r1[85,Do_nothing,Recover]
   + 2.0870346e+07 _r1[85,R1,Recover] + 1.4973644e+07 _r1[85,R2,Recover]
   + 9.116374e+06 _r1[85,R3,Recover]
   + 1.926495e+06 _r1[71,Do_nothing,Recover]
   + 1.92592e+06 _r1[71,R1,Recover] + 1.381771e+06 _r1[71,R2,Recover]
   + 841261 _r1[71,R3,Recover] + 1.3779239e+07 _r1[92,Do_nothing,Recover]
   + 1.3775128e+07 _r1[92,R1,Recover] + 9.883107e+06 _r1[92,R2,Recover]
   + 6.017113e+06 _r1[92,R3,Recover]
   + 1.6740245e+07 _r1[38,Do_nothing,Recover]
   + 1.6735251e+07 _r1[38,R1,Recover] + 1.2006878e+07 _r1[38,R2,Recover]
   + 7.310124e+06 _r1[38,R3,Recover]
   + 1.697441e+07 _r1[21,Do_nothing,Recover]
   + 1.6969346e+07 _r1[21,R1,Recover] + 1.2174832e+07 _r1[21,R2,Recover]
   + 7.412379e+06 _r1[21,R3,Recover]
   + 1.3483233e+07 _r1[41,Do_nothing,Recover]
   + 1.347921e+07 _r1[41,R1,Recover] + 9.670798e+06 _r1[41,R2,Recover]
   + 5.887853e+06 _r1[41,R3,Recover]
   + 1.5514733e+07 _r1[56,Do_nothing,Recover]
   + 1.5510104e+07 _r1[56,R1,Recover] + 1.1127884e+07 _r1[56,R2,Recover]
   + 6.774967e+06 _r1[56,R3,Recover]
   + 3.3460116e+07 _r1[80,Do_nothing,Recover]
   + 3.3450133e+07 _r1[80,R1,Recover] + 2.3999142e+07 _r1[80,R2,Recover]
   + 1.461135e+07 _r1[80,R3,Recover]
   + 1.3793197e+07 _r1[62,Do_nothing,Recover]
   + 1.3789358e+07 _r1[62,R1,Recover] + 1.0056913e+07 _r1[62,R2,Recover]
   + 6.381213e+06 _r1[62,R3,Recover]
   + 2.6317165e+07 _r1[20,Do_nothing,Recover]
   + 2.6309313e+07 _r1[20,R1,Recover] + 1.8875887e+07 _r1[20,R2,Recover]
   + 1.1492169e+07 _r1[20,R3,Recover]
   + 2.3539324e+07 _r1[10,Do_nothing,Recover]
   + 2.3532301e+07 _r1[10,R1,Recover] + 1.6883492e+07 _r1[10,R2,Recover]
   + 1.0279143e+07 _r1[10,R3,Recover]
   + 2.0671468e+07 _r1[55,Do_nothing,Recover]
   + 2.06653e+07 _r1[55,R1,Recover] + 1.4826533e+07 _r1[55,R2,Recover]
   + 9.026808e+06 _r1[55,R3,Recover]
   + 1.026872e+07 _r1[36,Do_nothing,Recover]
   + 1.0265657e+07 _r1[36,R1,Recover] + 7.365201e+06 _r1[36,R2,Recover]
   + 4.484141e+06 _r1[36,R3,Recover]
   + 2.556896e+06 _r1[44,Do_nothing,Recover]
   + 2.556133e+06 _r1[44,R1,Recover] + 1.833924e+06 _r1[44,R2,Recover]
   + 1.116544e+06 _r1[44,R3,Recover]
   + 1.677052e+07 _r1[48,Do_nothing,Recover]
   + 1.6765517e+07 _r1[48,R1,Recover] + 1.2028593e+07 _r1[48,R2,Recover]
   + 7.323344e+06 _r1[48,R3,Recover]
   + 1.6412617e+07 _r1[8,Do_nothing,Recover]
   + 1.640772e+07 _r1[8,R1,Recover] + 1.1771888e+07 _r1[8,R2,Recover]
   + 7.167055e+06 _r1[8,R3,Recover]
   + 2.4310998e+07 _r1[15,Do_nothing,Recover]
   + 2.4303745e+07 _r1[15,R1,Recover] + 1.7436972e+07 _r1[15,R2,Recover]
   + 1.0616117e+07 _r1[15,R3,Recover]
   + 1.664175e+06 _r1[61,Do_nothing,Recover]
   + 1.663678e+06 _r1[61,R1,Recover] + 1.193623e+06 _r1[61,R2,Recover]
   + 726711 _r1[61,R3,Recover] + 2.1085368e+07 _r1[78,Do_nothing,Recover]
   + 2.1079077e+07 _r1[78,R1,Recover] + 1.5123401e+07 _r1[78,R2,Recover]
   + 9.20755e+06 _r1[78,R3,Recover]
   + 5.4495756e+07 _r1[34,Do_nothing,Recover]
   + 5.4479497e+07 _r1[34,R1,Recover] + 3.9086876e+07 _r1[34,R2,Recover]
   + 2.3797185e+07 _r1[34,R3,Recover]
   + 1.3726875e+07 _r1[100,Do_nothing,Recover]
   + 1.3722779e+07 _r1[100,R1,Recover] + 9.845549e+06 _r1[100,R2,Recover]
   + 5.994246e+06 _r1[100,R3,Recover]
   + 3.5521762e+07 _r1[90,Do_nothing,Recover]
   + 3.5511165e+07 _r1[90,R1,Recover] + 2.547785e+07 _r1[90,R2,Recover]
   + 1.5511629e+07 _r1[90,R3,Recover]
   + 1.471966e+07 _r1[64,Do_nothing,Recover]
   + 1.4715269e+07 _r1[64,R1,Recover] + 1.0557621e+07 _r1[64,R2,Recover]
   + 6.427775e+06 _r1[64,R3,Recover]
   + 3.5959429e+07 _r1[79,Do_nothing,Recover]
   + 3.5948701e+07 _r1[79,R1,Recover] + 2.5791766e+07 _r1[79,R2,Recover]
   + 1.5702749e+07 _r1[79,R3,Recover]
   + 3.911387e+07 _r1[1,Do_nothing,Recover]
   + 3.9102201e+07 _r1[1,R1,Recover] + 2.8054276e+07 _r1[1,R2,Recover]
   + 1.7080229e+07 _r1[1,R3,Recover]
   + 1.2599232e+07 _r1[22,Do_nothing,Recover]
   + 1.2595474e+07 _r1[22,R1,Recover] + 9.036752e+06 _r1[22,R2,Recover]
   + 5.501828e+06 _r1[22,R3,Recover]
   + 2.6726068e+07 _r1[31,Do_nothing,Recover]
   + 2.6718095e+07 _r1[31,R1,Recover] + 1.9169172e+07 _r1[31,R2,Recover]
   + 1.1670729e+07 _r1[31,R3,Recover]
   + 1.9445107e+07 _r1[54,Do_nothing,Recover]
   + 1.9439306e+07 _r1[54,R1,Recover] + 1.394693e+07 _r1[54,R2,Recover]
   + 8.491282e+06 _r1[54,R3,Recover]
   + 1.3534085e+07 _r1[82,Do_nothing,Recover]
   + 1.3530047e+07 _r1[82,R1,Recover] + 9.707272e+06 _r1[82,R2,Recover]
   + 5.910059e+06 _r1[82,R3,Recover]
   + 1.1946244e+07 _r1[11,Do_nothing,Recover]
   + 1.194268e+07 _r1[11,R1,Recover] + 8.568399e+06 _r1[11,R2,Recover]
   + 5.216681e+06 _r1[11,R3,Recover]
   + 1.9233548e+07 _r1[75,Do_nothing,Recover]
   + 1.922781e+07 _r1[75,R1,Recover] + 1.379519e+07 _r1[75,R2,Recover]
   + 8.398898e+06 _r1[75,R3,Recover]
   + 7.9890053e+07 _r1[50,Do_nothing,Recover]
   + 7.9866218e+07 _r1[50,R1,Recover] + 5.730084e+07 _r1[50,R2,Recover]
   + 3.4886357e+07 _r1[50,R3,Recover]
   + 3.7772492e+07 _r1[19,Do_nothing,Recover]
   + 3.7761223e+07 _r1[19,R1,Recover] + 2.7092178e+07 _r1[19,R2,Recover]
   + 1.6494477e+07 _r1[19,R3,Recover]
   + 1.401132e+07 _r1[99,Do_nothing,Recover]
   + 1.400714e+07 _r1[99,R1,Recover] + 1.0049567e+07 _r1[99,R2,Recover]
   + 6.118458e+06 _r1[99,R3,Recover]
   + 2.3511427e+07 _r1[72,Do_nothing,Recover]
   + 2.3504412e+07 _r1[72,R1,Recover] + 1.6863483e+07 _r1[72,R2,Recover]
   + 1.0266961e+07 _r1[72,R3,Recover]
   + 1.0302942e+07 _r1[47,Do_nothing,Recover]
   + 1.0299868e+07 _r1[47,R1,Recover] + 7.389746e+06 _r1[47,R2,Recover]
   + 4.499084e+06 _r1[47,R3,Recover]
   + 3.5335076e+07 _r1[49,Do_nothing,Recover]
   + 3.5324534e+07 _r1[49,R1,Recover] + 2.534395e+07 _r1[49,R2,Recover]
   + 1.5430107e+07 _r1[49,R3,Recover]
   + 2.275954e+06 _r1[77,Do_nothing,Recover]
   + 2.275275e+06 _r1[77,R1,Recover] + 1.632419e+06 _r1[77,R2,Recover]
   + 993863 _r1[77,R3,Recover] + 3.0225358e+07 _r1[95,Do_nothing,Recover]
   + 3.0216341e+07 _r1[95,R1,Recover] + 2.1679025e+07 _r1[95,R2,Recover]
   + 1.3198797e+07 _r1[95,R3,Recover]
   + 8.36295e+06 _r1[98,Do_nothing,Recover]
   + 8.360455e+06 _r1[98,R1,Recover] + 5.998294e+06 _r1[98,R2,Recover]
   + 3.65193e+06 _r1[98,R3,Recover]
   + 2.399908e+07 _r1[53,Do_nothing,Recover]
   + 2.399192e+07 _r1[53,R1,Recover] + 1.721325e+07 _r1[53,R2,Recover]
   + 1.0479909e+07 _r1[53,R3,Recover]
   + 1.3037471e+07 _r1[65,Do_nothing,Recover]
   + 1.3033581e+07 _r1[65,R1,Recover] + 9.351077e+06 _r1[65,R2,Recover]
   + 5.693198e+06 _r1[65,R3,Recover]
   + 5.549173e+06 _r1[63,Do_nothing,Recover]
   + 5.547518e+06 _r1[63,R1,Recover] + 3.980124e+06 _r1[63,R2,Recover]
   + 2.423211e+06 _r1[63,R3,Recover]
   + 2.8131869e+07 _r1[14,Do_nothing,Recover]
   + 2.8123476e+07 _r1[14,R1,Recover] + 2.0177477e+07 _r1[14,R2,Recover]
   + 1.2284613e+07 _r1[14,R3,Recover]
   + 1.4888056e+07 _r1[23,Do_nothing,Recover]
   + 1.4883614e+07 _r1[23,R1,Recover] + 1.0678402e+07 _r1[23,R2,Recover]
   + 6.50131e+06 _r1[23,R3,Recover]
   + 1.4395332e+07 _r1[93,Do_nothing,Recover]
   + 1.4391038e+07 _r1[93,R1,Recover] + 1.0324998e+07 _r1[93,R2,Recover]
   + 6.286148e+06 _r1[93,R3,Recover]
   + 1.3562886e+07 _r1[28,Do_nothing,Recover]
   + 1.3558839e+07 _r1[28,R1,Recover] + 9.727929e+06 _r1[28,R2,Recover]
   + 5.922636e+06 _r1[28,R3,Recover]
   + 2.598607e+07 _r1[2,Do_nothing,Recover]
   + 2.5978609e+07 _r1[2,R1,Recover] + 1.8811219e+07 _r1[2,R2,Recover]
   + 1.1725292e+07 _r1[2,R3,Recover]
   + 1.1245096e+07 _r1[81,Do_nothing,Recover]
   + 1.1241741e+07 _r1[81,R1,Recover] + 8.065503e+06 _r1[81,R2,Recover]
   + 4.910504e+06 _r1[81,R3,Recover]
   + 2.0483269e+07 _r1[45,Do_nothing,Recover]
   + 2.0477158e+07 _r1[45,R1,Recover] + 1.4691548e+07 _r1[45,R2,Recover]
   + 8.944626e+06 _r1[45,R3,Recover]
   + 1.7600312e+07 _r1[88,Do_nothing,Recover]
   + 1.7595061e+07 _r1[88,R1,Recover] + 1.2623758e+07 _r1[88,R2,Recover]
   + 7.685697e+06 _r1[88,R3,Recover]
   + 7.828206e+06 _r1[89,Do_nothing,Recover]
   + 7.82587e+06 _r1[89,R1,Recover] + 5.614751e+06 _r1[89,R2,Recover]
   + 3.418418e+06 _r1[89,R3,Recover]
   + 5.800854e+06 _r1[25,Do_nothing,Recover]
   + 5.799124e+06 _r1[25,R1,Recover] + 4.160641e+06 _r1[25,R2,Recover]
   + 2.533115e+06 _r1[25,R3,Recover]
   + 2.0526875e+07 _r1[91,Do_nothing,Recover]
   + 2.0520751e+07 _r1[91,R1,Recover] + 1.4722824e+07 _r1[91,R2,Recover]
   + 8.963668e+06 _r1[91,R3,Recover]
   + 4.9197858e+07 _r1[43,Do_nothing,Recover]
   + 4.918318e+07 _r1[43,R1,Recover] + 3.5286979e+07 _r1[43,R2,Recover]
   + 2.1483701e+07 _r1[43,R3,Recover]
   + 1.09461187e+08 _r1[76,Do_nothing,Recover]
   + 1.0942853e+08 _r1[76,R1,Recover] + 7.8510625e+07 _r1[76,R2,Recover]
   + 4.7799468e+07 _r1[76,R3,Recover]
   + 3.3561304e+07 _r1[42,Do_nothing,Recover]
   + 3.3551292e+07 _r1[42,R1,Recover] + 2.4071719e+07 _r1[42,R2,Recover]
   + 1.4655537e+07 _r1[42,R3,Recover]
   + 1.7648288e+07 _r1[7,Do_nothing,Recover]
   + 1.7643023e+07 _r1[7,R1,Recover] + 1.2658169e+07 _r1[7,R2,Recover]
   + 7.706648e+06 _r1[7,R3,Recover]
   + 1.9181497e+07 _r1[17,Do_nothing,Recover]
   + 1.9175774e+07 _r1[17,R1,Recover] + 1.3757856e+07 _r1[17,R2,Recover]
   + 8.376168e+06 _r1[17,R3,Recover]
   + 1.0390541e+07 _r1[37,Do_nothing,Recover]
   + 1.0387441e+07 _r1[37,R1,Recover] + 7.452577e+06 _r1[37,R2,Recover]
   + 4.537337e+06 _r1[37,R3,Recover]
   + 2.0918355e+07 _r1[26,Do_nothing,Recover]
   + 2.0912114e+07 _r1[26,R1,Recover] + 1.5003611e+07 _r1[26,R2,Recover]
   + 9.134619e+06 _r1[26,R3,Recover]
   + 5.8801993e+07 _r1[9,Do_nothing,Recover]
   + 5.878445e+07 _r1[9,R1,Recover] + 4.2175508e+07 _r1[9,R2,Recover]
   + 2.5677631e+07 _r1[9,R3,Recover]
   + 1.5265061e+07 _r1[60,Do_nothing,Recover]
   + 1.5260507e+07 _r1[60,R1,Recover] + 1.0948808e+07 _r1[60,R2,Recover]
   + 6.665941e+06 _r1[60,R3,Recover]
   + 1.5422118e+07 _r1[16,Do_nothing,Recover]
   + 1.5417517e+07 _r1[16,R1,Recover] + 1.1061456e+07 _r1[16,R2,Recover]
   + 6.734524e+06 _r1[16,R3,Recover]
   + 1.3153079e+07 _r1[83,Do_nothing,Recover]
   + 1.3149154e+07 _r1[83,R1,Recover] + 9.433996e+06 _r1[83,R2,Recover]
   + 5.743681e+06 _r1[83,R3,Recover]
   + 7.2047115e+07 _r1[13,Do_nothing,Recover]
   + 7.202562e+07 _r1[13,R1,Recover] + 5.1675522e+07 _r1[13,R2,Recover]
   + 3.1461506e+07 _r1[13,R3,Recover]
   + 5.5945218e+07 _r1[67,Do_nothing,Recover]
   + 5.5931138e+07 _r1[67,R1,Recover] + 4.1676653e+07 _r1[67,R2,Recover]
   + 2.7818295e+07 _r1[67,R3,Recover]
   + 2.1147689e+07 _r1[59,Do_nothing,Recover]
   + 2.114138e+07 _r1[59,R1,Recover] + 1.51681e+07 _r1[59,R2,Recover]
   + 9.234765e+06 _r1[59,R3,Recover]
   + 1.4856103e+07 _r1[58,Do_nothing,Recover]
   + 1.485167e+07 _r1[58,R1,Recover] + 1.0655484e+07 _r1[58,R2,Recover]
   + 6.487357e+06 _r1[58,R3,Recover]
   + 2.2818833e+07 _r1[4,Do_nothing,Recover]
   + 2.2812026e+07 _r1[4,R1,Recover] + 1.6366723e+07 _r1[4,R2,Recover]
   + 9.964519e+06 _r1[4,R3,Recover]
   + 1.8591493e+07 _r1[68,Do_nothing,Recover]
   + 1.8585947e+07 _r1[68,R1,Recover] + 1.3334679e+07 _r1[68,R2,Recover]
   + 8.118526e+06 _r1[68,R3,Recover]
   + 4.873506e+06 _r1[87,Do_nothing,Recover]
   + 4.872052e+06 _r1[87,R1,Recover] + 3.495504e+06 _r1[87,R2,Recover]
   + 2.128161e+06 _r1[87,R3,Recover]
   + 2.1259203e+07 _r1[40,Do_nothing,Recover]
   + 2.1252861e+07 _r1[40,R1,Recover] + 1.5248084e+07 _r1[40,R2,Recover]
   + 9.28346e+06 _r1[40,R3,Recover]
   + 9.148287e+06 _r1[30,Do_nothing,Recover]
   + 9.145558e+06 _r1[30,R1,Recover] + 6.561575e+06 _r1[30,R2,Recover]
   + 3.994871e+06 _r1[30,R3,Recover]
   + 1.2238711e+07 _r1[27,Do_nothing,Recover]
   + 1.2235059e+07 _r1[27,R1,Recover] + 8.778169e+06 _r1[27,R2,Recover]
   + 5.344395e+06 _r1[27,R3,Recover]
   + 3.1009091e+07 _r1[73,Do_nothing,Recover]
   + 3.099984e+07 _r1[73,R1,Recover] + 2.2241154e+07 _r1[73,R2,Recover]
   + 1.3541038e+07 _r1[73,R3,Recover]
   + 2.2537319e+07 _r1[84,Do_nothing,Recover]
   + 2.2530595e+07 _r1[84,R1,Recover] + 1.6164808e+07 _r1[84,R2,Recover]
   + 9.841588e+06 _r1[84,R3,Recover]
   + 4.660656e+06 _r1[35,Do_nothing,Recover]
   + 4.659265e+06 _r1[35,R1,Recover] + 3.342838e+06 _r1[35,R2,Recover]
   + 2.035213e+06 _r1[35,R3,Recover]
   + 8.4516505e+07 _r1[46,Do_nothing,Recover]
   + 8.4491291e+07 _r1[46,R1,Recover] + 6.0619146e+07 _r1[46,R2,Recover]
   + 3.6906634e+07 _r1[46,R3,Recover]
   + 2.0923204e+07 _r1[39,Do_nothing,Recover]
   + 2.0916962e+07 _r1[39,R1,Recover] + 1.500709e+07 _r1[39,R2,Recover]
   + 9.136737e+06 _r1[39,R3,Recover]
   + 1.9955349e+08 _r1[51,Do_nothing,Recover]
   + 1.99493955e+08 _r1[51,R1,Recover] + 1.43128991e+08 _r1[51,R2,Recover]
   + 8.714094e+07 _r1[51,R3,Recover]
   + 7.9890532e+07 _r1[52,Do_nothing,Recover]
   + 7.9866697e+07 _r1[52,R1,Recover] + 5.7301184e+07 _r1[52,R2,Recover]
   + 3.4886566e+07 _r1[52,R3,Recover]
   + 1.6577583e+07 _r1[6,Do_nothing,Recover]
   + 1.6572637e+07 _r1[6,R1,Recover] + 1.1890209e+07 _r1[6,R2,Recover]
   + 7.239093e+06 _r1[6,R3,Recover]
   + 1.2827184e+07 _r1[96,Do_nothing,Recover]
   + 1.2823357e+07 _r1[96,R1,Recover] + 9.200249e+06 _r1[96,R2,Recover]
   + 5.60137e+06 _r1[96,R3,Recover]
   + 1.0275635e+07 _r1[74,Do_nothing,Recover]
   + 1.0272569e+07 _r1[74,R1,Recover] + 7.37016e+06 _r1[74,R2,Recover]
   + 4.48716e+06 _r1[74,R3,Recover] + 880257 _r1[86,Do_nothing,Recover]
   + 879995 _r1[86,R1,Recover] + 631361 _r1[86,R2,Recover]
   + 384390 _r1[86,R3,Recover] + 5.748969e+06 _r1[29,Do_nothing,Recover]
   + 5.747254e+06 _r1[29,R1,Recover] + 4.123427e+06 _r1[29,R2,Recover]
   + 2.510458e+06 _r1[29,R3,Recover]
   + 5.511117e+06 _r1[18,Do_nothing,Recover]
   + 5.509472e+06 _r1[18,R1,Recover] + 3.952828e+06 _r1[18,R2,Recover]
   + 2.406592e+06 _r1[18,R3,Recover]
   + 2.1123609e+07 _r1[24,Do_nothing,Recover]
   + 2.1117307e+07 _r1[24,R1,Recover] + 1.5150829e+07 _r1[24,R2,Recover]
   + 9.224249e+06 _r1[24,R3,Recover]
   + 4.315511e+06 _r1[94,Do_nothing,Recover]
   + 4.314224e+06 _r1[94,R1,Recover] + 3.095284e+06 _r1[94,R2,Recover]
   + 1.884496e+06 _r1[94,R3,Recover]
   + 5.1681149e+07 _r1[33,Do_nothing,Recover]
   + 5.166573e+07 _r1[33,R1,Recover] + 3.706811e+07 _r1[33,R2,Recover]
   + 2.2568104e+07 _r1[33,R3,Recover]
   + 1.378918e+06 _r1[12,Do_nothing,Recover]
   + 1.378506e+06 _r1[12,R1,Recover] + 989024 _r1[12,R2,Recover]
   + 602145 _r1[12,R3,Recover] + 4.2913776e+07 _r1[57,Do_nothing,Recover]
   + 4.2900973e+07 _r1[57,R1,Recover] + 3.0779745e+07 _r1[57,R2,Recover]
   + 1.8739571e+07 _r1[57,R3,Recover]
   + 2.092752e+06 _r1[97,Do_nothing,Recover]
   + 2.092128e+06 _r1[97,R1,Recover] + 1.501018e+06 _r1[97,R2,Recover]
   + 913862 _r1[97,R3,Recover] + 1.461207e+06 _r1[69,Do_nothing,Recover]
   + 1.460771e+06 _r1[69,R1,Recover] + 1.048045e+06 _r1[69,R2,Recover]
   + 638079 _r1[69,R3,Recover] + 1.8390093e+07 _r1[3,Do_nothing,Recover]
   + 1.8384607e+07 _r1[3,R1,Recover] + 1.3190225e+07 _r1[3,R2,Recover]
   + 8.030579e+06 _r1[3,R3,Recover] + 396329 _r1[66,Do_nothing,Recover]
   + 396210 _r1[66,R1,Recover] + 284265 _r1[66,R2,Recover]
   + 173069 _r1[66,R3,Recover] <= 1.5e+07
 R504: - _f[32,Do_nothing] + _r1[32,Do_nothing,Do_nothing]
   + _r1[32,Do_nothing,Recover] = 0
 R505: - _f[32,R1] + _r1[32,R1,Do_nothing] + _r1[32,R1,Recover] = 0
 R506: - _f[32,R2] + _r1[32,R2,Do_nothing] + _r1[32,R2,Recover] = 0
 R507: - _f[32,R3] + _r1[32,R3,Do_nothing] + _r1[32,R3,Recover] = 0
 R508: - _f[70,Do_nothing] + _r1[70,Do_nothing,Do_nothing]
   + _r1[70,Do_nothing,Recover] = 0
 R509: - _f[70,R1] + _r1[70,R1,Do_nothing] + _r1[70,R1,Recover] = 0
 R510: - _f[70,R2] + _r1[70,R2,Do_nothing] + _r1[70,R2,Recover] = 0
 R511: - _f[70,R3] + _r1[70,R3,Do_nothing] + _r1[70,R3,Recover] = 0
 R512: - _f[5,Do_nothing] + _r1[5,Do_nothing,Do_nothing]
   + _r1[5,Do_nothing,Recover] = 0
 R513: - _f[5,R1] + _r1[5,R1,Do_nothing] + _r1[5,R1,Recover] = 0
 R514: - _f[5,R2] + _r1[5,R2,Do_nothing] + _r1[5,R2,Recover] = 0
 R515: - _f[5,R3] + _r1[5,R3,Do_nothing] + _r1[5,R3,Recover] = 0
 R516: - _f[85,Do_nothing] + _r1[85,Do_nothing,Do_nothing]
   + _r1[85,Do_nothing,Recover] = 0
 R517: - _f[85,R1] + _r1[85,R1,Do_nothing] + _r1[85,R1,Recover] = 0
 R518: - _f[85,R2] + _r1[85,R2,Do_nothing] + _r1[85,R2,Recover] = 0
 R519: - _f[85,R3] + _r1[85,R3,Do_nothing] + _r1[85,R3,Recover] = 0
 R520: - _f[71,Do_nothing] + _r1[71,Do_nothing,Do_nothing]
   + _r1[71,Do_nothing,Recover] = 0
 R521: - _f[71,R1] + _r1[71,R1,Do_nothing] + _r1[71,R1,Recover] = 0
 R522: - _f[71,R2] + _r1[71,R2,Do_nothing] + _r1[71,R2,Recover] = 0
 R523: - _f[71,R3] + _r1[71,R3,Do_nothing] + _r1[71,R3,Recover] = 0
 R524: - _f[92,Do_nothing] + _r1[92,Do_nothing,Do_nothing]
   + _r1[92,Do_nothing,Recover] = 0
 R525: - _f[92,R1] + _r1[92,R1,Do_nothing] + _r1[92,R1,Recover] = 0
 R526: - _f[92,R2] + _r1[92,R2,Do_nothing] + _r1[92,R2,Recover] = 0
 R527: - _f[92,R3] + _r1[92,R3,Do_nothing] + _r1[92,R3,Recover] = 0
 R528: - _f[38,Do_nothing] + _r1[38,Do_nothing,Do_nothing]
   + _r1[38,Do_nothing,Recover] = 0
 R529: - _f[38,R1] + _r1[38,R1,Do_nothing] + _r1[38,R1,Recover] = 0
 R530: - _f[38,R2] + _r1[38,R2,Do_nothing] + _r1[38,R2,Recover] = 0
 R531: - _f[38,R3] + _r1[38,R3,Do_nothing] + _r1[38,R3,Recover] = 0
 R532: - _f[21,Do_nothing] + _r1[21,Do_nothing,Do_nothing]
   + _r1[21,Do_nothing,Recover] = 0
 R533: - _f[21,R1] + _r1[21,R1,Do_nothing] + _r1[21,R1,Recover] = 0
 R534: - _f[21,R2] + _r1[21,R2,Do_nothing] + _r1[21,R2,Recover] = 0
 R535: - _f[21,R3] + _r1[21,R3,Do_nothing] + _r1[21,R3,Recover] = 0
 R536: - _f[41,Do_nothing] + _r1[41,Do_nothing,Do_nothing]
   + _r1[41,Do_nothing,Recover] = 0
 R537: - _f[41,R1] + _r1[41,R1,Do_nothing] + _r1[41,R1,Recover] = 0
 R538: - _f[41,R2] + _r1[41,R2,Do_nothing] + _r1[41,R2,Recover] = 0
 R539: - _f[41,R3] + _r1[41,R3,Do_nothing] + _r1[41,R3,Recover] = 0
 R540: - _f[56,Do_nothing] + _r1[56,Do_nothing,Do_nothing]
   + _r1[56,Do_nothing,Recover] = 0
 R541: - _f[56,R1] + _r1[56,R1,Do_nothing] + _r1[56,R1,Recover] = 0
 R542: - _f[56,R2] + _r1[56,R2,Do_nothing] + _r1[56,R2,Recover] = 0
 R543: - _f[56,R3] + _r1[56,R3,Do_nothing] + _r1[56,R3,Recover] = 0
 R544: - _f[80,Do_nothing] + _r1[80,Do_nothing,Do_nothing]
   + _r1[80,Do_nothing,Recover] = 0
 R545: - _f[80,R1] + _r1[80,R1,Do_nothing] + _r1[80,R1,Recover] = 0
 R546: - _f[80,R2] + _r1[80,R2,Do_nothing] + _r1[80,R2,Recover] = 0
 R547: - _f[80,R3] + _r1[80,R3,Do_nothing] + _r1[80,R3,Recover] = 0
 R548: - _f[62,Do_nothing] + _r1[62,Do_nothing,Do_nothing]
   + _r1[62,Do_nothing,Recover] = 0
 R549: - _f[62,R1] + _r1[62,R1,Do_nothing] + _r1[62,R1,Recover] = 0
 R550: - _f[62,R2] + _r1[62,R2,Do_nothing] + _r1[62,R2,Recover] = 0
 R551: - _f[62,R3] + _r1[62,R3,Do_nothing] + _r1[62,R3,Recover] = 0
 R552: - _f[20,Do_nothing] + _r1[20,Do_nothing,Do_nothing]
   + _r1[20,Do_nothing,Recover] = 0
 R553: - _f[20,R1] + _r1[20,R1,Do_nothing] + _r1[20,R1,Recover] = 0
 R554: - _f[20,R2] + _r1[20,R2,Do_nothing] + _r1[20,R2,Recover] = 0
 R555: - _f[20,R3] + _r1[20,R3,Do_nothing] + _r1[20,R3,Recover] = 0
 R556: - _f[10,Do_nothing] + _r1[10,Do_nothing,Do_nothing]
   + _r1[10,Do_nothing,Recover] = 0
 R557: - _f[10,R1] + _r1[10,R1,Do_nothing] + _r1[10,R1,Recover] = 0
 R558: - _f[10,R2] + _r1[10,R2,Do_nothing] + _r1[10,R2,Recover] = 0
 R559: - _f[10,R3] + _r1[10,R3,Do_nothing] + _r1[10,R3,Recover] = 0
 R560: - _f[55,Do_nothing] + _r1[55,Do_nothing,Do_nothing]
   + _r1[55,Do_nothing,Recover] = 0
 R561: - _f[55,R1] + _r1[55,R1,Do_nothing] + _r1[55,R1,Recover] = 0
 R562: - _f[55,R2] + _r1[55,R2,Do_nothing] + _r1[55,R2,Recover] = 0
 R563: - _f[55,R3] + _r1[55,R3,Do_nothing] + _r1[55,R3,Recover] = 0
 R564: - _f[36,Do_nothing] + _r1[36,Do_nothing,Do_nothing]
   + _r1[36,Do_nothing,Recover] = 0
 R565: - _f[36,R1] + _r1[36,R1,Do_nothing] + _r1[36,R1,Recover] = 0
 R566: - _f[36,R2] + _r1[36,R2,Do_nothing] + _r1[36,R2,Recover] = 0
 R567: - _f[36,R3] + _r1[36,R3,Do_nothing] + _r1[36,R3,Recover] = 0
 R568: - _f[44,Do_nothing] + _r1[44,Do_nothing,Do_nothing]
   + _r1[44,Do_nothing,Recover] = 0
 R569: - _f[44,R1] + _r1[44,R1,Do_nothing] + _r1[44,R1,Recover] = 0
 R570: - _f[44,R2] + _r1[44,R2,Do_nothing] + _r1[44,R2,Recover] = 0
 R571: - _f[44,R3] + _r1[44,R3,Do_nothing] + _r1[44,R3,Recover] = 0
 R572: - _f[48,Do_nothing] + _r1[48,Do_nothing,Do_nothing]
   + _r1[48,Do_nothing,Recover] = 0
 R573: - _f[48,R1] + _r1[48,R1,Do_nothing] + _r1[48,R1,Recover] = 0
 R574: - _f[48,R2] + _r1[48,R2,Do_nothing] + _r1[48,R2,Recover] = 0
 R575: - _f[48,R3] + _r1[48,R3,Do_nothing] + _r1[48,R3,Recover] = 0
 R576: - _f[8,Do_nothing] + _r1[8,Do_nothing,Do_nothing]
   + _r1[8,Do_nothing,Recover] = 0
 R577: - _f[8,R1] + _r1[8,R1,Do_nothing] + _r1[8,R1,Recover] = 0
 R578: - _f[8,R2] + _r1[8,R2,Do_nothing] + _r1[8,R2,Recover] = 0
 R579: - _f[8,R3] + _r1[8,R3,Do_nothing] + _r1[8,R3,Recover] = 0
 R580: - _f[15,Do_nothing] + _r1[15,Do_nothing,Do_nothing]
   + _r1[15,Do_nothing,Recover] = 0
 R581: - _f[15,R1] + _r1[15,R1,Do_nothing] + _r1[15,R1,Recover] = 0
 R582: - _f[15,R2] + _r1[15,R2,Do_nothing] + _r1[15,R2,Recover] = 0
 R583: - _f[15,R3] + _r1[15,R3,Do_nothing] + _r1[15,R3,Recover] = 0
 R584: - _f[61,Do_nothing] + _r1[61,Do_nothing,Do_nothing]
   + _r1[61,Do_nothing,Recover] = 0
 R585: - _f[61,R1] + _r1[61,R1,Do_nothing] + _r1[61,R1,Recover] = 0
 R586: - _f[61,R2] + _r1[61,R2,Do_nothing] + _r1[61,R2,Recover] = 0
 R587: - _f[61,R3] + _r1[61,R3,Do_nothing] + _r1[61,R3,Recover] = 0
 R588: - _f[78,Do_nothing] + _r1[78,Do_nothing,Do_nothing]
   + _r1[78,Do_nothing,Recover] = 0
 R589: - _f[78,R1] + _r1[78,R1,Do_nothing] + _r1[78,R1,Recover] = 0
 R590: - _f[78,R2] + _r1[78,R2,Do_nothing] + _r1[78,R2,Recover] = 0
 R591: - _f[78,R3] + _r1[78,R3,Do_nothing] + _r1[78,R3,Recover] = 0
 R592: - _f[34,Do_nothing] + _r1[34,Do_nothing,Do_nothing]
   + _r1[34,Do_nothing,Recover] = 0
 R593: - _f[34,R1] + _r1[34,R1,Do_nothing] + _r1[34,R1,Recover] = 0
 R594: - _f[34,R2] + _r1[34,R2,Do_nothing] + _r1[34,R2,Recover] = 0
 R595: - _f[34,R3] + _r1[34,R3,Do_nothing] + _r1[34,R3,Recover] = 0
 R596: - _f[100,Do_nothing] + _r1[100,Do_nothing,Do_nothing]
   + _r1[100,Do_nothing,Recover] = 0
 R597: - _f[100,R1] + _r1[100,R1,Do_nothing] + _r1[100,R1,Recover] = 0
 R598: - _f[100,R2] + _r1[100,R2,Do_nothing] + _r1[100,R2,Recover] = 0
 R599: - _f[100,R3] + _r1[100,R3,Do_nothing] + _r1[100,R3,Recover] = 0
 R600: - _f[90,Do_nothing] + _r1[90,Do_nothing,Do_nothing]
   + _r1[90,Do_nothing,Recover] = 0
 R601: - _f[90,R1] + _r1[90,R1,Do_nothing] + _r1[90,R1,Recover] = 0
 R602: - _f[90,R2] + _r1[90,R2,Do_nothing] + _r1[90,R2,Recover] = 0
 R603: - _f[90,R3] + _r1[90,R3,Do_nothing] + _r1[90,R3,Recover] = 0
 R604: - _f[64,Do_nothing] + _r1[64,Do_nothing,Do_nothing]
   + _r1[64,Do_nothing,Recover] = 0
 R605: - _f[64,R1] + _r1[64,R1,Do_nothing] + _r1[64,R1,Recover] = 0
 R606: - _f[64,R2] + _r1[64,R2,Do_nothing] + _r1[64,R2,Recover] = 0
 R607: - _f[64,R3] + _r1[64,R3,Do_nothing] + _r1[64,R3,Recover] = 0
 R608: - _f[79,Do_nothing] + _r1[79,Do_nothing,Do_nothing]
   + _r1[79,Do_nothing,Recover] = 0
 R609: - _f[79,R1] + _r1[79,R1,Do_nothing] + _r1[79,R1,Recover] = 0
 R610: - _f[79,R2] + _r1[79,R2,Do_nothing] + _r1[79,R2,Recover] = 0
 R611: - _f[79,R3] + _r1[79,R3,Do_nothing] + _r1[79,R3,Recover] = 0
 R612: - _f[1,Do_nothing] + _r1[1,Do_nothing,Do_nothing]
   + _r1[1,Do_nothing,Recover] = 0
 R613: - _f[1,R1] + _r1[1,R1,Do_nothing] + _r1[1,R1,Recover] = 0
 R614: - _f[1,R2] + _r1[1,R2,Do_nothing] + _r1[1,R2,Recover] = 0
 R615: - _f[1,R3] + _r1[1,R3,Do_nothing] + _r1[1,R3,Recover] = 0
 R616: - _f[22,Do_nothing] + _r1[22,Do_nothing,Do_nothing]
   + _r1[22,Do_nothing,Recover] = 0
 R617: - _f[22,R1] + _r1[22,R1,Do_nothing] + _r1[22,R1,Recover] = 0
 R618: - _f[22,R2] + _r1[22,R2,Do_nothing] + _r1[22,R2,Recover] = 0
 R619: - _f[22,R3] + _r1[22,R3,Do_nothing] + _r1[22,R3,Recover] = 0
 R620: - _f[31,Do_nothing] + _r1[31,Do_nothing,Do_nothing]
   + _r1[31,Do_nothing,Recover] = 0
 R621: - _f[31,R1] + _r1[31,R1,Do_nothing] + _r1[31,R1,Recover] = 0
 R622: - _f[31,R2] + _r1[31,R2,Do_nothing] + _r1[31,R2,Recover] = 0
 R623: - _f[31,R3] + _r1[31,R3,Do_nothing] + _r1[31,R3,Recover] = 0
 R624: - _f[54,Do_nothing] + _r1[54,Do_nothing,Do_nothing]
   + _r1[54,Do_nothing,Recover] = 0
 R625: - _f[54,R1] + _r1[54,R1,Do_nothing] + _r1[54,R1,Recover] = 0
 R626: - _f[54,R2] + _r1[54,R2,Do_nothing] + _r1[54,R2,Recover] = 0
 R627: - _f[54,R3] + _r1[54,R3,Do_nothing] + _r1[54,R3,Recover] = 0
 R628: - _f[82,Do_nothing] + _r1[82,Do_nothing,Do_nothing]
   + _r1[82,Do_nothing,Recover] = 0
 R629: - _f[82,R1] + _r1[82,R1,Do_nothing] + _r1[82,R1,Recover] = 0
 R630: - _f[82,R2] + _r1[82,R2,Do_nothing] + _r1[82,R2,Recover] = 0
 R631: - _f[82,R3] + _r1[82,R3,Do_nothing] + _r1[82,R3,Recover] = 0
 R632: - _f[11,Do_nothing] + _r1[11,Do_nothing,Do_nothing]
   + _r1[11,Do_nothing,Recover] = 0
 R633: - _f[11,R1] + _r1[11,R1,Do_nothing] + _r1[11,R1,Recover] = 0
 R634: - _f[11,R2] + _r1[11,R2,Do_nothing] + _r1[11,R2,Recover] = 0
 R635: - _f[11,R3] + _r1[11,R3,Do_nothing] + _r1[11,R3,Recover] = 0
 R636: - _f[75,Do_nothing] + _r1[75,Do_nothing,Do_nothing]
   + _r1[75,Do_nothing,Recover] = 0
 R637: - _f[75,R1] + _r1[75,R1,Do_nothing] + _r1[75,R1,Recover] = 0
 R638: - _f[75,R2] + _r1[75,R2,Do_nothing] + _r1[75,R2,Recover] = 0
 R639: - _f[75,R3] + _r1[75,R3,Do_nothing] + _r1[75,R3,Recover] = 0
 R640: - _f[50,Do_nothing] + _r1[50,Do_nothing,Do_nothing]
   + _r1[50,Do_nothing,Recover] = 0
 R641: - _f[50,R1] + _r1[50,R1,Do_nothing] + _r1[50,R1,Recover] = 0
 R642: - _f[50,R2] + _r1[50,R2,Do_nothing] + _r1[50,R2,Recover] = 0
 R643: - _f[50,R3] + _r1[50,R3,Do_nothing] + _r1[50,R3,Recover] = 0
 R644: - _f[19,Do_nothing] + _r1[19,Do_nothing,Do_nothing]
   + _r1[19,Do_nothing,Recover] = 0
 R645: - _f[19,R1] + _r1[19,R1,Do_nothing] + _r1[19,R1,Recover] = 0
 R646: - _f[19,R2] + _r1[19,R2,Do_nothing] + _r1[19,R2,Recover] = 0
 R647: - _f[19,R3] + _r1[19,R3,Do_nothing] + _r1[19,R3,Recover] = 0
 R648: - _f[99,Do_nothing] + _r1[99,Do_nothing,Do_nothing]
   + _r1[99,Do_nothing,Recover] = 0
 R649: - _f[99,R1] + _r1[99,R1,Do_nothing] + _r1[99,R1,Recover] = 0
 R650: - _f[99,R2] + _r1[99,R2,Do_nothing] + _r1[99,R2,Recover] = 0
 R651: - _f[99,R3] + _r1[99,R3,Do_nothing] + _r1[99,R3,Recover] = 0
 R652: - _f[72,Do_nothing] + _r1[72,Do_nothing,Do_nothing]
   + _r1[72,Do_nothing,Recover] = 0
 R653: - _f[72,R1] + _r1[72,R1,Do_nothing] + _r1[72,R1,Recover] = 0
 R654: - _f[72,R2] + _r1[72,R2,Do_nothing] + _r1[72,R2,Recover] = 0
 R655: - _f[72,R3] + _r1[72,R3,Do_nothing] + _r1[72,R3,Recover] = 0
 R656: - _f[47,Do_nothing] + _r1[47,Do_nothing,Do_nothing]
   + _r1[47,Do_nothing,Recover] = 0
 R657: - _f[47,R1] + _r1[47,R1,Do_nothing] + _r1[47,R1,Recover] = 0
 R658: - _f[47,R2] + _r1[47,R2,Do_nothing] + _r1[47,R2,Recover] = 0
 R659: - _f[47,R3] + _r1[47,R3,Do_nothing] + _r1[47,R3,Recover] = 0
 R660: - _f[49,Do_nothing] + _r1[49,Do_nothing,Do_nothing]
   + _r1[49,Do_nothing,Recover] = 0
 R661: - _f[49,R1] + _r1[49,R1,Do_nothing] + _r1[49,R1,Recover] = 0
 R662: - _f[49,R2] + _r1[49,R2,Do_nothing] + _r1[49,R2,Recover] = 0
 R663: - _f[49,R3] + _r1[49,R3,Do_nothing] + _r1[49,R3,Recover] = 0
 R664: - _f[77,Do_nothing] + _r1[77,Do_nothing,Do_nothing]
   + _r1[77,Do_nothing,Recover] = 0
 R665: - _f[77,R1] + _r1[77,R1,Do_nothing] + _r1[77,R1,Recover] = 0
 R666: - _f[77,R2] + _r1[77,R2,Do_nothing] + _r1[77,R2,Recover] = 0
 R667: - _f[77,R3] + _r1[77,R3,Do_nothing] + _r1[77,R3,Recover] = 0
 R668: - _f[95,Do_nothing] + _r1[95,Do_nothing,Do_nothing]
   + _r1[95,Do_nothing,Recover] = 0
 R669: - _f[95,R1] + _r1[95,R1,Do_nothing] + _r1[95,R1,Recover] = 0
 R670: - _f[95,R2] + _r1[95,R2,Do_nothing] + _r1[95,R2,Recover] = 0
 R671: - _f[95,R3] + _r1[95,R3,Do_nothing] + _r1[95,R3,Recover] = 0
 R672: - _f[98,Do_nothing] + _r1[98,Do_nothing,Do_nothing]
   + _r1[98,Do_nothing,Recover] = 0
 R673: - _f[98,R1] + _r1[98,R1,Do_nothing] + _r1[98,R1,Recover] = 0
 R674: - _f[98,R2] + _r1[98,R2,Do_nothing] + _r1[98,R2,Recover] = 0
 R675: - _f[98,R3] + _r1[98,R3,Do_nothing] + _r1[98,R3,Recover] = 0
 R676: - _f[53,Do_nothing] + _r1[53,Do_nothing,Do_nothing]
   + _r1[53,Do_nothing,Recover] = 0
 R677: - _f[53,R1] + _r1[53,R1,Do_nothing] + _r1[53,R1,Recover] = 0
 R678: - _f[53,R2] + _r1[53,R2,Do_nothing] + _r1[53,R2,Recover] = 0
 R679: - _f[53,R3] + _r1[53,R3,Do_nothing] + _r1[53,R3,Recover] = 0
 R680: - _f[65,Do_nothing] + _r1[65,Do_nothing,Do_nothing]
   + _r1[65,Do_nothing,Recover] = 0
 R681: - _f[65,R1] + _r1[65,R1,Do_nothing] + _r1[65,R1,Recover] = 0
 R682: - _f[65,R2] + _r1[65,R2,Do_nothing] + _r1[65,R2,Recover] = 0
 R683: - _f[65,R3] + _r1[65,R3,Do_nothing] + _r1[65,R3,Recover] = 0
 R684: - _f[63,Do_nothing] + _r1[63,Do_nothing,Do_nothing]
   + _r1[63,Do_nothing,Recover] = 0
 R685: - _f[63,R1] + _r1[63,R1,Do_nothing] + _r1[63,R1,Recover] = 0
 R686: - _f[63,R2] + _r1[63,R2,Do_nothing] + _r1[63,R2,Recover] = 0
 R687: - _f[63,R3] + _r1[63,R3,Do_nothing] + _r1[63,R3,Recover] = 0
 R688: - _f[14,Do_nothing] + _r1[14,Do_nothing,Do_nothing]
   + _r1[14,Do_nothing,Recover] = 0
 R689: - _f[14,R1] + _r1[14,R1,Do_nothing] + _r1[14,R1,Recover] = 0
 R690: - _f[14,R2] + _r1[14,R2,Do_nothing] + _r1[14,R2,Recover] = 0
 R691: - _f[14,R3] + _r1[14,R3,Do_nothing] + _r1[14,R3,Recover] = 0
 R692: - _f[23,Do_nothing] + _r1[23,Do_nothing,Do_nothing]
   + _r1[23,Do_nothing,Recover] = 0
 R693: - _f[23,R1] + _r1[23,R1,Do_nothing] + _r1[23,R1,Recover] = 0
 R694: - _f[23,R2] + _r1[23,R2,Do_nothing] + _r1[23,R2,Recover] = 0
 R695: - _f[23,R3] + _r1[23,R3,Do_nothing] + _r1[23,R3,Recover] = 0
 R696: - _f[93,Do_nothing] + _r1[93,Do_nothing,Do_nothing]
   + _r1[93,Do_nothing,Recover] = 0
 R697: - _f[93,R1] + _r1[93,R1,Do_nothing] + _r1[93,R1,Recover] = 0
 R698: - _f[93,R2] + _r1[93,R2,Do_nothing] + _r1[93,R2,Recover] = 0
 R699: - _f[93,R3] + _r1[93,R3,Do_nothing] + _r1[93,R3,Recover] = 0
 R700: - _f[28,Do_nothing] + _r1[28,Do_nothing,Do_nothing]
   + _r1[28,Do_nothing,Recover] = 0
 R701: - _f[28,R1] + _r1[28,R1,Do_nothing] + _r1[28,R1,Recover] = 0
 R702: - _f[28,R2] + _r1[28,R2,Do_nothing] + _r1[28,R2,Recover] = 0
 R703: - _f[28,R3] + _r1[28,R3,Do_nothing] + _r1[28,R3,Recover] = 0
 R704: - _f[2,Do_nothing] + _r1[2,Do_nothing,Do_nothing]
   + _r1[2,Do_nothing,Recover] = 0
 R705: - _f[2,R1] + _r1[2,R1,Do_nothing] + _r1[2,R1,Recover] = 0
 R706: - _f[2,R2] + _r1[2,R2,Do_nothing] + _r1[2,R2,Recover] = 0
 R707: - _f[2,R3] + _r1[2,R3,Do_nothing] + _r1[2,R3,Recover] = 0
 R708: - _f[81,Do_nothing] + _r1[81,Do_nothing,Do_nothing]
   + _r1[81,Do_nothing,Recover] = 0
 R709: - _f[81,R1] + _r1[81,R1,Do_nothing] + _r1[81,R1,Recover] = 0
 R710: - _f[81,R2] + _r1[81,R2,Do_nothing] + _r1[81,R2,Recover] = 0
 R711: - _f[81,R3] + _r1[81,R3,Do_nothing] + _r1[81,R3,Recover] = 0
 R712: - _f[45,Do_nothing] + _r1[45,Do_nothing,Do_nothing]
   + _r1[45,Do_nothing,Recover] = 0
 R713: - _f[45,R1] + _r1[45,R1,Do_nothing] + _r1[45,R1,Recover] = 0
 R714: - _f[45,R2] + _r1[45,R2,Do_nothing] + _r1[45,R2,Recover] = 0
 R715: - _f[45,R3] + _r1[45,R3,Do_nothing] + _r1[45,R3,Recover] = 0
 R716: - _f[88,Do_nothing] + _r1[88,Do_nothing,Do_nothing]
   + _r1[88,Do_nothing,Recover] = 0
 R717: - _f[88,R1] + _r1[88,R1,Do_nothing] + _r1[88,R1,Recover] = 0
 R718: - _f[88,R2] + _r1[88,R2,Do_nothing] + _r1[88,R2,Recover] = 0
 R719: - _f[88,R3] + _r1[88,R3,Do_nothing] + _r1[88,R3,Recover] = 0
 R720: - _f[89,Do_nothing] + _r1[89,Do_nothing,Do_nothing]
   + _r1[89,Do_nothing,Recover] = 0
 R721: - _f[89,R1] + _r1[89,R1,Do_nothing] + _r1[89,R1,Recover] = 0
 R722: - _f[89,R2] + _r1[89,R2,Do_nothing] + _r1[89,R2,Recover] = 0
 R723: - _f[89,R3] + _r1[89,R3,Do_nothing] + _r1[89,R3,Recover] = 0
 R724: - _f[25,Do_nothing] + _r1[25,Do_nothing,Do_nothing]
   + _r1[25,Do_nothing,Recover] = 0
 R725: - _f[25,R1] + _r1[25,R1,Do_nothing] + _r1[25,R1,Recover] = 0
 R726: - _f[25,R2] + _r1[25,R2,Do_nothing] + _r1[25,R2,Recover] = 0
 R727: - _f[25,R3] + _r1[25,R3,Do_nothing] + _r1[25,R3,Recover] = 0
 R728: - _f[91,Do_nothing] + _r1[91,Do_nothing,Do_nothing]
   + _r1[91,Do_nothing,Recover] = 0
 R729: - _f[91,R1] + _r1[91,R1,Do_nothing] + _r1[91,R1,Recover] = 0
 R730: - _f[91,R2] + _r1[91,R2,Do_nothing] + _r1[91,R2,Recover] = 0
 R731: - _f[91,R3] + _r1[91,R3,Do_nothing] + _r1[91,R3,Recover] = 0
 R732: - _f[43,Do_nothing] + _r1[43,Do_nothing,Do_nothing]
   + _r1[43,Do_nothing,Recover] = 0
 R733: - _f[43,R1] + _r1[43,R1,Do_nothing] + _r1[43,R1,Recover] = 0
 R734: - _f[43,R2] + _r1[43,R2,Do_nothing] + _r1[43,R2,Recover] = 0
 R735: - _f[43,R3] + _r1[43,R3,Do_nothing] + _r1[43,R3,Recover] = 0
 R736: - _f[76,Do_nothing] + _r1[76,Do_nothing,Do_nothing]
   + _r1[76,Do_nothing,Recover] = 0
 R737: - _f[76,R1] + _r1[76,R1,Do_nothing] + _r1[76,R1,Recover] = 0
 R738: - _f[76,R2] + _r1[76,R2,Do_nothing] + _r1[76,R2,Recover] = 0
 R739: - _f[76,R3] + _r1[76,R3,Do_nothing] + _r1[76,R3,Recover] = 0
 R740: - _f[42,Do_nothing] + _r1[42,Do_nothing,Do_nothing]
   + _r1[42,Do_nothing,Recover] = 0
 R741: - _f[42,R1] + _r1[42,R1,Do_nothing] + _r1[42,R1,Recover] = 0
 R742: - _f[42,R2] + _r1[42,R2,Do_nothing] + _r1[42,R2,Recover] = 0
 R743: - _f[42,R3] + _r1[42,R3,Do_nothing] + _r1[42,R3,Recover] = 0
 R744: - _f[7,Do_nothing] + _r1[7,Do_nothing,Do_nothing]
   + _r1[7,Do_nothing,Recover] = 0
 R745: - _f[7,R1] + _r1[7,R1,Do_nothing] + _r1[7,R1,Recover] = 0
 R746: - _f[7,R2] + _r1[7,R2,Do_nothing] + _r1[7,R2,Recover] = 0
 R747: - _f[7,R3] + _r1[7,R3,Do_nothing] + _r1[7,R3,Recover] = 0
 R748: - _f[17,Do_nothing] + _r1[17,Do_nothing,Do_nothing]
   + _r1[17,Do_nothing,Recover] = 0
 R749: - _f[17,R1] + _r1[17,R1,Do_nothing] + _r1[17,R1,Recover] = 0
 R750: - _f[17,R2] + _r1[17,R2,Do_nothing] + _r1[17,R2,Recover] = 0
 R751: - _f[17,R3] + _r1[17,R3,Do_nothing] + _r1[17,R3,Recover] = 0
 R752: - _f[37,Do_nothing] + _r1[37,Do_nothing,Do_nothing]
   + _r1[37,Do_nothing,Recover] = 0
 R753: - _f[37,R1] + _r1[37,R1,Do_nothing] + _r1[37,R1,Recover] = 0
 R754: - _f[37,R2] + _r1[37,R2,Do_nothing] + _r1[37,R2,Recover] = 0
 R755: - _f[37,R3] + _r1[37,R3,Do_nothing] + _r1[37,R3,Recover] = 0
 R756: - _f[26,Do_nothing] + _r1[26,Do_nothing,Do_nothing]
   + _r1[26,Do_nothing,Recover] = 0
 R757: - _f[26,R1] + _r1[26,R1,Do_nothing] + _r1[26,R1,Recover] = 0
 R758: - _f[26,R2] + _r1[26,R2,Do_nothing] + _r1[26,R2,Recover] = 0
 R759: - _f[26,R3] + _r1[26,R3,Do_nothing] + _r1[26,R3,Recover] = 0
 R760: - _f[9,Do_nothing] + _r1[9,Do_nothing,Do_nothing]
   + _r1[9,Do_nothing,Recover] = 0
 R761: - _f[9,R1] + _r1[9,R1,Do_nothing] + _r1[9,R1,Recover] = 0
 R762: - _f[9,R2] + _r1[9,R2,Do_nothing] + _r1[9,R2,Recover] = 0
 R763: - _f[9,R3] + _r1[9,R3,Do_nothing] + _r1[9,R3,Recover] = 0
 R764: - _f[60,Do_nothing] + _r1[60,Do_nothing,Do_nothing]
   + _r1[60,Do_nothing,Recover] = 0
 R765: - _f[60,R1] + _r1[60,R1,Do_nothing] + _r1[60,R1,Recover] = 0
 R766: - _f[60,R2] + _r1[60,R2,Do_nothing] + _r1[60,R2,Recover] = 0
 R767: - _f[60,R3] + _r1[60,R3,Do_nothing] + _r1[60,R3,Recover] = 0
 R768: - _f[16,Do_nothing] + _r1[16,Do_nothing,Do_nothing]
   + _r1[16,Do_nothing,Recover] = 0
 R769: - _f[16,R1] + _r1[16,R1,Do_nothing] + _r1[16,R1,Recover] = 0
 R770: - _f[16,R2] + _r1[16,R2,Do_nothing] + _r1[16,R2,Recover] = 0
 R771: - _f[16,R3] + _r1[16,R3,Do_nothing] + _r1[16,R3,Recover] = 0
 R772: - _f[83,Do_nothing] + _r1[83,Do_nothing,Do_nothing]
   + _r1[83,Do_nothing,Recover] = 0
 R773: - _f[83,R1] + _r1[83,R1,Do_nothing] + _r1[83,R1,Recover] = 0
 R774: - _f[83,R2] + _r1[83,R2,Do_nothing] + _r1[83,R2,Recover] = 0
 R775: - _f[83,R3] + _r1[83,R3,Do_nothing] + _r1[83,R3,Recover] = 0
 R776: - _f[13,Do_nothing] + _r1[13,Do_nothing,Do_nothing]
   + _r1[13,Do_nothing,Recover] = 0
 R777: - _f[13,R1] + _r1[13,R1,Do_nothing] + _r1[13,R1,Recover] = 0
 R778: - _f[13,R2] + _r1[13,R2,Do_nothing] + _r1[13,R2,Recover] = 0
 R779: - _f[13,R3] + _r1[13,R3,Do_nothing] + _r1[13,R3,Recover] = 0
 R780: - _f[67,Do_nothing] + _r1[67,Do_nothing,Do_nothing]
   + _r1[67,Do_nothing,Recover] = 0
 R781: - _f[67,R1] + _r1[67,R1,Do_nothing] + _r1[67,R1,Recover] = 0
 R782: - _f[67,R2] + _r1[67,R2,Do_nothing] + _r1[67,R2,Recover] = 0
 R783: - _f[67,R3] + _r1[67,R3,Do_nothing] + _r1[67,R3,Recover] = 0
 R784: - _f[59,Do_nothing] + _r1[59,Do_nothing,Do_nothing]
   + _r1[59,Do_nothing,Recover] = 0
 R785: - _f[59,R1] + _r1[59,R1,Do_nothing] + _r1[59,R1,Recover] = 0
 R786: - _f[59,R2] + _r1[59,R2,Do_nothing] + _r1[59,R2,Recover] = 0
 R787: - _f[59,R3] + _r1[59,R3,Do_nothing] + _r1[59,R3,Recover] = 0
 R788: - _f[58,Do_nothing] + _r1[58,Do_nothing,Do_nothing]
   + _r1[58,Do_nothing,Recover] = 0
 R789: - _f[58,R1] + _r1[58,R1,Do_nothing] + _r1[58,R1,Recover] = 0
 R790: - _f[58,R2] + _r1[58,R2,Do_nothing] + _r1[58,R2,Recover] = 0
 R791: - _f[58,R3] + _r1[58,R3,Do_nothing] + _r1[58,R3,Recover] = 0
 R792: - _f[4,Do_nothing] + _r1[4,Do_nothing,Do_nothing]
   + _r1[4,Do_nothing,Recover] = 0
 R793: - _f[4,R1] + _r1[4,R1,Do_nothing] + _r1[4,R1,Recover] = 0
 R794: - _f[4,R2] + _r1[4,R2,Do_nothing] + _r1[4,R2,Recover] = 0
 R795: - _f[4,R3] + _r1[4,R3,Do_nothing] + _r1[4,R3,Recover] = 0
 R796: - _f[68,Do_nothing] + _r1[68,Do_nothing,Do_nothing]
   + _r1[68,Do_nothing,Recover] = 0
 R797: - _f[68,R1] + _r1[68,R1,Do_nothing] + _r1[68,R1,Recover] = 0
 R798: - _f[68,R2] + _r1[68,R2,Do_nothing] + _r1[68,R2,Recover] = 0
 R799: - _f[68,R3] + _r1[68,R3,Do_nothing] + _r1[68,R3,Recover] = 0
 R800: - _f[87,Do_nothing] + _r1[87,Do_nothing,Do_nothing]
   + _r1[87,Do_nothing,Recover] = 0
 R801: - _f[87,R1] + _r1[87,R1,Do_nothing] + _r1[87,R1,Recover] = 0
 R802: - _f[87,R2] + _r1[87,R2,Do_nothing] + _r1[87,R2,Recover] = 0
 R803: - _f[87,R3] + _r1[87,R3,Do_nothing] + _r1[87,R3,Recover] = 0
 R804: - _f[40,Do_nothing] + _r1[40,Do_nothing,Do_nothing]
   + _r1[40,Do_nothing,Recover] = 0
 R805: - _f[40,R1] + _r1[40,R1,Do_nothing] + _r1[40,R1,Recover] = 0
 R806: - _f[40,R2] + _r1[40,R2,Do_nothing] + _r1[40,R2,Recover] = 0
 R807: - _f[40,R3] + _r1[40,R3,Do_nothing] + _r1[40,R3,Recover] = 0
 R808: - _f[30,Do_nothing] + _r1[30,Do_nothing,Do_nothing]
   + _r1[30,Do_nothing,Recover] = 0
 R809: - _f[30,R1] + _r1[30,R1,Do_nothing] + _r1[30,R1,Recover] = 0
 R810: - _f[30,R2] + _r1[30,R2,Do_nothing] + _r1[30,R2,Recover] = 0
 R811: - _f[30,R3] + _r1[30,R3,Do_nothing] + _r1[30,R3,Recover] = 0
 R812: - _f[27,Do_nothing] + _r1[27,Do_nothing,Do_nothing]
   + _r1[27,Do_nothing,Recover] = 0
 R813: - _f[27,R1] + _r1[27,R1,Do_nothing] + _r1[27,R1,Recover] = 0
 R814: - _f[27,R2] + _r1[27,R2,Do_nothing] + _r1[27,R2,Recover] = 0
 R815: - _f[27,R3] + _r1[27,R3,Do_nothing] + _r1[27,R3,Recover] = 0
 R816: - _f[73,Do_nothing] + _r1[73,Do_nothing,Do_nothing]
   + _r1[73,Do_nothing,Recover] = 0
 R817: - _f[73,R1] + _r1[73,R1,Do_nothing] + _r1[73,R1,Recover] = 0
 R818: - _f[73,R2] + _r1[73,R2,Do_nothing] + _r1[73,R2,Recover] = 0
 R819: - _f[73,R3] + _r1[73,R3,Do_nothing] + _r1[73,R3,Recover] = 0
 R820: - _f[84,Do_nothing] + _r1[84,Do_nothing,Do_nothing]
   + _r1[84,Do_nothing,Recover] = 0
 R821: - _f[84,R1] + _r1[84,R1,Do_nothing] + _r1[84,R1,Recover] = 0
 R822: - _f[84,R2] + _r1[84,R2,Do_nothing] + _r1[84,R2,Recover] = 0
 R823: - _f[84,R3] + _r1[84,R3,Do_nothing] + _r1[84,R3,Recover] = 0
 R824: - _f[35,Do_nothing] + _r1[35,Do_nothing,Do_nothing]
   + _r1[35,Do_nothing,Recover] = 0
 R825: - _f[35,R1] + _r1[35,R1,Do_nothing] + _r1[35,R1,Recover] = 0
 R826: - _f[35,R2] + _r1[35,R2,Do_nothing] + _r1[35,R2,Recover] = 0
 R827: - _f[35,R3] + _r1[35,R3,Do_nothing] + _r1[35,R3,Recover] = 0
 R828: - _f[46,Do_nothing] + _r1[46,Do_nothing,Do_nothing]
   + _r1[46,Do_nothing,Recover] = 0
 R829: - _f[46,R1] + _r1[46,R1,Do_nothing] + _r1[46,R1,Recover] = 0
 R830: - _f[46,R2] + _r1[46,R2,Do_nothing] + _r1[46,R2,Recover] = 0
 R831: - _f[46,R3] + _r1[46,R3,Do_nothing] + _r1[46,R3,Recover] = 0
 R832: - _f[39,Do_nothing] + _r1[39,Do_nothing,Do_nothing]
   + _r1[39,Do_nothing,Recover] = 0
 R833: - _f[39,R1] + _r1[39,R1,Do_nothing] + _r1[39,R1,Recover] = 0
 R834: - _f[39,R2] + _r1[39,R2,Do_nothing] + _r1[39,R2,Recover] = 0
 R835: - _f[39,R3] + _r1[39,R3,Do_nothing] + _r1[39,R3,Recover] = 0
 R836: - _f[51,Do_nothing] + _r1[51,Do_nothing,Do_nothing]
   + _r1[51,Do_nothing,Recover] = 0
 R837: - _f[51,R1] + _r1[51,R1,Do_nothing] + _r1[51,R1,Recover] = 0
 R838: - _f[51,R2] + _r1[51,R2,Do_nothing] + _r1[51,R2,Recover] = 0
 R839: - _f[51,R3] + _r1[51,R3,Do_nothing] + _r1[51,R3,Recover] = 0
 R840: - _f[52,Do_nothing] + _r1[52,Do_nothing,Do_nothing]
   + _r1[52,Do_nothing,Recover] = 0
 R841: - _f[52,R1] + _r1[52,R1,Do_nothing] + _r1[52,R1,Recover] = 0
 R842: - _f[52,R2] + _r1[52,R2,Do_nothing] + _r1[52,R2,Recover] = 0
 R843: - _f[52,R3] + _r1[52,R3,Do_nothing] + _r1[52,R3,Recover] = 0
 R844: - _f[6,Do_nothing] + _r1[6,Do_nothing,Do_nothing]
   + _r1[6,Do_nothing,Recover] = 0
 R845: - _f[6,R1] + _r1[6,R1,Do_nothing] + _r1[6,R1,Recover] = 0
 R846: - _f[6,R2] + _r1[6,R2,Do_nothing] + _r1[6,R2,Recover] = 0
 R847: - _f[6,R3] + _r1[6,R3,Do_nothing] + _r1[6,R3,Recover] = 0
 R848: - _f[96,Do_nothing] + _r1[96,Do_nothing,Do_nothing]
   + _r1[96,Do_nothing,Recover] = 0
 R849: - _f[96,R1] + _r1[96,R1,Do_nothing] + _r1[96,R1,Recover] = 0
 R850: - _f[96,R2] + _r1[96,R2,Do_nothing] + _r1[96,R2,Recover] = 0
 R851: - _f[96,R3] + _r1[96,R3,Do_nothing] + _r1[96,R3,Recover] = 0
 R852: - _f[74,Do_nothing] + _r1[74,Do_nothing,Do_nothing]
   + _r1[74,Do_nothing,Recover] = 0
 R853: - _f[74,R1] + _r1[74,R1,Do_nothing] + _r1[74,R1,Recover] = 0
 R854: - _f[74,R2] + _r1[74,R2,Do_nothing] + _r1[74,R2,Recover] = 0
 R855: - _f[74,R3] + _r1[74,R3,Do_nothing] + _r1[74,R3,Recover] = 0
 R856: - _f[86,Do_nothing] + _r1[86,Do_nothing,Do_nothing]
   + _r1[86,Do_nothing,Recover] = 0
 R857: - _f[86,R1] + _r1[86,R1,Do_nothing] + _r1[86,R1,Recover] = 0
 R858: - _f[86,R2] + _r1[86,R2,Do_nothing] + _r1[86,R2,Recover] = 0
 R859: - _f[86,R3] + _r1[86,R3,Do_nothing] + _r1[86,R3,Recover] = 0
 R860: - _f[29,Do_nothing] + _r1[29,Do_nothing,Do_nothing]
   + _r1[29,Do_nothing,Recover] = 0
 R861: - _f[29,R1] + _r1[29,R1,Do_nothing] + _r1[29,R1,Recover] = 0
 R862: - _f[29,R2] + _r1[29,R2,Do_nothing] + _r1[29,R2,Recover] = 0
 R863: - _f[29,R3] + _r1[29,R3,Do_nothing] + _r1[29,R3,Recover] = 0
 R864: - _f[18,Do_nothing] + _r1[18,Do_nothing,Do_nothing]
   + _r1[18,Do_nothing,Recover] = 0
 R865: - _f[18,R1] + _r1[18,R1,Do_nothing] + _r1[18,R1,Recover] = 0
 R866: - _f[18,R2] + _r1[18,R2,Do_nothing] + _r1[18,R2,Recover] = 0
 R867: - _f[18,R3] + _r1[18,R3,Do_nothing] + _r1[18,R3,Recover] = 0
 R868: - _f[24,Do_nothing] + _r1[24,Do_nothing,Do_nothing]
   + _r1[24,Do_nothing,Recover] = 0
 R869: - _f[24,R1] + _r1[24,R1,Do_nothing] + _r1[24,R1,Recover] = 0
 R870: - _f[24,R2] + _r1[24,R2,Do_nothing] + _r1[24,R2,Recover] = 0
 R871: - _f[24,R3] + _r1[24,R3,Do_nothing] + _r1[24,R3,Recover] = 0
 R872: - _f[94,Do_nothing] + _r1[94,Do_nothing,Do_nothing]
   + _r1[94,Do_nothing,Recover] = 0
 R873: - _f[94,R1] + _r1[94,R1,Do_nothing] + _r1[94,R1,Recover] = 0
 R874: - _f[94,R2] + _r1[94,R2,Do_nothing] + _r1[94,R2,Recover] = 0
 R875: - _f[94,R3] + _r1[94,R3,Do_nothing] + _r1[94,R3,Recover] = 0
 R876: - _f[33,Do_nothing] + _r1[33,Do_nothing,Do_nothing]
   + _r1[33,Do_nothing,Recover] = 0
 R877: - _f[33,R1] + _r1[33,R1,Do_nothing] + _r1[33,R1,Recover] = 0
 R878: - _f[33,R2] + _r1[33,R2,Do_nothing] + _r1[33,R2,Recover] = 0
 R879: - _f[33,R3] + _r1[33,R3,Do_nothing] + _r1[33,R3,Recover] = 0
 R880: - _f[12,Do_nothing] + _r1[12,Do_nothing,Do_nothing]
   + _r1[12,Do_nothing,Recover] = 0
 R881: - _f[12,R1] + _r1[12,R1,Do_nothing] + _r1[12,R1,Recover] = 0
 R882: - _f[12,R2] + _r1[12,R2,Do_nothing] + _r1[12,R2,Recover] = 0
 R883: - _f[12,R3] + _r1[12,R3,Do_nothing] + _r1[12,R3,Recover] = 0
 R884: - _f[57,Do_nothing] + _r1[57,Do_nothing,Do_nothing]
   + _r1[57,Do_nothing,Recover] = 0
 R885: - _f[57,R1] + _r1[57,R1,Do_nothing] + _r1[57,R1,Recover] = 0
 R886: - _f[57,R2] + _r1[57,R2,Do_nothing] + _r1[57,R2,Recover] = 0
 R887: - _f[57,R3] + _r1[57,R3,Do_nothing] + _r1[57,R3,Recover] = 0
 R888: - _f[97,Do_nothing] + _r1[97,Do_nothing,Do_nothing]
   + _r1[97,Do_nothing,Recover] = 0
 R889: - _f[97,R1] + _r1[97,R1,Do_nothing] + _r1[97,R1,Recover] = 0
 R890: - _f[97,R2] + _r1[97,R2,Do_nothing] + _r1[97,R2,Recover] = 0
 R891: - _f[97,R3] + _r1[97,R3,Do_nothing] + _r1[97,R3,Recover] = 0
 R892: - _f[69,Do_nothing] + _r1[69,Do_nothing,Do_nothing]
   + _r1[69,Do_nothing,Recover] = 0
 R893: - _f[69,R1] + _r1[69,R1,Do_nothing] + _r1[69,R1,Recover] = 0
 R894: - _f[69,R2] + _r1[69,R2,Do_nothing] + _r1[69,R2,Recover] = 0
 R895: - _f[69,R3] + _r1[69,R3,Do_nothing] + _r1[69,R3,Recover] = 0
 R896: - _f[3,Do_nothing] + _r1[3,Do_nothing,Do_nothing]
   + _r1[3,Do_nothing,Recover] = 0
 R897: - _f[3,R1] + _r1[3,R1,Do_nothing] + _r1[3,R1,Recover] = 0
 R898: - _f[3,R2] + _r1[3,R2,Do_nothing] + _r1[3,R2,Recover] = 0
 R899: - _f[3,R3] + _r1[3,R3,Do_nothing] + _r1[3,R3,Recover] = 0
 R900: - _f[66,Do_nothing] + _r1[66,Do_nothing,Do_nothing]
   + _r1[66,Do_nothing,Recover] = 0
 R901: - _f[66,R1] + _r1[66,R1,Do_nothing] + _r1[66,R1,Recover] = 0
 R902: - _f[66,R2] + _r1[66,R2,Do_nothing] + _r1[66,R2,Recover] = 0
 R903: - _f[66,R3] + _r1[66,R3,Do_nothing] + _r1[66,R3,Recover] = 0
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
 _f[66,R1] _f[66,R2] _f[66,R3] _r0[32,Do_nothing,Do_nothing]
 _r0[32,Do_nothing,Recover] _r0[32,R1,Do_nothing] _r0[32,R1,Recover]
 _r0[32,R2,Do_nothing] _r0[32,R2,Recover] _r0[32,R3,Do_nothing]
 _r0[32,R3,Recover] _r0[70,Do_nothing,Do_nothing]
 _r0[70,Do_nothing,Recover] _r0[70,R1,Do_nothing] _r0[70,R1,Recover]
 _r0[70,R2,Do_nothing] _r0[70,R2,Recover] _r0[70,R3,Do_nothing]
 _r0[70,R3,Recover] _r0[5,Do_nothing,Do_nothing] _r0[5,Do_nothing,Recover]
 _r0[5,R1,Do_nothing] _r0[5,R1,Recover] _r0[5,R2,Do_nothing]
 _r0[5,R2,Recover] _r0[5,R3,Do_nothing] _r0[5,R3,Recover]
 _r0[85,Do_nothing,Do_nothing] _r0[85,Do_nothing,Recover]
 _r0[85,R1,Do_nothing] _r0[85,R1,Recover] _r0[85,R2,Do_nothing]
 _r0[85,R2,Recover] _r0[85,R3,Do_nothing] _r0[85,R3,Recover]
 _r0[71,Do_nothing,Do_nothing] _r0[71,Do_nothing,Recover]
 _r0[71,R1,Do_nothing] _r0[71,R1,Recover] _r0[71,R2,Do_nothing]
 _r0[71,R2,Recover] _r0[71,R3,Do_nothing] _r0[71,R3,Recover]
 _r0[92,Do_nothing,Do_nothing] _r0[92,Do_nothing,Recover]
 _r0[92,R1,Do_nothing] _r0[92,R1,Recover] _r0[92,R2,Do_nothing]
 _r0[92,R2,Recover] _r0[92,R3,Do_nothing] _r0[92,R3,Recover]
 _r0[38,Do_nothing,Do_nothing] _r0[38,Do_nothing,Recover]
 _r0[38,R1,Do_nothing] _r0[38,R1,Recover] _r0[38,R2,Do_nothing]
 _r0[38,R2,Recover] _r0[38,R3,Do_nothing] _r0[38,R3,Recover]
 _r0[21,Do_nothing,Do_nothing] _r0[21,Do_nothing,Recover]
 _r0[21,R1,Do_nothing] _r0[21,R1,Recover] _r0[21,R2,Do_nothing]
 _r0[21,R2,Recover] _r0[21,R3,Do_nothing] _r0[21,R3,Recover]
 _r0[41,Do_nothing,Do_nothing] _r0[41,Do_nothing,Recover]
 _r0[41,R1,Do_nothing] _r0[41,R1,Recover] _r0[41,R2,Do_nothing]
 _r0[41,R2,Recover] _r0[41,R3,Do_nothing] _r0[41,R3,Recover]
 _r0[56,Do_nothing,Do_nothing] _r0[56,Do_nothing,Recover]
 _r0[56,R1,Do_nothing] _r0[56,R1,Recover] _r0[56,R2,Do_nothing]
 _r0[56,R2,Recover] _r0[56,R3,Do_nothing] _r0[56,R3,Recover]
 _r0[80,Do_nothing,Do_nothing] _r0[80,Do_nothing,Recover]
 _r0[80,R1,Do_nothing] _r0[80,R1,Recover] _r0[80,R2,Do_nothing]
 _r0[80,R2,Recover] _r0[80,R3,Do_nothing] _r0[80,R3,Recover]
 _r0[62,Do_nothing,Do_nothing] _r0[62,Do_nothing,Recover]
 _r0[62,R1,Do_nothing] _r0[62,R1,Recover] _r0[62,R2,Do_nothing]
 _r0[62,R2,Recover] _r0[62,R3,Do_nothing] _r0[62,R3,Recover]
 _r0[20,Do_nothing,Do_nothing] _r0[20,Do_nothing,Recover]
 _r0[20,R1,Do_nothing] _r0[20,R1,Recover] _r0[20,R2,Do_nothing]
 _r0[20,R2,Recover] _r0[20,R3,Do_nothing] _r0[20,R3,Recover]
 _r0[10,Do_nothing,Do_nothing] _r0[10,Do_nothing,Recover]
 _r0[10,R1,Do_nothing] _r0[10,R1,Recover] _r0[10,R2,Do_nothing]
 _r0[10,R2,Recover] _r0[10,R3,Do_nothing] _r0[10,R3,Recover]
 _r0[55,Do_nothing,Do_nothing] _r0[55,Do_nothing,Recover]
 _r0[55,R1,Do_nothing] _r0[55,R1,Recover] _r0[55,R2,Do_nothing]
 _r0[55,R2,Recover] _r0[55,R3,Do_nothing] _r0[55,R3,Recover]
 _r0[36,Do_nothing,Do_nothing] _r0[36,Do_nothing,Recover]
 _r0[36,R1,Do_nothing] _r0[36,R1,Recover] _r0[36,R2,Do_nothing]
 _r0[36,R2,Recover] _r0[36,R3,Do_nothing] _r0[36,R3,Recover]
 _r0[44,Do_nothing,Do_nothing] _r0[44,Do_nothing,Recover]
 _r0[44,R1,Do_nothing] _r0[44,R1,Recover] _r0[44,R2,Do_nothing]
 _r0[44,R2,Recover] _r0[44,R3,Do_nothing] _r0[44,R3,Recover]
 _r0[48,Do_nothing,Do_nothing] _r0[48,Do_nothing,Recover]
 _r0[48,R1,Do_nothing] _r0[48,R1,Recover] _r0[48,R2,Do_nothing]
 _r0[48,R2,Recover] _r0[48,R3,Do_nothing] _r0[48,R3,Recover]
 _r0[8,Do_nothing,Do_nothing] _r0[8,Do_nothing,Recover]
 _r0[8,R1,Do_nothing] _r0[8,R1,Recover] _r0[8,R2,Do_nothing]
 _r0[8,R2,Recover] _r0[8,R3,Do_nothing] _r0[8,R3,Recover]
 _r0[15,Do_nothing,Do_nothing] _r0[15,Do_nothing,Recover]
 _r0[15,R1,Do_nothing] _r0[15,R1,Recover] _r0[15,R2,Do_nothing]
 _r0[15,R2,Recover] _r0[15,R3,Do_nothing] _r0[15,R3,Recover]
 _r0[61,Do_nothing,Do_nothing] _r0[61,Do_nothing,Recover]
 _r0[61,R1,Do_nothing] _r0[61,R1,Recover] _r0[61,R2,Do_nothing]
 _r0[61,R2,Recover] _r0[61,R3,Do_nothing] _r0[61,R3,Recover]
 _r0[78,Do_nothing,Do_nothing] _r0[78,Do_nothing,Recover]
 _r0[78,R1,Do_nothing] _r0[78,R1,Recover] _r0[78,R2,Do_nothing]
 _r0[78,R2,Recover] _r0[78,R3,Do_nothing] _r0[78,R3,Recover]
 _r0[34,Do_nothing,Do_nothing] _r0[34,Do_nothing,Recover]
 _r0[34,R1,Do_nothing] _r0[34,R1,Recover] _r0[34,R2,Do_nothing]
 _r0[34,R2,Recover] _r0[34,R3,Do_nothing] _r0[34,R3,Recover]
 _r0[100,Do_nothing,Do_nothing] _r0[100,Do_nothing,Recover]
 _r0[100,R1,Do_nothing] _r0[100,R1,Recover] _r0[100,R2,Do_nothing]
 _r0[100,R2,Recover] _r0[100,R3,Do_nothing] _r0[100,R3,Recover]
 _r0[90,Do_nothing,Do_nothing] _r0[90,Do_nothing,Recover]
 _r0[90,R1,Do_nothing] _r0[90,R1,Recover] _r0[90,R2,Do_nothing]
 _r0[90,R2,Recover] _r0[90,R3,Do_nothing] _r0[90,R3,Recover]
 _r0[64,Do_nothing,Do_nothing] _r0[64,Do_nothing,Recover]
 _r0[64,R1,Do_nothing] _r0[64,R1,Recover] _r0[64,R2,Do_nothing]
 _r0[64,R2,Recover] _r0[64,R3,Do_nothing] _r0[64,R3,Recover]
 _r0[79,Do_nothing,Do_nothing] _r0[79,Do_nothing,Recover]
 _r0[79,R1,Do_nothing] _r0[79,R1,Recover] _r0[79,R2,Do_nothing]
 _r0[79,R2,Recover] _r0[79,R3,Do_nothing] _r0[79,R3,Recover]
 _r0[1,Do_nothing,Do_nothing] _r0[1,Do_nothing,Recover]
 _r0[1,R1,Do_nothing] _r0[1,R1,Recover] _r0[1,R2,Do_nothing]
 _r0[1,R2,Recover] _r0[1,R3,Do_nothing] _r0[1,R3,Recover]
 _r0[22,Do_nothing,Do_nothing] _r0[22,Do_nothing,Recover]
 _r0[22,R1,Do_nothing] _r0[22,R1,Recover] _r0[22,R2,Do_nothing]
 _r0[22,R2,Recover] _r0[22,R3,Do_nothing] _r0[22,R3,Recover]
 _r0[31,Do_nothing,Do_nothing] _r0[31,Do_nothing,Recover]
 _r0[31,R1,Do_nothing] _r0[31,R1,Recover] _r0[31,R2,Do_nothing]
 _r0[31,R2,Recover] _r0[31,R3,Do_nothing] _r0[31,R3,Recover]
 _r0[54,Do_nothing,Do_nothing] _r0[54,Do_nothing,Recover]
 _r0[54,R1,Do_nothing] _r0[54,R1,Recover] _r0[54,R2,Do_nothing]
 _r0[54,R2,Recover] _r0[54,R3,Do_nothing] _r0[54,R3,Recover]
 _r0[82,Do_nothing,Do_nothing] _r0[82,Do_nothing,Recover]
 _r0[82,R1,Do_nothing] _r0[82,R1,Recover] _r0[82,R2,Do_nothing]
 _r0[82,R2,Recover] _r0[82,R3,Do_nothing] _r0[82,R3,Recover]
 _r0[11,Do_nothing,Do_nothing] _r0[11,Do_nothing,Recover]
 _r0[11,R1,Do_nothing] _r0[11,R1,Recover] _r0[11,R2,Do_nothing]
 _r0[11,R2,Recover] _r0[11,R3,Do_nothing] _r0[11,R3,Recover]
 _r0[75,Do_nothing,Do_nothing] _r0[75,Do_nothing,Recover]
 _r0[75,R1,Do_nothing] _r0[75,R1,Recover] _r0[75,R2,Do_nothing]
 _r0[75,R2,Recover] _r0[75,R3,Do_nothing] _r0[75,R3,Recover]
 _r0[50,Do_nothing,Do_nothing] _r0[50,Do_nothing,Recover]
 _r0[50,R1,Do_nothing] _r0[50,R1,Recover] _r0[50,R2,Do_nothing]
 _r0[50,R2,Recover] _r0[50,R3,Do_nothing] _r0[50,R3,Recover]
 _r0[19,Do_nothing,Do_nothing] _r0[19,Do_nothing,Recover]
 _r0[19,R1,Do_nothing] _r0[19,R1,Recover] _r0[19,R2,Do_nothing]
 _r0[19,R2,Recover] _r0[19,R3,Do_nothing] _r0[19,R3,Recover]
 _r0[99,Do_nothing,Do_nothing] _r0[99,Do_nothing,Recover]
 _r0[99,R1,Do_nothing] _r0[99,R1,Recover] _r0[99,R2,Do_nothing]
 _r0[99,R2,Recover] _r0[99,R3,Do_nothing] _r0[99,R3,Recover]
 _r0[72,Do_nothing,Do_nothing] _r0[72,Do_nothing,Recover]
 _r0[72,R1,Do_nothing] _r0[72,R1,Recover] _r0[72,R2,Do_nothing]
 _r0[72,R2,Recover] _r0[72,R3,Do_nothing] _r0[72,R3,Recover]
 _r0[47,Do_nothing,Do_nothing] _r0[47,Do_nothing,Recover]
 _r0[47,R1,Do_nothing] _r0[47,R1,Recover] _r0[47,R2,Do_nothing]
 _r0[47,R2,Recover] _r0[47,R3,Do_nothing] _r0[47,R3,Recover]
 _r0[49,Do_nothing,Do_nothing] _r0[49,Do_nothing,Recover]
 _r0[49,R1,Do_nothing] _r0[49,R1,Recover] _r0[49,R2,Do_nothing]
 _r0[49,R2,Recover] _r0[49,R3,Do_nothing] _r0[49,R3,Recover]
 _r0[77,Do_nothing,Do_nothing] _r0[77,Do_nothing,Recover]
 _r0[77,R1,Do_nothing] _r0[77,R1,Recover] _r0[77,R2,Do_nothing]
 _r0[77,R2,Recover] _r0[77,R3,Do_nothing] _r0[77,R3,Recover]
 _r0[95,Do_nothing,Do_nothing] _r0[95,Do_nothing,Recover]
 _r0[95,R1,Do_nothing] _r0[95,R1,Recover] _r0[95,R2,Do_nothing]
 _r0[95,R2,Recover] _r0[95,R3,Do_nothing] _r0[95,R3,Recover]
 _r0[98,Do_nothing,Do_nothing] _r0[98,Do_nothing,Recover]
 _r0[98,R1,Do_nothing] _r0[98,R1,Recover] _r0[98,R2,Do_nothing]
 _r0[98,R2,Recover] _r0[98,R3,Do_nothing] _r0[98,R3,Recover]
 _r0[53,Do_nothing,Do_nothing] _r0[53,Do_nothing,Recover]
 _r0[53,R1,Do_nothing] _r0[53,R1,Recover] _r0[53,R2,Do_nothing]
 _r0[53,R2,Recover] _r0[53,R3,Do_nothing] _r0[53,R3,Recover]
 _r0[65,Do_nothing,Do_nothing] _r0[65,Do_nothing,Recover]
 _r0[65,R1,Do_nothing] _r0[65,R1,Recover] _r0[65,R2,Do_nothing]
 _r0[65,R2,Recover] _r0[65,R3,Do_nothing] _r0[65,R3,Recover]
 _r0[63,Do_nothing,Do_nothing] _r0[63,Do_nothing,Recover]
 _r0[63,R1,Do_nothing] _r0[63,R1,Recover] _r0[63,R2,Do_nothing]
 _r0[63,R2,Recover] _r0[63,R3,Do_nothing] _r0[63,R3,Recover]
 _r0[14,Do_nothing,Do_nothing] _r0[14,Do_nothing,Recover]
 _r0[14,R1,Do_nothing] _r0[14,R1,Recover] _r0[14,R2,Do_nothing]
 _r0[14,R2,Recover] _r0[14,R3,Do_nothing] _r0[14,R3,Recover]
 _r0[23,Do_nothing,Do_nothing] _r0[23,Do_nothing,Recover]
 _r0[23,R1,Do_nothing] _r0[23,R1,Recover] _r0[23,R2,Do_nothing]
 _r0[23,R2,Recover] _r0[23,R3,Do_nothing] _r0[23,R3,Recover]
 _r0[93,Do_nothing,Do_nothing] _r0[93,Do_nothing,Recover]
 _r0[93,R1,Do_nothing] _r0[93,R1,Recover] _r0[93,R2,Do_nothing]
 _r0[93,R2,Recover] _r0[93,R3,Do_nothing] _r0[93,R3,Recover]
 _r0[28,Do_nothing,Do_nothing] _r0[28,Do_nothing,Recover]
 _r0[28,R1,Do_nothing] _r0[28,R1,Recover] _r0[28,R2,Do_nothing]
 _r0[28,R2,Recover] _r0[28,R3,Do_nothing] _r0[28,R3,Recover]
 _r0[2,Do_nothing,Do_nothing] _r0[2,Do_nothing,Recover]
 _r0[2,R1,Do_nothing] _r0[2,R1,Recover] _r0[2,R2,Do_nothing]
 _r0[2,R2,Recover] _r0[2,R3,Do_nothing] _r0[2,R3,Recover]
 _r0[81,Do_nothing,Do_nothing] _r0[81,Do_nothing,Recover]
 _r0[81,R1,Do_nothing] _r0[81,R1,Recover] _r0[81,R2,Do_nothing]
 _r0[81,R2,Recover] _r0[81,R3,Do_nothing] _r0[81,R3,Recover]
 _r0[45,Do_nothing,Do_nothing] _r0[45,Do_nothing,Recover]
 _r0[45,R1,Do_nothing] _r0[45,R1,Recover] _r0[45,R2,Do_nothing]
 _r0[45,R2,Recover] _r0[45,R3,Do_nothing] _r0[45,R3,Recover]
 _r0[88,Do_nothing,Do_nothing] _r0[88,Do_nothing,Recover]
 _r0[88,R1,Do_nothing] _r0[88,R1,Recover] _r0[88,R2,Do_nothing]
 _r0[88,R2,Recover] _r0[88,R3,Do_nothing] _r0[88,R3,Recover]
 _r0[89,Do_nothing,Do_nothing] _r0[89,Do_nothing,Recover]
 _r0[89,R1,Do_nothing] _r0[89,R1,Recover] _r0[89,R2,Do_nothing]
 _r0[89,R2,Recover] _r0[89,R3,Do_nothing] _r0[89,R3,Recover]
 _r0[25,Do_nothing,Do_nothing] _r0[25,Do_nothing,Recover]
 _r0[25,R1,Do_nothing] _r0[25,R1,Recover] _r0[25,R2,Do_nothing]
 _r0[25,R2,Recover] _r0[25,R3,Do_nothing] _r0[25,R3,Recover]
 _r0[91,Do_nothing,Do_nothing] _r0[91,Do_nothing,Recover]
 _r0[91,R1,Do_nothing] _r0[91,R1,Recover] _r0[91,R2,Do_nothing]
 _r0[91,R2,Recover] _r0[91,R3,Do_nothing] _r0[91,R3,Recover]
 _r0[43,Do_nothing,Do_nothing] _r0[43,Do_nothing,Recover]
 _r0[43,R1,Do_nothing] _r0[43,R1,Recover] _r0[43,R2,Do_nothing]
 _r0[43,R2,Recover] _r0[43,R3,Do_nothing] _r0[43,R3,Recover]
 _r0[76,Do_nothing,Do_nothing] _r0[76,Do_nothing,Recover]
 _r0[76,R1,Do_nothing] _r0[76,R1,Recover] _r0[76,R2,Do_nothing]
 _r0[76,R2,Recover] _r0[76,R3,Do_nothing] _r0[76,R3,Recover]
 _r0[42,Do_nothing,Do_nothing] _r0[42,Do_nothing,Recover]
 _r0[42,R1,Do_nothing] _r0[42,R1,Recover] _r0[42,R2,Do_nothing]
 _r0[42,R2,Recover] _r0[42,R3,Do_nothing] _r0[42,R3,Recover]
 _r0[7,Do_nothing,Do_nothing] _r0[7,Do_nothing,Recover]
 _r0[7,R1,Do_nothing] _r0[7,R1,Recover] _r0[7,R2,Do_nothing]
 _r0[7,R2,Recover] _r0[7,R3,Do_nothing] _r0[7,R3,Recover]
 _r0[17,Do_nothing,Do_nothing] _r0[17,Do_nothing,Recover]
 _r0[17,R1,Do_nothing] _r0[17,R1,Recover] _r0[17,R2,Do_nothing]
 _r0[17,R2,Recover] _r0[17,R3,Do_nothing] _r0[17,R3,Recover]
 _r0[37,Do_nothing,Do_nothing] _r0[37,Do_nothing,Recover]
 _r0[37,R1,Do_nothing] _r0[37,R1,Recover] _r0[37,R2,Do_nothing]
 _r0[37,R2,Recover] _r0[37,R3,Do_nothing] _r0[37,R3,Recover]
 _r0[26,Do_nothing,Do_nothing] _r0[26,Do_nothing,Recover]
 _r0[26,R1,Do_nothing] _r0[26,R1,Recover] _r0[26,R2,Do_nothing]
 _r0[26,R2,Recover] _r0[26,R3,Do_nothing] _r0[26,R3,Recover]
 _r0[9,Do_nothing,Do_nothing] _r0[9,Do_nothing,Recover]
 _r0[9,R1,Do_nothing] _r0[9,R1,Recover] _r0[9,R2,Do_nothing]
 _r0[9,R2,Recover] _r0[9,R3,Do_nothing] _r0[9,R3,Recover]
 _r0[60,Do_nothing,Do_nothing] _r0[60,Do_nothing,Recover]
 _r0[60,R1,Do_nothing] _r0[60,R1,Recover] _r0[60,R2,Do_nothing]
 _r0[60,R2,Recover] _r0[60,R3,Do_nothing] _r0[60,R3,Recover]
 _r0[16,Do_nothing,Do_nothing] _r0[16,Do_nothing,Recover]
 _r0[16,R1,Do_nothing] _r0[16,R1,Recover] _r0[16,R2,Do_nothing]
 _r0[16,R2,Recover] _r0[16,R3,Do_nothing] _r0[16,R3,Recover]
 _r0[83,Do_nothing,Do_nothing] _r0[83,Do_nothing,Recover]
 _r0[83,R1,Do_nothing] _r0[83,R1,Recover] _r0[83,R2,Do_nothing]
 _r0[83,R2,Recover] _r0[83,R3,Do_nothing] _r0[83,R3,Recover]
 _r0[13,Do_nothing,Do_nothing] _r0[13,Do_nothing,Recover]
 _r0[13,R1,Do_nothing] _r0[13,R1,Recover] _r0[13,R2,Do_nothing]
 _r0[13,R2,Recover] _r0[13,R3,Do_nothing] _r0[13,R3,Recover]
 _r0[67,Do_nothing,Do_nothing] _r0[67,Do_nothing,Recover]
 _r0[67,R1,Do_nothing] _r0[67,R1,Recover] _r0[67,R2,Do_nothing]
 _r0[67,R2,Recover] _r0[67,R3,Do_nothing] _r0[67,R3,Recover]
 _r0[59,Do_nothing,Do_nothing] _r0[59,Do_nothing,Recover]
 _r0[59,R1,Do_nothing] _r0[59,R1,Recover] _r0[59,R2,Do_nothing]
 _r0[59,R2,Recover] _r0[59,R3,Do_nothing] _r0[59,R3,Recover]
 _r0[58,Do_nothing,Do_nothing] _r0[58,Do_nothing,Recover]
 _r0[58,R1,Do_nothing] _r0[58,R1,Recover] _r0[58,R2,Do_nothing]
 _r0[58,R2,Recover] _r0[58,R3,Do_nothing] _r0[58,R3,Recover]
 _r0[4,Do_nothing,Do_nothing] _r0[4,Do_nothing,Recover]
 _r0[4,R1,Do_nothing] _r0[4,R1,Recover] _r0[4,R2,Do_nothing]
 _r0[4,R2,Recover] _r0[4,R3,Do_nothing] _r0[4,R3,Recover]
 _r0[68,Do_nothing,Do_nothing] _r0[68,Do_nothing,Recover]
 _r0[68,R1,Do_nothing] _r0[68,R1,Recover] _r0[68,R2,Do_nothing]
 _r0[68,R2,Recover] _r0[68,R3,Do_nothing] _r0[68,R3,Recover]
 _r0[87,Do_nothing,Do_nothing] _r0[87,Do_nothing,Recover]
 _r0[87,R1,Do_nothing] _r0[87,R1,Recover] _r0[87,R2,Do_nothing]
 _r0[87,R2,Recover] _r0[87,R3,Do_nothing] _r0[87,R3,Recover]
 _r0[40,Do_nothing,Do_nothing] _r0[40,Do_nothing,Recover]
 _r0[40,R1,Do_nothing] _r0[40,R1,Recover] _r0[40,R2,Do_nothing]
 _r0[40,R2,Recover] _r0[40,R3,Do_nothing] _r0[40,R3,Recover]
 _r0[30,Do_nothing,Do_nothing] _r0[30,Do_nothing,Recover]
 _r0[30,R1,Do_nothing] _r0[30,R1,Recover] _r0[30,R2,Do_nothing]
 _r0[30,R2,Recover] _r0[30,R3,Do_nothing] _r0[30,R3,Recover]
 _r0[27,Do_nothing,Do_nothing] _r0[27,Do_nothing,Recover]
 _r0[27,R1,Do_nothing] _r0[27,R1,Recover] _r0[27,R2,Do_nothing]
 _r0[27,R2,Recover] _r0[27,R3,Do_nothing] _r0[27,R3,Recover]
 _r0[73,Do_nothing,Do_nothing] _r0[73,Do_nothing,Recover]
 _r0[73,R1,Do_nothing] _r0[73,R1,Recover] _r0[73,R2,Do_nothing]
 _r0[73,R2,Recover] _r0[73,R3,Do_nothing] _r0[73,R3,Recover]
 _r0[84,Do_nothing,Do_nothing] _r0[84,Do_nothing,Recover]
 _r0[84,R1,Do_nothing] _r0[84,R1,Recover] _r0[84,R2,Do_nothing]
 _r0[84,R2,Recover] _r0[84,R3,Do_nothing] _r0[84,R3,Recover]
 _r0[35,Do_nothing,Do_nothing] _r0[35,Do_nothing,Recover]
 _r0[35,R1,Do_nothing] _r0[35,R1,Recover] _r0[35,R2,Do_nothing]
 _r0[35,R2,Recover] _r0[35,R3,Do_nothing] _r0[35,R3,Recover]
 _r0[46,Do_nothing,Do_nothing] _r0[46,Do_nothing,Recover]
 _r0[46,R1,Do_nothing] _r0[46,R1,Recover] _r0[46,R2,Do_nothing]
 _r0[46,R2,Recover] _r0[46,R3,Do_nothing] _r0[46,R3,Recover]
 _r0[39,Do_nothing,Do_nothing] _r0[39,Do_nothing,Recover]
 _r0[39,R1,Do_nothing] _r0[39,R1,Recover] _r0[39,R2,Do_nothing]
 _r0[39,R2,Recover] _r0[39,R3,Do_nothing] _r0[39,R3,Recover]
 _r0[51,Do_nothing,Do_nothing] _r0[51,Do_nothing,Recover]
 _r0[51,R1,Do_nothing] _r0[51,R1,Recover] _r0[51,R2,Do_nothing]
 _r0[51,R2,Recover] _r0[51,R3,Do_nothing] _r0[51,R3,Recover]
 _r0[52,Do_nothing,Do_nothing] _r0[52,Do_nothing,Recover]
 _r0[52,R1,Do_nothing] _r0[52,R1,Recover] _r0[52,R2,Do_nothing]
 _r0[52,R2,Recover] _r0[52,R3,Do_nothing] _r0[52,R3,Recover]
 _r0[6,Do_nothing,Do_nothing] _r0[6,Do_nothing,Recover]
 _r0[6,R1,Do_nothing] _r0[6,R1,Recover] _r0[6,R2,Do_nothing]
 _r0[6,R2,Recover] _r0[6,R3,Do_nothing] _r0[6,R3,Recover]
 _r0[96,Do_nothing,Do_nothing] _r0[96,Do_nothing,Recover]
 _r0[96,R1,Do_nothing] _r0[96,R1,Recover] _r0[96,R2,Do_nothing]
 _r0[96,R2,Recover] _r0[96,R3,Do_nothing] _r0[96,R3,Recover]
 _r0[74,Do_nothing,Do_nothing] _r0[74,Do_nothing,Recover]
 _r0[74,R1,Do_nothing] _r0[74,R1,Recover] _r0[74,R2,Do_nothing]
 _r0[74,R2,Recover] _r0[74,R3,Do_nothing] _r0[74,R3,Recover]
 _r0[86,Do_nothing,Do_nothing] _r0[86,Do_nothing,Recover]
 _r0[86,R1,Do_nothing] _r0[86,R1,Recover] _r0[86,R2,Do_nothing]
 _r0[86,R2,Recover] _r0[86,R3,Do_nothing] _r0[86,R3,Recover]
 _r0[29,Do_nothing,Do_nothing] _r0[29,Do_nothing,Recover]
 _r0[29,R1,Do_nothing] _r0[29,R1,Recover] _r0[29,R2,Do_nothing]
 _r0[29,R2,Recover] _r0[29,R3,Do_nothing] _r0[29,R3,Recover]
 _r0[18,Do_nothing,Do_nothing] _r0[18,Do_nothing,Recover]
 _r0[18,R1,Do_nothing] _r0[18,R1,Recover] _r0[18,R2,Do_nothing]
 _r0[18,R2,Recover] _r0[18,R3,Do_nothing] _r0[18,R3,Recover]
 _r0[24,Do_nothing,Do_nothing] _r0[24,Do_nothing,Recover]
 _r0[24,R1,Do_nothing] _r0[24,R1,Recover] _r0[24,R2,Do_nothing]
 _r0[24,R2,Recover] _r0[24,R3,Do_nothing] _r0[24,R3,Recover]
 _r0[94,Do_nothing,Do_nothing] _r0[94,Do_nothing,Recover]
 _r0[94,R1,Do_nothing] _r0[94,R1,Recover] _r0[94,R2,Do_nothing]
 _r0[94,R2,Recover] _r0[94,R3,Do_nothing] _r0[94,R3,Recover]
 _r0[33,Do_nothing,Do_nothing] _r0[33,Do_nothing,Recover]
 _r0[33,R1,Do_nothing] _r0[33,R1,Recover] _r0[33,R2,Do_nothing]
 _r0[33,R2,Recover] _r0[33,R3,Do_nothing] _r0[33,R3,Recover]
 _r0[12,Do_nothing,Do_nothing] _r0[12,Do_nothing,Recover]
 _r0[12,R1,Do_nothing] _r0[12,R1,Recover] _r0[12,R2,Do_nothing]
 _r0[12,R2,Recover] _r0[12,R3,Do_nothing] _r0[12,R3,Recover]
 _r0[57,Do_nothing,Do_nothing] _r0[57,Do_nothing,Recover]
 _r0[57,R1,Do_nothing] _r0[57,R1,Recover] _r0[57,R2,Do_nothing]
 _r0[57,R2,Recover] _r0[57,R3,Do_nothing] _r0[57,R3,Recover]
 _r0[97,Do_nothing,Do_nothing] _r0[97,Do_nothing,Recover]
 _r0[97,R1,Do_nothing] _r0[97,R1,Recover] _r0[97,R2,Do_nothing]
 _r0[97,R2,Recover] _r0[97,R3,Do_nothing] _r0[97,R3,Recover]
 _r0[69,Do_nothing,Do_nothing] _r0[69,Do_nothing,Recover]
 _r0[69,R1,Do_nothing] _r0[69,R1,Recover] _r0[69,R2,Do_nothing]
 _r0[69,R2,Recover] _r0[69,R3,Do_nothing] _r0[69,R3,Recover]
 _r0[3,Do_nothing,Do_nothing] _r0[3,Do_nothing,Recover]
 _r0[3,R1,Do_nothing] _r0[3,R1,Recover] _r0[3,R2,Do_nothing]
 _r0[3,R2,Recover] _r0[3,R3,Do_nothing] _r0[3,R3,Recover]
 _r0[66,Do_nothing,Do_nothing] _r0[66,Do_nothing,Recover]
 _r0[66,R1,Do_nothing] _r0[66,R1,Recover] _r0[66,R2,Do_nothing]
 _r0[66,R2,Recover] _r0[66,R3,Do_nothing] _r0[66,R3,Recover]
 _r1[32,Do_nothing,Do_nothing] _r1[32,Do_nothing,Recover]
 _r1[32,R1,Do_nothing] _r1[32,R1,Recover] _r1[32,R2,Do_nothing]
 _r1[32,R2,Recover] _r1[32,R3,Do_nothing] _r1[32,R3,Recover]
 _r1[70,Do_nothing,Do_nothing] _r1[70,Do_nothing,Recover]
 _r1[70,R1,Do_nothing] _r1[70,R1,Recover] _r1[70,R2,Do_nothing]
 _r1[70,R2,Recover] _r1[70,R3,Do_nothing] _r1[70,R3,Recover]
 _r1[5,Do_nothing,Do_nothing] _r1[5,Do_nothing,Recover]
 _r1[5,R1,Do_nothing] _r1[5,R1,Recover] _r1[5,R2,Do_nothing]
 _r1[5,R2,Recover] _r1[5,R3,Do_nothing] _r1[5,R3,Recover]
 _r1[85,Do_nothing,Do_nothing] _r1[85,Do_nothing,Recover]
 _r1[85,R1,Do_nothing] _r1[85,R1,Recover] _r1[85,R2,Do_nothing]
 _r1[85,R2,Recover] _r1[85,R3,Do_nothing] _r1[85,R3,Recover]
 _r1[71,Do_nothing,Do_nothing] _r1[71,Do_nothing,Recover]
 _r1[71,R1,Do_nothing] _r1[71,R1,Recover] _r1[71,R2,Do_nothing]
 _r1[71,R2,Recover] _r1[71,R3,Do_nothing] _r1[71,R3,Recover]
 _r1[92,Do_nothing,Do_nothing] _r1[92,Do_nothing,Recover]
 _r1[92,R1,Do_nothing] _r1[92,R1,Recover] _r1[92,R2,Do_nothing]
 _r1[92,R2,Recover] _r1[92,R3,Do_nothing] _r1[92,R3,Recover]
 _r1[38,Do_nothing,Do_nothing] _r1[38,Do_nothing,Recover]
 _r1[38,R1,Do_nothing] _r1[38,R1,Recover] _r1[38,R2,Do_nothing]
 _r1[38,R2,Recover] _r1[38,R3,Do_nothing] _r1[38,R3,Recover]
 _r1[21,Do_nothing,Do_nothing] _r1[21,Do_nothing,Recover]
 _r1[21,R1,Do_nothing] _r1[21,R1,Recover] _r1[21,R2,Do_nothing]
 _r1[21,R2,Recover] _r1[21,R3,Do_nothing] _r1[21,R3,Recover]
 _r1[41,Do_nothing,Do_nothing] _r1[41,Do_nothing,Recover]
 _r1[41,R1,Do_nothing] _r1[41,R1,Recover] _r1[41,R2,Do_nothing]
 _r1[41,R2,Recover] _r1[41,R3,Do_nothing] _r1[41,R3,Recover]
 _r1[56,Do_nothing,Do_nothing] _r1[56,Do_nothing,Recover]
 _r1[56,R1,Do_nothing] _r1[56,R1,Recover] _r1[56,R2,Do_nothing]
 _r1[56,R2,Recover] _r1[56,R3,Do_nothing] _r1[56,R3,Recover]
 _r1[80,Do_nothing,Do_nothing] _r1[80,Do_nothing,Recover]
 _r1[80,R1,Do_nothing] _r1[80,R1,Recover] _r1[80,R2,Do_nothing]
 _r1[80,R2,Recover] _r1[80,R3,Do_nothing] _r1[80,R3,Recover]
 _r1[62,Do_nothing,Do_nothing] _r1[62,Do_nothing,Recover]
 _r1[62,R1,Do_nothing] _r1[62,R1,Recover] _r1[62,R2,Do_nothing]
 _r1[62,R2,Recover] _r1[62,R3,Do_nothing] _r1[62,R3,Recover]
 _r1[20,Do_nothing,Do_nothing] _r1[20,Do_nothing,Recover]
 _r1[20,R1,Do_nothing] _r1[20,R1,Recover] _r1[20,R2,Do_nothing]
 _r1[20,R2,Recover] _r1[20,R3,Do_nothing] _r1[20,R3,Recover]
 _r1[10,Do_nothing,Do_nothing] _r1[10,Do_nothing,Recover]
 _r1[10,R1,Do_nothing] _r1[10,R1,Recover] _r1[10,R2,Do_nothing]
 _r1[10,R2,Recover] _r1[10,R3,Do_nothing] _r1[10,R3,Recover]
 _r1[55,Do_nothing,Do_nothing] _r1[55,Do_nothing,Recover]
 _r1[55,R1,Do_nothing] _r1[55,R1,Recover] _r1[55,R2,Do_nothing]
 _r1[55,R2,Recover] _r1[55,R3,Do_nothing] _r1[55,R3,Recover]
 _r1[36,Do_nothing,Do_nothing] _r1[36,Do_nothing,Recover]
 _r1[36,R1,Do_nothing] _r1[36,R1,Recover] _r1[36,R2,Do_nothing]
 _r1[36,R2,Recover] _r1[36,R3,Do_nothing] _r1[36,R3,Recover]
 _r1[44,Do_nothing,Do_nothing] _r1[44,Do_nothing,Recover]
 _r1[44,R1,Do_nothing] _r1[44,R1,Recover] _r1[44,R2,Do_nothing]
 _r1[44,R2,Recover] _r1[44,R3,Do_nothing] _r1[44,R3,Recover]
 _r1[48,Do_nothing,Do_nothing] _r1[48,Do_nothing,Recover]
 _r1[48,R1,Do_nothing] _r1[48,R1,Recover] _r1[48,R2,Do_nothing]
 _r1[48,R2,Recover] _r1[48,R3,Do_nothing] _r1[48,R3,Recover]
 _r1[8,Do_nothing,Do_nothing] _r1[8,Do_nothing,Recover]
 _r1[8,R1,Do_nothing] _r1[8,R1,Recover] _r1[8,R2,Do_nothing]
 _r1[8,R2,Recover] _r1[8,R3,Do_nothing] _r1[8,R3,Recover]
 _r1[15,Do_nothing,Do_nothing] _r1[15,Do_nothing,Recover]
 _r1[15,R1,Do_nothing] _r1[15,R1,Recover] _r1[15,R2,Do_nothing]
 _r1[15,R2,Recover] _r1[15,R3,Do_nothing] _r1[15,R3,Recover]
 _r1[61,Do_nothing,Do_nothing] _r1[61,Do_nothing,Recover]
 _r1[61,R1,Do_nothing] _r1[61,R1,Recover] _r1[61,R2,Do_nothing]
 _r1[61,R2,Recover] _r1[61,R3,Do_nothing] _r1[61,R3,Recover]
 _r1[78,Do_nothing,Do_nothing] _r1[78,Do_nothing,Recover]
 _r1[78,R1,Do_nothing] _r1[78,R1,Recover] _r1[78,R2,Do_nothing]
 _r1[78,R2,Recover] _r1[78,R3,Do_nothing] _r1[78,R3,Recover]
 _r1[34,Do_nothing,Do_nothing] _r1[34,Do_nothing,Recover]
 _r1[34,R1,Do_nothing] _r1[34,R1,Recover] _r1[34,R2,Do_nothing]
 _r1[34,R2,Recover] _r1[34,R3,Do_nothing] _r1[34,R3,Recover]
 _r1[100,Do_nothing,Do_nothing] _r1[100,Do_nothing,Recover]
 _r1[100,R1,Do_nothing] _r1[100,R1,Recover] _r1[100,R2,Do_nothing]
 _r1[100,R2,Recover] _r1[100,R3,Do_nothing] _r1[100,R3,Recover]
 _r1[90,Do_nothing,Do_nothing] _r1[90,Do_nothing,Recover]
 _r1[90,R1,Do_nothing] _r1[90,R1,Recover] _r1[90,R2,Do_nothing]
 _r1[90,R2,Recover] _r1[90,R3,Do_nothing] _r1[90,R3,Recover]
 _r1[64,Do_nothing,Do_nothing] _r1[64,Do_nothing,Recover]
 _r1[64,R1,Do_nothing] _r1[64,R1,Recover] _r1[64,R2,Do_nothing]
 _r1[64,R2,Recover] _r1[64,R3,Do_nothing] _r1[64,R3,Recover]
 _r1[79,Do_nothing,Do_nothing] _r1[79,Do_nothing,Recover]
 _r1[79,R1,Do_nothing] _r1[79,R1,Recover] _r1[79,R2,Do_nothing]
 _r1[79,R2,Recover] _r1[79,R3,Do_nothing] _r1[79,R3,Recover]
 _r1[1,Do_nothing,Do_nothing] _r1[1,Do_nothing,Recover]
 _r1[1,R1,Do_nothing] _r1[1,R1,Recover] _r1[1,R2,Do_nothing]
 _r1[1,R2,Recover] _r1[1,R3,Do_nothing] _r1[1,R3,Recover]
 _r1[22,Do_nothing,Do_nothing] _r1[22,Do_nothing,Recover]
 _r1[22,R1,Do_nothing] _r1[22,R1,Recover] _r1[22,R2,Do_nothing]
 _r1[22,R2,Recover] _r1[22,R3,Do_nothing] _r1[22,R3,Recover]
 _r1[31,Do_nothing,Do_nothing] _r1[31,Do_nothing,Recover]
 _r1[31,R1,Do_nothing] _r1[31,R1,Recover] _r1[31,R2,Do_nothing]
 _r1[31,R2,Recover] _r1[31,R3,Do_nothing] _r1[31,R3,Recover]
 _r1[54,Do_nothing,Do_nothing] _r1[54,Do_nothing,Recover]
 _r1[54,R1,Do_nothing] _r1[54,R1,Recover] _r1[54,R2,Do_nothing]
 _r1[54,R2,Recover] _r1[54,R3,Do_nothing] _r1[54,R3,Recover]
 _r1[82,Do_nothing,Do_nothing] _r1[82,Do_nothing,Recover]
 _r1[82,R1,Do_nothing] _r1[82,R1,Recover] _r1[82,R2,Do_nothing]
 _r1[82,R2,Recover] _r1[82,R3,Do_nothing] _r1[82,R3,Recover]
 _r1[11,Do_nothing,Do_nothing] _r1[11,Do_nothing,Recover]
 _r1[11,R1,Do_nothing] _r1[11,R1,Recover] _r1[11,R2,Do_nothing]
 _r1[11,R2,Recover] _r1[11,R3,Do_nothing] _r1[11,R3,Recover]
 _r1[75,Do_nothing,Do_nothing] _r1[75,Do_nothing,Recover]
 _r1[75,R1,Do_nothing] _r1[75,R1,Recover] _r1[75,R2,Do_nothing]
 _r1[75,R2,Recover] _r1[75,R3,Do_nothing] _r1[75,R3,Recover]
 _r1[50,Do_nothing,Do_nothing] _r1[50,Do_nothing,Recover]
 _r1[50,R1,Do_nothing] _r1[50,R1,Recover] _r1[50,R2,Do_nothing]
 _r1[50,R2,Recover] _r1[50,R3,Do_nothing] _r1[50,R3,Recover]
 _r1[19,Do_nothing,Do_nothing] _r1[19,Do_nothing,Recover]
 _r1[19,R1,Do_nothing] _r1[19,R1,Recover] _r1[19,R2,Do_nothing]
 _r1[19,R2,Recover] _r1[19,R3,Do_nothing] _r1[19,R3,Recover]
 _r1[99,Do_nothing,Do_nothing] _r1[99,Do_nothing,Recover]
 _r1[99,R1,Do_nothing] _r1[99,R1,Recover] _r1[99,R2,Do_nothing]
 _r1[99,R2,Recover] _r1[99,R3,Do_nothing] _r1[99,R3,Recover]
 _r1[72,Do_nothing,Do_nothing] _r1[72,Do_nothing,Recover]
 _r1[72,R1,Do_nothing] _r1[72,R1,Recover] _r1[72,R2,Do_nothing]
 _r1[72,R2,Recover] _r1[72,R3,Do_nothing] _r1[72,R3,Recover]
 _r1[47,Do_nothing,Do_nothing] _r1[47,Do_nothing,Recover]
 _r1[47,R1,Do_nothing] _r1[47,R1,Recover] _r1[47,R2,Do_nothing]
 _r1[47,R2,Recover] _r1[47,R3,Do_nothing] _r1[47,R3,Recover]
 _r1[49,Do_nothing,Do_nothing] _r1[49,Do_nothing,Recover]
 _r1[49,R1,Do_nothing] _r1[49,R1,Recover] _r1[49,R2,Do_nothing]
 _r1[49,R2,Recover] _r1[49,R3,Do_nothing] _r1[49,R3,Recover]
 _r1[77,Do_nothing,Do_nothing] _r1[77,Do_nothing,Recover]
 _r1[77,R1,Do_nothing] _r1[77,R1,Recover] _r1[77,R2,Do_nothing]
 _r1[77,R2,Recover] _r1[77,R3,Do_nothing] _r1[77,R3,Recover]
 _r1[95,Do_nothing,Do_nothing] _r1[95,Do_nothing,Recover]
 _r1[95,R1,Do_nothing] _r1[95,R1,Recover] _r1[95,R2,Do_nothing]
 _r1[95,R2,Recover] _r1[95,R3,Do_nothing] _r1[95,R3,Recover]
 _r1[98,Do_nothing,Do_nothing] _r1[98,Do_nothing,Recover]
 _r1[98,R1,Do_nothing] _r1[98,R1,Recover] _r1[98,R2,Do_nothing]
 _r1[98,R2,Recover] _r1[98,R3,Do_nothing] _r1[98,R3,Recover]
 _r1[53,Do_nothing,Do_nothing] _r1[53,Do_nothing,Recover]
 _r1[53,R1,Do_nothing] _r1[53,R1,Recover] _r1[53,R2,Do_nothing]
 _r1[53,R2,Recover] _r1[53,R3,Do_nothing] _r1[53,R3,Recover]
 _r1[65,Do_nothing,Do_nothing] _r1[65,Do_nothing,Recover]
 _r1[65,R1,Do_nothing] _r1[65,R1,Recover] _r1[65,R2,Do_nothing]
 _r1[65,R2,Recover] _r1[65,R3,Do_nothing] _r1[65,R3,Recover]
 _r1[63,Do_nothing,Do_nothing] _r1[63,Do_nothing,Recover]
 _r1[63,R1,Do_nothing] _r1[63,R1,Recover] _r1[63,R2,Do_nothing]
 _r1[63,R2,Recover] _r1[63,R3,Do_nothing] _r1[63,R3,Recover]
 _r1[14,Do_nothing,Do_nothing] _r1[14,Do_nothing,Recover]
 _r1[14,R1,Do_nothing] _r1[14,R1,Recover] _r1[14,R2,Do_nothing]
 _r1[14,R2,Recover] _r1[14,R3,Do_nothing] _r1[14,R3,Recover]
 _r1[23,Do_nothing,Do_nothing] _r1[23,Do_nothing,Recover]
 _r1[23,R1,Do_nothing] _r1[23,R1,Recover] _r1[23,R2,Do_nothing]
 _r1[23,R2,Recover] _r1[23,R3,Do_nothing] _r1[23,R3,Recover]
 _r1[93,Do_nothing,Do_nothing] _r1[93,Do_nothing,Recover]
 _r1[93,R1,Do_nothing] _r1[93,R1,Recover] _r1[93,R2,Do_nothing]
 _r1[93,R2,Recover] _r1[93,R3,Do_nothing] _r1[93,R3,Recover]
 _r1[28,Do_nothing,Do_nothing] _r1[28,Do_nothing,Recover]
 _r1[28,R1,Do_nothing] _r1[28,R1,Recover] _r1[28,R2,Do_nothing]
 _r1[28,R2,Recover] _r1[28,R3,Do_nothing] _r1[28,R3,Recover]
 _r1[2,Do_nothing,Do_nothing] _r1[2,Do_nothing,Recover]
 _r1[2,R1,Do_nothing] _r1[2,R1,Recover] _r1[2,R2,Do_nothing]
 _r1[2,R2,Recover] _r1[2,R3,Do_nothing] _r1[2,R3,Recover]
 _r1[81,Do_nothing,Do_nothing] _r1[81,Do_nothing,Recover]
 _r1[81,R1,Do_nothing] _r1[81,R1,Recover] _r1[81,R2,Do_nothing]
 _r1[81,R2,Recover] _r1[81,R3,Do_nothing] _r1[81,R3,Recover]
 _r1[45,Do_nothing,Do_nothing] _r1[45,Do_nothing,Recover]
 _r1[45,R1,Do_nothing] _r1[45,R1,Recover] _r1[45,R2,Do_nothing]
 _r1[45,R2,Recover] _r1[45,R3,Do_nothing] _r1[45,R3,Recover]
 _r1[88,Do_nothing,Do_nothing] _r1[88,Do_nothing,Recover]
 _r1[88,R1,Do_nothing] _r1[88,R1,Recover] _r1[88,R2,Do_nothing]
 _r1[88,R2,Recover] _r1[88,R3,Do_nothing] _r1[88,R3,Recover]
 _r1[89,Do_nothing,Do_nothing] _r1[89,Do_nothing,Recover]
 _r1[89,R1,Do_nothing] _r1[89,R1,Recover] _r1[89,R2,Do_nothing]
 _r1[89,R2,Recover] _r1[89,R3,Do_nothing] _r1[89,R3,Recover]
 _r1[25,Do_nothing,Do_nothing] _r1[25,Do_nothing,Recover]
 _r1[25,R1,Do_nothing] _r1[25,R1,Recover] _r1[25,R2,Do_nothing]
 _r1[25,R2,Recover] _r1[25,R3,Do_nothing] _r1[25,R3,Recover]
 _r1[91,Do_nothing,Do_nothing] _r1[91,Do_nothing,Recover]
 _r1[91,R1,Do_nothing] _r1[91,R1,Recover] _r1[91,R2,Do_nothing]
 _r1[91,R2,Recover] _r1[91,R3,Do_nothing] _r1[91,R3,Recover]
 _r1[43,Do_nothing,Do_nothing] _r1[43,Do_nothing,Recover]
 _r1[43,R1,Do_nothing] _r1[43,R1,Recover] _r1[43,R2,Do_nothing]
 _r1[43,R2,Recover] _r1[43,R3,Do_nothing] _r1[43,R3,Recover]
 _r1[76,Do_nothing,Do_nothing] _r1[76,Do_nothing,Recover]
 _r1[76,R1,Do_nothing] _r1[76,R1,Recover] _r1[76,R2,Do_nothing]
 _r1[76,R2,Recover] _r1[76,R3,Do_nothing] _r1[76,R3,Recover]
 _r1[42,Do_nothing,Do_nothing] _r1[42,Do_nothing,Recover]
 _r1[42,R1,Do_nothing] _r1[42,R1,Recover] _r1[42,R2,Do_nothing]
 _r1[42,R2,Recover] _r1[42,R3,Do_nothing] _r1[42,R3,Recover]
 _r1[7,Do_nothing,Do_nothing] _r1[7,Do_nothing,Recover]
 _r1[7,R1,Do_nothing] _r1[7,R1,Recover] _r1[7,R2,Do_nothing]
 _r1[7,R2,Recover] _r1[7,R3,Do_nothing] _r1[7,R3,Recover]
 _r1[17,Do_nothing,Do_nothing] _r1[17,Do_nothing,Recover]
 _r1[17,R1,Do_nothing] _r1[17,R1,Recover] _r1[17,R2,Do_nothing]
 _r1[17,R2,Recover] _r1[17,R3,Do_nothing] _r1[17,R3,Recover]
 _r1[37,Do_nothing,Do_nothing] _r1[37,Do_nothing,Recover]
 _r1[37,R1,Do_nothing] _r1[37,R1,Recover] _r1[37,R2,Do_nothing]
 _r1[37,R2,Recover] _r1[37,R3,Do_nothing] _r1[37,R3,Recover]
 _r1[26,Do_nothing,Do_nothing] _r1[26,Do_nothing,Recover]
 _r1[26,R1,Do_nothing] _r1[26,R1,Recover] _r1[26,R2,Do_nothing]
 _r1[26,R2,Recover] _r1[26,R3,Do_nothing] _r1[26,R3,Recover]
 _r1[9,Do_nothing,Do_nothing] _r1[9,Do_nothing,Recover]
 _r1[9,R1,Do_nothing] _r1[9,R1,Recover] _r1[9,R2,Do_nothing]
 _r1[9,R2,Recover] _r1[9,R3,Do_nothing] _r1[9,R3,Recover]
 _r1[60,Do_nothing,Do_nothing] _r1[60,Do_nothing,Recover]
 _r1[60,R1,Do_nothing] _r1[60,R1,Recover] _r1[60,R2,Do_nothing]
 _r1[60,R2,Recover] _r1[60,R3,Do_nothing] _r1[60,R3,Recover]
 _r1[16,Do_nothing,Do_nothing] _r1[16,Do_nothing,Recover]
 _r1[16,R1,Do_nothing] _r1[16,R1,Recover] _r1[16,R2,Do_nothing]
 _r1[16,R2,Recover] _r1[16,R3,Do_nothing] _r1[16,R3,Recover]
 _r1[83,Do_nothing,Do_nothing] _r1[83,Do_nothing,Recover]
 _r1[83,R1,Do_nothing] _r1[83,R1,Recover] _r1[83,R2,Do_nothing]
 _r1[83,R2,Recover] _r1[83,R3,Do_nothing] _r1[83,R3,Recover]
 _r1[13,Do_nothing,Do_nothing] _r1[13,Do_nothing,Recover]
 _r1[13,R1,Do_nothing] _r1[13,R1,Recover] _r1[13,R2,Do_nothing]
 _r1[13,R2,Recover] _r1[13,R3,Do_nothing] _r1[13,R3,Recover]
 _r1[67,Do_nothing,Do_nothing] _r1[67,Do_nothing,Recover]
 _r1[67,R1,Do_nothing] _r1[67,R1,Recover] _r1[67,R2,Do_nothing]
 _r1[67,R2,Recover] _r1[67,R3,Do_nothing] _r1[67,R3,Recover]
 _r1[59,Do_nothing,Do_nothing] _r1[59,Do_nothing,Recover]
 _r1[59,R1,Do_nothing] _r1[59,R1,Recover] _r1[59,R2,Do_nothing]
 _r1[59,R2,Recover] _r1[59,R3,Do_nothing] _r1[59,R3,Recover]
 _r1[58,Do_nothing,Do_nothing] _r1[58,Do_nothing,Recover]
 _r1[58,R1,Do_nothing] _r1[58,R1,Recover] _r1[58,R2,Do_nothing]
 _r1[58,R2,Recover] _r1[58,R3,Do_nothing] _r1[58,R3,Recover]
 _r1[4,Do_nothing,Do_nothing] _r1[4,Do_nothing,Recover]
 _r1[4,R1,Do_nothing] _r1[4,R1,Recover] _r1[4,R2,Do_nothing]
 _r1[4,R2,Recover] _r1[4,R3,Do_nothing] _r1[4,R3,Recover]
 _r1[68,Do_nothing,Do_nothing] _r1[68,Do_nothing,Recover]
 _r1[68,R1,Do_nothing] _r1[68,R1,Recover] _r1[68,R2,Do_nothing]
 _r1[68,R2,Recover] _r1[68,R3,Do_nothing] _r1[68,R3,Recover]
 _r1[87,Do_nothing,Do_nothing] _r1[87,Do_nothing,Recover]
 _r1[87,R1,Do_nothing] _r1[87,R1,Recover] _r1[87,R2,Do_nothing]
 _r1[87,R2,Recover] _r1[87,R3,Do_nothing] _r1[87,R3,Recover]
 _r1[40,Do_nothing,Do_nothing] _r1[40,Do_nothing,Recover]
 _r1[40,R1,Do_nothing] _r1[40,R1,Recover] _r1[40,R2,Do_nothing]
 _r1[40,R2,Recover] _r1[40,R3,Do_nothing] _r1[40,R3,Recover]
 _r1[30,Do_nothing,Do_nothing] _r1[30,Do_nothing,Recover]
 _r1[30,R1,Do_nothing] _r1[30,R1,Recover] _r1[30,R2,Do_nothing]
 _r1[30,R2,Recover] _r1[30,R3,Do_nothing] _r1[30,R3,Recover]
 _r1[27,Do_nothing,Do_nothing] _r1[27,Do_nothing,Recover]
 _r1[27,R1,Do_nothing] _r1[27,R1,Recover] _r1[27,R2,Do_nothing]
 _r1[27,R2,Recover] _r1[27,R3,Do_nothing] _r1[27,R3,Recover]
 _r1[73,Do_nothing,Do_nothing] _r1[73,Do_nothing,Recover]
 _r1[73,R1,Do_nothing] _r1[73,R1,Recover] _r1[73,R2,Do_nothing]
 _r1[73,R2,Recover] _r1[73,R3,Do_nothing] _r1[73,R3,Recover]
 _r1[84,Do_nothing,Do_nothing] _r1[84,Do_nothing,Recover]
 _r1[84,R1,Do_nothing] _r1[84,R1,Recover] _r1[84,R2,Do_nothing]
 _r1[84,R2,Recover] _r1[84,R3,Do_nothing] _r1[84,R3,Recover]
 _r1[35,Do_nothing,Do_nothing] _r1[35,Do_nothing,Recover]
 _r1[35,R1,Do_nothing] _r1[35,R1,Recover] _r1[35,R2,Do_nothing]
 _r1[35,R2,Recover] _r1[35,R3,Do_nothing] _r1[35,R3,Recover]
 _r1[46,Do_nothing,Do_nothing] _r1[46,Do_nothing,Recover]
 _r1[46,R1,Do_nothing] _r1[46,R1,Recover] _r1[46,R2,Do_nothing]
 _r1[46,R2,Recover] _r1[46,R3,Do_nothing] _r1[46,R3,Recover]
 _r1[39,Do_nothing,Do_nothing] _r1[39,Do_nothing,Recover]
 _r1[39,R1,Do_nothing] _r1[39,R1,Recover] _r1[39,R2,Do_nothing]
 _r1[39,R2,Recover] _r1[39,R3,Do_nothing] _r1[39,R3,Recover]
 _r1[51,Do_nothing,Do_nothing] _r1[51,Do_nothing,Recover]
 _r1[51,R1,Do_nothing] _r1[51,R1,Recover] _r1[51,R2,Do_nothing]
 _r1[51,R2,Recover] _r1[51,R3,Do_nothing] _r1[51,R3,Recover]
 _r1[52,Do_nothing,Do_nothing] _r1[52,Do_nothing,Recover]
 _r1[52,R1,Do_nothing] _r1[52,R1,Recover] _r1[52,R2,Do_nothing]
 _r1[52,R2,Recover] _r1[52,R3,Do_nothing] _r1[52,R3,Recover]
 _r1[6,Do_nothing,Do_nothing] _r1[6,Do_nothing,Recover]
 _r1[6,R1,Do_nothing] _r1[6,R1,Recover] _r1[6,R2,Do_nothing]
 _r1[6,R2,Recover] _r1[6,R3,Do_nothing] _r1[6,R3,Recover]
 _r1[96,Do_nothing,Do_nothing] _r1[96,Do_nothing,Recover]
 _r1[96,R1,Do_nothing] _r1[96,R1,Recover] _r1[96,R2,Do_nothing]
 _r1[96,R2,Recover] _r1[96,R3,Do_nothing] _r1[96,R3,Recover]
 _r1[74,Do_nothing,Do_nothing] _r1[74,Do_nothing,Recover]
 _r1[74,R1,Do_nothing] _r1[74,R1,Recover] _r1[74,R2,Do_nothing]
 _r1[74,R2,Recover] _r1[74,R3,Do_nothing] _r1[74,R3,Recover]
 _r1[86,Do_nothing,Do_nothing] _r1[86,Do_nothing,Recover]
 _r1[86,R1,Do_nothing] _r1[86,R1,Recover] _r1[86,R2,Do_nothing]
 _r1[86,R2,Recover] _r1[86,R3,Do_nothing] _r1[86,R3,Recover]
 _r1[29,Do_nothing,Do_nothing] _r1[29,Do_nothing,Recover]
 _r1[29,R1,Do_nothing] _r1[29,R1,Recover] _r1[29,R2,Do_nothing]
 _r1[29,R2,Recover] _r1[29,R3,Do_nothing] _r1[29,R3,Recover]
 _r1[18,Do_nothing,Do_nothing] _r1[18,Do_nothing,Recover]
 _r1[18,R1,Do_nothing] _r1[18,R1,Recover] _r1[18,R2,Do_nothing]
 _r1[18,R2,Recover] _r1[18,R3,Do_nothing] _r1[18,R3,Recover]
 _r1[24,Do_nothing,Do_nothing] _r1[24,Do_nothing,Recover]
 _r1[24,R1,Do_nothing] _r1[24,R1,Recover] _r1[24,R2,Do_nothing]
 _r1[24,R2,Recover] _r1[24,R3,Do_nothing] _r1[24,R3,Recover]
 _r1[94,Do_nothing,Do_nothing] _r1[94,Do_nothing,Recover]
 _r1[94,R1,Do_nothing] _r1[94,R1,Recover] _r1[94,R2,Do_nothing]
 _r1[94,R2,Recover] _r1[94,R3,Do_nothing] _r1[94,R3,Recover]
 _r1[33,Do_nothing,Do_nothing] _r1[33,Do_nothing,Recover]
 _r1[33,R1,Do_nothing] _r1[33,R1,Recover] _r1[33,R2,Do_nothing]
 _r1[33,R2,Recover] _r1[33,R3,Do_nothing] _r1[33,R3,Recover]
 _r1[12,Do_nothing,Do_nothing] _r1[12,Do_nothing,Recover]
 _r1[12,R1,Do_nothing] _r1[12,R1,Recover] _r1[12,R2,Do_nothing]
 _r1[12,R2,Recover] _r1[12,R3,Do_nothing] _r1[12,R3,Recover]
 _r1[57,Do_nothing,Do_nothing] _r1[57,Do_nothing,Recover]
 _r1[57,R1,Do_nothing] _r1[57,R1,Recover] _r1[57,R2,Do_nothing]
 _r1[57,R2,Recover] _r1[57,R3,Do_nothing] _r1[57,R3,Recover]
 _r1[97,Do_nothing,Do_nothing] _r1[97,Do_nothing,Recover]
 _r1[97,R1,Do_nothing] _r1[97,R1,Recover] _r1[97,R2,Do_nothing]
 _r1[97,R2,Recover] _r1[97,R3,Do_nothing] _r1[97,R3,Recover]
 _r1[69,Do_nothing,Do_nothing] _r1[69,Do_nothing,Recover]
 _r1[69,R1,Do_nothing] _r1[69,R1,Recover] _r1[69,R2,Do_nothing]
 _r1[69,R2,Recover] _r1[69,R3,Do_nothing] _r1[69,R3,Recover]
 _r1[3,Do_nothing,Do_nothing] _r1[3,Do_nothing,Recover]
 _r1[3,R1,Do_nothing] _r1[3,R1,Recover] _r1[3,R2,Do_nothing]
 _r1[3,R2,Recover] _r1[3,R3,Do_nothing] _r1[3,R3,Recover]
 _r1[66,Do_nothing,Do_nothing] _r1[66,Do_nothing,Recover]
 _r1[66,R1,Do_nothing] _r1[66,R1,Recover] _r1[66,R2,Do_nothing]
 _r1[66,R2,Recover] _r1[66,R3,Do_nothing] _r1[66,R3,Recover]
End
