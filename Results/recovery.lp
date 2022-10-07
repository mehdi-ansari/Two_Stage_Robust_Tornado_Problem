\ Model recoveryProblem_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  25 _r[290970107002023,Do_nothing,Do_nothing]
   + 25 _r[290970107002023,Do_nothing,Recover]
   + 25 _r[290970107002023,R1,Do_nothing]
   + 25 _r[290970107002023,R1,Recover]
   + 22 _r[290970107002023,R2,Do_nothing]
   + 19 _r[290970107002023,R2,Recover]
   + 20 _r[290970107002023,R3,Do_nothing]
   + 16 _r[290970107002023,R3,Recover]
   + 17 _r[290970107002022,Do_nothing,Do_nothing]
   + 17 _r[290970107002022,Do_nothing,Recover]
   + 17 _r[290970107002022,R1,Do_nothing]
   + 17 _r[290970107002022,R1,Recover]
   + 15 _r[290970107002022,R2,Do_nothing]
   + 13 _r[290970107002022,R2,Recover]
   + 14 _r[290970107002022,R3,Do_nothing]
   + 11 _r[290970107002022,R3,Recover]
   + 21 _r[291450205022018,Do_nothing,Do_nothing]
   + 21 _r[291450205022018,Do_nothing,Recover]
   + 21 _r[291450205022018,R1,Do_nothing]
   + 21 _r[291450205022018,R1,Recover]
   + 18 _r[291450205022018,R2,Do_nothing]
   + 16 _r[291450205022018,R2,Recover]
   + 17 _r[291450205022018,R3,Do_nothing]
   + 13 _r[291450205022018,R3,Recover]
   + 26 _r[290970106006020,Do_nothing,Do_nothing]
   + 26 _r[290970106006020,Do_nothing,Recover]
   + 26 _r[290970106006020,R1,Do_nothing]
   + 26 _r[290970106006020,R1,Recover]
   + 23 _r[290970106006020,R2,Do_nothing]
   + 19 _r[290970106006020,R2,Recover]
   + 21 _r[290970106006020,R3,Do_nothing]
   + 16 _r[290970106006020,R3,Recover]
Subject To
 R0: 273257 _r[290970112004035,Do_nothing,Recover]
   + 273175 _r[290970112004035,R1,Recover]
   + 195993 _r[290970112004035,R2,Recover]
   + 119326 _r[290970112004035,R3,Recover]
   + 248875 _r[290970112004038,Do_nothing,Recover]
   + 248801 _r[290970112004038,R1,Recover]
   + 178505 _r[290970112004038,R2,Recover]
   + 108679 _r[290970112004038,R3,Recover]
   + 466694 _r[290970107002023,Do_nothing,Recover]
   + 466555 _r[290970107002023,R1,Recover]
   + 334735 _r[290970107002023,R2,Recover]
   + 203796 _r[290970107002023,R3,Recover]
   + 459751 _r[290970107002022,Do_nothing,Recover]
   + 459614 _r[290970107002022,R1,Recover]
   + 329755 _r[290970107002022,R2,Recover]
   + 200764 _r[290970107002022,R3,Recover]
   + 2.807261e+06 _r[291450205022018,Do_nothing,Recover]
   + 2.806424e+06 _r[291450205022018,R1,Recover]
   + 2.013498e+06 _r[291450205022018,R2,Recover]
   + 1.225874e+06 _r[291450205022018,R3,Recover]
   + 228902 _r[290970101002004,Do_nothing,Recover]
   + 228834 _r[290970101002004,R1,Recover]
   + 164179 _r[290970101002004,R2,Recover]
   + 99957 _r[290970101002004,R3,Recover]
   + 500834 _r[290970103001034,Do_nothing,Recover]
   + 500684 _r[290970103001034,R1,Recover]
   + 359221 _r[290970103001034,R2,Recover]
   + 218704 _r[290970103001034,R3,Recover]
   + 86636 _r[290970101004032,Do_nothing,Recover]
   + 86610 _r[290970101004032,R1,Recover]
   + 62139 _r[290970101004032,R2,Recover]
   + 37832 _r[290970101004032,R3,Recover]
   + 485418 _r[290970106006020,Do_nothing,Recover]
   + 485273 _r[290970106006020,R1,Recover]
   + 348164 _r[290970106006020,R2,Recover]
   + 211972 _r[290970106006020,R3,Recover]
   + 291626 _r[290970104001042,Do_nothing,Recover]
   + 291539 _r[290970104001042,R1,Recover]
   + 209168 _r[290970104001042,R2,Recover]
   + 127347 _r[290970104001042,R3,Recover] <= 3e+07
 R1: _r[290970112004035,Do_nothing,Do_nothing]
   + _r[290970112004035,Do_nothing,Recover] = 1
 R2: _r[290970112004035,R1,Do_nothing] + _r[290970112004035,R1,Recover]
   = 0
 R3: _r[290970112004035,R2,Do_nothing] + _r[290970112004035,R2,Recover]
   = 0
 R4: _r[290970112004035,R3,Do_nothing] + _r[290970112004035,R3,Recover]
   = 0
 R5: _r[290970112004038,Do_nothing,Do_nothing]
   + _r[290970112004038,Do_nothing,Recover] = 1
 R6: _r[290970112004038,R1,Do_nothing] + _r[290970112004038,R1,Recover]
   = 0
 R7: _r[290970112004038,R2,Do_nothing] + _r[290970112004038,R2,Recover]
   = 0
 R8: _r[290970112004038,R3,Do_nothing] + _r[290970112004038,R3,Recover]
   = 0
 R9: _r[290970107002023,Do_nothing,Do_nothing]
   + _r[290970107002023,Do_nothing,Recover] = 1
 R10: _r[290970107002023,R1,Do_nothing] + _r[290970107002023,R1,Recover]
   = 0
 R11: _r[290970107002023,R2,Do_nothing] + _r[290970107002023,R2,Recover]
   = 0
 R12: _r[290970107002023,R3,Do_nothing] + _r[290970107002023,R3,Recover]
   = 0
 R13: _r[290970107002022,Do_nothing,Do_nothing]
   + _r[290970107002022,Do_nothing,Recover] = 1
 R14: _r[290970107002022,R1,Do_nothing] + _r[290970107002022,R1,Recover]
   = 0
 R15: _r[290970107002022,R2,Do_nothing] + _r[290970107002022,R2,Recover]
   = 0
 R16: _r[290970107002022,R3,Do_nothing] + _r[290970107002022,R3,Recover]
   = 0
 R17: _r[291450205022018,Do_nothing,Do_nothing]
   + _r[291450205022018,Do_nothing,Recover] = 1
 R18: _r[291450205022018,R1,Do_nothing] + _r[291450205022018,R1,Recover]
   = 0
 R19: _r[291450205022018,R2,Do_nothing] + _r[291450205022018,R2,Recover]
   = 0
 R20: _r[291450205022018,R3,Do_nothing] + _r[291450205022018,R3,Recover]
   = 0
 R21: _r[290970101002004,Do_nothing,Do_nothing]
   + _r[290970101002004,Do_nothing,Recover] = 1
 R22: _r[290970101002004,R1,Do_nothing] + _r[290970101002004,R1,Recover]
   = 0
 R23: _r[290970101002004,R2,Do_nothing] + _r[290970101002004,R2,Recover]
   = 0
 R24: _r[290970101002004,R3,Do_nothing] + _r[290970101002004,R3,Recover]
   = 0
 R25: _r[290970103001034,Do_nothing,Do_nothing]
   + _r[290970103001034,Do_nothing,Recover] = 1
 R26: _r[290970103001034,R1,Do_nothing] + _r[290970103001034,R1,Recover]
   = 0
 R27: _r[290970103001034,R2,Do_nothing] + _r[290970103001034,R2,Recover]
   = 0
 R28: _r[290970103001034,R3,Do_nothing] + _r[290970103001034,R3,Recover]
   = 0
 R29: _r[290970101004032,Do_nothing,Do_nothing]
   + _r[290970101004032,Do_nothing,Recover] = 1
 R30: _r[290970101004032,R1,Do_nothing] + _r[290970101004032,R1,Recover]
   = 0
 R31: _r[290970101004032,R2,Do_nothing] + _r[290970101004032,R2,Recover]
   = 0
 R32: _r[290970101004032,R3,Do_nothing] + _r[290970101004032,R3,Recover]
   = 0
 R33: _r[290970106006020,Do_nothing,Do_nothing]
   + _r[290970106006020,Do_nothing,Recover] = 1
 R34: _r[290970106006020,R1,Do_nothing] + _r[290970106006020,R1,Recover]
   = 0
 R35: _r[290970106006020,R2,Do_nothing] + _r[290970106006020,R2,Recover]
   = 0
 R36: _r[290970106006020,R3,Do_nothing] + _r[290970106006020,R3,Recover]
   = 0
 R37: _r[290970104001042,Do_nothing,Do_nothing]
   + _r[290970104001042,Do_nothing,Recover] = 1
 R38: _r[290970104001042,R1,Do_nothing] + _r[290970104001042,R1,Recover]
   = 0
 R39: _r[290970104001042,R2,Do_nothing] + _r[290970104001042,R2,Recover]
   = 0
 R40: _r[290970104001042,R3,Do_nothing] + _r[290970104001042,R3,Recover]
   = 0
Bounds
Binaries
 _r[290970112004035,Do_nothing,Do_nothing]
 _r[290970112004035,Do_nothing,Recover] _r[290970112004035,R1,Do_nothing]
 _r[290970112004035,R1,Recover] _r[290970112004035,R2,Do_nothing]
 _r[290970112004035,R2,Recover] _r[290970112004035,R3,Do_nothing]
 _r[290970112004035,R3,Recover] _r[290970112004038,Do_nothing,Do_nothing]
 _r[290970112004038,Do_nothing,Recover] _r[290970112004038,R1,Do_nothing]
 _r[290970112004038,R1,Recover] _r[290970112004038,R2,Do_nothing]
 _r[290970112004038,R2,Recover] _r[290970112004038,R3,Do_nothing]
 _r[290970112004038,R3,Recover] _r[290970107002023,Do_nothing,Do_nothing]
 _r[290970107002023,Do_nothing,Recover] _r[290970107002023,R1,Do_nothing]
 _r[290970107002023,R1,Recover] _r[290970107002023,R2,Do_nothing]
 _r[290970107002023,R2,Recover] _r[290970107002023,R3,Do_nothing]
 _r[290970107002023,R3,Recover] _r[290970107002022,Do_nothing,Do_nothing]
 _r[290970107002022,Do_nothing,Recover] _r[290970107002022,R1,Do_nothing]
 _r[290970107002022,R1,Recover] _r[290970107002022,R2,Do_nothing]
 _r[290970107002022,R2,Recover] _r[290970107002022,R3,Do_nothing]
 _r[290970107002022,R3,Recover] _r[291450205022018,Do_nothing,Do_nothing]
 _r[291450205022018,Do_nothing,Recover] _r[291450205022018,R1,Do_nothing]
 _r[291450205022018,R1,Recover] _r[291450205022018,R2,Do_nothing]
 _r[291450205022018,R2,Recover] _r[291450205022018,R3,Do_nothing]
 _r[291450205022018,R3,Recover] _r[290970101002004,Do_nothing,Do_nothing]
 _r[290970101002004,Do_nothing,Recover] _r[290970101002004,R1,Do_nothing]
 _r[290970101002004,R1,Recover] _r[290970101002004,R2,Do_nothing]
 _r[290970101002004,R2,Recover] _r[290970101002004,R3,Do_nothing]
 _r[290970101002004,R3,Recover] _r[290970103001034,Do_nothing,Do_nothing]
 _r[290970103001034,Do_nothing,Recover] _r[290970103001034,R1,Do_nothing]
 _r[290970103001034,R1,Recover] _r[290970103001034,R2,Do_nothing]
 _r[290970103001034,R2,Recover] _r[290970103001034,R3,Do_nothing]
 _r[290970103001034,R3,Recover] _r[290970101004032,Do_nothing,Do_nothing]
 _r[290970101004032,Do_nothing,Recover] _r[290970101004032,R1,Do_nothing]
 _r[290970101004032,R1,Recover] _r[290970101004032,R2,Do_nothing]
 _r[290970101004032,R2,Recover] _r[290970101004032,R3,Do_nothing]
 _r[290970101004032,R3,Recover] _r[290970106006020,Do_nothing,Do_nothing]
 _r[290970106006020,Do_nothing,Recover] _r[290970106006020,R1,Do_nothing]
 _r[290970106006020,R1,Recover] _r[290970106006020,R2,Do_nothing]
 _r[290970106006020,R2,Recover] _r[290970106006020,R3,Do_nothing]
 _r[290970106006020,R3,Recover] _r[290970104001042,Do_nothing,Do_nothing]
 _r[290970104001042,Do_nothing,Recover] _r[290970104001042,R1,Do_nothing]
 _r[290970104001042,R1,Recover] _r[290970104001042,R2,Do_nothing]
 _r[290970104001042,R2,Recover] _r[290970104001042,R3,Do_nothing]
 _r[290970104001042,R3,Recover]
End
