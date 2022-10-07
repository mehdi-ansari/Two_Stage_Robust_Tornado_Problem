\ Model recoveryProblem_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  32 _r[291450206011012,Do_nothing,Do_nothing]
   + 31 _r[291450206011012,Do_nothing,Recover]
   + 32 _r[291450206011012,R1,Do_nothing]
   + 31 _r[291450206011012,R1,Recover]
   + 28 _r[291450206011012,R2,Do_nothing]
   + 24 _r[291450206011012,R2,Recover]
   + 26 _r[291450206011012,R3,Do_nothing]
   + 20 _r[291450206011012,R3,Recover]
   + 19 _r[290970109002051,Do_nothing,Do_nothing]
   + 19 _r[290970109002051,Do_nothing,Recover]
   + 19 _r[290970109002051,R1,Do_nothing]
   + 19 _r[290970109002051,R1,Recover]
   + 17 _r[290970109002051,R2,Do_nothing]
   + 14 _r[290970109002051,R2,Recover]
   + 15 _r[290970109002051,R3,Do_nothing]
   + 12 _r[290970109002051,R3,Recover]
   + 30 _r[290970110005024,Do_nothing,Do_nothing]
   + 29 _r[290970110005024,Do_nothing,Recover]
   + 30 _r[290970110005024,R1,Do_nothing]
   + 29 _r[290970110005024,R1,Recover]
   + 26 _r[290970110005024,R2,Do_nothing]
   + 22 _r[290970110005024,R2,Recover]
   + 24 _r[290970110005024,R3,Do_nothing]
   + 19 _r[290970110005024,R3,Recover]
   + 29 _r[290970108003014,Do_nothing,Do_nothing]
   + 28 _r[290970108003014,Do_nothing,Recover]
   + 29 _r[290970108003014,R1,Do_nothing]
   + 28 _r[290970108003014,R1,Recover]
   + 25 _r[290970108003014,R2,Do_nothing]
   + 22 _r[290970108003014,R2,Recover]
   + 24 _r[290970108003014,R3,Do_nothing]
   + 18 _r[290970108003014,R3,Recover]
Subject To
 R0: 1.340499e+06 _r[291450206011012,Do_nothing,Recover]
   + 1.340099e+06 _r[291450206011012,R1,Recover]
   + 961468 _r[291450206011012,R2,Recover]
   + 585369 _r[291450206011012,R3,Recover]
   + 261190 _r[290970112004069,Do_nothing,Recover]
   + 261112 _r[290970112004069,R1,Recover]
   + 187337 _r[290970112004069,R2,Recover]
   + 114056 _r[290970112004069,R3,Recover]
   + 684624 _r[290970109002051,Do_nothing,Recover]
   + 684420 _r[290970109002051,R1,Recover]
   + 491044 _r[290970109002051,R2,Recover]
   + 298961 _r[290970109002051,R3,Recover]
   + 654504 _r[290970110005024,Do_nothing,Recover]
   + 654309 _r[290970110005024,R1,Recover]
   + 469441 _r[290970110005024,R2,Recover]
   + 285809 _r[290970110005024,R3,Recover]
   + 701538 _r[290970108003014,Do_nothing,Recover]
   + 701329 _r[290970108003014,R1,Recover]
   + 503176 _r[290970108003014,R2,Recover]
   + 306347 _r[290970108003014,R3,Recover]
   + 400526 _r[290970101001022,Do_nothing,Recover]
   + 400407 _r[290970101001022,R1,Recover]
   + 287276 _r[290970101001022,R2,Recover]
   + 174902 _r[290970101001022,R3,Recover]
   + 590995 _r[290970106004001,Do_nothing,Recover]
   + 590818 _r[290970106004001,R1,Recover]
   + 423889 _r[290970106004001,R2,Recover]
   + 258075 _r[290970106004001,R3,Recover]
   + 285251 _r[290970106003025,Do_nothing,Recover]
   + 285166 _r[290970106003025,R1,Recover]
   + 204596 _r[290970106003025,R2,Recover]
   + 124563 _r[290970106003025,R3,Recover]
   + 1.632938e+06 _r[290970102001050,Do_nothing,Recover]
   + 1.632451e+06 _r[290970102001050,R1,Recover]
   + 1.171219e+06 _r[290970102001050,R2,Recover]
   + 713071 _r[290970102001050,R3,Recover]
   + 1.202945e+06 _r[290970104003071,Do_nothing,Recover]
   + 1.202586e+06 _r[290970104003071,R1,Recover]
   + 862808 _r[290970104003071,R2,Recover]
   + 525302 _r[290970104003071,R3,Recover] <= 3e+07
 R1: _r[291450206011012,Do_nothing,Do_nothing]
   + _r[291450206011012,Do_nothing,Recover] = 1
 R2: _r[291450206011012,R1,Do_nothing] + _r[291450206011012,R1,Recover]
   = 0
 R3: _r[291450206011012,R2,Do_nothing] + _r[291450206011012,R2,Recover]
   = 0
 R4: _r[291450206011012,R3,Do_nothing] + _r[291450206011012,R3,Recover]
   = 0
 R5: _r[290970112004069,Do_nothing,Do_nothing]
   + _r[290970112004069,Do_nothing,Recover] = 1
 R6: _r[290970112004069,R1,Do_nothing] + _r[290970112004069,R1,Recover]
   = 0
 R7: _r[290970112004069,R2,Do_nothing] + _r[290970112004069,R2,Recover]
   = 0
 R8: _r[290970112004069,R3,Do_nothing] + _r[290970112004069,R3,Recover]
   = 0
 R9: _r[290970109002051,Do_nothing,Do_nothing]
   + _r[290970109002051,Do_nothing,Recover] = 1
 R10: _r[290970109002051,R1,Do_nothing] + _r[290970109002051,R1,Recover]
   = 0
 R11: _r[290970109002051,R2,Do_nothing] + _r[290970109002051,R2,Recover]
   = 0
 R12: _r[290970109002051,R3,Do_nothing] + _r[290970109002051,R3,Recover]
   = 0
 R13: _r[290970110005024,Do_nothing,Do_nothing]
   + _r[290970110005024,Do_nothing,Recover] = 1
 R14: _r[290970110005024,R1,Do_nothing] + _r[290970110005024,R1,Recover]
   = 0
 R15: _r[290970110005024,R2,Do_nothing] + _r[290970110005024,R2,Recover]
   = 0
 R16: _r[290970110005024,R3,Do_nothing] + _r[290970110005024,R3,Recover]
   = 0
 R17: _r[290970108003014,Do_nothing,Do_nothing]
   + _r[290970108003014,Do_nothing,Recover] = 1
 R18: _r[290970108003014,R1,Do_nothing] + _r[290970108003014,R1,Recover]
   = 0
 R19: _r[290970108003014,R2,Do_nothing] + _r[290970108003014,R2,Recover]
   = 0
 R20: _r[290970108003014,R3,Do_nothing] + _r[290970108003014,R3,Recover]
   = 0
 R21: _r[290970101001022,Do_nothing,Do_nothing]
   + _r[290970101001022,Do_nothing,Recover] = 1
 R22: _r[290970101001022,R1,Do_nothing] + _r[290970101001022,R1,Recover]
   = 0
 R23: _r[290970101001022,R2,Do_nothing] + _r[290970101001022,R2,Recover]
   = 0
 R24: _r[290970101001022,R3,Do_nothing] + _r[290970101001022,R3,Recover]
   = 0
 R25: _r[290970106004001,Do_nothing,Do_nothing]
   + _r[290970106004001,Do_nothing,Recover] = 1
 R26: _r[290970106004001,R1,Do_nothing] + _r[290970106004001,R1,Recover]
   = 0
 R27: _r[290970106004001,R2,Do_nothing] + _r[290970106004001,R2,Recover]
   = 0
 R28: _r[290970106004001,R3,Do_nothing] + _r[290970106004001,R3,Recover]
   = 0
 R29: _r[290970106003025,Do_nothing,Do_nothing]
   + _r[290970106003025,Do_nothing,Recover] = 1
 R30: _r[290970106003025,R1,Do_nothing] + _r[290970106003025,R1,Recover]
   = 0
 R31: _r[290970106003025,R2,Do_nothing] + _r[290970106003025,R2,Recover]
   = 0
 R32: _r[290970106003025,R3,Do_nothing] + _r[290970106003025,R3,Recover]
   = 0
 R33: _r[290970102001050,Do_nothing,Do_nothing]
   + _r[290970102001050,Do_nothing,Recover] = 1
 R34: _r[290970102001050,R1,Do_nothing] + _r[290970102001050,R1,Recover]
   = 0
 R35: _r[290970102001050,R2,Do_nothing] + _r[290970102001050,R2,Recover]
   = 0
 R36: _r[290970102001050,R3,Do_nothing] + _r[290970102001050,R3,Recover]
   = 0
 R37: _r[290970104003071,Do_nothing,Do_nothing]
   + _r[290970104003071,Do_nothing,Recover] = 1
 R38: _r[290970104003071,R1,Do_nothing] + _r[290970104003071,R1,Recover]
   = 0
 R39: _r[290970104003071,R2,Do_nothing] + _r[290970104003071,R2,Recover]
   = 0
 R40: _r[290970104003071,R3,Do_nothing] + _r[290970104003071,R3,Recover]
   = 0
Bounds
Binaries
 _r[291450206011012,Do_nothing,Do_nothing]
 _r[291450206011012,Do_nothing,Recover] _r[291450206011012,R1,Do_nothing]
 _r[291450206011012,R1,Recover] _r[291450206011012,R2,Do_nothing]
 _r[291450206011012,R2,Recover] _r[291450206011012,R3,Do_nothing]
 _r[291450206011012,R3,Recover] _r[290970112004069,Do_nothing,Do_nothing]
 _r[290970112004069,Do_nothing,Recover] _r[290970112004069,R1,Do_nothing]
 _r[290970112004069,R1,Recover] _r[290970112004069,R2,Do_nothing]
 _r[290970112004069,R2,Recover] _r[290970112004069,R3,Do_nothing]
 _r[290970112004069,R3,Recover] _r[290970109002051,Do_nothing,Do_nothing]
 _r[290970109002051,Do_nothing,Recover] _r[290970109002051,R1,Do_nothing]
 _r[290970109002051,R1,Recover] _r[290970109002051,R2,Do_nothing]
 _r[290970109002051,R2,Recover] _r[290970109002051,R3,Do_nothing]
 _r[290970109002051,R3,Recover] _r[290970110005024,Do_nothing,Do_nothing]
 _r[290970110005024,Do_nothing,Recover] _r[290970110005024,R1,Do_nothing]
 _r[290970110005024,R1,Recover] _r[290970110005024,R2,Do_nothing]
 _r[290970110005024,R2,Recover] _r[290970110005024,R3,Do_nothing]
 _r[290970110005024,R3,Recover] _r[290970108003014,Do_nothing,Do_nothing]
 _r[290970108003014,Do_nothing,Recover] _r[290970108003014,R1,Do_nothing]
 _r[290970108003014,R1,Recover] _r[290970108003014,R2,Do_nothing]
 _r[290970108003014,R2,Recover] _r[290970108003014,R3,Do_nothing]
 _r[290970108003014,R3,Recover] _r[290970101001022,Do_nothing,Do_nothing]
 _r[290970101001022,Do_nothing,Recover] _r[290970101001022,R1,Do_nothing]
 _r[290970101001022,R1,Recover] _r[290970101001022,R2,Do_nothing]
 _r[290970101001022,R2,Recover] _r[290970101001022,R3,Do_nothing]
 _r[290970101001022,R3,Recover] _r[290970106004001,Do_nothing,Do_nothing]
 _r[290970106004001,Do_nothing,Recover] _r[290970106004001,R1,Do_nothing]
 _r[290970106004001,R1,Recover] _r[290970106004001,R2,Do_nothing]
 _r[290970106004001,R2,Recover] _r[290970106004001,R3,Do_nothing]
 _r[290970106004001,R3,Recover] _r[290970106003025,Do_nothing,Do_nothing]
 _r[290970106003025,Do_nothing,Recover] _r[290970106003025,R1,Do_nothing]
 _r[290970106003025,R1,Recover] _r[290970106003025,R2,Do_nothing]
 _r[290970106003025,R2,Recover] _r[290970106003025,R3,Do_nothing]
 _r[290970106003025,R3,Recover] _r[290970102001050,Do_nothing,Do_nothing]
 _r[290970102001050,Do_nothing,Recover] _r[290970102001050,R1,Do_nothing]
 _r[290970102001050,R1,Recover] _r[290970102001050,R2,Do_nothing]
 _r[290970102001050,R2,Recover] _r[290970102001050,R3,Do_nothing]
 _r[290970102001050,R3,Recover] _r[290970104003071,Do_nothing,Do_nothing]
 _r[290970104003071,Do_nothing,Recover] _r[290970104003071,R1,Do_nothing]
 _r[290970104003071,R1,Recover] _r[290970104003071,R2,Do_nothing]
 _r[290970104003071,R2,Recover] _r[290970104003071,R3,Do_nothing]
 _r[290970104003071,R3,Recover]
End
