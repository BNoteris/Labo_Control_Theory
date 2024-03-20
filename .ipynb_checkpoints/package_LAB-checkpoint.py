import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

   
def LL_RT(MV,Kp,Ts,Tle,Tla,PV,PVInit=0,method='EBD'):
    
    """
    The function "FO_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tla: lag time constant [s]
    :Ts: sampling period [s]
    :Tle: lead time constant [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "FO_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
        

    if (Tla != 0):
        K = Ts/Tla
        if len(PV) == 0:
            PV.append(PVInit)
        else:  # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1 / (1 + K)) * PV[-1] + (K * Kp / (1 + K)) * (MV[-1] * (1 + Tle / Ts) - MV[-2] * Tle / Ts))
            elif method == 'EFD':
                PV.append((1 - K) * PV[-1] + K * Kp * (MV[-1] * Tle / Ts + (1 - Tle / Ts) * MV[-2]))
            elif method == 'TRAP':
                # Assuming T was mistakenly used instead of Tla for consistency
                PV.append((1 / (2 * Tla + Ts)) * ((2 * Tla - Ts) * PV[-1] + Kp * Ts * (MV[-1] + MV[-2])))
            else:
                # This else block seems to repeat the EBD formula; might not be necessary unless it's a fallback for an unspecified method
                PV.append((1 / (1 + K)) * PV[-1] + (K * Kp / (1 + K)) * (MV[-1] * (1 + Tle / Ts) - MV[-2] * Tle / Ts))
    else:
        PV.append(Kp * MV[-1])
