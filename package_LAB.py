import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

   
def LL_RT(MV,Kp,Ts,Tle,Tla,PV,PVInit=0,method='EBD'):
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    
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
    
    The function "LL_RT" appends a value to the output vector "PV".
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



def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD-EBD'):
    '''
The function "PID_RT" needs to be included in a "for or while loop". 

:SP: SP (or SetPoint) vector 
:PV: PV (or Process Value) vector 
:Man: Man (or Manual controller mode) vector (True or False) 
:MVMan: MVMan (or Manual value for MV) vector 
:MVFF: MVFF (or Feedforward) vector 

:Kc: controller gain 
:Ti: integral time constant [s] 
:Td: derivative time constant [s] 
:alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s] 
:Ts: sampling period [s] 

:MVMin: minimum value for MV (used for saturation and anti wind-up) 
:MVMax: maximum value for MV (used for saturation and anti wind-up) 

:MV: MV (or Manipulated Value) vector 
:MVP: MVP (or Propotional part of MV) vector 
:MVI: MVI (or Integral part of MV) vector 
:MVD: MVD (or Derivative part of MV) vector 
:E: E (or control Error) vector 

:ManFF: Activated FF in manual mode (optional: default boolean value is False) 
:PVInit: Initial value for PV (optional: default value is 0): used if PID_RT is ran first in the squence and no value of PV is available yet. 

:method: discretisation method (optional: default value is 'EBD') 
    EBD-EBD: EBD for integral action and EBD for derivative action 
    EBD-TRAP: EBD for integral action and TRAP for derivative action 
    TRAP-EBD: TRAP for integral action and EBD for derivative action 
    TRAP-TRAP: TRAP for integral action and TRAP for derivative action 

The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD". The appended values are based on the PID algorithm, the controller mode, and feedforward. Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up. 
    '''
 if len(PV) == 0: #Si pas de valeur initiale
        E.append(SP[-1] - PVInit)  # 
    else:
        E.append(SP[-1] - PV[-1]) 

    methodI, methodD = method.split('-')

    ''' Initialisation de MVI'''
    if len(MVI) == 0:
        MVI.append((Kc*Ts/Ti)*E[-1])
    else:
        if methodI == 'TRAP':
            MVI.append(MVI[-1] + (0.5*Kc*Ts/Ti)*(E[-1]+E[-2]))
        else:
            MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])

    '''Initialisation de MVD : voir slide 193  '''
    Tfd = alpha * Td
    if Td > 0:
        if len(MVD) != 0:
            if len(E) == 1:
                MVD.append((Tfd / (Tfd + Ts)) *
                           MVD[-1] + ((Kc * Td) / (Tfd + Ts)) * (E[-1]))
            else:
                MVD.append((Tfd / (Tfd + Ts)) *
                           MVD[-1] + ((Kc * Td) / (Tfd + Ts)) * (E[-1] - E[-2]))
        else:
            if len(E) == 1:
                MVD.append((Kc * Td) / (Tfd + Ts) * (E[-1]))
            else:
                MVD.append((Kc * Td) / (Tfd + Ts) * (E[-1] - E[-2]))

    '''Actualisation de MVP'''
    MVP.append(E[-1] * Kc)

    '''Activation Feedforward'''
    if ManFF:
        MVFFI = MVFF[-1]
    else:
        MVFFI = 0
    '''Mode manuel et anti-wind-up'''

    if Man[-1]:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFFI

    '''Limitation de MV'''

    MV_TEMP = MVP[-1] + MVI[-1] + MVD[-1] + MVFFI

    if MV_TEMP >= MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFFI
        MV_TEMP = MVMax

    if MV_TEMP <= MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFFI
        MV_TEMP = MVMin

    MV.append(MV_TEMP)

