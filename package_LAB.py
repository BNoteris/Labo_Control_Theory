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



# PID functions

def Proportional_action(MVP, Kc, E):
    """
    Action proportionnel du PID
    MVP
    """
    MVP.append(Kc*E[-1])
    return MVP

def Intergral_action(MVI, Kc, Ts, Ti, E, methodI):
    """
    Action intégrale du PID
    MVI
    """
    # MV[k] is MV[-1] and MV[k-1] is MV[-2]
    if len(MVI)==0:
        MVI.append(((Kc*Ts)/Ti)*(E[-1]))
    else:
        if methodI == 'EBD':
            MVI.append(MVI[-1] + ((Kc*Ts)/Ti)*(E[-1]))
        elif methodI == 'TRAP':
            MVI.append(MVI[-1] + ((Kc*Ts)/(2*Ti))*(E[-1] + E[-2]))        
        else: # EBD
            MVI.append(MVI[-1] + ((Kc*Ts)/Ti)*(E[-1]))
    return MVI

def Derivative_action(MVD,Tfd, Ts, Kc, Td, E, methodD):
    """
    Action derivative du PID
    MVD
    """
    if len(MVD)==0:
        MVD.append(((Kc*Td)/(Tfd+Ts))*(E[-1]))
    else:
        if methodD == 'EBD':
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1] + ((Kc*Td)/(Tfd+Ts))*(E[-1] - E[-2])) 
        elif methodD =='TRAP':
            MVD.append((((Tfd-(Ts/2))/(Tfd+(Ts/2)))*MVD[-1] + ((Kc*Td)/(Tfd+(Ts/2)))*(E[-1] - E[-2])))     
        else:  # EBD
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1] + ((Kc*Td)/(Tfd+Ts))*(E[-1] - E[-2]))
    return MVD

def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVmin, MVmax, MV, MVP, MVI, MVD, E, ManFF=False, PVinit=0, method="EBD-EBD"):
    """
    SP = setpoint vector
    PV = process value vector
    Man = manual controller mode vector : bool
    MVMan = Manual value for MV vector
    MVFF = feedforward vector
    Kc = controller gain
    Ti= integral time constant
    Td = derivative time constant
    alpha = Tfd = alpha*Td = where Tfd is derivative filter time constant

    MVMin = min value for MV
    MVMax = max value for MV

    MV = Maniplated value vector
    MVP = proportional part of MV vector
    MVI =  integral part of MV vector
    MVD = derivative part of MV vector
    E = control error vector

    ManFF = activated FF in manuel mode
    PVInit = initial value of PV
    Method : discreditisation value for PV
        EBD-EBD: Euler Backward difference
        TRAP-TRAP: Trapezoïdal method

    The function PID_RT appends new values to the vector MV, MVP, MVI, MVD based on PID algorithm, controller mode and FF
    """

    # Initialisation 
    methodI, methodD = method.split('-')
    Tfd = alpha*Td
    if len(PV) == 0:
        E.append(SP[-1] - PVinit)
    else: 
        E.append(SP[-1] - PV[-1])
    
    # Calcul des 3 parties
    MVP = Proportional_action(MVP, Kc, E)
    MVI = Intergral_action(MVI, Kc, Ts, Ti, E, methodI)
    MVD = Derivative_action(MVD,Tfd, Ts, Kc, Td, E, methodD)

    # Mode manuel + anti-wind up (integrator help)
    if Man[-1] == True:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] 
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]
            
    # Saturation
    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) > MVmax:
        MVI[-1] = MVmax - MVP[-1] - MVD[-1] - MVFF[-1]
    elif (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) < MVmin:
        MVI[-1] = MVmin - MVP[-1] - MVD[-1] - MVFF[-1]

    # Ajout sur MV
    MV.append(MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1])



def IMCTuning(K, Tlag1, Tlag2=0, theta=0, gamma = 0.5, process="FOPDT-PI"):
    """
    IMCTuning computes the IMC PID tuning parameters for FOPDT and SOPDT processes.
    K: process gain (Kp)
    Tlag1: first (main) lag time constant [s]
    Tlag2: second lag time constant [s]
    theta: delay [s]
    gamma: used to computed the desired closed loop time constant Tclp [s] (range [0.2 -> 0.9])
    process:
        FOPDT-PI: First Order Plus Dead Time for P-I control (IMC tuning case G)
        FOPDT-PID: First Order Plus Dead Time for P-I-D control (IMC tuning case H)
        SOPDT :Second Order Plus Dead Time for P-I-D control (IMC tuning case I)
        
    return : PID controller parameters Kc, Ti and Td
    """
    Tclp = gamma*Tlag1 
    if process=="FOPDT-PI":
        Kc = (Tlag1/(Tclp+theta))/K
        Ti = Tlag1
        Td = 0
    elif process=="FOPDT-PID":
        Kc= ((Tlag1 + theta/2)/(Tclp + theta/2))/K
        Ti = Tlag1 + theta/2
        Td = (Tlag1*theta)/(2*Tlag1+theta)
    elif process=="SOPDT": 
        Kc = ((Tlag1 + Tlag2)/(Tclp + theta))/K
        Ti = (Tlag1 +Tlag2)
        Td = ((Tlag1*Tlag2))/(Tlag1+Tlag2)
    else:
        Kc = (Tlag1/(Tclp+theta))/K
        Ti = Tlag1
        Td = 0
    return (Kc, Ti, Td)
    


class PID:
    def __init__(self, parameters):

        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 0.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 0.0

def Margin(Ps,C,omega,Show=True):
    """
    Calculate the gain margin and phase margin. They allow us to analyze the robustness of the PID.
    :Ps : Process
    :C: Controller Transfer Function
    :omega : frequency vector
    :show : allows graphical display

    """
    # Initialisation des paramètres
    s = 1j*omega
    Kc = C.parameters['Kc']
    Ti = C.parameters['Ti']
    Td = C.parameters['Td']
    Tfd = C.parameters['Tfd']
    
    # Calcul du Controller 
    Cs = Kc*(1 + 1/(Ti*s)+ (Td*s)/(Tfd*s +1))
    
    # Loop gain L(s) = P(s)C(s)
    Ls =Cs*Ps 

    # Plot de L(s)
    if Show == True:
        fig, (ax_freq, ax_time) = plt.subplots(2, 1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Amplitude
        ax_freq.semilogx(omega, 20*np.log10(np.abs(Ls)), label='L(s)')
        gain_min = np.min(20*np.log10(np.abs(Ls)/5))
        gain_max = np.max(20*np.log10(np.abs(Ls)*5))
        ax_freq.set_xlim([np.min(omega), np.max(omega)])
        ax_freq.set_ylim([gain_min, gain_max])
        ax_freq.set_ylabel('Amplitude |P| [db]')
        ax_freq.set_title('Bode plot of P')
        ax_freq.legend(loc='best')

        # Find crossover frequency where amplitude is approximately 0 dB (gain is 1)
        crossover_freq = omega[np.argmin(np.abs(20*np.log10(np.abs(Ls)) - 0))]
        ax_freq.axvline(x=crossover_freq, color='red', linestyle='--', linewidth=1)  

        # Phase
        ax_time.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)), label='L(s)')
        ax_time.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_time.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_time.set_ylabel(r'Phase $\angle P$ [°]')
        ax_time.legend(loc='best')
        ax_freq.axhline(y=0, color='black')
        ax_time.axhline(y=-180, color='black')

        # Find crossover frequency for phase
        crossover_phase = omega[np.argmin(np.abs((180/np.pi)*np.unwrap(np.angle(Ls)) + 180))]
        ax_time.axvline(x=crossover_phase, color='red', linestyle='--', linewidth=1)  

    # Crossover frequency
    i = 0
    for value in Ls:   # slide 69
        i+=1
        dB = 20*np.log10(np.abs(value))
        if dB < 0.05 and dB > -0.05:
            OmegaC =  omega[i-1]
            PhaseC = np.angle(value,deg=True)
            break        

    # Ultimate Frequency
    n = 0
    for value in Ls:
        n+=1
        deg = np.angle(value,deg=True)
        if deg < -179.5 and deg > -180.5:
            OmegaU = omega[n-1]
            u_freq = 20*np.log10(np.abs(value))
            break
    
    # Affichage graphique
    if Show ==True:
        ax_freq.plot([OmegaU, OmegaU], [0, u_freq], color='red', linewidth=5)
        ax_freq.plot([OmegaU, OmegaU], [ph_min, ph_max], linestyle='--', color='red')
        ax_time.plot([OmegaC,OmegaC],[PhaseC,-180], color='red', linewidth=5)
        ax_time.plot([OmegaC, OmegaC], [ph_min, ph_max], linestyle='--', color='red')
    print('Gain margin :',-u_freq,'dB at the ultimate frequency :',OmegaU,'rad/s')
    print('Phase margin : ',PhaseC +180,'° at the crossover frequency :',OmegaC,'rad/s')        