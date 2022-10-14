import math
import sympy
from scipy.optimize import fsolve

from tavernmag import core, calculate_classes, sample_class

# Functions meant to be used internally by TAVERNmag

def fugCO2(KCO2, fCO, fO2):
    """
    Parameters
    ----------
    KCO2: float
        equilibrium constant for CO2 as calculated by calculate_classes
    
    fCO: float
        fugacity of CO
    
    fO2: float
        fugacity of O2
    
    Returns
    -------
    float
        Fugacity of CO2
    """
    return KCO2 * fCO * math.sqrt(fO2)

def fugSO2(KSO2, fS2, fO2):
    """
    Parameters
    ----------
    KSO2: float
        equilibrium constant for SO2 as calculated by calculate_classes
    
    fS2: float
        fugacity of S2
    
    fO2: float
        fugacity of O2
    
    Returns
    -------
    float
        Fugacity of SO2
    """
    return KSO2 * math.sqrt(fS2) * fO2

def fugH2S(KH2S, fH2, fS2):
    """
    Parameters
    ----------
    KH2S: float
        equilibrium constant for H2S as calculated by calculate_classes
    
    fH2: float
        fugacity of H2
    
    fS2: float
        fugacity of S2
    
    Returns
    -------
    float
        Fugacity of H2S
    """
    return KH2S * fH2 * math.sqrt(fS2)

def fugH2O(KH2O, fO2, fH2):
    """
    Parameters
    ----------
    KH2O: float
        equilibrium constant for H2O as calculated by calculate_classes
    fO2: float
        fugacity of O2
    fH2: float
        fugacity of H2
    
    Returns
    -------
    float
        Fugacity of H2O
    """
    return KH2O * math.sqrt(fO2) * fH2

def fugCO(KCO2, fO2, XCtot, pressure, gammaCO, gammaCO2):
    """
    Parameters
    ----------
    KCO2: float
        equilibrium constant for CO2 as calculated by calculate_classes
    fO2: float
        fugacity of O2
    XCtot: float
        Mole fraction of total C in the system.
    pressure: float
        Pressure in bars.
    gammaCO: float
        Fugacity coefficient of CO
    
    Returns
    -------
    float
        Fugacity of CO
    """
    if XCtot == 0:
        fCO = 0
    else:
        fCO = sympy.symbols('fCO') #for sympy

        equation = (((KCO2 * fCO * math.sqrt(fO2)) /
                    (sympy.Rational(3.0) * gammaCO2 * pressure)) +
                    ((fCO)/(sympy.Rational(2.0) * gammaCO * pressure))  -
                    XCtot)
        fCO = sympy.solve(equation, fCO)[0] #newly implemented sympy way
    
    return fCO


def fugs_H2_S2(opt, gammas, K_vals, pressure, XHtot, XStot, fO2):
    """
    Uses one of three user-chosen optimization functions to solve simultaneously for fH2 and fS2.

    Parameters
    ----------
    opt: str
        Optional. Default value is 'gekko' in which case the GEKKO library will be used for
        optimization of fS2 and fH2. Other options are 'leastsq' (scipy.optimize.leastsq) and
        'least_squares' (scipy.optimize.least_squares). If the GEKKO method fails to converge on
        a solution or hits the maximum iteration number of 250, an internal exception will be
        thrown and the 'leastsq' method will be attempted instead.
    gammas: dict
        Optional. Fugacity coefficients.
        If gamma values are not passed, they will be calculated.
        Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    K_vals: dict
        Optional. Equilibrium constants of reaction.
        If K values are not passed, they will be calculated.
        Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    pressure:   float
        Pressure of the system in bars.
    XHtot:  float
        Mole fraction of total H in the system.
    XStot:  float
        Mole fraction of total S in the system.
    fO2:    float
        Oxygen fugacity.

    Returns
    -------
    floats
        Fugacities of H2 and S2 as fH2, fS2
    """
    def optimize_leastsq(gammas, K_vals, pressure, XHtot, XStot, fO2):
        if XStot == 0:
            fS2 = 0
        if XHtot == 0:
            fH2 = 0
        def equations(p):
            fH2, fS2 = p
            return  [
                        ((fH2/(gammas['H2'] * pressure)) +
                            ((sympy.Rational(2.0) * K_vals['H2O'] * fH2 * math.sqrt(fO2)) /
                                (sympy.Rational(3.0) * gammas['H2O'] * pressure))   +
                            ((sympy.Rational(2.0) * K_vals['H2S'] * fH2 * math.sqrt(abs(fS2))) /
                                (sympy.Rational(3.0) * gammas['H2S'] * pressure)) -
                            XHtot), 
                        ((fS2/(gammas['S2'] * pressure)) +
                            ((K_vals['H2S'] * fH2 * math.sqrt(abs(fS2))) /
                                (sympy.Rational(3.0) * gammas['H2S'] * pressure)) + 
                            ((K_vals['SO2'] * fO2 * math.sqrt(abs(fS2))) /
                                (sympy.Rational(3.0) * gammas['SO2'] * pressure)) -
                            XStot)
                    ]

        # fH2_a, fS2_a = fsolve(equations, (1, 0))
        # leastsq outperforms fsolve, particularly at low fO2 conditions
        # where H2 and S2 (co-solved for above) are abundant
        from scipy.optimize import leastsq 
        fH2_a, fS2_a = leastsq(equations, (0,0))[0]

        if XHtot == 0:
            fH2 = 0
        else:
            fH2 = abs(fH2_a)

        if XStot  == 0:
            fS2 = 0
        else:
            fS2 = abs(fS2_a)
        
        return fH2, fS2

    def optimize_gekko(gammas, K_vals, pressure, XHtot, XStot, fO2):
        # gekko non-linear eq solver
        from gekko import GEKKO

        m = GEKKO() # create GEKKO model
        fH2 = m.Var(value=0) # define new variable, initial value=0
        fS2 = m.Var(value=0) # define new variable, initial value=0
        fH2.lower = 0 # set lower bound to 0
        fS2.lower = 0 # set lower bound to 0
        m.Equations([
                            ((fH2/(gammas['H2'] * pressure)) +
                                ((sympy.Rational(2.0) * K_vals['H2O'] * fH2 * m.sqrt(fO2)) /
                                    (sympy.Rational(3.0) * gammas['H2O'] * pressure))   +
                                ((sympy.Rational(2.0) * K_vals['H2S'] * fH2 * m.sqrt(fS2)) /
                                    (sympy.Rational(3.0) * gammas['H2S'] * pressure)) -
                                XHtot==0), 
                            ((fS2/(gammas['S2'] * pressure)) +
                                ((K_vals['H2S'] * fH2 * m.sqrt(fS2)) /
                                    (sympy.Rational(3.0) * gammas['H2S'] * pressure)) + 
                                ((K_vals['SO2'] * fO2 * m.sqrt(fS2)) /
                                    (sympy.Rational(3.0) * gammas['SO2'] * pressure)) -
                                XStot==0)
                        ])
        m.solve(disp=False)
        
        return fH2.value[0], fS2.value[0]

    def optimize_least_squares(gammas, K_vals, pressure, XHtot, XStot, fO2):
        # scipy optimize least_squares
        if XStot == 0:
            fS2 = 0
        if XHtot == 0:
            fH2 = 0
        def equations(p):
            vals = np.zeros(p.size)
            vals[0] = ((p[0]/(gammas['H2'] * pressure)) +
                            ((sympy.Rational(2.0) * K_vals['H2O'] * p[0] * math.sqrt(fO2)) /
                                (sympy.Rational(3.0) * gammas['H2O'] * pressure))   +
                            ((sympy.Rational(2.0) * K_vals['H2S'] * p[0] * math.sqrt(p[1])) /
                                (sympy.Rational(3.0) * gammas['H2S'] * pressure)) -
                            XHtot)
            vals[1] = ((p[1]/(gammas['S2'] * pressure)) +
                            ((K_vals['H2S'] * p[0] * math.sqrt(p[1])) /
                                (sympy.Rational(3.0) * gammas['H2S'] * pressure)) + 
                            ((K_vals['SO2'] * fO2 * math.sqrt(p[1])) /
                                (sympy.Rational(3.0) * gammas['SO2'] * pressure)) -
                            XStot)
            return vals

        # fH2_a, fS2_a = fsolve(equations, (1, 0))
        # leastsq outperforms fsolve, particularly at low fO2 conditions
        # where H2 and S2 (co-solved for above) are abundant
        from scipy.optimize import least_squares 
        fH2_a, fS2_a = least_squares(equations, (0,0), bounds=([0,0],np.inf)).x

        if XHtot == 0:
            fH2 = 0
        else:
            fH2 = fH2_a

        if XStot  == 0:
            fS2 = 0
        else:
            fS2 = fS2_a
        
        return fH2, fS2
    
    if opt == 'leastsq':
        fH2, fS2 = optimize_leastsq(gammas=gammas, K_vals=K_vals,
                                            pressure=pressure,
                                            XHtot=XHtot, XStot=XStot, fO2=fO2)
    
    if opt == 'gekko':
        try:
            fH2, fS2 = optimize_gekko(gammas=gammas, K_vals=K_vals,
                                            pressure=pressure, XHtot=XHtot, XStot=XStot, fO2=fO2)
        except:
            fH2, fS2 = optimize_leastsq(gammas=gammas, K_vals=K_vals,
                                            pressure=pressure, XHtot=XHtot, XStot=XStot,
                                            fO2=fO2)
    
    if opt == 'least_squares':
        fH2, fS2 = optimize_least_squares(gammas=gammas, K_vals=K_vals,
                                                pressure=pressure, XHtot=XHtot, XStot=XStot,
                                                fO2=fO2)

    if XHtot == 0:
        fH2 = 0
    else:
        fH2 = abs(fH2)

    if XStot  == 0:
        fS2 = 0
    else:
        fS2 = abs(fS2)
    
    return fH2, fS2