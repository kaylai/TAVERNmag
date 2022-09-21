from abc import abstractmethod
import math
import sympy
import numpy as np
from scipy.optimize import fsolve
import warnings as w

from tavern import core, fO2_buffers, sample_class, sums

from copy import deepcopy, copy

class Calculate(object):
    """ The Calculate object is a template for implementing user-friendly methods for running
    calculations on MagmaticFluid or SilicateMelt objects.
    """
    def __init__(self, **kwargs):
        """
        Initializes the calculation.

        Parameters
        ----------
        sample: MagmaticFluid or SilicateMelt class
            The composition of a magmatic fluid or silicate melt as a tavern class.
        """
        self.result = self.calculate(**kwargs)

    @abstractmethod
    def calculate(self):
        """ """

class calculate_equilibrium_constants(Calculate):
    """Returns log of the equilibrium constants from lookup table.

    Parameters
    ----------
    temperature: float
        Temperature in degrees C.

    species: str
        Choose which species to calculate.
        Options are: CO2, H2O, H2S, SO2, or all.
        If all is passed, a dictionary of values is returned.

    return_as: str
        Default value is 'standard', which returns equilibrium constant K.
        Optional value is 'log', which returns the logK.

    Returns
    -------
    Calculate object
        Calculate object. Access results by fetching the result property. Equilibrium constant for
        chosen species. If one species is passed, returns float. If "all" is passed, returns dict
        with keys CO2, H2O, H2S,
    """
    def calculate(self, temperature, species='all', return_as='standard', **kwargs):
        tempK = temperature + 273.15

        # CO + 1/2O2 = CO2
        """From Wagman et al (1945). Text table was fit with sixth-order polynomial to
        arrive at the following coefficients."""
        CO2_logK =       (9.11899*10.0**(-17.0) *   tempK**6.0
                    +   -5.65762*10.0**(-13.0)  *   tempK**5.0
                    +    1.44706*10.0**(-9.0)   *   tempK**4.0
                    +   -1.96925*10.0**(-6.0)   *   tempK**3.0
                    +    1.53277*10.0**(-3.0)   *   tempK**2.0
                    +   -6.79138*10.0**(-1.0)   *   tempK
                    +    1.53288*10.0**(2.0))

        # H2 + 1/2O2 = H2O
        """From Robie and Hemingway (1995). Text table was fit with sixth-order polynomial to
        arrive at the following coefficients."""
        H2O_logK =      (3.3426*10.0**(-17.0)   *   tempK**6.0
                    +   -2.40813*10.0**(-13.0)  *   tempK**5.0
                    +   7.10612*10.0**(-10.0)   *   tempK**4.0
                    +   -1.10671*10.0**(-6.0)   *   tempK**3.0
                    +   9.76038*10.0**(-4.0)    *   tempK**2.0
                    +   -4.84445*10.0**(-1.0)   *   tempK
                    +   1.21953*10.0**(2.0))

        # H2 + 1/2S2 = H2S
        """From Robie and Hemingway (1995). Text table was fit with sixth-order polynomial to
        arrive at the following coefficients."""
        H2S_logK =      (1.882737*10.0**(-18.0) * tempK**6.0
                    +   -1.779266*10.0**(-14.0) * tempK**5.0
                    +    6.319209*10.0**(-11.0) * tempK**4.0
                    +   -1.092048*10.0**(-7.0) * tempK**3.0
                    +    9.824774*10.0**(-5.0) * tempK**2.0
                    +   -4.805344*10.0**(-2.0) * tempK
                    +    1.389857*10.0**(1.0))

        # 1/2S2 + O2 = SO2
        """From Robie and Hemingway (1995). Text table was fit with sixth-order polynomial to
        arrive at the following coefficients."""
        SO2_logK =       (4.01465*10.0**(-17.0) * tempK**6.0
                    +   -2.93845*10.0**(-13.0)  * tempK**5.0
                    +    8.78818*10.0**(-10.0)  * tempK**4.0
                    +   -1.38066*10.0**(-6.0)   * tempK**3.0
                    +    1.21978*10.0**(-3.0)   * tempK**2.0
                    +   -6.03476*10.0**(-1.0)   * tempK
                    +    1.54350*10.0**(2.0))

        if return_as == 'standard':
            if species == 'CO2':
                return 10.0**CO2_logK
            if species == 'H2O':
                return 10.0**H2O_logK
            if species == 'H2S':
                return 10.0**H2S_logK
            if species == 'SO2':
                return 10.0**SO2_logK
            if species == 'all':
                return {'CO2': 10.0**CO2_logK, 'H2O': 10.0**H2O_logK, 'H2S': 10.0**H2S_logK,
                        'SO2': 10.0**SO2_logK}
        if return_as == 'log':
            if species == 'CO2':
                return CO2_logK
            if species == 'H2O':
                return H2O_logK
            if species == 'H2S':
                return H2S_logK
            if species == 'SO2':
                return SO2_logK
            if species == 'all':
                return {'CO2': CO2_logK, 'H2O': H2O_logK, 'H2S': H2S_logK, 'SO2': SO2_logK}
        #TODO raise exception if something other than 'standard' or 'log' is passed?


class calculate_fugacity_coefficients(Calculate):
    """Returns fugacity coefficients calculated using the Redlich Kwong Equation of State.
    Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 
    30 October 2003.

    Parameters
    ----------
    temperature: float
        Temperature in degrees C.

    pressure: float
        Pressure in bars.

    species: str
        Choose which species to calculate. Options are: 'CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2',
        'S2', 'SO2', or all. If all is passed, a dictionary of values is returned. Default value
        is 'all'.

    Returns
    -------
    Calculate object
        Calculate object. Access results by fetching the result property. Fugacity coefficient for
        passed species. If single species is passed, float. If "all" is passed, dictionary with keys
        'CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2', 'S2', 'SO2'.
    """
    def calculate(self, temperature, pressure, species='all', **kwargs):
        tempK = temperature + 273.15
        R = 8.3145

        gamma_dict = {}

        for species in core.fluid_species_names:
            #Calculate a and b parameters (depend only on critical parameters)...
            a = (0.42748 * R**2.0 * core.critical_params[species]["cT"]**(2.5) / 
                 (core.critical_params[species]["cP"] * 10.0**5))
            b = (0.08664 * R * core.critical_params[species]["cT"] / 
                 (core.critical_params[species]["cP"] * 10.0**5))
            kappa = 0.0

            #Calculate coefficients in the cubic equation of state...
            #coeffs: (C0, C1, C2, A, B)
            A = a * pressure * 10.0**5 / (math.sqrt(tempK) * (R * tempK)**2.0)
            B = b * pressure * 10.0**5 / (R * tempK)
            C2 = -1.0
            C1 = A - B - B * B
            C0 = -A * B

            #Solve the cubic equation for Z0 - Z2, D...
            Q1 = C2 * C1 / 6.0 - C0 / 2.0 - C2**3.0 / 27.0
            P1 = C2**2.0 / 9.0 - C1 / 3.0
            D = Q1**2.0 - P1**3.0

            if D >= 0:
                kOneThird = 1.0 / 3.0

                absQ1PsqrtD = math.fabs(Q1 + math.sqrt(D))
                temp1 = absQ1PsqrtD**kOneThird
                temp1 *= (Q1 + math.sqrt(D)) / absQ1PsqrtD

                absQ1MsqrtD = math.fabs(Q1 - math.sqrt(D))
                temp2 = absQ1MsqrtD**kOneThird
                temp2 *= (Q1 - math.sqrt(D)) / absQ1MsqrtD

                Z0 = temp1 + temp2 - C2 / 3.0
            else:
                temp1 = Q1**2.0 / (P1**3.0)
                temp2 = math.sqrt(1.0 - temp1) / math.sqrt(temp1)
                temp2 *= Q1 / math.fabs(Q1)

                gamma = math.atan(temp2)

                if gamma < 0:
                    gamma = gamma + math.pi 

                Z0 = 2.0 * math.sqrt(P1) * math.cos(gamma/3.0) - C2 / 3.0
                Z1 = 2.0 * math.sqrt(P1) * math.cos((gamma + 2.0 * math.pi) / 3.0) - C2/3.0
                Z2 = 2.0 * math.sqrt(P1) * math.cos((gamma + 4.0 * math.pi) / 3.0) - C2/3.0

                if Z0 < Z1:
                    temp0 = Z0
                    Z0 = Z1
                    Z1 = temp0

                if Z1 < Z2:
                    temp0 = Z1
                    Z1 = Z2
                    Z2 = temp0

                if Z0 < Z1:
                    temp0 = Z0
                    Z0 = Z1
                    Z1 = temp0

            #Determine the fugacity coefficient of first root and departure functions...
            #calcdepfns(coeffs[3],  coeffs[4],  paramsab[0],    Z[0])
            #calcdepfns(A,          B,          kappa,          Z)

            #Calculate Departure Functions
            gamma = math.exp(Z0 - 1.0 - math.log(Z0-B) - A * math.log(1.0+B/Z0)/B)
            Hdep = R * tempK * (Z0 - 1.0 - 1.5*A*math.log(1.0+B/Z0)/B)
            Sdep = R * (math.log(Z0-B) - 0.5*A*math.log(1.0+B/Z0)/B)

            gamma_dict[species] = gamma

            #gamma_tuple = tuple(gamma_dict.values())

        if species == 'CO':
            return gamma_dict["CO"]
        if species == 'CO2':
            return gamma_dict["CO2"]
        if species == 'H2':
            return gamma_dict["H2"]
        if species == 'H2O':
            return gamma_dict["H2O"]
        if species == 'H2S':
            return gamma_dict["H2S"]
        if species == 'O2':
            return gamma_dict["O2"]
        if species == 'S2':
            return gamma_dict["S2"]
        if species == 'SO2':
            return gamma_dict
        if species == 'all':
            return gamma_dict


class calculate_fugacities(Calculate):
    # TODO check that fugacities sum to total P?
    """Returns fugacity values for each species.

    Parameters
    ----------
    sample: MagmaticFluid class
        Composition of the magmatic fluid as MagmaticFluid object. Can pass simple (H2O, CO2, S) or
        complex (any species covering H, C, and S budgets) fluid.
    pressure:   float
        Pressure of the system in bars.
    temperature:    float
        Temperature of the system in degrees C.
    fO2_buffer: str
        Name of buffer for which to calculate fO2. Can be one of QFM.
    fO2_delta: float
        Deviation from named buffer in log units. For example, for QFM+2, enter 'QFM' as the
        fO2_buffer and 2 as the fO2_delta.
    gammas: dict
        Optional. Fugacity coefficients.
        If gamma values are not passed, they will be calculated.
        Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    K_vals: dict
        Optional. Equilibrium constants of reaction.
        If K values are not passed, they will be calculated.
        Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    opt: str
        Optional. Default value is 'gekko' in which case the GEKKO library will be used for
        optimization of fS2 and fH2. Other options are 'leastsq' (scipy.optimize.leastsq) and
        'least_squares' (scipy.optimize.least_squares).

    Returns
    -------
    Calculate object
        Calculate object. Access results by fetching the result property. Fugacities of all species.

    """
    def calculate(self, sample, pressure, temperature, fO2_buffer, fO2_delta, gammas='calculate',
                  K_vals='calculate', opt='gekko', **kwargs):
        # Get fO2 value, gammas, and K_values
        logfO2 = fO2_buffers.calc_logfO2_from_buffer(pressure=pressure, temperature=temperature,
                                                     buffer=fO2_buffer, delta=fO2_delta)
        fO2 = 10.0**logfO2

        if gammas == 'calculate':
            gammas = calculate_fugacity_coefficients(temperature=temperature, pressure=pressure,
                                                     species="all").result

        if K_vals == 'calculate':
            K_vals = calculate_equilibrium_constants(temperature=temperature, species="all").result           

        # check if user is passing a simplified composition
        required_species = {"H2O": "H2O", "CO2": "CO2", "S": "S"}
        if set(sample.get_composition().keys()) <= set(required_species.keys()):
            # user has passed a simple composition
            # composition_molfrac = sample.get_composition(units='molfrac', normalization='standard') 

            # XH2Otot = composition_molfrac["H2O"]
            # XCO2tot = composition_molfrac["CO2"]
            # XStot = composition_molfrac["S"]

            # XHtot_unnorm = XH2Otot * 2
            # XStot_unnorm = XStot
            # XCtot_unnorm = XCO2tot

            # # normalize H, S, C tots
            # HCS_tots = sample_class.MagmaticFluid({"H": XHtot_unnorm,
            #                                        "S": XStot_unnorm,
            #                                        "C": XCtot_unnorm},
            #                                        units='molfrac',
            #                                        normalization='standard')

            HCS_tots = sample.get_simplified_fluid_composition(H_species="H",
                                                               C_species="C",
                                                               S_species="S",
                                                               units='molfrac',
                                                               asSampleClass=True)
            XHtot = HCS_tots.get_composition(species="H", units='molfrac')
            XStot = HCS_tots.get_composition(species="S", units='molfrac')
            XCtot = HCS_tots.get_composition(species="C", units='molfrac')

        else: #user passed complex composition
            composition_molfrac = sample.get_simplified_fluid_composition(units='molfrac',
                                                                          H_species='H',
                                                                          C_species='C',
                                                                          S_species='S')

            XHtot = composition_molfrac["H"]
            XStot = composition_molfrac["S"]
            XCtot = composition_molfrac["C"]

        #FIRST calculate fH2 and fS2 using fsolve, two eqns; two unknowns (eqn 9 in Iacovino, 2015)
        # scipy optimize process
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

        # gekko non-linear eq solver
        from gekko import GEKKO
        def optimize_gekko(gammas, K_vals, pressure, XHtot, XStot, fO2):
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

        # scipy optimize least_squares
        def optimize_least_squares(gammas, K_vals, pressure, XHtot, XStot, fO2):
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

        # choose optimization routine based on user input
        if opt == 'leastsq':
            fH2, fS2 = optimize_leastsq(gammas, K_vals, pressure, XHtot, XStot, fO2)
        
        if opt == 'gekko':
            fH2, fS2 = optimize_gekko(gammas, K_vals, pressure, XHtot, XStot, fO2)
        
        if opt == 'least_squares':
            fH2, fS2 = optimize_least_squares(gammas, K_vals, pressure, XHtot, XStot, fO2)
        # end optimization options

        if XHtot == 0:
            fH2 = 0
        else:
            fH2 = abs(fH2)

        if XStot  == 0:
            fS2 = 0
        else:
            fS2 = abs(fS2)

        # SECOND calculate fCO (eqn 10 in Iacovino, 2015) using sympy
        if XCtot  == 0:
            fCO = 0
        else:
            fCO = sympy.symbols('fCO') #for sympy

            equation = (((K_vals['CO2'] * fCO * math.sqrt(fO2)) /
                        (sympy.Rational(3.0) * gammas['CO2'] * pressure)) +
                        ((fCO)/(sympy.Rational(2.0) * gammas['CO'] * pressure))  -
                        XCtot)
            fCO = sympy.solve(equation, fCO)[0] #newly implemented sympy way

        # THIRD calculate fCO2 using calc'd fCO and known fO2 value
        fCO2 = K_vals['CO2'] * fCO * math.sqrt(fO2)

        # FOURTH calcualte fSO2 using calc'd fS2 and known fO2 value
        fSO2 = K_vals['SO2'] * math.sqrt(fS2) * fO2

        # FIFTH calculate fH2S using calc'd fH2 and fS2 values
        fH2S = K_vals['H2S'] * fH2 * math.sqrt(fS2)

        # SIXTH calculate fH2O using calc'd fH2 and knwn fO2 value
        fH2O = K_vals['H2O'] * math.sqrt(fO2) * fH2

        # TODO raise exception if a fugacity is negative or zero.

        return_dict = {'CO': fCO, 'CO2': fCO2, 'H2': fH2, 'H2O': fH2O, 'H2S': fH2S, 'O2': fO2,
                       'S2': fS2, 'SO2': fSO2}

        # perform checks/debugging
        f_error = False

        # check fH2O/fH2 ratio
        if XHtot == 0:
            fH_ratio_constraint = np.nan 
        else:
            fH_ratio_constraint = float(K_vals['H2O'] * math.sqrt(fO2))
            fH_ratio_calculated = float(fH2O/fH2)
        if XHtot != 0:
            if round(fH_ratio_constraint,5) != round(fH_ratio_calculated,5):
                f_error = True
                print("fH ratio error")
                print("fH2O/fH2 should be: " + str(fH_ratio_constraint))
                print("fH2O/fH2 ratio is: " + str(fH_ratio_calculated))

        # check fCO2/fCO ratio
        if XCtot == 0:
            fC_ratio_constraint = np.nan 
        else:
            fC_ratio_constraint = float(K_vals['CO2'] * math.sqrt(fO2))
            fC_ratio_calculated = float(fCO2/fCO)
        if XCtot != 0:
            if round(fC_ratio_constraint,5) != round(fC_ratio_calculated,5):
                f_error = True
                print("fC ratio error")
                print("fCO2/fCO should be: " + str(fC_ratio_constraint))
                print("fCO2/fCO ratio is: " + str(fC_ratio_calculated))

        # check fSO2/math.sqrt(fS2) ratio
        if XStot == 0:
            fS_ratio_constraint = np.nan 
        else:
            fS_ratio_constraint = float(K_vals['SO2'] * fO2)
            fS_ratio_calculated = float(fSO2/math.sqrt(fS2))
        if XStot != 0:
            if round(fS_ratio_constraint,5) != round(fS_ratio_calculated,5):
                f_error = True
                print("fS ratio error")
                print("fSO2/math.sqrt(fS2) should be: " + str(fS_ratio_constraint))
                print("fSO2/math.sqrt(fS2) ratio is: " + str(fS_ratio_calculated))

        # check fH2S/math.sqrt(fS2) ratio
        if XStot == 0:
            fHS_ratio_constraint = np.nan 
        else:
            fHS_ratio_constraint = float(K_vals['H2S'] * fH2)
            fHS_ratio_calculated = float(fH2S/math.sqrt(fS2))
        if XStot != 0:
            if round(fHS_ratio_constraint,5) != round(fHS_ratio_calculated,5):
                f_error = True
                print("fHS ratio error")
                print("fH2S/math.sqrt(fS2) should be: " + str(fHS_ratio_constraint))
                print("fH2S/math.sqrt(fS2) ratio is: " + str(fHS_ratio_calculated))

        if f_error is True:
            # print out all fugacities to check they are sensible
            for k, v in return_dict.items():
                print("f" + str(k) + ": " + str(v))

        return return_dict


class calculate_speciation(Calculate):
    """Speciates a fluid given bulk composition.

    Parameters
    ----------
    sample: MagmaticFluid class
        Composition of the magmatic fluid as MagmaticFluid object.
    pressure: float
        Pressure of the system in bars.
    temperature: float
        Temperature of the system in degrees C.
    fO2_buffer: str
        Name of buffer for which to calculate fO2. Can be one of QFM.
    fO2_delta: float
        Deviation from named buffer in log units. For example, for QFM+2, enter 'QFM' as the
        fO2_buffer and 2 as the fO2_delta.
    gammas: dict
        Optional. Fugacity coefficients. If gamma values are not passed, they will be
        calculated. Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    K_vals: dict
        Optional. Equilibrium constants of reaction. If K values are not passed, they will be 
        calculated. Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    fugacities: dict
        Optional. Fugacity values for each species. If fugacity values are not passed, they will
        be calculated. Format: {'CO': value, 'CO2': value, 'H2': value, 'H2O': value,
        'H2S': value, 'O2': value, 'S2': value, 'SO2': value}

    Returns
    -------
    Calculate object
        Calculate object. Access results by fetching the result property. Speciated fluid
        composition in default units set for MagmaticFluid object.
    """
    def normalize_fluid_FixedOxygen(self, composition, units='wtpercent', **kwargs):
        """
        Normalizes a fluid composition to 100%, including O2. The O2 wt%
        will remain fixed, whilst the other species are reduced proprotionally
        so that the total is 100 wt%

        Parameters
        ----------
        composition: dict
            fluid composition

        units:  str
            The units of composition. Should be one of:
            - wtpercent (default)
            - molpercent
            - molfrac

        Returns
        -------
        dict
            Normalized fluid composition.
        """
        normalized = {}
        O2_sum = composition['O2']

        for spec in composition.keys():
            if spec != 'O2':
                normalized[spec] = composition[spec]
        normsum = sum(normalized.values())

        if units == 'wtpercent' or units == 'molpercent':
            normalized = {k: v/normsum*(100.-O2_sum) for k, v in normalized.items()}
        elif units == 'molfrac':
            normalized = {k: v/normsum*(1.-O2_sum) for k, v in normalized.items()}
        else:
            raise Exception("Units must be one of 'wtpercent', "
                            "'molpercent', or 'molfrac'.")

        normalized['O2'] = composition['O2']

        return normalized

    def return_default_units(self, sample, calc_result, units, **kwargs):
        """ Checkes the default units set for the sample_class.MagmaticFluid
        object and returns the result of a calculation in those units.

        Parameters
        ----------
        sample: MagmaticFluid class
            The fluid composition as a MagmaticFluid object.

        calc_result: dict
            Result of a calculate_speciation() calculation on a sample.

        units:  str
            Units of calc_result dict. Should be one of:
            - wtpercent (default)
            - molpercent
            - molfrac
        """
        default_units = sample.default_units

        bulk_comp = deepcopy(calc_result)
        bulk_comp = sample_class.MagmaticFluid(bulk_comp, units=units)

        return dict(bulk_comp.get_composition(units=default_units))


    def calculate(self, sample, pressure, temperature, fO2_buffer, fO2_delta, gammas='calculate',
                  K_vals='calculate', fugacities='calculate', **kwargs):
        if gammas == 'calculate':
            gammas = calculate_fugacity_coefficients(temperature=temperature, pressure=pressure,
                                                     species="all").result

        if K_vals == 'calculate':
            K_vals = calculate_equilibrium_constants(temperature=temperature, species="all").result

        if fugacities == 'calculate':
            fugacities = calculate_fugacities(sample=sample, pressure=pressure,
                                              temperature=temperature, fO2_buffer=fO2_buffer,
                                              fO2_delta=fO2_delta, gammas=gammas,
                                              K_vals=K_vals, **kwargs).result

        X_dict = {}
        for species in core.fluid_species_names:
            X = fugacities[species] / (gammas[species] * pressure)
            X_dict[species] = X

        # Normalize and fix O2 since fO2 is input and should not be adjusted during norm
        norm = self.normalize_fluid_FixedOxygen(X_dict, units='molfrac')

        speciated_default_units = self.return_default_units(sample, norm, units='molfrac', **kwargs)

        return_sample = sample_class.MagmaticFluid(speciated_default_units,
                                                   units=sample.default_units,
                                                   default_units=sample.default_units)

        return return_sample


class respeciate_fluid(Calculate):
    """Takes in a speciated fluid of given composition in terms of CO, CO2, H2, H2O, H2S, O2, S2,
    and SO2 and returns that fluid respeciated at given conditions.

    Parameters
    ----------
    sample: MagmaticFluid class
        Composition of the magmatic fluid as MagmaticFluid object.

    pressure: float
        Pressure at which to respeciate the fluid, in bars.

    temperature: float
        Temperature at which to respeciate the fluid, in degrees C.

    fO2_buffer: str
        Name of new buffer for which to calculate fO2. Can be one of 'QFM' or 'NNO'.

    fO2_delta: float
        Deviation from named buffer in log units. For example, for QFM+2, enter 'QFM' as the
        fO2_buffer and 2 as the fO2_delta.

    Returns
    -------
    MagmaticFluid object
        Fluid composition after respeciation.
    """

    def calculate(self, sample, pressure, temperature, fO2_buffer, fO2_delta):
        # Calculate new fugacity coefficients and equilibrium constants at given conditions.
        gammas = calculate_fugacity_coefficients(pressure=pressure, temperature=temperature).result
        Ks = calculate_equilibrium_constants(temperature=temperature).result

        #Recalculate fO2 at new pressure
        logfO2_res = fO2_buffers.calc_logfO2_from_buffer(temperature=temperature, buffer=fO2_buffer,
                                             delta=fO2_delta, pressure=pressure)
        fO2_res = 10**(logfO2_res)

        # get composition in mole fraction
        fluid_comp = sample.get_composition(units='molfrac')

        # Name some local variables for clarity in below equations. If a particular species is not
        # passed, it will be assigned a value of 0.0 
        if 'CO' in fluid_comp:
            CO = fluid_comp["CO"]
        else:
            CO = 0.0
        if 'CO2' in fluid_comp:
            CO2 = fluid_comp["CO2"]
        else:
            CO2 = 0.0
        if 'H2' in fluid_comp:
            H2 = fluid_comp["H2"]
        else:
            H2 = 0.0
        if 'H2O' in fluid_comp:
            H2O = fluid_comp["H2O"]
        else:
            H2O = 0.0
        if 'H2S' in fluid_comp:
            H2S = fluid_comp["H2S"]
        else:
            H2S = 0.0
        if 'O2' in fluid_comp:
            O2 = fluid_comp["O2"]
        else:
            O2 = 0.0
        if 'S2' in fluid_comp:
            S2 = fluid_comp["S2"]
        else:
            S2 = 0.0
        if 'SO2' in fluid_comp:
            SO2 = fluid_comp["SO2"]
        else:
            SO2 = 0.0

        XHtot = H2 + (0.666666)*H2O + (0.666666)*H2S
        XStot = S2 + (0.333333)*H2S + (0.333333)*SO2
        XCtot = (0.333333)*CO2 + (0.5)*CO
        XOtot = O2 + (0.666666)*CO2 + (0.5)*CO + (0.333333)*H2O + (0.666666)*SO2

        #FIRST calculate fH2 and fS2 using fsolve, two eqns; two unknowns (eqn 9 in Iacovino, 2015)

        #TODO - figure out how to do these two equations with two unknowns and set bounds that roots must be >0
        #Used least_squares to do this in a similar implimentation in this very script, but I can't
        #figureo out how to get it to work with two equations instead of just one.
        #THIS IS IMPORTANT since right now it's a total non-mathematical flub.

        def equations(p):
            fH2, fS2 = p
            return  (
                     ((fH2/(gammas["H2"]*pressure))    +
                        ((sympy.Rational(2.0) * Ks["H2O"] * fH2 * sympy.sqrt(fO2_res))/
                            (sympy.Rational(3.0) * gammas["H2O"] * pressure)) +
                        ((sympy.Rational(2.0) * Ks["H2S"] * fH2 * sympy.sqrt(abs(fS2)))/
                            (sympy.Rational(3.0) * gammas["H2S"] * pressure)) -
                        XHtot), 
                     ((fS2/(gammas["S2"] * pressure))  +
                        ((Ks["H2S"] * fH2 * sympy.sqrt(abs(fS2)))/
                            (sympy.Rational(3.0) * gammas["H2S"] * pressure)) +
                        ((Ks["SO2"] * fO2_res * sympy.sqrt(abs(fS2)))/
                            (sympy.Rational(3.0) * gammas["SO2"] * pressure)) -
                        XStot)
                    )

        fH2_a, fS2_a = fsolve(equations, (1, 1))
        fH2 = abs(fH2_a)
        fS2 = abs(fS2_a)

        #SECOND calculate fCO (eqn 10 in Iacovino, 2015)
        fCO = sympy.symbols('fCO') #for sympy
        equation = (((Ks["CO2"] * fCO * sympy.sqrt(fO2_res))/(3.0 * gammas["CO2"] * pressure)) +
                    ((fCO)/(2.0 * gammas["CO"] * pressure)) -
                    XCtot)
        fCO = sympy.solve(equation, fCO)[0] #newly implemented sympy way

        #THIRD calculate fCO2 using calc'd fCO and known fO2 value
        fCO2 = Ks["CO2"] * fCO * math.sqrt(fO2_res)

        #FOURTH calcualte fSO2 using calc'd fS2 and known fO2 value
        fSO2 = Ks["SO2"] * math.sqrt(fS2) * fO2_res

        #FIFTH calculate fH2S using calc'd fH2 and fS2 values
        fH2S = Ks["H2S"] * fH2 * math.sqrt(fS2)

        #SIXTH calculate fH2O using calc'd fH2 and knwn fO2 value
        fH2O = Ks["H2O"] * math.sqrt(fO2_res) * fH2 

        new_fugacities = {"CO": fCO,
                            "CO2": fCO2,
                            "H2": fH2,
                            "H2O": fH2O,
                            "H2S": fH2S,
                            "O2": fO2_res,
                            "S2": fS2,
                            "SO2": fSO2}

        X_dict = {}
        for species in core.fluid_species_names:
            X = new_fugacities[species] / (gammas[species] * pressure)
            X_dict[species] = X

        X_dict = {key: value/sum(X_dict.values()) for key,value in X_dict.items()}

        return_sample = sample_class.MagmaticFluid(X_dict, units="molfrac")
        
        return return_sample


class calculate_fH2O_from_melt(Calculate):
    """ Calculates the H2O fugacity using the model of Moore (1998) based on major element silicate
    melt composition, including H2O concentration. NOTE: Only valid up to 3,000 bars.

    Parameters
    ----------
    sample: SilicateMelt class
        The silicate melt composition, including H2O, as a SilicateMelt object.
    pressure: float
        Pressure of the system in bars.
    temperature: float
        Temperature of the system in degrees C.

    Returns
    -------
    Calculate object. Access results by fetching the result property.
    """

    def calculate(self, sample, pressure, temperature, **kwargs):
        #define constants
        a = 2565.0
        a_err = 362.0 #standard error
        b_Al2O3 = -1.997
        b_Al2O3_err = 0.706
        b_FeOt = -0.9275
        b_FeOt_err = 0.394
        b_Na2O = 2.736
        b_Na2O_err = 0.871
        c = 1.171
        c_err = 0.069
        d = -14.21
        d_err = 0.54
        XH2O_err = 0.148 #reported 1 sigma error

        melt_comp_X = sample.get_composition(units="molfrac")
        tempk = temperature + 273.15

        #get FeOt from iron input as FeO and Fe2O3
        FeOt = melt_comp_X["FeO"] + 0.8998*melt_comp_X["Fe2O3"]

        #TODO assert melt_comp_X["H2O"] > 0

        #note math.log with one argument passed returns natural log
        ln_fH2O = (((2.0 * math.log(melt_comp_X["H2O"])) - (a/tempk) - 
                    ((b_Al2O3 * melt_comp_X["Al2O3"] * (pressure/tempk)) + 
                     (b_FeOt * FeOt * (pressure/tempk)) + 
                     (b_Na2O * melt_comp_X["Na2O"] * (pressure/tempk))) - 
                     d)/c)

        return math.exp(ln_fH2O)


class calculate_degassed_fluid_composition(Calculate):
    """ Calculates the composition of fluid degassed from a melt as a magma moves from a deeper
    storage retion (magma a) to a more shallow storage (magma b).

    Parameters
    ----------
    deep_volatiles: tuple
        Tuple of volatile concentration in wt percent in magma a (H2O, CO2, S)
    shallow_volatiles: tuple
        Tuple of volatile concentration in wt percent in magma b (H2O, CO2, S)
    F:  float
        Optional. If no value is passed, F is set to 1.0.
        Value of F, or 100 minus the amount of crystallization differentiation required 
        to produce daughter magma b from parent magma a. This is equivalent to the 
        percentage residual melt that magma b represents. Use F if magma a and magma b are
        related by fractionalal crystallization. This adjusts appropriately for the difference in
        mass of the two magmatic bodies.
        
    Returns
    -------
    MagmaticFluid object
        Bulk fluid composition in wt percent with keys "H2O", "CO2", "S". NOTE: Fluid is not
        speciated as fO2 is not considerred.
    """

    def calculate(self, deep_volatiles, shallow_volatiles, F=1.0, **kwargs):
        deep_vols = deep_volatiles.get_composition(units='wtpercent')
        shallow_vols = shallow_volatiles.get_composition(units='wtpercent')

        delta_H2O = (deep_vols["H2O"] / F) - shallow_vols["H2O"]
        delta_CO2 = (deep_vols["CO2"] / F) - shallow_vols["CO2"]
        delta_S = (deep_vols["S"] / F) - shallow_vols["S"]

        sum_vols = delta_H2O + delta_CO2 + delta_S

        norm_H2O = 100.0 * delta_H2O / sum_vols
        norm_CO2 = 100.0 * delta_CO2 / sum_vols
        norm_S = 100.0 * delta_S / sum_vols

        result = {"H2O":norm_H2O, "CO2":norm_CO2, "S":norm_S}
        result_sample = sample_class.MagmaticFluid(result, units="wtpercent")

        return result_sample


class match(Calculate):
    """Runs the simple matching model and returns all possible combinations of sub_gases, scaled
    each at 0-100%, that could produce the surface_gas composition to within a user defined error
    (called "threshold").

    Parameters
    ----------
    sub_gases: dict
        Dictionary of dictionaries of gas compositions where keys are user-defined names for each
        gas and values are MagmaticFluid objects. Gas compositions will be simplified for matching.

    surface_gas: MagmaticFluid
        MagmaticFluid object with gas composition of surface gas. Gas compositions will be
        simplified for matching.

    threshold: float
        Optional. Default value is 0.10 (10%). The threshold defines how closely a combination of
        sub_gases must match the surface_gas, expressed as relative percent of the surface_gas
        value.

    Returns
    -------
    ???
    """
    def sums(self, length, total_sum, step):
        """Returns a list of all possible arrays of integers, where the sum of all array elements is
        'total_sum', where each array is 'length' values long.

        NOTE: This function can take an extremely long amount of time. Each +1 increase to the
        "length" value results in an exponential increase in computation time. Timing was tested on
        a 2014 MacBook Pro with 3 GHz processor, 16 GB memory, and solid state drive, using one
        core. Timing for length=5 ~2 seconds. Timing for length=6 ~230 seconds.

        Parameters
        ----------
        length: int
            Length of generated arrays. Example: length=3 gives a list of arrays (value, value,
            value)

        total_sum: int
            All values within each array must sum to total_sum. Example: length=3, total_sum=2
            gives: (2, 0, 0), (1, 0, 1), (1, 1, 0), (0, 0, 2), (0, 2, 0), (0, 1, 1)

        step: float
            Step size between each value of the array

        Returns
        -------
        generator object
            Generator object of all generated arrays. To return a list of arrays, pass
            list(sums(length, total_sum))

        """
        if length == 1:
            yield (total_sum,)
        else:
            for value in range(total_sum + 1):
                for permutation in self.sums(length - 1, total_sum - value, step):
                    remainder = value % step
                    if remainder == 0:
                        yield (value,) + permutation

    def calculate(self, sub_gases, surface_gas, threshold=0.1, **kwargs):
        # ensure only H2O, CO2, SO2 in surface_gas
        required_species = {"H2O": "H2O", "CO2": "CO2", "SO2": "SO2"}
        if set(surface_gas.get_composition().keys()) <= set(required_species.keys()):
            # get dict of composition from surface_gas
            surface_gas_dict = surface_gas.get_composition(units='wtpercent')
        else:
            surface_gas_dict = surface_gas.get_simplified_fluid_composition(units='wtpercent',
                                                                            S_species="SO2",
                                                                            warnings=False)

        # if user passing alread speciated composition, simplify it
        sub_gases_simpl = {}
        for gasname, gascomp in sub_gases.items():
            if set(gascomp.get_composition().keys()) <= set(required_species.keys()):
                composition_wtper = gascomp.get_composition(units='wtpercent') 
            else:
                composition_wtper = gascomp.get_simplified_fluid_composition(units='wtpercent',
                                                                             S_species="SO2",
                                                                             warnings=False)                
            sub_gases_simpl[gasname] = composition_wtper

        # split up gases into CO2, H2O, and SO2 groups
        sub_CO2 = [gas_dict["CO2"] for gas_name, gas_dict in sub_gases_simpl.items()]
        sub_H2O = [gas_dict["H2O"] for gas_name, gas_dict in sub_gases_simpl.items()]
        sub_SO2 = [gas_dict["SO2"] for gas_name, gas_dict in sub_gases_simpl.items()]

        # get number of sub_gases passed
        number_of_subgases = len(list(sub_gases.keys()))

        # generate all possible gas combinations
        print("Generating all possible gas combinations...")
        if number_of_subgases < 6:
            pos = list(self.sums(number_of_subgases, 100, 1))
        elif number_of_subgases == 6:
            pos = sums.sums_6_100_1() #TODO pulling in this file isn't working...
        # elif number_of_subgases == 7: #THIS DIDN'T EVEN FINISH ON MY COMPUTER FOR HOURS!
        #     pos = sums.sums_7_100_1()
        else:
            w.warn("This is a very large number of sub_gases and will generate a very large set of "
                   "possible gases. This may take hours...",
                   RuntimeWarning, stacklevel=2)
            pos = list(self.sums(number_of_subgases, 100, 1))

        print(pos)

        result_list = []
        sum_CO2 = 0.0
        sum_H2O = 0.0
        sum_SO2 = 0.0
        
        iterno = 0
        for i in range(len(sub_gases_simpl)):
            for combo in pos:
                iterno += 1
                percent = iterno/len(sub_gases_simpl)
                core.status_bar.status_bar(percent, list(sub_gases_simpl.keys())[i])
                sum_CO2 += sub_CO2[i] * (combo[i]/100.0)
                sum_H2O += sub_H2O[i] * (combo[i]/100.0)
                sum_SO2 += sub_SO2[i] * (combo[i]/100.0)

                if (sum_CO2 < (surface_gas_dict["CO2"] + surface_gas_dict["CO2"]*threshold) and
                    sum_CO2 > (surface_gas_dict["CO2"] - surface_gas_dict["CO2"]*threshold)):
                    if (sum_H2O < (surface_gas_dict["H2O"] + surface_gas_dict["H2O"]*threshold) and
                        sum_H2O > (surface_gas_dict["H2O"] - surface_gas_dict["H2O"]*threshold)):
                        if (sum_SO2 < (surface_gas_dict["SO2"] + surface_gas_dict["SO2"]*threshold) and
                            sum_SO2 > (surface_gas_dict["SO2"] - surface_gas_dict["SO2"]*threshold)):
                            result_list.append(combo)

        return result_list
