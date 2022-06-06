from abc import abstractmethod
import math
import sympy
import numpy as np

from tavern import core, fO2_buffers, sample_class

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
        Composition of the magmatic fluid as MagmaticFluid object.
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

    Returns
    -------
    Calculate object
        Calculate object. Access results by fetching the result property. Fugacities of all species.

    """
    def calculate(self, sample, pressure, temperature, fO2_buffer, fO2_delta, gammas='calculate',
                  K_vals='calculate', **kwargs):
        composition_molfrac = sample.get_composition(units='molfrac')
        logfO2 = fO2_buffers.calc_logfO2_from_buffer(pressure=pressure, temperature=temperature,
                                                     buffer=fO2_buffer, delta=fO2_delta)
        fO2 = 10.0**logfO2

        if gammas == 'calculate':
            gammas = calculate_fugacity_coefficients(temperature=temperature, pressure=pressure,
                                                     species="all").result

        if K_vals == 'calculate':
            K_vals = calculate_equilibrium_constants(temperature=temperature, species="all").result

        XH2Otot = composition_molfrac["H2O"]
        XCO2tot = composition_molfrac["CO2"]
        XStot = composition_molfrac["S"]

        XHtot = XH2Otot * (0.6666)
        XStot = XStot
        XCtot = XCO2tot * (0.3333)

        #FIRST calculate fH2 and fS2 using fsolve, two eqns; two unknowns (eqn 9 in Iacovino, 2015)
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

    def return_default_units(self, sample, calc_result, units='wtpercent', **kwargs):
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
        bulk_comp = sample_class.Sample(bulk_comp, units=units)

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
                                              K_vals=K_vals).result

        X_dict = {}
        for species in core.fluid_species_names:
            X = fugacities[species] / (gammas[species] * pressure)
            X_dict[species] = X

        # Normalize and fix O2 since fO2 is input and should not be adjusted during norm
        norm = self.normalize_fluid_FixedOxygen(X_dict, units='molfrac')

        speciated_default_units = self.return_default_units(sample, norm, units='molfrac', **kwargs)

        return speciated_default_units
