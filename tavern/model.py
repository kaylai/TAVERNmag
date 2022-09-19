from tavern import core, fO2_buffers, sample_class, calculate_classes

import sympy
import math


class Model(object):
    """An object with arguments describing all needed thermodynamic parameters:

    Attributes
    ----------
        pressure:  float
            Pressure in bars.
        temperature:   float
            Temperature in degrees C.
        fO2_buffer: str
            Name of buffer for which to calculate fO2. Can be one of QFM.
        fO2_delta: float
            Deviation from named buffer in log units. For example, for QFM+2, enter 'QFM' as the
            fO2_buffer and 2 as the fO2_delta.
        H2O_param: float
            The H2O parameter: either fugacity, partial pressure, or mole fraction in the fluid.            
        H2O_param_type: str
            Possible param_type strings are: fugacity, ppressure, molfrac.
        SC_param: dict
            Either the sulfur fugacity or any CO2 parameter: fugacity, partial pressure, or
            mole fraction in the fluid.
        SC_param_type: str
            Possible param_type strings are: Sfugacity, CO2fugacity, CO2ppressure, CO2molfrac.
    """

    def __init__(self, pressure, temperature, fO2_buffer, fO2_delta, H2O_param, H2O_param_type,
                 SC_param, SC_param_type):
        """Return a Model object whose parameters are defined here."""
        self.pressure = pressure
        self.temperature = temperature
        self.fO2_buffer = fO2_buffer
        self.fO2_delta = fO2_delta
        self.H2O_param = H2O_param
        self.H2O_param_type = H2O_param_type
        self.SC_param = SC_param
        self.SC_param_type = SC_param_type
        self.logfO2 = fO2_buffers.calc_logfO2_from_buffer(pressure=pressure,
                                                          temperature=temperature,
                                                          buffer=fO2_buffer,
                                                          delta=fO2_delta)

        self.tempK = self.temperature + 273.15
        self.fO2 = 10.0 **self.logfO2

        #calculate gammas and Ks
        self.gammas = calculate_classes.calculate_fugacity_coefficients(temperature=self.temperature, 
                                                                        pressure=self.pressure,
                                                                        species="all").result
        self.Ks = calculate_classes.calculate_equilibrium_constants(temperature=self.temperature,
                                                                    species="all").result

        #TODO write if-then statements to handle any param_type of H2O_param and SC_param that might get passed.
        #For now, just write the case for testing with Augustine (pass fH2O, XCO2fluid)
        #H2O param: Possible param_type strings are: fugacity, ppressure, molfrac.
        #SC param: Possible param_type strings are: Sfugacity, CO2fugacity, CO2ppressure, CO2molfrac.

        #Parse H2O param type:
        if self.H2O_param_type == "fugacity":
            self.fH2O = self.H2O_param
            self.PH2O = self.fH2O / self.gammas["H2O"]

        if self.H2O_param_type == "ppressure":
            pass

        if self.H2O_param_type == "molfrac":
            pass

        #Parse SC param type:
        if SC_param_type == "Sfugacity":
            pass

        if SC_param_type == "CO2fugacity":
            pass

        if SC_param_type == "CO2ppressure":
            pass

        if SC_param_type == "CO2molfrac":
            self.XCO2 = self.SC_param
            self.fCO2 = self.pressure * self.gammas["CO2"] * self.XCO2
            self.PCO2 = self.fCO2 / self.gammas["CO2"]

    def get_fugacities(self):
        """Returns fugacities of all fluid species, if they have been calculated.

        Returns
        -------
        dict
            Fugacity of each fluid species.
        """
        if self.fugacities:
            return self.fugacities
        else:
            raise core.Error("Fugacities have not been calculated for this sample.") 

    def set_fugacities(self, fugacities):
        """Sets self.fugacities value for MagmaticFluid object so that they can be retrieved using
        the get_fugacities() method.

        Parameters
        ----------
        fugacities: dict
            Fugacities of all species.

        Returns
        -------
        Pass - sets self.fugacities as dict
        """
        self.fugacities = fugacities 

        pass

    def get_partial_pressures(self):
        """Returns partial pressures of all fluid species, if they have been calculated.

        Returns
        -------
        dict
            Partial pressure of each fluid species.
        """
        return self.partial_pressures

    def get_fugacity_coefficients(self):
        """Returns fugacity coefficients of all fluid species.

        Returns
        -------
        dict
            Fugacity coefficient of each fluid species.
        """
        return self.gammas

    def get_equilibrium_constants(self):
        """Returns equilibrium constants (K values) of all fluid reactions.

        Returns
        -------
        dict
            Equilibrium constants of each fluid reaction.
        """
        return self.Ks

    def equilibrium_fluid(self, return_units='molfrac'):
        """Return the composition of an equilibrium fluid.

        Parameters
        ----------
        self, inherited from Class

        return_units: str
            Optional. Returns composition of fluid in form passed by the user. Default
            value is 'molfrac', which returns as mole fraction. Other options are 'wtpercent',
            and 'molpercent'.

        Returns
        -------
        MagmaticFluid object
            TAVERN MagmaticFluid object with composition of equilibrium fluid as: CO, CO2, H2, H2O,
            H2S, O2, S2, and SO2. 
        """

        fO2 = self.fO2
        logfO2 = self.logfO2
        fH2O = self.fH2O
        fCO2 = self.fCO2
        pressure = self.pressure
        temperature = self.temperature
        Ks = self.Ks
        gammas = self.gammas

        PH2O = self.PH2O
        PCO2 = self.PCO2 #TODO edit this to be in an if statement when you write options for passing other SC_param types
        PO2 = fO2 / gammas["O2"]
        #All P's = CO2, H2O, O2 | All f's = CO2, H2O, O2

        #Calculate fH2 from fH2O and fO2, equation 2 Iacovino (2015) EPSL
        fH2 = fH2O / (Ks["H2O"] * math.sqrt(fO2))
        PH2 = fH2 / gammas["H2"]
        #All P's  = CO2, H2, H2O, O2 | All f's = CO2, H2, H2O, O2

        #calculate fCO from fCO2 and fO2
        fCO = self.fCO2 / (Ks["CO2"] * math.sqrt(fO2))
        PCO = fCO / gammas["CO"]
        #All P's  = CO, CO2, H2, H2O, O2 | All f's = CO, CO2, H2, H2O, O2

        #calculate PStot by difference
        PStot = pressure - (PCO + PCO2 + PH2 + PH2O + PO2)

        #Calculate fS2 with sympy, equation 7
        fS2 = sympy.symbols('fS2') #for sympy

        gS2 = gammas["S2"]
        KSO2 = Ks["SO2"]
        gSO2 = gammas["SO2"]
        KH2S = Ks["H2S"]
        gH2S = gammas["H2S"]
        
        equation = ((fS2 / gS2) + (KSO2 * sympy.sqrt(fS2) * fO2 / gSO2) + 
                    (KH2S * fH2 * sympy.sqrt(fS2) / gH2S) - PStot)

        fS2 = sympy.solve(equation, fS2)[0] #newly implemented sympy way
        
        PS2 = fS2 / gammas["S2"]
        #All P's  = CO, CO2, H2, H2O, O2, S2 | All f's = CO, CO2, H2, H2O, O2, S2

        fSO2 = Ks["SO2"] * math.sqrt(fS2) * fO2
        PSO2 = fSO2 / gammas["SO2"]
        #All P's  = CO, CO2, H2, H2O, O2, S2, SO2 | All f's = CO, CO2, H2, H2O, O2, S2, SO2

        fH2S = Ks["H2S"] * fH2 * math.sqrt(fS2)
        PH2S = fH2S / gammas["H2S"]
        #All P's  = CO, CO2, H2, H2O, H2S, O2, S2, SO2 | All f's = CO, H2, H2O, H2S, O2, S2

        partial_pressures = { "CO": PCO,
                              "CO2": PCO2,
                              "H2": PH2,
                              "H2O": PH2O,
                              "H2S": PH2S,
                              "O2": PO2,
                              "S2": PS2,
                              "SO2": PSO2 }
        self.partial_pressures = partial_pressures
        self.recalc_pressure = sum(partial_pressures.values()) #this is a sanity check, mostly, but keeps things internally consistent.

        fugacities = {    "CO": fCO,
                          "CO2": fCO2,
                          "H2": fH2,
                          "H2O": fH2O,
                          "H2S": fH2S,
                          "O2": fO2,
                          "S2": fS2,
                          "SO2": fSO2 }
        self.fugacities = fugacities


        mol_fractions = {}
        for key, value in partial_pressures.items():
            mol_fractions[key] = fugacities[key] / (gammas[key] * self.recalc_pressure)

        #For debugging:
        P_err = self.recalc_pressure - pressure

        if P_err != 0.0:
            print("Calculated pressure error is " + str(P_err) + " bars")

        retSamp = sample_class.MagmaticFluid(mol_fractions, units='molfrac', default_units='molfrac')
        if return_units == 'molfrac':
            return retSamp

        if return_units == 'wtpercent':
            return retSamp.get_composition(units='wtpercent', asSampleClass=True)

        if return_units == 'molpercent':
            return retSamp.get_composition(units='molpercent', asSampleClass=True)

    def respeciate(self, sample, pressure=None, temperature=None, fO2_buffer=None,
                   fO2_delta=None):
            """Takes in a fluid of given composition and returns that fluid respeciated at a given
            pressure.

            Parameters
            ----------
            sample: MagmaticFluid class
                Composition of the magmatic fluid as MagmaticFluid object.

            pressure: float
                Pressure at which to respeciate the fluid, in bars. If None, previously defined
                value will be used.

            temperature: float
                Temperature at which to respeciate the fluid, in degrees C. If None, previously
                defined value will be used.

            fO2_buffer: str
                Name of new buffer for which to calculate fO2. Can be one of 'QFM' or 'NNO'.  If
                None, previously defined value will be used.
        
            fO2_delta: float
                Deviation from named buffer in log units. For example, for QFM+2, enter 'QFM' as the
                fO2_buffer and 2 as the fO2_delta.  If None, previously defined value will be used.

            Returns
            -------
            MagmaticFluid object
                Fluid composition after respeciation
            """

            # see what parameters have changed and pass as necessary to respeciate_fluid() class
            if pressure is None:
                pressure = self.pressure

            if temperature is None:
                temperature = self.temperature

            if fO2_buffer is None:
                fO2_buffer = self.fO2_buffer

            if fO2_delta is None:
                fO2_delta = self.fO2_delta

            result = calculate_classes.respeciate_fluid(sample=sample, pressure=pressure,
                                                        temperature=temperature,
                                                        fO2_buffer=fO2_buffer,
                                                        fO2_delta=fO2_delta).result
            return result
