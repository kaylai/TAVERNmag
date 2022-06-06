"""
TAVERN

A thermodynamic model code for magmatic volatiles.
"""

__version__ = "0.1.0"
__author__ = "Kayla Iacovino"

# -------------------- IMPORTS -------------------- #
import tavern.core
import tavern.calculate_classes
import tavern.fO2_buffers
import tavern.sample_class

import pandas as pd

from copy import deepcopy, copy

# ------------ CALCULATION DEFINITIONS ------------ #
class calculate_fugacities(tavern.calculate_classes.calculate_fugacities):
	pass

class calculate_speciation(tavern.calculate_classes.calculate_speciation):
	pass

class calculate_fugacity_coefficients(tavern.calculate_classes.calculate_fugacity_coefficients):
	pass

class calculate_equilibrium_constants(tavern.calculate_classes.calculate_equilibrium_constants):
	pass

# ------------ SAMPLE CLASS DEFINITIONS ------------ #
class MagmaticFluid(sample_class.Sample):
    """ A fluid derived from a magma. A fluid must have the following properties:

    Attributes
    ----------
        composition: dict
            Dict of volatile concentration in wt percent, mol percent, or mol fraction 
            in fluid with keys H2O, CO2, S.
        units: str
            String defining whether fluid_comp is input as wt percent ("wtpercent"), 
            mole percent ("molpercent"), or mole fraction ("molfrac"). Default is "wtpercent".
        default_normalization:     None or str
            The type of normalization to apply to the data by default. One of:
                - None (no normalization)
                - 'standard' (default): Normalizes an input composition to 100%.
        default_units:     str
            The type of composition to return by default, one of:
            - wtpercent (default)
            - molpercent
            - molfrac
    """

    def __init__(self, composition, units='wtpercent', default_normalization='none',
                 default_units='wtpercent', **kwargs):
        """Return a MagmaticFluid object whose parameters are defined here."""
        composition = deepcopy(composition)
        
        if isinstance(composition, dict):
            composition = pd.Series(composition, dtype='float64')
        elif isinstance(composition, pd.Series) is False:
            raise core.InputError("The composition must be given as either a dictionary or a "
                                  "pandas Series.")
        
        if units == "wtpercent":
            self._composition = composition
        elif units == "molpercent":
            self._composition = self._moles_to_wtpercent(composition) 
        elif units == "molfrac":
            self._composition = self._moles_to_wtpercent(composition) 
        else:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        if default_normalization not in ['none', 'standard']:
            raise core.InputError("For a MagmaticFluid object, normalization can be one of "
                                  "'none' or 'standard'.")

        if default_units not in ['wtpercent', 'molpercent', 'molfrac']:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        composition = deepcopy(composition)
        #self._composition = composition
        self.units = units
        self.default_normalization = default_normalization
        self.default_units = default_units
        self.sample_type = 'MagmaticFluid'