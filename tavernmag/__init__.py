"""
TAVERNmag

A thermodynamic model code for magmatic volatiles.
"""

__version__ = "0.2.0"
__author__ = "Kayla Iacovino"

# -------------------- IMPORTS -------------------- #
import tavernmag.core
import tavernmag.calculate_classes
import tavernmag.fO2_buffers
import tavernmag.sample_class
import tavernmag.model

import pandas as pd

from copy import deepcopy, copy

# ------------ CALCULATION DEFINITIONS ------------ #
class calculate_fugacities(tavernmag.calculate_classes.calculate_fugacities):
	pass


class calculate_speciation(tavernmag.calculate_classes.calculate_speciation):
	pass


class calculate_fugacity_coefficients(tavernmag.calculate_classes.calculate_fugacity_coefficients):
	pass


class calculate_equilibrium_constants(tavernmag.calculate_classes.calculate_equilibrium_constants):
	pass


class calculate_degassed_fluid_composition(tavernmag.calculate_classes.calculate_degassed_fluid_composition):
    pass


class calculate_fH2O_from_melt(tavernmag.calculate_classes.calculate_fH2O_from_melt):
    pass

class match(tavernmag.calculate_classes.match):
    pass

# ------------ MODEL CLASS DEFINITIONS ------------ #
class Model(tavernmag.model.Model):
    pass

# ------------ SAMPLE CLASS DEFINITIONS ------------ #
class MagmaticFluid(sample_class.MagmaticFluid):
    """ A generic class describing a fluid derived from a magma. A fluid must have the following
    properties:

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
    pass


class SilicateMelt(sample_class.SilicateMelt):
    """ A silicate melt major oxide composition. A melt must have the following properties:

    Attributes
    ----------
        composition: dict
            Dict of major oxide concentrations in wt percent, mol percent, or mol fraction.
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
    pass
