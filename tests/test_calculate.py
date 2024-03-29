import unittest
import math
import tavernmag as tv
import pandas as pd

class TestEquilibriumConstants(unittest.TestCase):
    def setUp(self):
        # parameters for the calculation
        self.temperature = 1000.0

        # Equilibrium constants calculated with TAVERN rounded to 2 dp
        self.eq_constants_known = {'CO2': 11609339.69,
                                   'H2O': 17101426.28,
                                   'H2S': 13.72,
                                   'SO2': 98548280692.69}
        self.eq_constants_known_log = {'CO2': math.log10(11609339.69),
                                       'H2O': math.log10(17101426.28),
                                       'H2S': math.log10(13.72),
                                       'SO2': math.log10(98548280692.69)}

    def test_calculate_equilibrium_constants_return_std(self):
        calcd_result = tv.calculate_equilibrium_constants(temperature=self.temperature,
                                                               return_as='standard').result 
        known_result = self.eq_constants_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=2)

    def test_calculate_equilibrium_constants_return_log(self):
        calcd_result = tv.calculate_equilibrium_constants(temperature=self.temperature,
                                                               return_as='log').result 
        known_result = self.eq_constants_known_log
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=2)


class TestFugacityCoefficients(unittest.TestCase):
    def setUp(self):
        # parameters for the calculation
        self.temperature = 1000.0
        self.pressure = 1000.0

        # Fugacity coefficients calculated with TAVERN rounded to 4 dp
        self.fugacity_coefficients_known = {'CO':   1.2557,
                                            'CO2':  1.1767,
                                            'H2':   1.1872,
                                            'H2O':  0.9014,
                                            'H2S':  1.1270,
                                            'O2':   1.1900,
                                            'S2':   1.1552,
                                            'SO2':  1.1376}

    def test_calculate_fugacity_coefficients(self):
        calcd_result = tv.calculate_fugacity_coefficients(temperature=self.temperature,
                                                          pressure=self.pressure).result
        known_result = self.fugacity_coefficients_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=4)

class TestFugacities(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'H2O':    47.00,
                                 'CO2':    35.00,
                                 'S':      16.00,
                                })
        self.majors_molpercent = pd.Series({'H2O':    66.841024,
                                            'CO2':    20.374933,
                                            'S':      12.784044,
                                           })
        self.majors_molfrac = pd.Series({'H2O':    0.668410,
                                         'CO2':    0.203749,
                                         'S':      0.127840,
                                        })

        self.sample = tv.MagmaticFluid(self.majors)
        self.sample_unitsMolPercent = tv.MagmaticFluid(self.majors_molpercent, units='molpercent')
        self.sample_unitsMolFrac = tv.MagmaticFluid(self.majors_molfrac, units='molfrac')

        # parameters for the calculation
        self.pressure = 1000.0
        self.temperature = 926.0
        self.fO2_buffer = 'QFM'
        self.fO2_delta = 1

        # fugacities calculated with TAVERN rounded to 4 dp
        self.fugacities_known = {'CO':  4.26919,
                                 'CO2': 713.5232351,
                                 'H2':  7.29685,
                                 'H2O': 1423.6482,
                                 'H2S': 400.9292817,
                                 'O2':  7.91101e-12,
                                 'S2':  4.97656,
                                 'SO2': 14.22078616}

    def test_calculate_fugacities_wtpercent(self):
        calcd_result = tv.calculate_fugacities(sample=self.sample, pressure=self.pressure,
                                               temperature=self.temperature,
                                               fO2_buffer=self.fO2_buffer,
                                               fO2_delta=self.fO2_delta).result
        known_result = self.fugacities_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=3)

    def test_calculate_fugacities_unitsMolPercent(self):
        calcd_result = tv.calculate_fugacities(sample=self.sample_unitsMolPercent, 
                                               pressure=self.pressure,
                                               temperature=self.temperature,
                                               fO2_buffer=self.fO2_buffer,
                                               fO2_delta=self.fO2_delta).result
        known_result = self.fugacities_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=3)

    def test_calculate_fugacities_unitsMolFrac(self):
        calcd_result = tv.calculate_fugacities(sample=self.sample_unitsMolFrac, 
                                               pressure=self.pressure,
                                               temperature=self.temperature,
                                               fO2_buffer=self.fO2_buffer,
                                               fO2_delta=self.fO2_delta).result
        known_result = self.fugacities_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=2)


class TestSpeciation(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'H2O':    47.00,
                                 'CO2':    35.00,
                                 'S':      16.00,
                                })
        self.majors_molpercent = pd.Series({'H2O':    66.841024,
                                            'CO2':    20.374933,
                                            'S':      12.784044,
                                           })
        self.majors_molfrac = pd.Series({'H2O':    0.668410,
                                         'CO2':    0.203749,
                                         'S':      0.127840,
                                        })

        self.sample = tv.MagmaticFluid(self.majors)
        self.sample_unitsMolPercent = tv.MagmaticFluid(self.majors_molpercent, units='molpercent')
        self.sample_unitsMolFrac = tv.MagmaticFluid(self.majors_molfrac, units='molfrac')

        # parameters for the calculation
        self.pressure = 1000.0
        self.temperature = 926.0
        self.fO2_buffer = 'QFM'
        self.fO2_delta = 1

        # speciation calculated with TAVERN to 4 dp
        self.speciation_known = {'CO':  0.096679984,
                                 'CO2': 38.35293618,
                                 'H2':  0.017629767,
                                 'H2O': 42.42249584,
                                 'H2S': 17.54972334,
                                 'S2':  0.395331538,
                                 'SO2': 1.165203359,
                                 'O2':  7.975e-13}

    def test_calculate_speciation_wtpercent(self):
        calcd_result = tv.calculate_speciation(sample=self.sample, pressure=self.pressure,
                                               temperature=self.temperature,
                                               fO2_buffer=self.fO2_buffer,
                                               fO2_delta=self.fO2_delta).result
        known_result = self.speciation_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=3)

    def test_calculate_speciation_unitsMolPercent(self):
        calcd_result = tv.calculate_speciation(sample=self.sample_unitsMolPercent,
                                               pressure=self.pressure,
                                               temperature=self.temperature,
                                               fO2_buffer=self.fO2_buffer,
                                               fO2_delta=self.fO2_delta).result
        known_result = self.speciation_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=3)

    def test_calculate_speciation_unitsMolFrac(self):
        calcd_result = tv.calculate_speciation(sample=self.sample_unitsMolFrac,
                                               pressure=self.pressure,
                                               temperature=self.temperature,
                                               fO2_buffer=self.fO2_buffer,
                                               fO2_delta=self.fO2_delta).result
        known_result = self.speciation_known
        for k in calcd_result.keys():
            self.assertAlmostEqual(calcd_result[k], known_result[k], places=3)

