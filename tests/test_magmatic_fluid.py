import unittest
import tavern as tv
import numpy as np
import pandas as pd

class TestCreateSample(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'H2O':    47.00,
                                 'CO2':    35.00,
                                 'S':      16.00,
                                })

        # Standard normalization calculated externally
        self.majors_normed = pd.Series({'H2O':    47.959,
                                        'CO2':    35.714,
                                        'S':      16.327,
                                       })

        # Mol fractions calculated externally
        self.majors_molfrac = pd.Series({'H2O':   0.66841,
                                        'CO2':    0.20375,
                                        'S':      0.12784,
                                       })
        self.majors_molpercent = self.majors_molfrac*100

        # self.majors_all_species = pd.Series({'SiO2':    47.95,
        #                                      'TiO2':    1.67,
        #                                      'Al2O3':   17.32,
        #                                      'FeO':     10.24,
        #                                      'Fe2O3':   0.1,
        #                                      'MgO':     5.76,
        #                                      'CaO':     10.93,
        #                                      'Na2O':    3.45,
        #                                      'K2O':     1.99,
        #                                      'P2O5':    0.51,
        #                                      'MnO':     0.1,
        #                                      'Cr2O3':   0.05,
        #                                      'NiO':     0.04,
        #                                      'CoO':     0.02,
        #                                      'H2O':     2.00,
        #                                      'CO2':     0.12,
        #                                      'F2O':     0.05
        #                                         })

        self.sample = tv.MagmaticFluid(self.majors, units='wtpercent')

    def test_createSample(self):
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox],self.majors[ox])

    def test_setdefault_noargs(self):
        self.assertEqual(self.sample.default_normalization,'none')
        self.assertEqual(self.sample.default_units,'wtpercent')

    def test_setdefaults_none_wtpt(self):
        sample = tv.MagmaticFluid(self.majors, default_normalization='none',
                                  default_units='wtpercent')
        self.assertEqual(sample.default_normalization,'none')
        self.assertEqual(sample.default_units,'wtpercent')

    def test_setdefaults_standard_molfrac(self):
        sample = tv.MagmaticFluid(self.majors, default_normalization='standard',
                                  default_units='molfrac')
        self.assertEqual(sample.default_normalization,'standard')
        self.assertEqual(sample.default_units,'molfrac')

    def test_setdefaults_garbageNorm(self):
        with self.assertRaises(tv.core.InputError):
            tv.MagmaticFluid(composition=self.majors, default_normalization='garbage')

    def test_setdefaults_garbageType(self):
        with self.assertRaises(tv.core.InputError):
            tv.MagmaticFluid(composition=self.majors, default_units='garbage')

    def test_type_garbage(self):
        with self.assertRaises(tv.core.InputError):
            tv.MagmaticFluid(composition=self.majors, units='garbage')

    def test_type_wtptoxides(self):
        sample = tv.MagmaticFluid(self.majors,units='wtpercent')
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox], self.majors[ox])

    def test_type_molfrac(self):
        sample = tv.MagmaticFluid(self.majors_molfrac, units='molfrac')
        for ox in self.majors.index:
            self.assertEqual(np.round(sample._composition[ox],2),np.round(self.majors_normed[ox],2))

    def test_type_molpercent(self):
        sample = tv.MagmaticFluid(self.majors_molpercent,units='molpercent')
        for ox in self.majors.index:
            self.assertEqual(np.round(sample._composition[ox],2),
                             np.round(self.majors_normed[ox],2))


class TestGetComposition(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'H2O':    47.00,
                                 'CO2':    35.00,
                                 'S':      16.00,
                                })

        # Standard normalization calculated externally
        self.majors_normed = pd.Series({'H2O':    47.959,
                                        'CO2':    35.714,
                                        'S':      16.327,
                                       })

        # Mol fractions calculated externally
        self.majors_molfrac = pd.Series({'H2O':   0.66841,
                                        'CO2':    0.20375,
                                        'S':      0.12784,
                                       })
        self.majors_molpercent = self.majors_molfrac*100

        # self.majors_all_species = pd.Series({'SiO2':    47.95,
        #                                      'TiO2':    1.67,
        #                                      'Al2O3':   17.32,
        #                                      'FeO':     10.24,
        #                                      'Fe2O3':   0.1,
        #                                      'MgO':     5.76,
        #                                      'CaO':     10.93,
        #                                      'Na2O':    3.45,
        #                                      'K2O':     1.99,
        #                                      'P2O5':    0.51,
        #                                      'MnO':     0.1,
        #                                      'Cr2O3':   0.05,
        #                                      'NiO':     0.04,
        #                                      'CoO':     0.02,
        #                                      'H2O':     2.00,
        #                                      'CO2':     0.12,
        #                                      'F2O':     0.05
        #                                         })

        self.sample = tv.MagmaticFluid(self.majors)

    def test_default(self):
        composition = self.sample.get_composition()
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majors[ox])

    def test_wtptoxides_none(self):
        composition = self.sample.get_composition(normalization='none')
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majors[ox])

    def test_wtptoxides_std(self):
        composition = self.sample.get_composition(normalization='standard')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],2),np.round(self.majors_normed[ox],2))

    def test_molpercent(self):
        composition = self.sample.get_composition(units='molpercent')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_molpercent[ox],3))

    def test_molfrac(self):
        composition = self.sample.get_composition(units='molfrac')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_molfrac[ox],3))
