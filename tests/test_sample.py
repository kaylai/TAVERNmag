import unittest
import tavern as tv
import numpy as np

class TestCreateMagmaticFluid(unittest.TestCase):
    def setUp(self):
        self.comp = {"H2O": 79.0,
                     "CO2": 15.0,
                     "S":   6.0}
        
        self.comp_molpt = {"H2O": 89.2543955,
                           "CO2": 6.93707961,
                           "S": 3.80852486}
        
        self.comp_molfrac = {"H2O": 0.892543955,
                             "CO2": 0.0693707961,
                             "S": 0.0380852486}

        self.sample = tv.MagmaticFluid(self.comp)

    def test_createSample(self):
        for ox in self.comp.keys():
            self.assertEqual(self.sample._composition[ox],self.comp[ox])
    
    def test_setdefault_noargs(self):
        self.assertEqual(self.sample.default_normalization,'none')
        self.assertEqual(self.sample.default_units,'wtpercent')

    def test_setdefaults_none_wtpt(self):
        sample = tv.MagmaticFluid(self.comp, default_normalization='none',default_units='wtpercent')
        self.assertEqual(sample.default_normalization,'none')
        self.assertEqual(sample.default_units,'wtpercent')

    def test_setdefaults_standard_molpt(self):
        sample = tv.MagmaticFluid(self.comp, default_normalization='standard',
                                  default_units='molpercent')
        self.assertEqual(sample.default_normalization,'standard')
        self.assertEqual(sample.default_units,'molpercent')

    def test_setdefaults_garbageNorm(self):
        with self.assertRaises(tv.core.InputError):
            tv.MagmaticFluid(composition=self.comp,default_normalization='garbage')

    def test_setdefaults_garbageType(self):
        with self.assertRaises(tv.core.InputError):
            tv.MagmaticFluid(composition=self.comp,default_units='garbage')

    def test_units_garbage(self):
        with self.assertRaises(tv.core.InputError):
            tv.MagmaticFluid(composition=self.comp,units='garbage')

    def test_units_wtpt(self):
        sample = tv.MagmaticFluid(self.comp,units='wtpercent')
        for ox in self.comp.keys():
            self.assertEqual(sample._composition[ox],self.comp[ox])

    def test_units_molpt(self):
        sample = tv.MagmaticFluid(self.comp_molpt,units='molpercent')
        for ox in self.comp.keys():
            self.assertEqual(np.round(sample._composition[ox],2),np.round(self.comp[ox],2))

    def test_type_molfrac(self):
        sample = tv.MagmaticFluid(self.comp_molfrac,units='molfrac')
        for ox in self.comp.keys():
            self.assertEqual(np.round(sample._composition[ox],2),np.round(self.comp[ox],2))


class TestGetComposition(unittest.TestCase):
    def setUp(self):
        # Create generic samples to test Sample class
        self.generic = {"H2O":  20,
                        "SiO2": 50,
                        "FeO":  10,
                        "CaO":  10,
                        "MgO":  10}
        
        self.generic_molpt = {"H2O":    44.265803,
                              "SiO2":   33.1811178,
                              "FeO":    5.5498611,
                              "CaO":    7.11029871,
                              "MgO":    9.89291933}

        self.generic_molfrac = {"H2O":    0.44265803,
                                "SiO2":   0.331811178,
                                "FeO":    0.055498611,
                                "CaO":    0.0711029871,
                                "MgO":    0.0989291933}

        # Create fluid samples to test MagmaticFluid class
        self.fluid = {"H2O": 79.0,
                      "CO2": 15.0,
                      "S":   6.0}
        
        self.fluid_molpt = {"H2O": 89.2543955,
                            "CO2": 6.93707961,
                            "S": 3.80852486}
        
        self.fluid_molfrac = {"H2O": 0.892543955,
                              "CO2": 0.0693707961,
                              "S": 0.0380852486}
        
        # create silicate samples to test SilicateMelt class
        self.silicate = {'SiO2':    47.94,
                         'TiO2':    1.67,
                         'Al2O3':   17.32,
                         'FeO':     10.24,
                         'Fe2O3':   0.1,
                         'MgO':     5.76,
                         'CaO':     10.93,
                         'Na2O':    3.45,
                         'K2O':     1.99,
                         'P2O5':    0.51,
                         'MnO':     0.1}
        
        self.silicate_normed = {'SiO2':    47.94,
                                'TiO2':    1.67,
                                'Al2O3':   17.32,
                                'FeO':     10.24,
                                'Fe2O3':   0.1,
                                'MgO':     5.76,
                                'CaO':     10.93,
                                'Na2O':    3.45,
                                'K2O':     1.99,
                                'P2O5':    0.51,
                                'MnO':     0.1}
        
        self.silicate_molpt = {'SiO2':  51.429,
                               'TiO2':  1.348,
                               'Al2O3': 10.949,
                               'FeO':   9.187,
                               'Fe2O3': 0.04,
                               'MnO':   0.091,
                               'MgO':   9.212,
                               'CaO':   12.563,
                               'Na2O':  3.588,
                               'K2O':   1.362,
                               'P2O5': 0.232}
        
        self.silicate_molfrac = {'SiO2':    0.51429,
                                 'TiO2':    0.01348,
                                 'Al2O3':   0.10949,
                                 'FeO':     0.09187,
                                 'Fe2O3':   0.0004,
                                 'MnO':     0.00091,
                                 'MgO':     0.09212,
                                 'CaO':     0.12563,
                                 'Na2O':    0.03588,
                                 'K2O':     0.01362,
                                 'P2O5':    0.00232}
        
        self.silicatev = {'SiO2':   47.94,
                          'TiO2':   1.67,
                          'Al2O3':  17.32,
                          'FeO':    10.24,
                          'Fe2O3':  0.1,
                          'MgO':    5.76,
                          'CaO':    10.93,
                          'Na2O':   3.45,
                          'K2O':    1.99,
                          'P2O5':   0.51,
                          'MnO':    0.1,
                          'CO2':    0.08,
                          'H2O':    4.0}
        
        self.silicatev_normed = {'SiO2':  46.06,
                                'TiO2':   1.60,
                                'Al2O3':  16.64,
                                'FeO':    9.84,
                                'Fe2O3':  0.1,
                                'MgO':    5.53,
                                'CaO':    10.50,
                                'Na2O':   3.31,
                                'K2O':    1.91,
                                'P2O5':   0.49,
                                'MnO':    0.1,
                                'CO2':    0.08,
                                'H2O':    3.84}
        
        self.silicatev_molpt = {'SiO2':  44.944,
                                'TiO2':  1.178,
                                'Al2O3': 9.568,
                                'FeO':   8.029,
                                'Fe2O3': 0.035,
                                'MnO':   0.079,
                                'MgO':   8.050,
                                'CaO':   10.979,
                                'Na2O':  3.135,
                                'K2O':   1.190,
                                'P2O5':  0.202,
                                'H2O':   12.507,
                                'CO2':   0.102}
        
        self.silicatev_molfrac = {'SiO2':   0.45,
                                'TiO2':     0.01178,
                                'Al2O3':    0.09568,
                                'FeO':      0.08028,
                                'Fe2O3':    0.00035,
                                'MnO':      0.00079,
                                'MgO':      0.08050,
                                'CaO':      0.10978,
                                'Na2O':     0.03135,
                                'K2O':      0.01190,
                                'P2O5':     0.00202,
                                'H2O':      0.12503,
                                'CO2':      0.00102}
        
        # create Sample type objects
        self.genericSample = tv.sample_class.Sample(self.generic)
        self.fluidSample = tv.MagmaticFluid(self.fluid)
        self.silicateSample = tv.SilicateMelt(self.silicate)
        self.silicatevSample = tv.SilicateMelt(self.silicatev)
    
    def test_default(self):
        composition_gen = self.genericSample.get_composition()
        composition_fl = self.fluidSample.get_composition()
        composition_sil = self.silicateSample.get_composition()
        composition_silv = self.silicatevSample.get_composition()
        for ox in composition_gen.keys():
            self.assertEqual(composition_gen[ox],self.generic[ox])
        for ox in composition_fl.keys():
            self.assertEqual(composition_fl[ox],self.fluid[ox])
        for ox in composition_sil.keys():
            self.assertEqual(composition_sil[ox],self.silicate[ox])
        for ox in composition_silv.keys():
            self.assertEqual(composition_silv[ox],self.silicatev[ox])

    def test_wtptoxides_none(self):
        composition_gen = self.genericSample.get_composition(normalization='none')
        composition_fl = self.fluidSample.get_composition(normalization='none')
        composition_sil = self.silicateSample.get_composition(normalization='none')
        composition_silv = self.silicatevSample.get_composition(normalization='none')
        for ox in composition_gen.keys():
            self.assertEqual(composition_gen[ox],self.generic[ox])
        for ox in composition_fl.keys():
            self.assertEqual(composition_fl[ox],self.fluid[ox])
        for ox in composition_sil.keys():
            self.assertEqual(composition_sil[ox],self.silicate[ox])
        for ox in composition_silv.keys():
            self.assertEqual(composition_silv[ox],self.silicatev[ox])

    def test_wtptoxides_none_exclV(self):
        composition_sil = self.silicateSample.get_composition(normalization='none',
                                                              exclude_volatiles=True)
        composition_silv = self.silicatevSample.get_composition(normalization='none',
                                                                exclude_volatiles=True)
        for ox in composition_sil.keys():
            self.assertEqual(composition_sil[ox], self.silicate[ox])
        for ox in composition_silv.keys():
            self.assertEqual(composition_silv[ox], self.silicatev[ox])

    def test_wtptoxides_std(self):
        composition_silv = self.silicatevSample.get_composition(normalization='standard')
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicatev_normed[ox],2))

    def test_wtptoxides_std_exclV(self):
        composition_silv = self.silicatevSample.get_composition(normalization='standard',
                                                                exclude_volatiles=True)
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicate_normed[ox],2))

    def test_wtptoxides_fixedVolatiles(self):
        composition_silv = self.silicatevSample.get_composition(normalization='fixedvolatiles')
        for ox in composition_silv.keys():
            if ox not in ['CO2','H2O']:
                self.assertEqual(np.round(composition_silv[ox],2),
                                 np.round(self.silicate_normed[ox]*
                                          (100-self.silicatev['CO2']-self.silicatev['H2O'])/100,2))
            else:
                self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicatev[ox],2))

    def test_wtptoxides_fixedVolatiles_exclV(self):
        composition_silv = self.silicatevSample.get_composition(normalization='fixedvolatiles',
                                                                exclude_volatiles=True)
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicate_normed[ox],2))

    def test_wtptoxides_additionalVolatiles(self):
        composition_silv = self.silicatevSample.get_composition(normalization='additionalvolatiles')
        for ox in composition_silv.keys():
            if ox not in ['CO2','H2O']:
                self.assertEqual(np.round(composition_silv[ox],2),
                                 np.round(self.silicate_normed[ox],2))
            else:
                self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicatev[ox],2))

    def test_wtptoxides_additionalVolatiles_exclV(self):
        composition_silv = self.silicatevSample.get_composition(normalization='additionalvolatiles',
                                                                exclude_volatiles=True)
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicate_normed[ox],2))

    def test_molpt(self):
        composition_gen = self.genericSample.get_composition(units='molpercent')
        composition_fl = self.fluidSample.get_composition(units='molpercent')
        composition_sil = self.silicateSample.get_composition(units='molpercent')
        composition_silv = self.silicatevSample.get_composition(units='molpercent')
        for ox in composition_gen.keys():
            self.assertEqual(np.round(composition_gen[ox],3),np.round(self.generic_molpt[ox],3))
        for ox in composition_fl.keys():
            self.assertEqual(np.round(composition_fl[ox],3),np.round(self.fluid_molpt[ox],3))
        for ox in composition_sil.keys():
            self.assertEqual(np.round(composition_sil[ox],3),np.round(self.silicate_molpt[ox],3))
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],2),np.round(self.silicatev_molpt[ox],2))

    def test_molpt_exclV(self):
        composition_silv = self.silicatevSample.get_composition(units='molpercent',
                                                                exclude_volatiles=True)
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],3),np.round(self.silicate_molpt[ox],3))

    def test_molfrac(self):
        composition_silv = self.silicatevSample.get_composition(units='molfrac')
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],2),
                             np.round(self.silicatev_molfrac[ox],2))

    def test_molfrac_exclV(self):
        composition_silv = self.silicatevSample.get_composition(units='molfrac',
                                                                exclude_volatiles=True)
        for ox in composition_silv.keys():
            self.assertEqual(np.round(composition_silv[ox],3),np.round(self.silicate_molfrac[ox],3))