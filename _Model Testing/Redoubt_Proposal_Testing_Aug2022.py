from DegassThermoModel import *
from math import *
import xlsxwriter

##----USER INPUTS----##
#Define regions: shallow_magma, int_magma, deep_magma

##NOTES ON MAGMA PARAMETERS FROM THE LITERATURE##
"""
REDOUBT COOMBS ET AL., 2013
---------------------------
SAMPLES
For the 2009 eruption, three magma types were involved: Low silica andesite (LSA), intermediate
silica andesite (ISA) and high silica andesite (HSA). All are thought to have staged in the upper
crust prior to eruption. The LSA is the most likely candidate for having ascended from much
greater depth prior to eruption.

HIGH-SILICA ANDESITE
**Shallow Storage**
glass SiO2 > 76 wt%
Effusive-phase ISA and HSA
T = 800 C
fO2 = NNO+1.44
P = 1500 bars??? (nd)

INTERMEDIATE-SILICA ANDESITE
**Shallow Storage**
Does this reflect mixing between LSA and HSA???
T = 927 C
fO2 = NNO+1.44
P = 1500 bars

LOW-SILICA ANDESITE
**Shallow Storage**
This is the LSA with conditions during shallow staging/mixing
Glass <69 wt% SiO2
T = 950 C
fO2 = NNO+1.35
P = 1500 bars

LSA DEEP MAGMA
**Deep Storage**
This is the LSA that may have ascended from a greater depth
T = 980?
fO2 = ?
P = >2000 bars


OXYGEN FUGACITY
Coexisting iron-titanium oxides.

SULFUR

"""

# MELT COMPOSITIONS
# High-silica andesite glass from Coombs et al., 2013
# Average of table 3 samples >76 wt% SiO2
HSA_melt = {"SiO2": 78.38, 
					"TiO2": 0.28,
					"Al2O3": 11.71,
					"FeO": 0.98,
					"Fe2O3": 0.,
					"MnO": 0.04,
					"MgO": 0.12,
					"CaO": 0.58,
					"Na2O": 3.43,
					"K2O": 4.38,
					"P2O5": 0.03,
					"H2O": 4.5, # table 4 ???
					"CO2": 0.8, # ????
					"S": 0.02 # ????
					}

HSA_vols = (HSA_melt["H2O"],
			HSA_melt["CO2"],
			HSA_melt["S"])

# int-silica andesite
# ISA glass from Coombs et al., 2013
# Average of table 3 samples 72-76 wt% SiO2
ISA_melt = {"SiO2": 73.76, #Values in wt% 
					"TiO2": 0.45, 
					"Al2O3": 13.82,
					"FeO": 2.04, 
					"Fe2O3": 0.,
					"MnO": 0.06,
					"MgO": 0.48,
					"CaO": 1.74,
					"Na2O": 4.25,
					"K2O": 3.18,
					"P2O5": 0.09,
					"H2O": 4.5, # table 4 
					"CO2": 0.011, # ????
					"S": 0.02 # ????
					}

ISA_vols = (ISA_melt["H2O"],
			ISA_melt["CO2"],
			ISA_melt["S"])

# Low-silica andesite glass from Coombs et al., 2013
# Average of table 3 samples <69 wt% SiO2
LSA_shallow_melt = {"SiO2": 67.53, 
					"TiO2": 0.63,
					"Al2O3": 15.55,
					"FeO": 3.96,
					"Fe2O3": 0.,
					"MnO": 0.12,
					"MgO": 1.96,
					"CaO": 3.4,
					"Na2O": 4.67,
					"K2O": 2.59,
					"P2O5": 0.24,
					"H2O": 4.5, # section 6.3
					"CO2": 0.8, # section 6.3
					"S": 0.02 # ????
					}

LSA_shallow_vols = (LSA_shallow_melt["H2O"],
					LSA_shallow_melt["CO2"],
					LSA_shallow_melt["S"])

# Low-silica andesite glass from Coombs et al., 2013
# Average of table 3 samples <69 wt% SiO2
LSA_deep_melt = {"SiO2": 67.53, 
					"TiO2": 0.63,
					"Al2O3": 15.55,
					"FeO": 3.96,
					"Fe2O3": 0.,
					"MnO": 0.12,
					"MgO": 1.96,
					"CaO": 3.4,
					"Na2O": 4.67,
					"K2O": 2.59,
					"P2O5": 0.24,
					"H2O": 4.5, # ????
					"CO2": 0.8, # ????
					"S": 0.02 # ????
					}

LSA_deep_vols = (LSA_deep_melt["H2O"],
				 LSA_deep_melt["CO2"],
				 LSA_deep_melt["S"])

##############################
# calculate equilibrium fluids
# HSA
HSA_melt = SilicateMelt(HSA_melt)
HSA_XCO2 = 0.10 # ????
HSA_temp = 800.0
HSA_press = 1500.0
HSA_fH2O = HSA_melt.calc_fH2O(HSA_temp, HSA_press)
HSA_obj = Model(temp=HSA_temp, press=HSA_press, fO2_buffer='NNO', fO2_delta=1.44, H2O_param=HSA_fH2O, 
				H2O_param_type="fugacity", SC_param=HSA_XCO2, SC_param_type="CO2molfrac")
HSA_eq_highP = HSA_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
HSA_eq = HSA_obj.respeciate(HSA_eq_highP, newP=1, input_type='molfrac')


# ISA
ISA_melt = SilicateMelt(ISA_melt)
ISA_XCO2 = 0.10 # ????
ISA_temp = 927.0
ISA_press = 1500.0
ISA_fH2O = int_melt.calc_fH2O(ISA_temp, ISA_press)
ISA_obj = Model(temp=ISA_temp, press=ISA_press, fO2_buffer='NNO', fO2_delta=1.44, H2O_param=ISA_fH2O, 
				H2O_param_type="fugacity", SC_param=ISA_XCO2, SC_param_type="CO2molfrac")
ISA_eq_highP = ISA_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
ISA_eq = ISA_obj.respeciate(ISA_eq_highP, newP=1, input_type='molfrac')

# LSA SHALLOW
LSA_shallow_melt = SilicateMelt(LSA_shallow_melt)
LSA_shallow_XCO2 = 0.10 # ????
LSA_shallow_temp = 950.0
LSA_shallow_press = 1500.0
LSA_shallow_fH2O = LSA_shallow_melt.calc_fH2O(LSA_shallow_temp, LSA_shallow_press)
LSA_shallow_obj = Model(temp=LSA_shallow_temp, press=LSA_shallow_press, fO2_buffer='NNO',
						fO2_delta=1.35, H2O_param=LSA_shallow_fH2O, H2O_param_type="fugacity",
						SC_param=LSA_shallow_XCO2, SC_param_type="CO2molfrac")
LSA_shallow_eq_highP = LSA_shallow_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
LSA_shallow_eq = LSA_shallow_obj.respeciate(LSA_shallow_eq_highP, newP=1, input_type='molfrac')

# LSA DEEP
LSA_deep_melt = SilicateMelt(LSA_deep_melt)
LSA_deep_XCO2 = 0.10 # ????
LSA_deep_temp = 980.0
LSA_deep_press = 3000.0
LSA_deep_fH2O = LSA_deep_melt.calc_fH2O(LSA_deep_temp, LSA_deep_press)
LSA_deep_obj = Model(temp=LSA_deep_temp, press=LSA_deep_press, fO2_buffer='NNO',
						fO2_delta=1.35, H2O_param=LSA_deep_fH2O, H2O_param_type="fugacity",
						SC_param=LSA_deep_XCO2, SC_param_type="CO2molfrac")
LSA_deep_eq_highP = LSA_deep_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
LSA_deep_eq = LSA_deep_obj.respeciate(LSA_deep_eq_highP, newP=1, input_type='molfrac')

###########################
# calculate degassing fluids
# HSA DEGASSING (shallow to surface)
HSA_degassing_un = Degas_MI(HSA_vols, (0,0,0))
HSA_degassing_obj = MagmaticFluid(HSA_degassing_un, press=1, temp=HSA_temp, fO2_buffer='NNO',
								  fO2_delta=1.44)
HSA_degassing = HSA_degassing_obj.speciate(return_as='molfrac')

# ISA DEGASSING (shallow to surface)
ISA_degassing_un = Degas_MI(ISA_vols, (0,0,0))
ISA_degassing_obj = MagmaticFluid(ISA_degassing_un, press=1, temp=ISA_temp, fO2_buffer='NNO',
								  fO2_delta=1.44)
ISA_degassing = ISA_degassing_obj.speciate(return_as='molfrac')

# LSA SHALLOW DEGASSING (shallow to surface)
LSA_shallow_degassing_un = Degas_MI(LSA_shallow_vols, (0,0,0))
LSA_shallow_degassing_obj = MagmaticFluid(LSA_shallow_degassing_un, press=1, temp=LSA_shallow_temp,
										  fO2_buffer='NNO', fO2_delta=1.35)
LSA_shallow_degassing = LSA_shallow_degassing_obj.speciate(return_as='molfrac')

# LSA DEEP DEGASSING (deep to shallow)
LSA_deep_degassing_un = Degas_MI(LSA_deep_vols, LSA_shallow_vols)
LSA_deep_degassing_obj = MagmaticFluid(LSA_deep_degassing_un, press=1, temp=LSA_deep_temp,
										  fO2_buffer='NNO', fO2_delta=1.35)
LSA_deep_degassing = LSA_deep_degassing_obj.speciate(return_as='molfrac')

# sub_gases = {#"Shallow_EQ": shallow_eq,
# 			 "Int_EQ": int_eq,
# 			 "Deep_EQ": deep_eq,
# 			 "Shallow_Degassing": shallow_degassing,
# 			 "Int_Degassing": int_degassing,
# 			 "Deep_Degassing": deep_degassing}


# surface_gas = {"H2O": 0.9238, "CO2": 0.0330, "SO2": 0.0432}

# mixer = Match(sub_gases, surface_gas, input_type='molfrac')

# mixer.matchmodel(threshold=0.5)

print("HSA EQ Fluid = ")
print(HSA_eq)
print(sum(HSA_eq.values()))
print("ISA EQ Fluid = ")
print(ISA_eq)
print(sum(ISA_eq.values()))
print("LSA Shallow EQ Fluid = ")
print(LSA_shallow_eq)
print(sum(LSA_shallow_eq.values()))
print("LSA Deep EQ Fluid = ")
print(LSA_deep_eq)
print(sum(LSA_deep_eq.values()))
print("HSA degassing = ")
print(HSA_degassing)
print(sum(HSA_degassing.values()))
print("ISA degassing = ")
print(ISA_degassing)
print(sum(ISA_degassing.values()))
print("LSA Shallow degassing = ")
print(LSA_shallow_degassing)
print(sum(LSA_shallow_degassing.values()))
print("LSA Deep degassing = ")
print(LSA_deep_degassing)
print(sum(LSA_deep_degassing.values()))

workbook = xlsxwriter.Workbook('Redoubt_output.xlsx')
worksheet = workbook.add_worksheet()

worksheet.write(1,0, "HSA EQ Fluid")
worksheet.write(3,0, "ISA EQ Fluid")
worksheet.write(5,0, "LSA Shallow EQ Fluid")
worksheet.write(7,0, "LSA Deep EQ Fluid")
worksheet.write(9,0, "HSA degassing")
worksheet.write(11,0, "ISA degassing")
worksheet.write(13,0, "LSA Shallow degassing")
worksheet.write(15,0, "LSA Deep degassing")

row = 0
col = 1
for key in HSA_eq.keys():
	worksheet.write(row, col, key)
	col += 1
col = 1
row += 2

for key in ISA_eq.keys():
	worksheet.write(row, col, key)
	col += 1
col = 1
row += 2

for key in LSA_shallow_eq.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in LSA_deep_eq.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in HSA_degassing.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in ISA_degassing.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in LSA_shallow_degassing.keys():
	worksheet.write(row, col, key)
	col += 1

for key in LSA_deep_degassing.keys():
	worksheet.write(row, col, key)
	col += 1

row = 1
col = 1
for key, value in HSA_eq.items():
	worksheet.write(row, col, value)
	col += 1

row = 3
col = 1
for key, value in ISA_eq.items():
	worksheet.write(row, col, value)
	col += 1

row = 5
col = 1
for key, value in LSA_shallow_eq.items():
	worksheet.write(row, col, value)
	col += 1

row = 7
col = 1
for key, value in LSA_deep_eq.items():
	worksheet.write(row, col, value)
	col += 1

row = 9
col = 1
for key, value in HSA_degassing.items():
	worksheet.write(row, col, value)
	col += 1

row = 11
col = 1
for key, value in ISA_degassing.items():
	worksheet.write(row, col, value)
	col += 1

row = 13
col = 1
for key, value in LSA_shallow_degassing.items():
	worksheet.write(row, col, value)
	col += 1

row = 15
col = 1
for key, value in LSA_deep_degassing.items():
	worksheet.write(row, col, value)
	col += 1

workbook.close()








