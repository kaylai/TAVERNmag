from DegassThermoModel import *
from math import *
import xlsxwriter

##----USER INPUTS----##
#Define regions: shallow_magma, int_magma, deep_magma

##NOTES ON MAGMA PARAMETERS FROM THE LITERATURE##
"""
AUGUSTINE WEBSTER ET AL., 2010
------------------------------
SAMPLES
For the 2006 eruption we refer to the shallow storage as the Low Silica Andesite (LSA), 
the intermediate storage as High Silica Andesite (HSA) and the deep primitive source is 
inferred but not well constrained.

LOW-SILICA ANDESITE
**Shallow Storage**
T = 900-970 deg C (970 might be the primitive recharge)
fO2 = 
P = 100-200 MPa (4-8 km)

HIGH-SILICA ANDESITE
**Intermediate Storage**
T = 850-880 deg C
fO2 = 
P = 130-160 MPa

DEEP PRIMITIVE MAGMA
**Deep Storage**
T = 970-? deg C
fO2 = 
P = 

GENERAL
The magmatic fluids were relatively oxidizing and included H2O-enriched 
S-, S2-, and SO2+/- CO2-bearing vapors; hydrosa-line aqueous 
liquids largely enriched in Cl-, SO42-, alkalis, and H2O; and moderately 
saline, H2O-poor liquids containing Cl-, SO42-, and alkali elements.

OXYGEN FUGACITY
Coexisting iron-titanium
oxides of 2006 rock samples, which are generally consistent
with those of prior eruptive materials, indicate fO2
values of approximately NNO+1.5 to NNO+2.5 and oxide crystalliza-
tion temperatures of 835 to 1,052 degC.

SULFUR
As sulfides are rare and only present in Augustine rocks
at abundance levels well below that required to remove
significant S from the residual melts and because magma
mixing is not solely responsible for the strong reduction in
S concentrations with melt evolution, the involvement of
a fluid or fluid phases in S sequestration is required.
"""

#MELT COMPOSITIONS
#Low-silica andesite groundmass from Larsen et al. 2010
shallow_melt = {	"SiO2": 71.16, 
					"TiO2": 0.97, # 06AUMC005c
					"Al2O3": 13.72,
					"FeO": 4.03,
					"Fe2O3": 0.,
					"MnO": 0.12,
					"MgO": 0.90,
					"CaO": 3.13,
					"Na2O": 3.37,
					"K2O": 2.06,
					"P2O5": 0.23,
					"H2O": 2.00, #TODO need better vals
					"CO2": 0.002,
					"S": 0.02
					}

shallow_vols = (shallow_melt["H2O"],
				shallow_melt["CO2"],
				shallow_melt["S"])

#high-silica andesite
#HSA matrix glass from Larsen et al. 2010
int_melt = {	"SiO2": 75.86, #Values in wt% 
					"TiO2": 0.44, 
					"Al2O3": 12.81,
					"FeO": 2.08, 
					"Fe2O3": 0.,
					"MnO": 0.08,
					"MgO": 0.41,
					"CaO": 2.09,
					"Na2O": 3.78,
					"K2O": 0.99,
					"P2O5": 0.14,
					"H2O": 4.10, #Highest reported evolved augustine MI in Webster et al 2010
					"CO2": 0.011, #Highest reported evolved augustine MI in Webster et al 2010
					"S": 0.02
					}

int_vols = (int_melt["H2O"],
			int_melt["CO2"],
			int_melt["S"])


#primitive mafic recharge
deep_melt = {	"SiO2": 54.92, #Values in wt% #Augustine basalt sample from Webster et al 2010 (cpx-hosted MI)
					"TiO2": 0.62, #RB-W91A137 cpx MI
					"Al2O3": 14.62,
					"FeO": 6.86, 
					"Fe2O3": 0.,
					"MnO": 0.16,
					"MgO": 5.89,
					"CaO": 10.33,
					"Na2O": 2.04,
					"K2O": 0.61,
					"P2O5": 0.15,
					"H2O": 4.00, #Water by difference from reheated MI, max value
					"CO2": 0.02, #TODO... this is a WAG
					"S": 0.23
					}

deep_vols = (deep_melt["H2O"],
			 deep_melt["CO2"],
			 deep_melt["S"])

#calculate equilibrium fluids
#SHALLOW
shallow_melt = SilicateMelt(shallow_melt)
shallow_XCO2 = 0.10
shallow_temp = 900.0
shallow_press = 1000.0
shallow_fH2O = shallow_melt.calc_fH2O(shallow_temp, shallow_press)
shallow_obj = Model(temp=shallow_temp, press=shallow_press, fO2_buffer='QFM', fO2_delta=1.0, H2O_param=shallow_fH2O, 
				H2O_param_type="fugacity", SC_param=shallow_XCO2, SC_param_type="CO2molfrac")
shallow_eq_highP = shallow_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
shallow_eq = shallow_obj.respeciate(shallow_eq_highP, newP=1, input_type='molfrac')


#INTERMEDIATE
int_melt = SilicateMelt(int_melt)
int_XCO2 = 0.10
int_temp = 850.0
int_press = 1500.0
int_fH2O = int_melt.calc_fH2O(int_temp, int_press)
int_obj = Model(temp=int_temp, press=int_press, fO2_buffer='QFM', fO2_delta=1.0, H2O_param=int_fH2O, 
				H2O_param_type="fugacity", SC_param=int_XCO2, SC_param_type="CO2molfrac")
int_eq_highP = int_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
int_eq = int_obj.respeciate(int_eq_highP, newP=1, input_type='molfrac')

#DEEP
deep_melt = SilicateMelt(deep_melt)
deep_XCO2 = 0.10
deep_temp = 970.0
deep_press = 3000.0
deep_fH2O = deep_melt.calc_fH2O(deep_temp, deep_press)
deep_obj = Model(temp=deep_temp, press=deep_press, fO2_buffer='QFM', fO2_delta=1.0, H2O_param=deep_fH2O, 
				H2O_param_type="fugacity", SC_param=deep_XCO2, SC_param_type="CO2molfrac")
deep_eq_highP = deep_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
deep_eq = deep_obj.respeciate(deep_eq_highP, newP=1, input_type='molfrac')

#calculate degassing fluids
#SHALLOW DEGASSING (shallow to surface)
shallow_degassing_un = Degas_MI(shallow_vols, (0,0,0)) #TODO get matrix glass values for erupted LSA
shallow_degassing_obj = MagmaticFluid(shallow_degassing_un, press=1, temp=shallow_temp, fO2_buffer='QFM', fO2_delta=1.0)
shallow_degassing = shallow_degassing_obj.speciate(return_as='molfrac')

#INT DEGASSING (int to shallow)
int_degassing_un = Degas_MI(int_vols, shallow_vols)
int_degassing_obj = MagmaticFluid(int_degassing_un, press=1, temp=int_temp, fO2_buffer='QFM', fO2_delta=1.0)
int_degassing = int_degassing_obj.speciate(return_as='molfrac')

#DEEP DEGASSING (deep to int)
#TODO, figure out what to degas this to, since it's not liquid line of descent to INT magma, necessarily
#TODO cont... model sat'n volatile comp at lower P?
deep_degassing_un = Degas_MI(deep_vols, (0,0,0))
deep_degassing_obj = MagmaticFluid(deep_degassing_un, press=1, temp=deep_temp, fO2_buffer='QFM', fO2_delta=1.0)
deep_degassing = deep_degassing_obj.speciate(return_as='molfrac')

# sub_gases = {#"Shallow_EQ": shallow_eq,
# 			 "Int_EQ": int_eq,
# 			 "Deep_EQ": deep_eq,
# 			 "Shallow_Degassing": shallow_degassing,
# 			 "Int_Degassing": int_degassing,
# 			 "Deep_Degassing": deep_degassing}


# surface_gas = {"H2O": 0.9238, "CO2": 0.0330, "SO2": 0.0432}

# mixer = Match(sub_gases, surface_gas, input_type='molfrac')

# mixer.matchmodel(threshold=0.5)

print "Shallow EQ Fluid = "
print shallow_eq
print sum(shallow_eq.values())
print "Int EQ Fluid = "
print int_eq
print sum(int_eq.values())
print "Deep EQ Fluid = "
print deep_eq
print sum(deep_eq.values())
print "Shallow degassing = "
print shallow_degassing
print sum(shallow_degassing.values())
print "Int degassing = "
print int_degassing
print sum(int_degassing.values())
print "Deep degassing = "
print deep_degassing
print sum(deep_degassing.values())

workbook = xlsxwriter.Workbook('Augustine_output.xlsx')
worksheet = workbook.add_worksheet()

worksheet.write(1,0, "Shallow EQ Fluid")
worksheet.write(3,0, "Int EQ Fluid")
worksheet.write(5,0, "Deep EQ Fluid")
worksheet.write(7,0, "Shallow degassing")
worksheet.write(9,0, "Int degassing")
worksheet.write(11,0, "Deep degassing")

row = 0
col = 1
for key in shallow_eq.keys():
	worksheet.write(row, col, key)
	col += 1
col = 1
row += 2

for key in int_eq.keys():
	worksheet.write(row, col, key)
	col += 1
col = 1
row += 2

for key in deep_eq.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in shallow_degassing.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in int_degassing.keys():
	worksheet.write(row, col, key)
	col += 1
row += 2
col = 1

for key in deep_degassing.keys():
	worksheet.write(row, col, key)
	col += 1

row = 1
col = 1
for key, value in shallow_eq.iteritems():
	worksheet.write(row, col, value)
	col += 1

row = 3
col = 1
for key, value in int_eq.iteritems():
	worksheet.write(row, col, value)
	col += 1

row = 5
col = 1
for key, value in deep_eq.iteritems():
	worksheet.write(row, col, value)
	col += 1

row = 7
col = 1
for key, value in shallow_degassing.iteritems():
	worksheet.write(row, col, value)
	col += 1

row = 9
col = 1
for key, value in int_degassing.iteritems():
	worksheet.write(row, col, value)
	col += 1

row = 11
col = 1
for key, value in deep_degassing.iteritems():
	worksheet.write(row, col, value)
	col += 1

workbook.close()








