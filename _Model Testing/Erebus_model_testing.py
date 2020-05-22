from DegassThermoModel import *
from math import *

#MELT COMPOSITIONS
#EA1 average
# shallow_melt = {	"SiO2": 55.096, 
# 					"TiO2": 0.955, 
# 					"Al2O3": 19.919,
# 					"FeO": 4.997,
# 					"Fe2O3": 0.,
# 					"MnO": 0.24,
# 					"MgO": 0.82,
# 					"CaO": 1.833,
# 					"Na2O": 8.842,
# 					"K2O": 6.302,
# 					"P2O5": 0.305,
# 					"H2O": 0.167, 
# 					"CO2": 0.0691,
# 					"S": 0.03805
# 					}

# shallow_vols = (shallow_melt["H2O"],
# 				shallow_melt["CO2"],
# 				shallow_melt["S"])

# #Tephriphonolite 97009, 97010, 97011, and AW-82033 (same avg used in paper)
# int_melt = {	"SiO2": 49.83, #Values in wt% 
# 					"TiO2": 2.38, 
# 					"Al2O3": 18.67,
# 					"FeO": 7.33, 
# 					"Fe2O3": 1.65,
# 					"MnO": 0.22,
# 					"MgO": 1.99,
# 					"CaO": 5.87,
# 					"Na2O": 6.97,
# 					"K2O": 4.17,
# 					"P2O5": 0.91,
# 					"H2O": 0.28, 
# 					"CO2": 0.1638, 
# 					"S": 0.0918
# 					}

# int_vols = (int_melt["H2O"],
# 			int_melt["CO2"],
# 			int_melt["S"])

#Basanite DVDP3-295
deep_melt = {	"SiO2": 41.70, #Values in wt% 
					"TiO2": 4.21, 
					"Al2O3": 14.80,
					"FeO": 8.05, 
					"Fe2O3": 5.72,
					"MnO": 0.16,
					"MgO": 5.95,
					"CaO": 13.07,
					"Na2O": 3.74,
					"K2O": 1.65,
					"P2O5": 0.94,
					"H2O": 1.50, 
					"CO2": 0.5537, 
					"S": 0.2166
					}

deep_vols = (deep_melt["H2O"],
			 deep_melt["CO2"],
			 deep_melt["S"])

# #calculate equilibrium fluids
# #SHALLOW
# shallow_melt = SilicateMelt(shallow_melt)
# shallow_logfS2 = -2.7
# shallow_fS2 = 10**shallow_logfS2
# shallow_logfO2 = -11.56 #QFM+1
# shallow_temp = 900.0
# shallow_press = 1000.0
# shallow_fH2O = shallow_melt.calc_fH2O(shallow_temp, shallow_press)
# shallow_obj = Model(temp=shallow_temp, press=shallow_press, logfO2=shallow_logfO2, H2O_param=shallow_fH2O, 
# 				H2O_param_type="fugacity", SC_param=shallow_XCO2, SC_param_type="CO2molfrac")
# shallow_eq_highP = shallow_obj.eq_fluid(return_as='molfrac')
# #respecieate at 1 bar
# shallow_eq = shallow_obj.respeciate(shallow_eq_highP, newP=1, input_type='molfrac')


# #INTERMEDIATE
# int_melt = SilicateMelt(int_melt)
# int_XCO2 = 0.10
# int_logfO2 = -12.46 #QFM+1
# int_temp = 850.0
# int_press = 1500.0
# int_fH2O = int_melt.calc_fH2O(int_temp, int_press)
# int_obj = Model(temp=int_temp, press=int_press, logfO2=int_logfO2, H2O_param=int_fH2O, 
# 				H2O_param_type="fugacity", SC_param=int_XCO2, SC_param_type="CO2molfrac")
# int_eq_highP = int_obj.eq_fluid(return_as='molfrac')
# #respecieate at 1 bar
# int_eq = int_obj.respeciate(int_eq_highP, newP=1, input_type='molfrac')


#DEEP
deep_melt = SilicateMelt(deep_melt)
deep_XCO2 = 0.93
deep_logfO2 = -7.63
deep_temp = 1100.0
deep_press = 4445.0
deep_fH2O = 333.9
deep_obj = Model(temp=deep_temp, press=deep_press, logfO2=deep_logfO2, H2O_param=deep_fH2O, 
				H2O_param_type="fugacity", SC_param=deep_XCO2, SC_param_type="CO2molfrac")
deep_eq_highP = deep_obj.eq_fluid(return_as='molfrac')
#respecieate at 1 bar
deep_eq = deep_obj.respeciate(deep_eq_highP, newP=1, input_type='molfrac')





