from DegassThermoModel import *

my_fluid = {"H2O": 1.498, "CO2":0.5537, "S":0.2166} #DVDP3-295 average delta MI
my_fluid = normalize(my_fluid)

my_fluid = MagmaticFluid(my_fluid, 5000, 1100, -7.63)




#######silicate melt comp
#Augustine shallow melt comp
my_melt = {"SiO2":75.86, 
		 "TiO2":0.43,
		 "Al2O3":13.10, 
		 "FeO":1.52, 
		 "Fe2O3":0.0, 
		 "MgO":0.25, 
		 "MnO":0.05, 
		 "CaO":1.71, 
		 "Na2O":4.28, 
		 "K2O":2.72, 
		 "P2O5":0.07, 
		 # "Cr2O3":comp["Cr2O3"], 
		 # "NiO":comp["NiO"], 
		 # "F":comp["F"], 
		 # "Cl":comp["Cl"], 
		 "H2O":0.01, 
		 "CO2":0.003, 
		 "S":0.01}

my_melt = SilicateMelt(my_melt)

#Erebus EA1-a for testing
my_melt = {"SiO2":54.99, 
		 "TiO2":1.015,
		 "Al2O3":20.049, 
		 "FeO":4.87, 
		 "Fe2O3":0.0, 
		 "MgO":0.846, 
		 "MnO":0.269, 
		 "CaO":1.881, 
		 "Na2O":9.099, 
		 "K2O":6.092, 
		 "P2O5":0.329, 
		 "H2O":0.031}

from DegassThermoModel import *

#Erebus DVDP3-295a for testing
my_melt = {"SiO2":41.9445, 
		 "TiO2":4.115,
		 "Al2O3":15.166, 
		 "FeO":10.190, 
		 "Fe2O3":0.0, 
		 "MgO":5.992, 
		 "MnO":0.151, 
		 "CaO":12.345, 
		 "Na2O":4.215, 
		 "K2O":1.807, 
		 "P2O5":0.923,
		 "H2O":1.66}

my_melt = SilicateMelt(my_melt)
fH2O = my_melt.calc_fH2O(1100.0, 5000.0)
XCO2 = 0.917
model_obj = Model(temp=1100.0, press=4920.0, logfO2=-7.63, H2O_param=399.16, H2O_param_type="fugacity", SC_param=XCO2, SC_param_type="CO2molfrac")
model_eq = model_obj.eq_fluid()

#Test depressurization
model_lowP = model_obj.respeciate(model_eq, newP=1, input_type='molfrac')




