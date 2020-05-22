#Python library written by Kayla Iacovino after Iacovino (2015)
#Defines methods for the thermodynamic model

# VERSION 0.1 - January 2020

from math import *
from scipy.optimize import fsolve
from scipy.optimize import least_squares
from scipy.optimize import root
import numpy as np
from sympy import *

#TODO - better to return dicts instead of tuples
#TODO - write exceptions
#TODO - test ALL methods and calcs!

##-----------Define Global Variables------------##
#Molecular weights of components in grams per mol
oxideMass = {'SiO2':28.085+32.0,'MgO':24.305+16.0,'FeO':55.845+16.0,'CaO':40.078+16.0,'Al2O3':2.0*26.982+16.0*3.0,
			 'Na2O':22.99*2.0+16.0, 'K2O':39.098*2.0+16.0,'MnO':54.938+16.0,'TiO2':47.867+32.0,
			 'P2O5':2.0*30.974+5.0*16.0,'Cr2O3':51.996*2.0+3.0*16.0, 'NiO':58.693+16.0,'F':18.998,
			 'Cl':35.45,'Fe2O3':55.845*2.0+16.0*3.0,'CO': 28.01,'CO2':44.01,'H2': 2.01588,'H2O':18.02,
			 'H2S':34.1, 'O':15.999, 'O2':15.999*2.0, 'S':32.065, 'Stot': 32.065, 'S2':32.065*2.0, 'SO2':32.065+32.0}

fluid_species_names = ['CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2', 'S2', 'SO2']

#Critical parameters cP, cT, o for relevant species
#dict of dicts
critical_params = {'CO':{	"cT": 	133.15,
							"cP": 	34.9571,
							"o": 	0.049
						},
				   'CO2':{	"cT": 	304.15,
							"cP": 	73.8659,
							"o": 	0.225
						},
				   'H2':{	"cT": 	33.25,
							"cP": 	12.9696,
							"o": 	-0.218
						},
				   'H2O':{	"cT": 	647.25,
							"cP": 	221.1925,
							"o": 	0.334
						},
				   'H2S':{	"cT": 	373.55,
							"cP": 	90.0779,
							"o": 	0.081
						},
				   'O2':{	"cT": 	154.75,
							"cP": 	50.7638,
							"o": 	0.021
						},
				   'S2':{	"cT": 	208.15,
							"cP": 	72.954,
							"o": 	0.0 #need omega value for S2
						},
				   'SO2':{	"cT": 	430.95,
							"cP": 	77.87295,
							"o": 	0.256
						}
					}

def normalize(composition):
	"""Normalizes an input composition to 100%

	Parameters
	----------
	composition: dict
		Dictionary of a composition with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

	Returns
	-------
	dict
		Dictionary of a composition, normalized to 100%, with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
	"""
	return {k: 100.0 * v / sum(composition.values()) for k, v in composition.iteritems()}

def molfrac_to_wtpercent(molfrac):
	"""Converts composition in mole fraction to wt percent.
	NOTE: this also works to conver mole percent to wt percent.

	Parameters
	----------
	molfrac: dict
		Dictionary of a composition in mole fraction with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

	Returns
	-------
	dict
		Dictionary of a composition in wt percent with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
	"""

	MPO_dict = {}
	for key, value in molfrac.iteritems():
		MPO_dict[key] = value * oxideMass[key]

	MPO_sum = sum(MPO_dict.values())

	wtpercent_dict = {}
	for key, value in MPO_dict.iteritems():
		wtpercent_dict[key] = 100.0 * value / MPO_sum

	return wtpercent_dict #TODO works in terminal, need to spot check that values are correct.

def wtpercent_to_molfrac(wtpercent):
	"""Converts composition in wt percent to mole fraction.

	Parameters
	----------
	wtpercent: dict
		Dictionary of a composition in wt percent with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

	Returns
	-------
	dict
		Dictionary of a composition in mole fraction with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
	"""

	MPO_dict = {}
	for key, value in wtpercent.iteritems():
		MPO_dict[key] = value / oxideMass[key]

	MPO_sum = sum(MPO_dict.values())

	molfrac_dict = {}
	for key, value in MPO_dict.iteritems():
		molfrac_dict[key] = value / MPO_sum

	return molfrac_dict #TODO need to spot check that values are correct.

def wtpercent_to_molpercent(wtpercent):
	"""Converts composition in wt percent to mole percent.

	Parameters
	----------
	wtpercent: dict
		Dictionary of a composition in wt percent with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

	Returns
	-------
	dict
		Dictionary of a composition in mole percent with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
	"""

	MPO_dict = {}
	for key, value in wtpercent.iteritems():
		MPO_dict[key] = value / oxideMass[key]

	MPO_sum = sum(MPO_dict.values())

	molfrac_dict = {}
	for key, value in MPO_dict.iteritems():
		molfrac_dict[key] = 100.0 * value / MPO_sum

	return molfrac_dict #TODO need to spot check that values are correct.


def get_iron_speciation(meltcomp, logfO2):
	#TODO write this... should it go in a theoretical SilicateMelt Class?

	"""Speciates FeO and Fe2O3 at a given fO2 after Kress and Carmichael

	Parameters
	----------
	meltcomp: dict
		Dictionary of silicate melt composition in wt percent with possible keys:
		SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
		H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
	logfO2: float
		Log of the oxygen fugacity of the fluid. For example, QFM at 900 degrees C 
	    and 1 bar is input as -12.657235.

	Returns
	-------
	dict
		Dictionary with original melt composition for non-iron species and speciated FeO and
		Fe2O3 for iron species.
	"""

	#TODO CURRENTLY DOES NOTHING - THIS IS A PLACEHOLDER

def get_fO2_buffer(logfO2, temp, buffer='buffer', press='default'):
	#TODO write this... prob already written in another script, just move it here.

	"""Returns oxygen fugacity in terms of the difference from common buffers

	Parameters
	----------
	logfO2: float
		Log of the oxygen fugacity of the fluid. For example, QFM at 900 degrees C 
	    and 1 bar is input as -12.657235.
	temp: float
		Temperature of the melt, in degrees C.
	buffer: str
		Required. User uses this arg to set which buffer to calculate. Options are: QFM, NNO,
		IW, ...#TODO add more buffers
	press: float
		Optional. Pressure of the melt, in bars. Default is 1 bar.

	Returns
	-------
	float
		Value of log units away from chosen buffer. For example:
		get_fO2_buffer(-13.657235, 900, buffer='QFM') will return -1
	"""

	if isinstance(press,str):
		press = 1.0

	if buffer == 'QFM':
		pass

	if buffer == 'NNO':
		pass

	if buffer == 'IW':
		pass

    #TODO CURRENTLY DOES NOTHING - THIS IS A PLACEHOLDER

def calc_logfO2_from_buffer(temp, buffer, delta, press='default'):
	"""Returns logfO2 value given value relative to common solid buffers

	Parameters
	----------
	temp: float
		Temperature in degrees C.
	buffer: str
		Name of solid buffer referenced. Possible strings are: QFM... will be
		adding more soon! #TODO add more buffers!
	delta: float
		Number of log units away from given buffer. Example: for QFM+2, one would enter
		buffer='QFM' and delta=2
	press: float
		Optional. Pressure in bars. If no value is passed, default value of
		1.0 bars will be used

	Returns
	-------
	float
		Value of logfO2
	"""

	tempK = temp + 273.15

	if isinstance(press,str):
		press = 1.0

	if buffer == 'QFM':
		a = -25096.3
		b = 8.735
		c = 0.110
		pass

	if buffer == 'NNO':
		#TODO
		pass

	if buffer == 'IW':
		#TODO
		pass

	return a/tempK + b + c * ((press - 1.0) / tempK) + delta


def Degas_MI(vols_a, vols_b, F='default'): 
    """Composition of fluids from degassing

    Calculates the composition of a fluid phase produced by degassing as magma moves 
    from a deeper storage region (magma a) to a more shallow storage region (magma b).

    Parameters
    ----------
    vols_a:	tuple
        Tuple of volatile concentration in wt percent in magma a (H2O, CO2, S)
    vols_b:	tuple
        Tuple of volatile concentration in wt percent in magma b (H2O, CO2, S)
    F:	float
    	Optional. If no value is passed, F is set to 1.0.
    	Value of F, or 100 minus the amount of crystallization differentiation required 
    	to produce daughter magma b from parent magma a. This is equivalent to the 
    	percentage residual melt that magma b represents.
        
    Returns
    -------
    dictionary
    	Bulk fluid composition in wt percent with keys "H2O", "CO2", "S".
    	NOTE: Fluid is not speciated as fO2 is not considerred."""

    if isinstance(F,str):
    	F = 1.0

    delta_H2O = (vols_a[0] / F) - vols_b[0]
    delta_CO2 = (vols_a[1] / F) - vols_b[1]
    delta_S = (vols_a[2] / F) - vols_b[2]

    sum_vols = delta_H2O + delta_CO2 + delta_S

    norm_H2O = 100.0 * delta_H2O / sum_vols
    norm_CO2 = 100.0 * delta_CO2 / sum_vols
    norm_S = 100.0 * delta_S / sum_vols

    return {"H2O":norm_H2O, "CO2":norm_CO2, "S":norm_S} #TODO spot check that this returns correct values


def calc_Ks(temp, species='all', return_as='standard'):
		"""Returns log of the equilibrium constants from lookup table.

		Parameters
		----------
		temp: float
			Temperature in degrees C.

		species: str
			Choose which species to calculate.
			Options are: CO2, H2O, H2S, SO2, or all.
			If all is passed, a dictionary of values is returned.

		return_as: str
			Default value is 'standard', which returns equilibrium constant K.
			Optional value is 'log', which returns the logK.

		Returns
		-------
		float or dict
			Equilibrium constant for chosen species. If one species is passed,
			returns float. If "all" is passed, returns dict with keys CO2, H2O, H2S,
		"""

		tempK = temp + 273.15

		# CO + 1/2O2 = CO2
		"""From Wagman et al (1945). Text table was fit with sixth-order polynomial to
		arrive at the following coefficients."""
		CO2_logK =		 (9.11899*10.0**(-17.0)	* 	tempK**6.0
					+	-5.65762*10.0**(-13.0) 	* 	tempK**5.0
					+	 1.44706*10.0**(-9.0) 	*	tempK**4.0
					+	-1.96925*10.0**(-6.0) 	*	tempK**3.0
					+	 1.53277*10.0**(-3.0) 	*	tempK**2.0
					+	-6.79138*10.0**(-1.0) 	*	tempK
					+	 1.53288*10.0**(2.0))

		# H2 + 1/2O2 = H2O
		"""From Robie and Hemingway (1995). Text table was fit with sixth-order polynomial to
		arrive at the following coefficients."""
		H2O_logK = 		(3.3426*10.0**(-17.0) 	* 	tempK**6.0
					+	-2.40813*10.0**(-13.0)	*	tempK**5.0
					+	7.10612*10.0**(-10.0)	*	tempK**4.0
					+	-1.10671*10.0**(-6.0)	*	tempK**3.0
					+	9.76038*10.0**(-4.0)	*	tempK**2.0
					+	-4.84445*10.0**(-1.0)	*	tempK
					+	1.21953*10.0**(2.0))

		# H2 + 1/2S2 = H2S
		"""From Robie and Hemingway (1995). Text table was fit with sixth-order polynomial to
		arrive at the following coefficients."""
		H2S_logK = 		(1.882737*10.0**(-18.0) * tempK**6.0
					+	-1.779266*10.0**(-14.0) * tempK**5.0
					+	 6.319209*10.0**(-11.0) * tempK**4.0
					+	-1.092048*10.0**(-7.0) * tempK**3.0
					+	 9.824774*10.0**(-5.0) * tempK**2.0
					+	-4.805344*10.0**(-2.0) * tempK
					+	 1.389857*10.0**(1.0))

		# S2 + 1/2O2 = SO2
		"""From Robie and Hemingway (1995). Text table was fit with sixth-order polynomial to
		arrive at the following coefficients."""
		SO2_logK =		 (4.01465*10.0**(-17.0) * tempK**6.0
					+	-2.93845*10.0**(-13.0) 	* tempK**5.0
					+	 8.78818*10.0**(-10.0) 	* tempK**4.0
					+	-1.38066*10.0**(-6.0) 	* tempK**3.0
					+	 1.21978*10.0**(-3.0) 	* tempK**2.0
					+	-6.03476*10.0**(-1.0) 	* tempK
					+	 1.54350*10.0**(2.0))

		if return_as == 'standard':
			if species == 'CO2':
				return 10.0**CO2_logK
			if species == 'H2O':
				return 10.0**H2O_logK
			if species == 'H2S':
				return 10.0**H2S_logK
			if species == 'SO2':
				return 10.0**SO2_logK
			if species == 'all':
				return {'CO2': 10.0**CO2_logK, 'H2O': 10.0**H2O_logK, 'H2S': 10.0**H2S_logK, 'SO2': 10.0**SO2_logK}
		if return_as == 'log':
			if species == 'CO2':
				return CO2_logK
			if species == 'H2O':
				return H2O_logK
			if species == 'H2S':
				return H2S_logK
			if species == 'SO2':
				return SO2_logK
			if species == 'all':
				return {'CO2': CO2_logK, 'H2O': H2O_logK, 'H2S': H2S_logK, 'SO2': SO2_logK}
		#TODO raise exception if something other than 'standard' or 'log' is passed?

def calc_gammas(temp, press, species='all'):
		"""Returns fugacity coefficients calculated using the Redlich Kwong Equation of State.
		Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 
		30 October 2003.

		Parameters
		----------
		temp: float
			Temperature in degrees C.

		press: float
			Pressure in bars.

		species: str
			Choose which species to calculate.
			Options are: 'CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2', 'S2', 'SO2', or all.
			If all is passed, a dictionary of values is returned. Default value is 'all'.

		Returns
		-------
		float or dict
			Fugacity coefficient for passed species.
			If single species is passed, float.
			If "all" is passed, dictionary with keys 'CO', 'CO2', 'H2', 'H2O', 'H2S', 'O2', 'S2', 'SO2'
		"""

		tempK = temp + 273.15
		R = 8.3145

		gamma_dict = {}

		for species in fluid_species_names:
			#Calculate a and b parameters (depend only on critical parameters)...
			a = 0.42748 * R**2.0 * critical_params[species]["cT"]**(2.5) / (critical_params[species]["cP"] * 10.0**5)
			b = 0.08664 * R * critical_params[species]["cT"] / (critical_params[species]["cP"] * 10.0**5)
			kappa = 0.0

			#Calculate coefficients in the cubic equation of state...
			#coeffs: (C0, C1, C2, A, B)
			A = a * press * 10.0**5 / (sqrt(tempK) * (R * tempK)**2.0)
			B = b * press * 10.0**5 / (R * tempK)
			C2 = -1.0
			C1 = A - B - B * B
			C0 = -A * B

			#Solve the cubic equation for Z0 - Z2, D...
			Q1 = C2 * C1 / 6.0 - C0 / 2.0 - C2**3.0 / 27.0
			P1 = C2**2.0 / 9.0 - C1 / 3.0
			D = Q1**2.0 - P1**3.0

			if D >= 0:
				kOneThird = 1.0 / 3.0

				absQ1PSqrtD = fabs(Q1 + sqrt(D))
				temp1 = absQ1PSqrtD**kOneThird
				temp1 *= (Q1 + sqrt(D)) / absQ1PSqrtD

				absQ1MSqrtD = fabs(Q1 - sqrt(D))
				temp2 = absQ1MSqrtD**kOneThird
				temp2 *= (Q1 - sqrt(D)) / absQ1MSqrtD

				Z0 = temp1 + temp2 - C2 / 3.0
			else:
				temp1 = Q1**2.0 / (P1**3.0)
				temp2 = sqrt(1.0 - temp1) / sqrt(temp1)
				temp2 *= Q1 / fabs(Q1)

				gamma = atan(temp2)

				if gamma < 0:
					gamma = gamma + pi 

				Z0 = 2.0 * sqrt(P1) * cos(gamma/3.0) - C2 / 3.0
				Z1 = 2.0 * sqrt(P1) * cos((gamma + 2.0 * pi) / 3.0) - C2/3.0
				Z2 = 2.0 * sqrt(P1) * cos((gamma + 4.0 * pi) / 3.0) - C2/3.0

				if Z0 < Z1:
					temp0 = Z0
					Z0 = Z1
					Z1 = temp0

				if Z1 < Z2:
					temp0 = Z1
					Z1 = Z2
					Z2 = temp0

				if Z0 < Z1:
					temp0 = Z0
					Z0 = Z1
					Z1 = temp0

			#Determine the fugacity coefficient of first root and departure functions...
			#calcdepfns(coeffs[3], 	coeffs[4], 	paramsab[0], 	Z[0])
			#calcdepfns(A, 			B, 			kappa, 			Z)

			#Calculate Departure Functions
			gamma = exp(Z0 - 1.0 - log(Z0-B) - A * log(1.0+B/Z0)/B)
			Hdep = R * tempK * (Z0 - 1.0 - 1.5*A*log(1.0+B/Z0)/B)
			Sdep = R * (log(Z0-B) - 0.5*A*log(1.0+B/Z0)/B)

			gamma_dict[species] = gamma

			#gamma_tuple = tuple(gamma_dict.values())

		if species == 'CO':
			return gamma_dict["CO"]
		if species == 'CO2':
			return gamma_dict["CO2"]
		if species == 'H2':
			return gamma_dict["H2"]
		if species == 'H2O':
			return gamma_dict["H2O"]
		if species == 'H2S':
			return gamma_dict["H2S"]
		if species == 'O2':
			return gamma_dict["O2"]
		if species == 'S2':
			return gamma_dict["S2"]
		if species == 'SO2':
			return gamma_dict
		if species == 'all':
			return gamma_dict


#TODO write code to take in mol% and mol frac. Currently only takes in wt%
class SilicateMelt(object):
	"""A silicate melt with the following properties:

	Attributes
	----------
		comp: dict
			Dict of melt composition in weight percent, mole percent, or mole fraction with possible
			keys SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
			H2O, CO2, S. If any key is not passed, it will inherit the default value of 0.0.
	"""

	def __init__(self, comp):
		"""Return a SilicateMelt object whose parameters are defined here."""
		#set default values for oxides to 0.0
		self.comp = {"SiO2":0.0, 
					 "TiO2":0.0,
					 "Al2O3":0.0, 
					 "FeO":0.0, 
					 "Fe2O3":0.0, 
					 "MgO":0.0, 
					 "MnO":0.0, 
					 "CaO":0.0, 
					 "Na2O":0.0, 
					 "K2O":0.0, 
					 "P2O5":0.0, 
					 "Cr2O3":0.0, 
					 "NiO":0.0, 
					 "F":0.0, 
					 "Cl":0.0, 
					 "H2O":0.0, 
					 "CO2":0.0, 
					 "S":0.0}

		#for loop that only appends an oxide value if the user passes one
		for key, value in comp.iteritems():
			if key in self.comp:
				self.comp[key] = value
			else:
				pass

	#TODO don't let the user go above 3 kbar!!! GMM: "You have been warned."
	#TODO - testing done for EA1 erebus comp
	def calc_fH2O(self, temp, press, model="Moore1998"):
		"""Returns the fugacity of H2O in the fluid given melt composition, including dissolved H2O

		Parameters
		----------
		self, inherited from Class
		temp: float
			temperature of the melt in degrees C
		press: float
			pressure of the melt in bars
		model: str
			Currently, this argument does nothing. It is here to notify the user they are using the 
			fugacity model of Moore et al. (1998). In future, this method may includ the option to
			use other models.

		Returns
		-------
		float
			H2O fugacity in the fluid.
		"""

		#define constants
		a = 2565.0
		a_err = 362.0 #standard error
		b_Al2O3 = -1.997
		b_Al2O3_err = 0.706
		b_FeOt = -0.9275
		b_FeOt_err = 0.394
		b_Na2O = 2.736
		b_Na2O_err = 0.871
		c = 1.171
		c_err = 0.069
		d = -14.21
		d_err = 0.54
		XH2O_err = 0.148 #reported 1 sigma error

		melt_comp_wt = self.comp
		tempk = temp + 273.15
		melt_comp_X = wtpercent_to_molfrac(melt_comp_wt)

		#get FeOt from iron input as FeO and Fe2O3
		FeOt = melt_comp_X["FeO"] + 0.8998*melt_comp_X["Fe2O3"]

		#TODO assert melt_comp_X["H2O"] > 0

		#note math.log with one argument passed returns natural log
		ln_fH2O = (((2.0 * log(melt_comp_X["H2O"])) - (a/tempk) - 
					((b_Al2O3 * melt_comp_X["Al2O3"] * (press/tempk)) + 
					 (b_FeOt * FeOt * (press/tempk)) + 
					 (b_Na2O * melt_comp_X["Na2O"] * (press/tempk))) - 
					 d)/c)

		return exp(ln_fH2O)


class MagmaticFluid(object):
	"""A fluid derived from a magma. A fluid must have the following properties:

	Attributes
	----------
		fluid_comp:	dict
	        Dict of volatile concentration in wt percent, mol percent, or mol fraction 
	        in fluid with keys H2O, CO2, S.
	    press:	float
	        Pressure at which to speciate fluid, in bars.
	    temp:	float
	    	Temperature at which to speciate fluid, in degrees C.
	    logfO2:	float
	    	Log of the oxygen fugacity of the fluid. For example, QFM at 900 degrees C 
	    	and 1 bar is input as -12.657235.
	    input_type: str
	    	String defining whether fluid_comp is input as wt percent ("wtpercent"), 
	    	mole percent ("molpercent"), or mole fraction ("molfrac"). Default is "wtpercent".
	"""

	def __init__(self, fluid_comp, press, temp, fO2_buffer, fO2_delta, input_type='wtpercent'):
		"""Return a MagmaticFluid object whose parameters are defined here."""
		self.fluid_comp = {"H2O":fluid_comp["H2O"], "CO2":fluid_comp["CO2"], "S":fluid_comp["S"]}
		self.press = press
		self.temp = temp
		self.fO2_buffer = fO2_buffer
		self.fO2_delta = fO2_delta
		self.input_type = input_type

		self.tempK = self.temp + 273.15
		self.logfO2 = calc_logfO2_from_buffer(press=press, temp=temp, buffer=fO2_buffer, delta=fO2_delta)
		self.fO2 = 10.0**self.logfO2

		#TODO test all input types produce correct values
		if input_type == "wtpercent":
			self.fluid_comp_wt = self.fluid_comp
			self.fluid_comp_molpercent = wtpercent_to_molpercent(self.fluid_comp)
			self.fluid_comp_molfrac = wtpercent_to_molfrac(self.fluid_comp)

		if input_type == "molpercent":
			self.fluid_comp_wt = molfrac_to_wtpercent(self.fluid_comp)
			self.fluid_comp_molpercent = self.fluid_comp
			self.fluid_comp_molfrac = {k: v / 100.0 for k, v in self.fluid_comp.iteritems()}

		if input_type == "molfrac":
			self.flud_comp_wt = molfrac_to_wtpercent(self.fluid_comp)
			self.flud_comp_molpercent = {k: v * 100.0 for k, v in self.fluid_comp.iteritems()}
			self.fliud_comp_molfrac = self.fluid_comp

	def fugacities(self, gammas='calculate', K_vals='calculate'):
		"""Returns fugacity values for each species.

		Parameters
    	----------
    	self, inherited from Class
		Notably, uses fluid_comp, temp, press, fO2

    	gammas: dict
    		Optional. Fugacity coefficients.
    		If gamma values are not passed, they will be calculated.
    		Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    	K_vals: dict
    		Optional. Equilibrium constants of reaction.
    		If K values are not passed, they will be calculated.
    		Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}

		Returns
		-------
		dict
			Fugacities of all species

		"""

		if isinstance(gammas,str):
			gammas = calc_gammas(self.temp, self.press, species="all")
			self.gammas = gammas

		if isinstance(K_vals,str):
			K_vals = calc_Ks(self.temp, species="all")

		XH2Otot = self.fluid_comp_molfrac["H2O"]
		XCO2tot = self.fluid_comp_molfrac["CO2"]
		XStot = self.fluid_comp_molfrac["S"]

		XHtot = XH2Otot * (0.6666)
		XStot = XStot
		XCtot = XCO2tot * (0.3333)
		XOtot = XH2Otot * (0.3333) + XCO2tot * (0.6666)

		B = gammas['H2']
		P = self.press
		C = K_vals['H2O']
		D = self.fO2
		sD = sqrt(D)
		E = gammas['H2O']
		F = K_vals['H2S']
		G = gammas['H2S']
		J = gammas['S2']
		K = K_vals['SO2']
		L = gammas['SO2']
		M = K_vals['CO2']
		N = gammas['CO2']
		Q = gammas['CO']

		#FIRST calculate fH2 and fS2 using fsolve, two eqns; two unknowns (eqn 9 in Iacovino, 2015)
		def equations(p):
			fH2, fS2 = p
			return 	(
						( (fH2/(B*P)) 	+ ((Rational(2.0) * C * fH2 * sD)/(Rational(3.0) * E * P))		+ ((Rational(2.0) * F * fH2 * sqrt(abs(fS2)))/(Rational(3.0) * G * P)) 	- XHtot), 
						( (fS2/(J * P)) 	+ ((F * fH2 * sqrt(abs(fS2)))/(Rational(3.0) * G * P))	+ ((K * D * sqrt(abs(fS2)))/(Rational(3.0) * L * P))		- XStot)
					)

		fH2_a, fS2_a = fsolve(equations, (1, 1))
		fH2 = abs(fH2_a)
		fS2 = abs(fS2_a)
		print fS2

		#SECOND calculate fCO (eqn 10 in Iacovino, 2015) using sympy
		fCO = symbols('fCO') #for sympy

		equation = (((M * fCO * sD)/(Rational(3.0) * N * P)) + ((fCO)/(Rational(2.0) * Q * P))	- XCtot)
		fCO = solve(equation, fCO)[0] #newly implemented sympy way

		#THIRD calculate fCO2 using calc'd fCO and known fO2 value
		fCO2 = M * fCO * sD

		#FOURTH calcualte fSO2 using calc'd fS2 and known fO2 value
		fSO2 = K * sqrt(fS2) * D

		#FIFTH calculate fH2S using calc'd fH2 and fS2 values
		fH2S = F * fH2 * sqrt(fS2)

		#SIXTH calculate fH2O using calc'd fH2 and knwn fO2 value
		fH2O = C * sD * fH2

		#TODO raise exception if a fugacity is negative or zero.

		return {'CO': fCO, 'CO2': fCO2, 'H2': fH2, 'H2O': fH2O, 'H2S': fH2S, 'O2': D, 'S2': fS2, 'SO2': fSO2}


	def speciate(self, return_as='wtpercent', gammas='calculate', K_vals='calculate', fugacities='calculate'):
		"""Speciates a fluid given bulk composition.

		Parameters
    	----------
    	self, inherited from Class
		Notably, passes fluid_comp, temp, press, fO2 to other methods

    	return_as: str
    		Default value is 'wtpercent', which returns composition in wt percent.
			Optional values are 'molpercent', which returns composition as mol percent and
			'molfrac', which returns composition as mole fraction. 
    	gammas: dict
    		Optional. Fugacity coefficients.
    		If gamma values are not passed, they will be calculated.
    		Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    	K_vals: dict
    		Optional. Equilibrium constants of reaction.
    		If K values are not passed, they will be calculated.
    		Format: {'H2O': value, 'CO2': value, 'SO2': value, 'H2S': value}
    	fugacities: dict
    		Optional. Fugacity values for each species.
    		If fugacity values are not passed, they will be calculated.
    		Format: {'CO': value, 'CO2': value, 'H2': value, 'H2O': value, 'H2S': value, 
    		'O2': value, 'S2': value, 'SO2': value}

		Returns
		-------
		dict
			Speciated fluid composition in wt percent, mol percent, or mol fraction.
		"""

		press = self.press 

		if isinstance(gammas,str):
			gammas = calc_gammas(self.temp, self.press, species="all")

		if isinstance(K_vals,str):
			K_vals = calc_Ks(self.temp, species="all")

		if isinstance(fugacities,str):
			fugacities = self.fugacities()

		X_dict = {}
		for species in fluid_species_names:
			X = fugacities[species] / (gammas[species] * press)
			X_dict[species] = X

		#TODO test all these types...
		if return_as == 'wtpercent':
			return molfrac_to_wtpercent(X_dict)

		if return_as == 'molpercent':
			molpercent_dict = {}
			for key, value in X_dict.iteritems():
				molpercent_dict[key] = value*100.0
			return molpercent_dict

		if return_as == 'molfrac':
			return {key: value/sum(X_dict.values()) for key,value in X_dict.iteritems()}

		#TODO raise exception if something other than 'wtpercent', 'molpercent' or 'molfrac' is passed?


class Model(object):
	"""An object with arguments describing all needed thermodynamic parameters:

	Attributes
	----------
		press:	float
	        Pressure in bars.
	    temp:	float
	    	Temperature in degrees C.
	    logfO2:	float
	    	Log of the oxygen fugacity. For example, QFM at 900 degrees C 
	    	and 1 bar is input as -12.657235.
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

	def __init__(self, press, temp, fO2_buffer, fO2_delta, H2O_param, H2O_param_type, SC_param, SC_param_type):
		"""Return a Model object whose parameters are defined here."""
		self.press = press
		self.temp = temp
		self.fO2_buffer = fO2_buffer
		self.fO2_delta = fO2_delta
		self.H2O_param = H2O_param
		self.H2O_param_type = H2O_param_type
		self.SC_param = SC_param
		self.SC_param_type = SC_param_type
		self.logfO2 = calc_logfO2_from_buffer(press=press, temp=temp, buffer=fO2_buffer, delta=fO2_delta)

		self.tempK = self.temp + 273.15
		self.fO2 = 10.0 **self.logfO2

		#calculate gammas and Ks
		self.gammas = calc_gammas(self.temp, self.press, species="all")
		self.Ks = calc_Ks(self.temp, species="all")

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
			self.fCO2 = self.press * self.gammas["CO2"] * self.XCO2
			self.PCO2 = self.fCO2 / self.gammas["CO2"]

	#Don't need these because, since they are calc'd in init, you can just call object.gammas or object.Ks
	# def get_gammas(self):
	# 	"""Returns gamma values for all species as a dict"""
	# 	gammas = self.gammas 
	# 	return gammas

	# def get_Ks(self):
	# 	"""Returns equilibrium constant values for all species as dict"""
	# 	Ks = self.Ks
	# 	return Ks


	def eq_fluid(self, return_as='molfrac'):
		"""Return the composition of an equilibrium fluid.

		Parameters
		----------
		self, inherited from Class

		return_as: str
			Optional. Returns composition of fluid in form passed by the user. Default
			value is 'molfrac', which returns as mole fraction. Other options are 'wtpercent',
			'molpercent', 'ppress' (partial pressures), and 'fugacities'.

		#TODO add argument "return_as" and allow return as mole frac, mol%, or wt%

		Returns
		-------
		dict
			Dictionary with composition of equilibrium fluid as:
			{'CO': value, 'CO2': value, 'H2': value, 'H2O': value, 'H2S': value, 
    		'O2': value, 'S2': value, 'SO2': value}
		"""

		fO2 = self.fO2
		logfO2 = self.logfO2
		fH2O = self.fH2O
		fCO2 = self.fCO2
		press = self.press
		temp = self.temp
		Ks = self.Ks
		gammas = self.gammas

		PH2O = self.PH2O
		PCO2 = self.PCO2 #TODO edit this to be in an if statement when you write options for passing other SC_param types
		PO2 = fO2 / gammas["O2"]
		#All P's = CO2, H2O, O2 | All f's = CO2, H2O, O2

		#Calculate fH2 from fH2O and fO2, equation 2 Iacovino (2015) EPSL
		fH2 = fH2O / (Ks["H2O"] * sqrt(fO2))
		PH2 = fH2 / gammas["H2"]
		#All P's  = CO2, H2, H2O, O2 | All f's = CO2, H2, H2O, O2

		#calculate fCO from fCO2 and fO2
		fCO = self.fCO2 / (Ks["CO2"] * sqrt(fO2))
		PCO = fCO / gammas["CO"]
		#All P's  = CO, CO2, H2, H2O, O2 | All f's = CO, CO2, H2, H2O, O2

		#calculate PStot by difference
		PStot = press - (PCO + PCO2 + PH2 + PH2O + PO2)

		#Calculate fS2 with sympy, equation 7
		fS2 = symbols('fS2') #for sympy

		gS2 = gammas["S2"]
		KSO2 = Ks["SO2"]
		gSO2 = gammas["SO2"]
		KH2S = Ks["H2S"]
		gH2S = gammas["H2S"]
		
		equation = ((fS2 / gS2) + (KSO2 * sqrt(fS2) * fO2 / gSO2) + (KH2S * fH2 * sqrt(fS2) / gH2S) - PStot)

		# G = gammas["S2"]
		# K = Ks["SO2"]
		# F = fO2
		# H = gammas["SO2"]
		# J = Ks["H2S"]
		# X = fH2
		# L = gammas["H2S"]
		# fS2 = (0.5) * (-sqrt(G**2*(-(F**2*G*K**2 / H**2) - (2.0*F*G*J*K*X / (H * L)) - (G*J**2*X**2 / L**2) - 2.0*PStot)**2 
		# 					- 4.0*G**2*PStot**2) - G*(-F**2*G*K**2/H**2 - 2.0*F*G*J*K*X/(H*L) - G*J**2*X**2/(L**2) - 2.0*PStot))

		# fS2_array = fsolve(equation, 0.0001)
		#fS2_array = least_squares(equation, 0.0001, bounds=(0, np.inf))
		#fS2 = fS2_array.x[0]
		fS2 = solve(equation, fS2)[0] #newly implemented sympy way
		
		PS2 = fS2 / gammas["S2"]
		#All P's  = CO, CO2, H2, H2O, O2, S2 | All f's = CO, CO2, H2, H2O, O2, S2

		fSO2 = Ks["SO2"] * sqrt(fS2) * fO2
		PSO2 = fSO2 / gammas["SO2"]
		#All P's  = CO, CO2, H2, H2O, O2, S2, SO2 | All f's = CO, CO2, H2, H2O, O2, S2, SO2

		fH2S = Ks["H2S"] * fH2 * sqrt(fS2)
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

		fugacities = { 	  "CO": fCO,
						  "CO2": fCO2,
						  "H2": fH2,
						  "H2O": fH2O,
						  "H2S": fH2S,
						  "O2": fO2,
						  "S2": fS2,
						  "SO2": fSO2 }
		self.fugacities = fugacities


		mol_fractions = {}
		for key, value in partial_pressures.iteritems():
			mol_fractions[key] = fugacities[key] / (gammas[key] * self.recalc_pressure)

		#For debugging:
		P_err = self.recalc_pressure - press

		if P_err != 0.0:
			print "Calculated pressure error is " + str(P_err) + " bars"

		if return_as == 'molfrac':
			return mol_fractions

		if return_as == 'wtpercent':
			return molfrac_to_wtpercent(mol_fractions)

		if return_as == 'molpercent':
			molpercent_dict = {}
			for key, value in mol_fractions.iteritems():
				molpercent_dict[key] = value*100.0
			return molpercent_dict

		if return_as == 'fugacities':
			return fugacities

		if return_as == 'ppress':
			return partial_pressures

	def respeciate(self, fluid_comp, newP, input_type='wtpercent'):
		"""Takes in a fluid of given composition and returns that fluid respeciated at a given pressure.

		Parameters
		----------
		self, inherited from Class

		fluid_comp: dict
			Dictionary of fluid composition with possible keys:
			CO, CO2, H2, H2O, H2S, O2, S2, SO2

		newP: float
			Pressure at which to respeciate the fluid, in bars.

		input_type: str
			String defining whether fluid_comp is input as wt percent ("wtpercent"), 
		    mole percent ("molpercent"), or mole fraction ("molfrac"). Default is "wtpercent".

		Returns
		-------
		dict
			Dictionary of fluid composition in mole fraction after respeciation.
			#TODO give other return type options.
		"""

		#Get some needed Class variables
		fugacities = self.fugacities
		gammas = calc_gammas(press=newP, temp=self.temp)
		Ks = self.Ks

		#Recalculate fO2 at new pressure
		logfO2_res = calc_logfO2_from_buffer(press=newP, temp=self.temp, buffer=self.fO2_buffer, delta=self.fO2_delta)
		self.fO2_res = 10**(logfO2_res)

		#Take any input type and convert to mole fraction
		if input_type == "wtpercent":
			fluid_comp = wtpercent_to_molfrac(fluid_comp)

		if input_type == "molpercent":
			fluid_comp = {k: v / 100.0 for k, v in fluid_comp.iteritems()}

		if input_type == "molfrac":
			pass

		#Name some local variables for clarity in below equations. If a particular species is not passed,
		#it will be assigned a value of 0.0	
		if 'CO' in fluid_comp:
			CO = fluid_comp["CO"]
		else:
			CO = 0.0
		if 'CO2' in fluid_comp:
			CO2 = fluid_comp["CO2"]
		else:
			CO2 = 0.0
		if 'H2' in fluid_comp:
			H2 = fluid_comp["H2"]
		else:
			H2 = 0.0
		if 'H2O' in fluid_comp:
			H2O = fluid_comp["H2O"]
		else:
			H2O = 0.0
		if 'H2S' in fluid_comp:
			H2S = fluid_comp["H2S"]
		else:
			H2S = 0.0
		if 'O2' in fluid_comp:
			O2 = fluid_comp["O2"]
		else:
			O2 = 0.0
		if 'S2' in fluid_comp:
			S2 = fluid_comp["S2"]
		else:
			S2 = 0.0
		if 'SO2' in fluid_comp:
			SO2 = fluid_comp["SO2"]
		else:
			SO2 = 0.0

		XHtot = H2 + (0.666666)*H2O + (0.666666)*H2S
		XStot = S2 + (0.333333)*H2S + (0.333333)*SO2
		XCtot = (0.333333)*CO2 + (0.5)*CO
		XOtot = O2 + (0.666666)*CO2 + (0.5)*CO + (0.333333)*H2O + (0.666666)*SO2

		#FIRST calculate fH2 and fS2 using fsolve, two eqns; two unknowns (eqn 9 in Iacovino, 2015)
		#TODO - User might want the fO2 to be set by a buffer. If so, the fO2 used in the equation below must
		#be recalculated at 1 bar!!!

		#TODO - figure out how to do these two equations with two unknowns and set bounds that roots must be >0
		#Used least_squares to do this in a similar implimentation in this very script, but I can't
		#figureo out how to get it to work with two equations instead of just one.
		#THIS IS IMPORTANT since right now it's a total non-mathematical flub.

		def equations(p):
			fH2, fS2 = p
			return 	(
						( (fH2/(self.gammas["H2"]*newP)) 	+ ((Rational(2.0) * self.Ks["H2O"] * fH2 * sqrt(self.fO2_res))/(Rational(3.0) * self.gammas["H2O"] * newP)) + ((Rational(2.0) * self.Ks["H2S"] * fH2 * sqrt(abs(fS2)))/(Rational(3.0) * self.gammas["H2S"] * newP)) - XHtot), 
						( (fS2/(self.gammas["S2"] * newP)) 	+ ((self.Ks["H2S"] * fH2 * sqrt(abs(fS2)))/(Rational(3.0) * self.gammas["H2S"] * newP))	+ ((self.Ks["SO2"] * self.fO2_res * sqrt(abs(fS2)))/(Rational(3.0) * self.gammas["SO2"] * newP))		- XStot)
					)

		fH2_a, fS2_a = fsolve(equations, (1, 1))
		fH2 = abs(fH2_a)
		fS2 = abs(fS2_a)
		#res = least_squares(equations, (0.0001, 0.0001), bounds=((0.0, 0.0), (10.0, 10.0)))
		#fH2 = res.x[0]
		#fS2 = res.x[1]

		#SECOND calculate fCO (eqn 10 in Iacovino, 2015)
		fCO = symbols('fCO') #for sympy
		equation = (((Ks["CO2"] * fCO * sqrt(self.fO2_res))/(3.0 * gammas["CO2"] * newP)) + ((fCO)/(2.0 * gammas["CO"] * newP))	- XCtot)
		fCO = solve(equation, fCO)[0] #newly implemented sympy way

		# def fCO_func(fCO):
		# 	return (((Ks["CO2"] * fCO * sqrt(fugacities["O2"]))/(3 * gammas["CO2"] * newP)) + ((fCO)/(2 * gammas["CO"] * newP))	- XCtot)
		# fCO_array = least_squares(fCO_func, 0.001, bounds=(0, np.inf))
		# fCO = fCO_array.x[0]

		#THIRD calculate fCO2 using calc'd fCO and known fO2 value
		fCO2 = self.Ks["CO2"] * fCO * sqrt(self.fO2_res)

		#FOURTH calcualte fSO2 using calc'd fS2 and known fO2 value
		fSO2 = self.Ks["SO2"] * sqrt(fS2) * self.fO2_res

		#FIFTH calculate fH2S using calc'd fH2 and fS2 values
		fH2S = self.Ks["H2S"] * fH2 * sqrt(fS2)

		#SIXTH calculate fH2O using calc'd fH2 and knwn fO2 value
		fH2O = self.Ks["H2O"] * sqrt(self.fO2_res) * fH2 

		new_fugacities = {"CO": fCO,
							"CO2": fCO2,
							"H2": fH2,
							"H2O": fH2O,
							"H2S": fH2S,
							"O2": self.fO2_res,
							"S2": fS2,
							"SO2": fSO2}
		self.new_fugacities = new_fugacities

		X_dict = {}
		for species in fluid_species_names:
			X = new_fugacities[species] / (gammas[species] * newP)
			X_dict[species] = X

		return {key: value/sum(X_dict.values()) for key,value in X_dict.iteritems()}

#TODO test this entire class and init
#TODO rewrite this to use lists, not dicts, since dicts can be passed in whatever order!
class Match(object):
	"""A Match object that takes in gas compositions for use in the mixing model:

	Attributes
	----------
		sub_gases: dict
			Dictionary of dictionaries of gas compositions
			Example argument to pass: {"Gas1": {'H2O'=val, 'CO2'=val, 'SO2'=val}, "Gas2": {'H2O'=val, 'CO2'=val}}

		surface_gas: dict
			Dictionary of surface gas composition with possible keys H2O, CO2, SO2, H2S, Stot

		input_type: str
			Optional. Default value is 'wtpercent'. String defining the units of input data. 
			Possible strings are 'wtpercent', 'molpercent', and 'molfrac'.

	"""

	def __init__(self, sub_gases, surface_gas, input_type='wtpercent'):
		"""Return a Match object whose parameters are defined here."""
		self.input_type = input_type
		self.surface_gas = surface_gas

		# #set default values for surface gas comp to 0.0
		# self.surface_gas = {"H2O":0.0, 
		# 					 "CO2":0.0,
		# 					 "SO2":0.0, 
		# 					 "H2S":0.0, 
		# 					 "Stot":0.0}

		# #for loop that only appends a gas species value if the user passes one
		# for species, value in surface_gas.iteritems():
		# 	if species in self.surface_gas:
		# 			self.surface_gas[species] = value
		# 	else:
		# 		pass

		#set default species values for subsurface gas comps to 0.0
		#TODO write to take whatever gas species are passed (e.g., optionally include CO, H2S, etc...)
		blank_sub_gases = {}
		gas_name_list = list(sub_gases.keys())
		species_list = ["H2O", "CO2", "SO2"]
		for gas_name in gas_name_list:
			blank_sub_gases[gas_name] = {"H2O":0.0, 
									 "CO2":0.0,
									 "SO2":0.0}

		#for loop that only appends a gas species value to each sub_gas if the user passes one and only if that
		#gas is also in the surface_gas dictionary.
		for gas_name, gas_dict in sub_gases.iteritems():
			gas_dict["SO2"] = gas_dict["SO2"] + gas_dict["H2S"]
			for key in gas_dict.keys():
				if key not in species_list:
					del gas_dict[key]
			for species in species_list:
				if species not in gas_dict:
					gas_dict[species] = 0.0	


		#normalize sub_gases
		self.sub_gases = {gas_name: {species: (100.0 * value/sum(gas_dict.values())) for species, value in gas_dict.iteritems()} for gas_name, gas_dict in sub_gases.iteritems()}

		#TODO test input types
		if input_type == "wtpercent":
			self.sub_gases_wt = self.sub_gases
			self.sub_gases_molpercent = {gas_name: wtpercent_to_molfrac(gas_dict) for gas_name, gas_dict in self.sub_gases_wt.iteritems()}
			self.sub_gases_molfrac = {gas_name: {species: value / 100.0 for species, value in gas_dict.iteritems()} for gas_name, gas_dict in self.sub_gases_molpercent.iteritems()}

			self.surface_gas_wt = self.surface_gas 
			self.surface_gas_molpercent = wtpercent_to_molpercent(self.surface_gas)
			self.surface_gas_molfrac = wtpercent_to_molfrac(self.surface_gas)

		if input_type == "molpercent":
			self.sub_gases_molpercent = self.sub_gases
			self.sub_gases_molfrac = {gas_name: {species: value / 100.0 for species, value in gas_dict.iteritems()} for gas_name, gas_dict in self.sub_gases_molpercent.iteritems()}
			self.sub_gases_wt = {gas_name: molfrac_to_wtpercent(gas_dict) for gas_name, gas_dict in self.sub_gases_molfrac.iteritems()}

			self.surface_gas_molpercent = self.surface_gas
			self.surface_gas_molfrac = {k: v / 100.0 for k, v in self.surface_gas.iteritems()}
			self.surface_gas_wt = molfrac_to_wtpercent(self.surface_gas_molfrac)

		if input_type == "molfrac":
			self.sub_gases_molfrac = self.sub_gases
			self.sub_gases_molpercent = {gas_name: {species: value * 100.0 for species, value in gas_dict.iteritems()} for gas_name, gas_dict in self.sub_gases_molfrac.iteritems()}
			self.sub_gases_wt = {gas_name: molfrac_to_wtpercent(gas_dict) for gas_name, gas_dict in self.sub_gases_molfrac.iteritems()}

			self.surface_gas_molfrac = self.surface_gas
			self.surface_gas_molpercent = {k: v * 100.0 for k, v in self.surface_gas.iteritems()}
			self.surface_gas_wt = molfrac_to_wtpercent(self.surface_gas)

	def matchmodel(self, threshold=0.10):
		"""Runs the simple matching model and returns all possible combinations of sub_gases, scaled each at 0-100%, that could produce
		the surface_gas composition to within a user defined error (called "threshold").

		Parameters
		----------
		self, inherited from Class

		threshold: float
			Optional. Default value is 0.10 (10%). The threshold defines how closely a combination of sub_gases must match the 
			surface_gas, expressed as relative percent of the surface_gas value.


		Returns
		-------
		dict???
			???
		"""
		sub_gases = self.sub_gases_molpercent
		surface_gas = self.surface_gas_molpercent
		self.threshold = threshold

		sub_CO2 = [gas_dict["CO2"] for gas_name, gas_dict in sub_gases.iteritems()]
		sub_H2O = [gas_dict["H2O"] for gas_name, gas_dict in sub_gases.iteritems()]
		#TODO make SO2 actual SO2, not Stot
		sub_SO2 = [gas_dict["SO2"] for gas_name, gas_dict in sub_gases.iteritems()]
		#sub_H2S = [gas_dict["H2S"] for gas_name, gas_dict in sub_gases.iteritems()]
		#TODO write this to be able to take in whatever keys are passed
		#sub_Stot = [gas_dict["Stot"] for gas_name, gas_dict in sub_gases.iteritems()]

		def sums(length, total_sum):
			"""Returns a list of all possible arrays of integers, where the sum of all array elements is 'total_sum", 
			where each array is 'length' values long.

			NOTE: This function can take an extremely long amount of time. Each +1 increase to the "length" value 
			results in an exponential increase in computation time. Timing was tested on a 2014 MacBook Pro with 3 GHz 
			processor, 16 GB memory, and solid state drive, using one core. Timing for length=5 ~2 seconds. Timing for
			length=6 ~230 seconds.

			Parameters
			----------
			length: int
				Length of generated arrays. Example: length=3 gives a list of arrays (value, value, value)

			total_sum: int
				All values within each array must sum to total_sum. Example: length=3, total_sum=2 gives:
				(2, 0, 0), (1, 0, 1), (1, 1, 0), (0, 0, 2), (0, 2, 0), (0, 1, 1)

			Returns
			-------
			generator object
				Generator object of all generated arrays. To return a list of arrays, pass list(sums(length, total_sum))

			"""
			if length == 1:
				yield (total_sum,)
			else:
				for value in range(total_sum + 1):
					for permutation in sums(length - 1, total_sum - value):
						yield (value,) + permutation

		pos = list(sums(len(sub_gases), 100))

		result_list = []
		sum_CO2 = 0.0
		sum_H2O = 0.0
		sum_SO2 = 0.0
		#sum_H2S = 0.0
		#sum_Stot = 0.0
		for combo in pos:
			for i in range(len(sub_gases)):
				sum_CO2 += sub_CO2[i] * (combo[i]/100.0)
				sum_H2O += sub_H2O[i] * (combo[i]/100.0)
				sum_SO2 += sub_SO2[i] * (combo[i]/100.0)
				#sum_H2S += sub_H2S[i] * (combo[i]/100.0)
				#sum_Stot += sub_Stot[i] * (combo[i]/100.0)

				if sum_CO2 < (surface_gas["CO2"] + surface_gas["CO2"]*threshold) and sum_CO2 > (surface_gas["CO2"] - surface_gas["CO2"]*threshold):
					if sum_H2O < (surface_gas["H2O"] + surface_gas["H2O"]*threshold) and sum_H2O > (surface_gas["H2O"] - surface_gas["H2O"]*threshold):
						if sum_SO2 < (surface_gas["SO2"] + surface_gas["SO2"]*threshold) and sum_SO2 > (surface_gas["SO2"] - surface_gas["SO2"]*threshold):
							#if sum_H2S < (surface_gas["H2S"] + surface_gas["H2S"]*threshold) and sum_H2S > (surface_gas["H2S"] - surface_gas["H2S"]*threshold):
								#if sum_Stot < (surface_gas["Stot"] + surface_gas["Stot"]*threshold) and sum_Stot > (surface_gas["Stot"] - surface_gas["Stot"]*threshold):
							result_list.append(combo)

		print result_list




















		






