from math import *
from scipy.optimize import fsolve

def equation(p):
			fS2 = p
			gS2 = 
			KS2 = 
			fO2 =
			gSO2 =  
			KHS2 = 
			fH2 = 
			gH2S = 
			PStot = 
			return 	( (fS2 / gS2) + (KSO2 * sqrt(fS2) * fO2 / gSO2) + (KH2S * fH2 * sqrt(fS2) / gH2S) - PStot)

		fS2 = fsolve(equation, (0.0001))