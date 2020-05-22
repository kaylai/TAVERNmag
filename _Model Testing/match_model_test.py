from math import *
import time
import csv

surface_gas_molpercent = {"H2O": 93.32,
						  "CO2": 3.33,
						  "SO2": 3.34}

shallow_eq = {"H2O": 71.41,
			  "CO2": 19.93,
			  "SO2": 8.67}

int_eq = {"H2O": 24.43,
			  "CO2": 3.83,
			  "SO2": 71.74}

deep_eq = {"H2O": 35.03,
			  "CO2": 4.68,
			  "SO2": 60.29}

shallow_degassing = {"H2O": 98.36,
			  "CO2": 0.04,
			  "SO2": 1.60}

int_degassing = {"H2O": 27.43,
			  "CO2": 0.06,
			  "SO2": 72.51}

deep_degassing = {"H2O": 90.97,
			  "CO2": 0.19,
			  "SO2": 8.84}

sub_gases = {"Shallow_EQ_Fluid": shallow_eq,
			 "Int_EQ_Fluid": int_eq,
			 "Deep_EQ_Fluid": deep_eq,
			 "Shallow_Degassing": shallow_degassing,
			 "Int_Degassing": int_degassing,
			 "Deep_Degassing": deep_degassing}

threshold = 1.0 #in absolute mol% value

sub_CO2 = {gas_name: gas_dict["CO2"] for gas_name, gas_dict in sub_gases.iteritems()}
sub_H2O = {gas_name: gas_dict["H2O"] for gas_name, gas_dict in sub_gases.iteritems()}
sub_SO2 = {gas_name: gas_dict["SO2"] for gas_name, gas_dict in sub_gases.iteritems()}

start_time = time.time()
def sums(length, total_sum):
    if length == 1:
        yield (total_sum,)
    else:
        for value in range(total_sum + 1):
            for permutation in sums(length - 1, total_sum - value):
                yield (value,) + permutation

final_list = []
for l in list(sums(6,100)):
	sum_CO2 = ( l[0]/100. * sub_CO2["Shallow_EQ_Fluid"] + 
				l[1]/100. * sub_CO2["Int_EQ_Fluid"] + 
				l[2]/100. * sub_CO2["Deep_EQ_Fluid"] +
				l[3]/100. * sub_CO2["Shallow_Degassing"] +
				l[4]/100. * sub_CO2["Int_Degassing"] +
				l[5]/100. * sub_CO2["Deep_Degassing"] )

	if sum_CO2 < surface_gas_molpercent["CO2"] + threshold and sum_CO2 > surface_gas_molpercent["CO2"] - threshold:
		sum_H2O = ( l[0]/100. * sub_H2O["Shallow_EQ_Fluid"] + 
					l[1]/100. * sub_H2O["Int_EQ_Fluid"] + 
					l[2]/100. * sub_CO2["Deep_EQ_Fluid"] +
					l[3]/100. * sub_H2O["Shallow_Degassing"] +
					l[4]/100. * sub_H2O["Int_Degassing"] +
					l[5]/100. * sub_H2O["Deep_Degassing"] )

		if sum_H2O < surface_gas_molpercent["H2O"] + threshold and sum_H2O > surface_gas_molpercent["H2O"] - threshold:
			sum_SO2 = ( l[0]/100. * sub_SO2["Shallow_EQ_Fluid"] + 
						l[1]/100. * sub_SO2["Int_EQ_Fluid"] + 
						l[2]/100. * sub_CO2["Deep_EQ_Fluid"] +
						l[3]/100. * sub_SO2["Shallow_Degassing"] +
						l[4]/100. * sub_SO2["Int_Degassing"] +
						l[5]/100. * sub_SO2["Deep_Degassing"] )

			if sum_SO2 < surface_gas_molpercent["SO2"] + threshold and sum_SO2 > surface_gas_molpercent["SO2"] - threshold:
				final_list.append(l)

end_time = time.time()
print ("--- %s seconds ---" % (end_time - start_time))

print final_list

with open("match_model_test_output.csv", "wb") as f:
	writer = csv.writer(f)
	writer.writerows(final_list)


#L = list(sums(6,100))
# print L[-100:]




# start_time = time.time()
# possible_gas_mixture = []
# for a in range(101):
# 	for b in range(101):
# 		if a+b > 100:
# 			break
# 		for c in range(101):
# 			if a+b+c > 100:
# 				break
# 			for d in range(101):
# 				if a+b+c+d > 100:
# 					break
# 				for e in range(101):
# 					if a+b+c+d+e > 100:
# 						break
# 					for f in range(101):
# 						if a+b+c+d+e+f == 100:
# 							possible_gas_mixture.append((a, b, c, d, e, f))

# print possible_gas_mixture
# print ("--- %s seconds ---" % (time.time() - start_time))



# start_time = time.time()
# import itertools
# numbers = [1, 2, 3, 7, 7, 9, 10]
# result = [seq for i in range(len(numbers), 0, -1) for seq in itertools.combinations(numbers, i) if sum(seq) == 10]
# print result
# print ("--- %s seconds ---" % (time.time() - start_time))

# result = itertools.takewhile(a+b+c+d <= 100, [a, ])
# print result













