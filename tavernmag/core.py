import sys

# -----------Define Global Variables ----------- #
#Molecular weights of components in grams per mol
oxideMass = {'SiO2':  60.083,
             'MgO':   40.304,
             'FeO':   71.844,
             'CaO':   56.077,
             'Al2O3': 101.961,
             'Na2O':  61.979,
             'K2O':   94.195,
             'MnO':   70.937,
             'TiO2':  79.867,
             'P2O5':  141.943,
             'Cr2O3': 151.992,
             'NiO':   74.692,
             'CoO':   44.01,
             'Fe2O3': 159.687,
             'H': 	  1.01,
             'H2':	  2.016,
             'H2O':   18.015,
             'C': 	  12.0107,
             'CO2':   44.010,
             'F2O':   37.997,
             'F':	  18.992,
             'Cl':	  35.453,
             'CO': 	  20.010,
             'O':	  15.999,
             'O2':	  15.999*2.0,
             'S':	  32.065,
             'Stot':  32.065,
             'S2':	  32.065*2.0,
             'SO2':	  64.066,
             'H2S':	  34.100}

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

class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class status_bar(object):
    """Various styles of status bars that display the progress of a calculation
    within a loop
    """
    def __init__():
        pass

    def status_bar(percent, sample_name=None, btext=None, barLen=20):
        """
        Prints an updating status bar to the terminal or jupyter notebook.

        Parameters
        ----------
        percent: float
            Percent value of progress from 0 to 1

        sample_name: string
            Name of the current sample being calculated

        btext: string
            Any extra text to display next to status bar

        barLen: int
            Length of bar to print
        """
        sys.stdout.write("\r")
        sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent),
                                                   barLen, percent * 100))

        sample_string = str(sample_name)
        # Set max number of characters in sample name
        max_name_length = 25
        if len(str(sample_name)) >= max_name_length:
            sample_string = str(sample_name)[0:max_name_length-1] + "..."

        # Write out sample name and trailing spaces to cover contents of
        # previous sample names left over on line
        if sample_name is not None:
            sys.stdout.write("  Working on sample " + sample_string +
                             "                            ")
        if btext is not None:
            sys.stdout.write(" " + str(btext))
        if percent == 1.0:
            sys.stdout.write("\n")
        sys.stdout.flush()
