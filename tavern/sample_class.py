import pandas as pd

from tavern import core

from copy import deepcopy, copy

class Sample(object):
    """ Generic sample that could be a MagmaticFluid or SilicateMelt object.

    Attributes
    ----------
    composition:    dict
        Dictionary of major element or volatile concentrations
    """

    def __init__(self, composition, units='wtpercent', default_normalization='none',
                 default_units='wtpercent', **kwargs):

        composition = deepcopy(composition)
        
        if isinstance(composition, dict):
            composition = pd.Series(composition, dtype='float64')
        elif isinstance(composition, pd.Series) is False:
            raise core.InputError("The composition must be given as either a dictionary or a "
                                  "pandas Series.")
        
        if units == "wtpercent":
            self._composition = composition
        elif units == "molpercent":
            self._composition = self._moles_to_wtpercent(composition) 
        elif units == "molfrac":
            self._composition = self._moles_to_wtpercent(composition) 
        else:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        self.set_default_normalization(default_normalization)
        self.set_default_units(default_units)

        # handle possibly passed FeOT values, convert to FeO and warn user
        possible_FeOT_names = ['FeOT', 'FeO*', 'FeOtot', 'FeOt', 'FeOtotal',
                               'FeOstar']

        for name in possible_FeOT_names:
            if name in composition:
                if 'FeO' in composition:
                    w.warn("FeO and " + str(name) + " oxides passed. Discarding " + str(name) +
                           " oxide.", RuntimeWarning, stacklevel=2)
                    if isinstance(composition, dict):
                        del composition[name]
                    if isinstance(composition, pd.Series):
                        composition.drop(labels=[name], inplace=True)
                else:
                    w.warn(str(name) + " oxide found. Using " + str(name) + " for FeO value.",
                           RuntimeWarning, stacklevel=2)
                    composition['FeO'] = composition[name]
                    if isinstance(composition, dict):
                        del composition[name]
                    if isinstance(composition, pd.Series):
                        composition.drop(labels=[name], inplace=True)

    def set_default_normalization(self, default_normalization):
        """ Set the default type of normalization to use with the get_composition() method.

        Parameters
        ----------
        default_normalization:    str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            # - 'fixedoxygen': For MagmaticFluid objects only. Normalizes a fluid composition to 100%,
            #   including O2. The O2 wt% will remain fixed, whilst the other species are reduced
            #   proprotionally so that the total is 100 wt%
            - 'fixedvolatiles': For SilicateMelt objects only. Normalizes major element oxides to
              100 wt%, including volatiles. The volatile wt% will remain fixed, whilst the other
              major element oxides are reduced proportionally so that the total is 100 wt%.
            - 'additionalvolatiles': For SilicateMelt objects only. Normalises major element oxide
              wt% to 100%, assuming it is volatile-free. If H2O or CO2 are passed to the function,
              their un-normalized values will be retained in addition to the normalized non-volatile
              oxides, summing to >100%.
        """
        if default_normalization in ['none', 'standard', 'fixedvolatiles', 'additionalvolatiles']:
            self.default_normalization = default_normalization
        else:
            raise core.InputError("The normalization method must be one of 'none', 'standard', "
                                  "'fixedvolatiles', or 'additionalvolatiles'.")

    def set_default_units(self, default_units):
        """ Set the default units of composition to return when using the get_composition() method.

        Parameters
        ----------
        default_units     str
            The type of composition to return, one of:
            - wtpercent (default)
            - molpercent
            - molfrac
        """
        if default_units in ['wtpercent', 'molpercent', 'molfrac']:
            self.default_units = default_units
        else:
            raise core.InputError("The units must be one of 'wtpercent', 'molpercent', "
                                  "or 'molfrac'.")

    def get_composition(self, species=None, normalization=None, units=None, exclude_volatiles=False,
                        asSampleClass=False, oxide_masses={}):
        """ Returns the composition in the format requested, normalized as requested.

        Parameters
        ----------
        species:    NoneType or str
            The name of the oxide or cation to return the concentration of. If NoneType (default)
            the whole composition will be returned as a pandas.Series. If an oxide is passed, the
            value in wtpt will be returned unless units is set to 'molpercent', even if the
            default units for the sample object are molpercent. If an element is passed, the
            concentration will be returned as molfrac, unless 'mol_singleO' is specified as
            units, even if the default units for the sample object are mol_singleO. Unless
            normalization is specified in the method call, none will be applied.

        normalization:     NoneType or str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            If NoneType is passed the default normalization option will be used
            (self.default_normalization).

        units:     NoneType or str
            The units of composition to return, one of:
            - wtpercent (default)
            - molpercent
            - molfrac
            If NoneType is passed the default units option will be used (self.default_type).

        exclude_volatiles   bool
            If True, volatiles H2O, CO2, and S will be excluded from the returned composition, prior
            to normalization and conversion. Can only be set to true for SilicateMelt objects.

        asSampleClass:  bool
            If True, the sample composition will be returned as a sample class, with default
            options. In this case any normalization instructions will be ignored.

        oxide_masses:  dict
            Specify here any oxide masses that should be changed from the tavern default. This
            might be useful for recreating other implementations of models that use slightly
            different molecular masses. The default values in tavern are given to 3 dp.

        Returns
        -------
        pandas.Series, float, or Sample class
            The sample composition, as specified.
        """
        # Process the oxide_masses variable, if necessary:
        oxideMass = copy(core.oxideMass)
        for ox in oxide_masses:
            if ox in oxideMass:
                oxideMass[ox] = oxide_masses[ox]
            else:
                raise core.InputError("The oxide name provided in oxide_masses is not recognised.")

        # Fetch the default return types if not specified in function call
        if normalization is None and species is None:
            normalization = self.default_normalization
        if units is None and species is None:
            units = self.default_units

        # Check whether to exclude volatiles
        # note that here composition is gotten as wtpercent
        if exclude_volatiles:
            composition = self._composition.copy()
            if 'H2O' in composition.index:
                composition = composition.drop(index='H2O')
            if 'CO2' in composition.index:
                composition = composition.drop(index='CO2')
            if 'S' in composition.index:
                composition = composition.drop(index='S')
        else:
            composition = self._composition.copy()

        # Check for a species being provided, if so, work out which units to return.
        if isinstance(species, str):
            if species in composition.index:  # if the requested species has a value, proceed
                if species in core.oxides:
                    if units in ['molpercent, molfrac'] or units is None:
                        units = 'wtpercent'
                else:
                    raise core.InputError(species + " was not recognised, check spelling, " +
                                          "capitalization and stoichiometry.")
                if normalization is None:
                    normalization = 'none'
            else:
                return 0.0  # if the requested species has no set value, return a float of 0.0
        elif species is not None:
            raise core.InputError("Species must be either a string or a NoneType.")

        # Get the requested type of composition
        if units == 'wtpercent':
            converted = composition
        elif units == 'molpercent':
            converted = self._wtpercent_to_molpercent(composition, oxideMass=oxideMass)
        elif units == 'molfrac':
            converted = self._wtpercent_to_molfrac(composition, oxideMass=oxideMass)
        else:
            raise core.InputError("The units must be one of 'wtpercent', 'molpercent', "
                                  "or 'molfrac'.")

        # Do requested normalization
        if normalization == 'none':
            final = converted
        elif normalization == 'standard':
            final = self._normalize_Standard(converted, units=units)
        else:
            raise core.InputError("The normalization method must be one of 'none', 'standard', "
                                  "'fixedvolatiles', or 'additionalvolatiles'.")

        if species is None:
            if asSampleClass is False:
                return final
            else:
                return Sample(final)
        elif isinstance(species, str):
            if asSampleClass:
                w.warn("Cannot return single species as Sample class. Returning as float.",
                       RuntimeWarning, stacklevel=2)
            return final[species]

    def change_composition(self, new_composition, units='wtpercent', inplace=True):
        """
        Change the concentration of some component of the composition.

        If the units are moles, they are read as moles relative to the present composition,
        i.e. if you wish to double the moles of MgO, if the present content is 0.1 moles,
        you should provide {'MgO':0.2}. The composition will then be re-normalized. If the
        original composition was provided in un-normalized wt%, the unnormalized total will
        be lost.

        Parameters
        ----------
        new_composition:    dict or pandas.Series
            The components to be updated.
        units:      str
            The units of new_composition. Should be one of:
            - wtpercent (default)
            - molpercent
            - molfrac
        inplace:    bool
            If True the object will be modified in place. If False, a copy of the Sample
            object will be created, modified, and then returned.

        Returns
        -------
        Sample class
            Modified Sample class.
        """

        # if new_composition is pandas.Series, convert to dict
        if isinstance(new_composition, pd.Series):
            new_composition = dict(new_composition)

        if inplace is False:
            newsample = deepcopy(self)
            return newsample.change_composition(new_composition, units=units)

        if units == 'wtpercent':
            for ox in new_composition:
                self._composition[ox] = new_composition[ox]

        elif units == 'molpercent':
            _comp = self.get_composition(units='molpercent')
            for ox in new_composition:
                _comp[ox] = new_composition[ox]
            self._composition = self._moles_to_wtpercent(_comp)

        elif units == 'molfrac':
            _comp = self.get_composition(units='molfrac')
            for el in new_composition:
                _comp[el] = new_composition[el]
            self._composition = self._moles_to_wtpercent(_comp)

        else:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        return self

    def delete_oxide(self, oxide, inplace=True):
        """ Allows user to remove a given oxide from the Sample composition

        Parameters
        ----------
        oxide:  str or list
            Name or names of the oxide(s) to remove.

        inplace:    bool
            If True the object will be modified in place. If False, a copy of the Sample
            object will be created, modified, and then returned.

        Returns
        -------
        Sample class
            Modified Sample class.
        """
        # if new_composition is pandas.Series, convert to dict
        if isinstance(oxide, str):
            oxide = [oxide]

        if inplace is False:
            newsample = deepcopy(self)
            return newsample.delete_oxide(oxide)

        self._composition.drop(index=oxide, inplace=True)

        return self

    def check_oxide(self, oxide):
        """
        Check whether the sample composition contains the given oxide.

        Parameters
        ----------
        oxide:  str
            Oxide name to check composition for.

        Returns
        -------
        bool
            Whether the composition contains the given oxide, or not.
        """

        if oxide not in core.oxides:
            w.warn("Oxide name not recognised. If it is in your sample, unexpected behaviour "
                   "might occur!",
                   RuntimeWarning, stacklevel=2)
        return oxide in self._composition

    def _normalize_Standard(self, composition, units='wtpercent'):
        """
        Normalizes the given composition to 100 wt%, including volatiles. This method
        is intended only to be called by the get_composition() method.

        Parameters
        ----------
        composition:     pandas.Series
            A rock composition with oxide names as keys and concentrations as values.

        units:      str
            The units of composition. Should be one of:
            - wtpercent (default)
            - molpercent
            - molfrac

        Returns
        -------
        pandas.Series
            Normalized oxides in wt%.
        """
        comp = composition.copy()
        comp = dict(comp)

        if units == 'wtpercent':
            normed = pd.Series({k: 100.0 * v / sum(comp.values()) for k, v in comp.items()})
        elif units == 'molpercent' or units == 'molfrac':
            normed = pd.Series({k: v / sum(comp.values()) for k, v in comp.items()})
        else:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        return normed

    # def _normalize_fixedoxygen(composition, units='wtpercent'):
    # """
    # Normalizes a fluid composition to 100%, including O2. The O2 wt%
    # will remain fixed, whilst the other species are reduced proprotionally
    # so that the total is 100 wt%

    # Parameters
    # ----------
    # composition: pandas.Series
    #     fluid composition

    # units:  str
    #     The units of composition. Should be one of:
    #     - wtpercent (default)
    #     - molpercent
    #     - molfrac

    # Returns
    # -------
    # pandas.Series
    #     Normalized fluid composition.
    # """
    # composition = composition.copy()
    # composition = dict(composition)

    # normalized = {}
    # O2_sum = composition['O2']

    # for spec in composition.keys():
    #     if spec != 'O2':
    #         normalized[spec] = composition[spec]
    # normsum = sum(normalized.values())

    # if units == 'wtpercent' or units == 'molpercent':
    #     normalized = {k: v/normsum*(100.-O2_sum) for k, v in normalized.items()}
    # elif units == 'molfrac':
    #     normalized = {k: v/normsum*(1.-O2_sum) for k, v in normalized.items()}
    # else:
    #     raise Exception("Units must be one of 'wtpercent', "
    #                     "'molpercent', or 'molfrac'.")

    # normalized['O2'] = composition['O2']
    # normalized = pd.Series(normalized)

    # return normalized

    def _normalize_FixedVolatiles(self, composition, units='wtpercent'):
        """
        Normalizes major element oxides to 100 wt%, including volatiles. The volatile wt% will
        remain fixed, whilst the other major element oxides are reduced proportionally so that the
        total is 100 wt%.

        Intended to be called only by the get_composition() method.

        Parameters
        ----------
        composition:     pandas Series
            Major element composition

        units:      str
            The units of composition. Should be one of:
            - wtpercent (default)
            - molpercent
            - molfrac

        Returns
        -------
        pandas Series
            Normalized major element oxides.
        """
        comp = composition.copy()
        normalized = pd.Series({}, dtype=float)
        volatiles = 0
        if 'CO2' in list(comp.index):
            volatiles += comp['CO2']
        if 'H2O' in list(comp.index):
            volatiles += comp['H2O']
        if 'S' in list(comp.index):
            volatiles += comp['S']

        for ox in list(comp.index):
            if ox != 'H2O' and ox != 'CO2' and ox != 'S':
                normalized[ox] = comp[ox]

        if units == 'wtpercent':
            normalized = normalized/np.sum(normalized)*(100-volatiles)
        elif units == 'molpercent' or units == 'molfrac':
            normalized = normalized/np.sum(normalized)*(1-volatiles)
        else:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        if 'CO2' in list(comp.index):
            normalized['CO2'] = comp['CO2']
        if 'H2O' in list(comp.index):
            normalized['H2O'] = comp['H2O']
        if 'S' in list(comp.index):
            normalized['S'] = comp['S']

        return normalized

    def _normalize_AdditionalVolatiles(self, composition, units='wtpercent'):
        """
        Normalises major element oxide wt% to 100%, assuming it is volatile-free. If H2O or CO2
        are passed to the function, their un-normalized values will be retained in addition to the
        normalized non-volatile oxides, summing to >100%.

        Intended to be called only by the get_composition() method.

        Parameters
        ----------
        composition:     pandas.Series
            Major element composition

        units:      str
            The units of composition. Should be one of:
            - wtpercent (default)
            - molpercent
            - molfrac

        Returns
        -------
        pandas.Series
            Normalized major element oxides.
        """
        comp = composition.copy()
        normalized = pd.Series({}, dtype=float)
        for ox in list(comp.index):
            if ox != 'H2O' and ox != 'CO2' and ox != 'S':
                normalized[ox] = comp[ox]

        if units == 'wtpercent':
            normalized = normalized/np.sum(normalized)*100
        elif units == 'molpercent' or units == 'molfrac':
            normalized = normalized/np.sum(normalized)
        else:
            raise core.InputError("Units must be one of 'wtpercent', 'molpercent', or "
                                  "'molfrac'.")

        if 'H2O' in comp.index:
            normalized['H2O'] = comp['H2O']
        if 'CO2' in comp.index:
            normalized['CO2'] = comp['CO2']
        if 'S' in comp.index:
            normalized['S'] = comp['S']

        return normalized

    def _moles_to_wtpercent(self, composition, oxideMass=core.oxideMass):
        """Converts composition in mole fraction or mole percent to wt percent.

        Parameters
        ----------
        composition: pandas.Series
            Composition in mole fraction with possible keys:
            SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
            H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

        Returns
        -------
        pandas.Series
            Composition in wt percent with possible keys:
            SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
            H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
        """
        moles = composition.copy()

        MPO_dict = {}
        for key, value in moles.items():
            MPO_dict[key] = value * oxideMass[key]

        MPO_sum = sum(MPO_dict.values())

        wtpercent_dict = {}
        for key, value in MPO_dict.items():
            wtpercent_dict[key] = 100.0 * value / MPO_sum

        wtpercent = pd.Series(wtpercent_dict)

        return wtpercent

    def _wtpercent_to_molfrac(self, composition, oxideMass=core.oxideMass):
        """Converts composition in wt percent to mole fraction.

        Parameters
        ----------
        composition: pandas.Series
            Composition in wt percent with possible keys:
            SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
            H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

        Returns
        -------
        pandas.Series
            Composition in mole fraction with possible keys:
            SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
            H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
        """
        wtpercent = composition.copy()

        MPO_dict = {}
        for key, value in wtpercent.items():
            MPO_dict[key] = value / oxideMass[key]

        MPO_sum = sum(MPO_dict.values())

        molfrac_dict = {}
        for key, value in MPO_dict.items():
            molfrac_dict[key] = value / MPO_sum

        molfrac = pd.Series(molfrac_dict)

        return molfrac

    def _wtpercent_to_molpercent(self, composition, oxideMass=core.oxideMass):
        """Converts composition in wt percent to mole percent.

        Parameters
        ----------
        wtpercent: pandas.Series
            Composition in wt percent with possible keys:
            SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
            H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2

        Returns
        -------
        pandas.Series
            Composition in mole percent with possible keys:
            SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, Cr2O3, NiO, F, Cl, 
            H2O, CO2, S, S2, CO, H2, H2S, SO2, O, O2
        """
        wtpercent = composition.copy()

        MPO_dict = {}
        for key, value in wtpercent.items():
            MPO_dict[key] = value / oxideMass[key]

        MPO_sum = sum(MPO_dict.values())

        molfrac_dict = {}
        for key, value in MPO_dict.items():
            molfrac_dict[key] = 100.0 * value / MPO_sum

        molpercent = pd.Series(molfrac_dict)

        return molpercent
