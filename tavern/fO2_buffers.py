def calc_logfO2_from_buffer(temperature, buffer, delta, pressure=1.0):
    """Returns logfO2 value given value relative to common solid buffers

    Parameters
    ----------
    temperature: float
        Temperature in degrees C.
    pressure: float
        Optional. Pressure in bars. If no value is passed, default value of
        1.0 bars will be used
    buffer: str
        Name of solid buffer referenced. Possible strings are: QFM... will be
        adding more soon! #TODO add more buffers!
    delta: float
        Number of log units away from given buffer. Example: for QFM+2, one would enter
        buffer='QFM' and delta=2

    Returns
    -------
    float
        Value of logfO2
    """
    def calc_QFM(P, T):
        """
        Quartz-Fayalite-Magnetite (QFM)
        ===============================
        Define QFM buffer value at P

        Parameters
        ----------
        P: float
            Pressure in bars

        T: float
            Temperature in degrees C

        Returns
        -------
        float or numpy array
            fO2

        References
        ----------
        B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
        """
        tempK = temperature + 273.15

        if tempK<573:
            fO2 = (-26445.3/tempK) + 10.344 + 0.092 * (P-1)/tempK
        if tempK>=573:
            fO2 = (-25096.3/tempK) + 8.735 + 0.11 * (P-1)/tempK

        return fO2

    if buffer == 'QFM':
        fO2 = calc_QFM(pressure, temperature)

    if buffer == 'NNO':
        #TODO
        pass

    if buffer == 'IW':
        #TODO
        pass

    return fO2 + delta