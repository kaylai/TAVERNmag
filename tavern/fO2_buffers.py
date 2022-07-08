""" Generic functions to calculate logfO2 values of common buffers at P, T """

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
    float
        logfO2

    References
    ----------
    B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
    """
    tempK = T + 273.15

    if tempK<573:
        log_fO2 = (-26445.3/tempK) + 10.344 + 0.092 * (P-1)/tempK
    if tempK>=573:
        log_fO2 = (-25096.3/tempK) + 8.735 + 0.11 * (P-1)/tempK

    return log_fO2

def calc_NNO(P, T):
    """ 
    Ni-NiO (NNO)
    ============
    Define NNO buffer value at P

    References
    ----------
    Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide
    oxygen fugacity buffers

    Parameters
    ----------
    P: float
        Pressure in bars

    T: float or numpy array
        Temperature in degrees C

    Returns
    -------
    float
        logfO2

    Polynomial coefficients
    -----------------------  
    log fO2  =  (a0+a1*P+a2*P^2+a3*P^3+a4*P^4) + (b0+b1*P+b2*P^2+b3*P^3)/T
    a0: 8.699
    a1: 0.01642
    a2: -0.0002755
    a3: 0.000002683
    a4: -1.015E-08
    b0: -24205
    b1: 444.73
    b2: -0.59288
    b3:0.0015292                            
    """
    P_GPa = P/10000
    T_K = T + 273.15

    log_fO2 = ((8.699 + 0.01642*P_GPa - 0.0003*P_GPa**2 + (2.7*10**(-6))*P_GPa**3 -
                (10**(-8))*P_GPa**4) +
               (-24205 + 444.73*P_GPa - 0.5929*P_GPa**2 + 0.00153*P_GPa**3)/T_K)

    return log_fO2

def calc_IW(P, T):
    """
    Fe-FeO (Iron-Wustite)
    =====================
    Define IW buffer value at P
    
    References
    ----------
    Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers

    Parameters
    ----------
    P: float
        Pressure in bars

    T: float or numpy array
        Temperature in degrees C

    Returns
    -------
    float or numpy array
        log_fO2

    Polynomial coefficients
    -----------------------
    log fO2  =  (a0+a1*P) + (b0+b1*P+b2*P^2+b3*P^3)/T
    a0: 6.54106 
    a1: 0.0012324
    b0: -28163.6
    b1: 546.32
    b2: -1.13412
    b3: 0.0019274               
    """
    P_GPa = P/10000
    T_K = T + 273.15

    log_fO2 = ((6.54106+0.0012324*P_GPa) + 
               (-28163.6+546.32*P_GPa-1.13412*P_GPa**2+0.0019274*P_GPa**3)/T_K)

    return log_fO2

""" User function to claculate the logfO2 value from a buffer and delta """
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
    if buffer == 'QFM':
        log_fO2 = calc_QFM(pressure, temperature)

    if buffer == 'NNO':
        log_fO2 = calc_NNO(pressure, temperature)
        pass

    if buffer == 'IW':
        log_fO2 = calc_IW(pressure, temperature)
        pass

    return log_fO2 + delta
