import numpy as np

def std_atmosphere_SI(ha_meters):
    # Convert input altitude from meters to feet for internal calculations
    ha = ha_meters / 0.3048  # meters to feet
    
    # constants
    Ra = 1716.49   # specific gas constant for air (ft lb / slug R)
    g = 32.174049  # acceleration due to gravity (ft/s^2)
    gamma = 1.4
    Su = 198.72    # Sutherland's constant (R)
    bu = 2.2697E-8 # (slug/ft-s-sqrt(R))

    # SI constants
    Ra_SI = 287.058    # specific gas constant for air (J/kg·K)
    g_SI = 9.80665     # acceleration due to gravity (m/s^2)
    Su_SI = 110.4      # Sutherland's constant (K)
    bu_SI = 1.458e-6   # (kg/m-s-sqrt(K))

    # Conversion factors
    FT_TO_M = 0.3048           # feet to meters
    PSF_TO_PA = 47.880258      # pounds per sq foot to Pascal
    SLUG_TO_KG = 14.593903     # slug to kilogram
    RANKINE_TO_KELVIN = 5/9    # Rankine to Kelvin
    SLUGFT_TO_KGMS = 1.488164  # slug/ft-s to kg/m-s

    ### parameters that define the standard atmosphere...
    # altitudes (ft)
    H1 = [0, 36089, 65617, 104987, 154199, 167323, 232940, 278385, 298556, 360892, 393701]
    
    # lapse rates (R/ft)
    A1 = [-3.566E-3, 0, 5.486E-4, 1.540E-3, 0, -1.540E-3, -1.1004E-3, 0, -137.382, 6.584E-3, 0]
    
    # base layer Temp definitions (R)
    T1 = [518.67, 389.97, 389.97, 411.57, 487.17, 487.17, 386.37, 336.361, 473.743, 432]
    T1.append(T1[9] + A1[9]*(H1[10]-H1[9]))  # temp at 120 km
    
    # pressure definitions (psf)
    P1 = [2116.21662, 472.688, 114.345, 18.129, 2.31634, 1.39805, 0.082632, 0.0077986]
    P1.append(P1[7]*np.exp(-1*(g/Ra/T1[7])*(H1[8]-H1[7])))
    P1.append(P1[8]*np.exp(-1*(g/Ra/T1[8])*(H1[9]-H1[8])))
    P1.append(P1[9]*np.exp(-1*(g/Ra/T1[9])*(H1[10]-H1[9])))
    
    Rho1 = np.divide(np.array(P1[0:10]), np.array(T1[0:10]))/Ra

    idx_0 = list(map(lambda i: i > ha, H1)).index(True) - 1
    a0 = A1[idx_0]    # lapse rate for the layer
    T0 = T1[idx_0]    # temperature at the base of the layer
    H0 = H1[idx_0]    # altitude at the base of the layer
    P0 = P1[idx_0]    # pressure at the base of the layer
    Rho0 = Rho1[idx_0]  # density at the base of the layer
    
    if idx_0 != 8:  # the linear and zero lapse layers (all layers except 91-120km)
        Ta = T0 + a0 * (ha - H0)
        if a0 != 0:
            Pa = P0 * (Ta/T0)**(-g/a0/Ra)
            Rhoa = Rho0 * (Ta/T0)**(-1*(g/a0/Ra+1))
        else:
            Pa = P0 * np.exp(-1*(g/Ra/Ta) * (ha-H0))    
            Rhoa = Rho0 * np.exp(-1*(g/Ra/Ta) * (ha-H0))  
    else:  # 91-120 km elliptic temperature profile layer
        Ta = T0 + a0 * np.sqrt(1 - ((ha - H0) / (-19.9429*3281))**2)
        Pa = P0 * (Ta/T0)**(-g/a0/Ra)
        Rhoa = Rho0 * (Ta/T0)**(-1*(g/a0/Ra+1))
    
    # Calculate Imperial values
    a = np.sqrt(gamma * Ra * Ta)
    mu = bu*Ta**1.5/(Ta+Su)
    delta = Pa/P1[0]
    theta = Ta/T1[0]
    sigma = Rhoa/Rho1[0]
    
    # Calculate SI values
    ha_SI = ha * FT_TO_M
    Ta_SI = Ta * RANKINE_TO_KELVIN
    Pa_SI = Pa * PSF_TO_PA
    Rhoa_SI = Rhoa * SLUG_TO_KG / (FT_TO_M**3)
    a_SI = a * FT_TO_M
    mu_SI = mu * SLUGFT_TO_KGMS
        
    return (delta, theta, sigma, 
            Pa, Ta, Rhoa, a, mu,  # Imperial
            Pa_SI, Ta_SI, Rhoa_SI, a_SI, mu_SI)  # SI

loop = True
while loop:
    print('')
    # user input
    h_in = input('altitude (m): '.ljust(16))
    
    h_in = float(h_in)
    max_alt_meters = 120000  # approximately 393,701 ft
    if not (0 <= h_in <= max_alt_meters):
        raise ValueError(f'altitude must be between 0 and {max_alt_meters} meters')
    
    # define the standard atmosphere
    delta, theta, sigma, Ph, Th, Rhoh, ah, muh, Ph_SI, Th_SI, Rhoh_SI, ah_SI, muh_SI = std_atmosphere_SI(h_in)
    
    # output with both unit systems
    print('\nInput altitude:', f'{h_in:.1f} m', f'({h_in/0.3048:.1f} ft)')
    
    print('\nNon-dimensional parameters:')
    print('   delta:'.ljust(21) + '{:0.4f}'.format(delta))
    print('   theta:'.ljust(21) + '{:0.4f}'.format(theta))
    print('   sigma:'.ljust(21) + '{:0.4f}'.format(sigma))
    
    print('\nSI Units:')
    print('   P (Pa):'.ljust(21) + '{:0.2f}'.format(Ph_SI))
    print('   T (K):'.ljust(21) + '{:0.2f}'.format(Th_SI))
    print('   rho (kg/m^3):'.ljust(21) + '{:0.3e}'.format(Rhoh_SI))
    print('   a (m/s):'.ljust(21) + '{:0.1f}'.format(ah_SI))
    print('   mu (kg/m-s):'.ljust(21) + '{:0.3e}'.format(muh_SI))
    
    print('\nImperial Units:')
    print('   P (psf):'.ljust(21) + '{:0.2f}'.format(Ph))
    print('   T (R):'.ljust(21) + '{:0.2f}'.format(Th))
    print('   rho (slug/ft^3):'.ljust(21) + '{:0.3e}'.format(Rhoh))
    print('   a (fps):'.ljust(21) + '{:0.1f}'.format(ah))
    print('   mu (slug/ft-s):'.ljust(21) + '{:0.3e}'.format(muh))