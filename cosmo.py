from astropy.cosmology import FlatLambdaCDM

def flagship_cosmology():
    # Flagship cosmological parameters
    H0 = 67 # Hubble parameter (today) [km/s/Mpc]
    Om0 = 0.319 # Total matter density parameter (today) [dimensionless]
    Ob0 = 0.049 # Baryon density parameter (today) [dimensionless]
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    return cosmo
