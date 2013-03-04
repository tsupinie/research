import numpy as np

def computeReflectivity(**kwargs):
    if 'method' in kwargs:
        if kwargs['method'] == 'zhang':
            return computeReflectivityZhang(**kwargs)
        elif kwargs['method'] == 'jzx':
            return computeReflectivityJZX(**kwargs)
    else:
        return computeReflectivityZhang(**kwargs)

def computeReflectivityZhang(**kwargs):
    N0_rain = 4e5
    N0_snow = 3e6
    N0_hail = 4e4

    density_rain = 1000.
    density_snow = 100.
    density_hail = 913.
    density_ice = 917.

    dielectric_water = 0.93
    dielectric_ice = 0.176

    if 'pt' in kwargs:
        temperature = kwargs['pt'] * (kwargs['p'] / 100000.) ** (2. / 7.)
    elif 't' in kwargs:
        temperature = kwargs['t']

    density = kwargs['p'] / (287. * temperature)

    is_wet = np.where(temperature > 273.15, 1, 0)

    Zr_factor = 720 * 1e18 / ((np.pi * density_rain) ** 1.75 * N0_rain ** 0.75)
    Zs_wet_factor = 720 * 1e18 / ((np.pi * density_snow) ** 1.75 * N0_snow ** 0.75)
    Zs_dry_factor = 720 * 1e18 * dielectric_ice * density_snow ** 0.25 / (np.pi ** 1.75 * dielectric_water * N0_snow ** 0.75 * density_ice ** 2)
    Zh_factor = 720 * 1e18 / ((np.pi * density_hail) ** 1.75 * N0_hail ** 0.75)

    Zr = Zr_factor * (density * kwargs['qr']) ** 1.75
    Zs = Zs_wet_factor * is_wet * (density * kwargs['qs']) ** 1.75 + Zs_dry_factor * (1 - is_wet) * (density * kwargs['qs']) ** 1.75
    Zh = (Zh_factor * (density * kwargs['qh']) ** 1.75) ** 0.95

    return 10 * np.log10(Zr + Zs + Zh)

def computeReflectivityJZX(**kwargs):
    N0_rain = 4e5
    N0_snow = 3e6
    N0_hail = 4e4

    density_rain = 1000.
    density_snow = 100.
    density_hail = 913.
    density_ice = 917.

    dielectric_water = 0.93
    dielectric_ice = 0.176

    frac_max = 0.5

    if 'pt' in kwargs:
        temperature = kwargs['pt'] * (kwargs['p'] / 100000.) ** (2. / 7.)
    elif 't' in kwargs:
        temperature = kwargs['t']

    density = kwargs['p'] / (287. * temperature)

    melt_frac = frac_max * np.min(kwargs['qs'] / kwargs['qr'], kwargs['qr'] / kwargs['qs']) ** 0.3
    qr_unmix = (1 - melt_frac) * kwargs['qr']
    qs_unmix = (1 - melt_frac) * kwargs['qs']
    water_frac = kwargs['qr'] / (kwargs['qr'] + kwargs['qs'])

    density_mix = density_snow * (1 - water_frac ** 2) + density_rain * water_frac ** 2

    alpha_ra = alpha_rb = 4.28e-4
    beta_ra = 3.04
    beta_rb = 2.77

    alpha_sa = (0.194 + 7.094 * water_frac + 2.135 * water_frac ** 2 - 5.225 * water_frac ** 3) * 10 ** -4
    alpha_sb = (0.191 + 6.916 * water_frac - 2.841 * water_frac ** 2 - 1.160 * water_frac ** 3) * 10 ** -4
    alpha_ha = (0.191 + 2.390 * water_frac - 12.57 * water_frac ** 2 + 38.71 * water_frac ** 3 - 65.53 * water_frac ** 4 + 56.16 * water_frac ** 5 - 18.98 * water_frac ** 6) * 10 ** -3
    alpha_hb = (0.165 + 1.720 * water_frac - 9.920 * water_frac ** 2 + 32.15 * water_frac ** 3 - 56.00 * water_frac ** 4 + 48.83 * water_frac ** 5 - 16.69 * water_frac ** 6) * 10 ** -3
    beta_sa = beta_sb = beta_ha = beta_hb = 3

    A_coeff_r, A_coeff_s, A_coeff_h = 0., 1., 1.
    B_coeff_r, B_coeff_s, B_coeff_h = 0., 0., 0.
    C_coeff_r, C_coeff_s, C_coeff_h = 0., 0., 0.

    return

def computeVorticity(**kwargs):
    if 'u' in kwargs: u = kwargs['u']
    if 'v' in kwargs: v = kwargs['v']

    if 'grid_spacing' in kwargs:
        dx = kwargs['grid_spacing']
        dy = kwargs['grid_spacing']
    elif 'dx' in kwargs and 'dy' in kwargs:
        dx = kwargs['dx']
        dy = kwargs['dy']

    ndim = len(u.shape)
    if ndim < 2: return

    spacings = [ 1 ] * (ndim - 2)
    spacings.extend([dy, dx])

    dvdx = np.gradient(v, *spacings)[ndim - 1]
    dudy = np.gradient(u, *spacings)[ndim - 2]

    return dvdx - dudy

toRecArray = lambda **d: np.array(zip(*[ d[k].ravel() for k in sorted(d.keys()) ]), dtype=[ (k, float) for k in sorted(d.keys()) ]).reshape(d[d.keys()[0]].shape)

def theta2Temperature(**kwargs):
    pres = kwargs['p']
    pt = kwargs['pt']

    return pt * (pres / 100000.) ** (2. / 7.)

def qv2Dewpoint(**kwargs):
    pres = kwargs['p']
#   pt = kwargs['pt']
    qv = kwargs['qv']

    vapor_pres = (pres * qv) / (qv + 0.622)
    return 1. / (1. / 273.15 - 461.5 / 2.5e6 * np.log(vapor_pres / 611.))

def vectorTuple(**kwargs):
    vectors = np.empty(kwargs['u'].shape, dtype=[('u', np.float32), ('v', np.float32)])
    vectors['u'] = kwargs['u']
    vectors['v'] = kwargs['v']
    return vectors

def computePMSL(**kwargs):
    p = kwargs['p']
    z = kwargs['z']
    T = theta2Temperature(**kwargs)
    qv = kwargs['qv']

    epsilon = 0.622
    g = 9.806
    R_d = 287.

    virtual_temp_const = 1. / epsilon - 1
    T_v = (1 + virtual_temp_const * qv / (1 + qv)) * T

    return p * np.exp( ( g * z ) / ( R_d * T_v ))
