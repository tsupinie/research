import numpy as np
from scipy.integrate import trapz
from scipy.stats import nanmean

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

def coordDerivative(data, coords, axis=0):
    def cta(coord):
        return tuple(([ slice(None) ] * axis) + [ coord ])

    l_slc = slice(0, -2)
    c_slc = slice(1, -1)
    r_slc = slice(2,  None)

    derivative = np.zeros(data.shape)
    alpha = (coords[cta(r_slc)] - coords[cta(c_slc)]) / (coords[cta(c_slc)] - coords[cta(l_slc)])

    derivative[cta(c_slc)] = (data[cta(r_slc)] + (alpha - 1) * data[cta(c_slc)] - alpha * data[cta(l_slc)]) / (2 * alpha * (coords[cta(c_slc)] - coords[cta(l_slc)]))
    derivative[cta(0)] = (data[cta(1)] - data[cta(0)]) / (coords[cta(1)] - coords[cta(0)])
    derivative[cta(-1)] = (data[cta(-1)] - data[cta(-2)]) / (coords[cta(-1)] - coords[cta(-2)])

    return derivative

def compute3DVorticity(**kwargs):
    if 'vg_tensor' not in kwargs:
        vg_tensor = computeVGTensor(**kwargs)
    else:
        vg_tensor = kwargs['vg_tensor']

    vort3d = np.empty(vg_tensor.shape, dtype=[('zvort', np.float32), ('yvort', np.float32), ('xvort', np.float32)])

    vort3d['zvort'] = vg_tensor['dvdx'] - vg_tensor['dudy']
    vort3d['yvort'] = vg_tensor['dudz'] - vg_tensor['dwdx']
    vort3d['xvort'] = vg_tensor['dwdy'] - vg_tensor['dvdz']
    return vort3d

def computeVGTensor(**kwargs):
    dtype = []
    for wc in ['u', 'v', 'w']:
        for ax in ['x', 'y', 'z']:
            dtype.append(("d%sd%s" % (wc, ax), np.float32))

    vg_tensor = np.empty(kwargs['u'].shape, dtype=dtype)

    nz, ny, nx = kwargs['z'].shape
    kwargs['x'] = np.atleast_3d(kwargs['x']).reshape((1, 1, nx))
    kwargs['y'] = np.atleast_3d(kwargs['y']).reshape((1, ny, 1))

    for wc in ['u', 'v', 'w']:
        for ax_no, ax in enumerate(['z', 'y', 'x']):
            vg_tensor["d%sd%s" % (wc, ax)] = coordDerivative(kwargs[wc], kwargs[ax], axis=ax_no)
    return vg_tensor

def computeBuoyancy(**kwargs):
    hydro_mass = kwargs['qc'] + kwargs['qi'] + kwargs['qr'] + kwargs['qs'] + kwargs['qh']
    theta_rho = kwargs['pt'] * (1 + kwargs['qv'] / 0.622) / (1 + hydro_mass)
    theta_rho_bar = nanmean(nanmean(theta_rho, axis=-1), axis=-1)[..., np.newaxis, np.newaxis]
    p_bar = nanmean(nanmean(kwargs['p'], axis=-1), axis=-1)[..., np.newaxis, np.newaxis]    

    return 9.806 * ((theta_rho - theta_rho_bar) / theta_rho_bar + (2. / 7. - 1) * (kwargs['p'] - p_bar) / p_bar)

def computeBaroclinicVortGen(**kwargs):
#   baroc_vort = np.empty(kwargs['pt'].shape, dtype=[('bvgx', np.float32), ('bvgy', np.float32)])

    buoyancy = computeBuoyancy(**kwargs)

    bvgx =  coordDerivative(buoyancy, kwargs['y'][:, np.newaxis], axis=0)
    bvgy = -coordDerivative(buoyancy, kwargs['x'][np.newaxis, :], axis=1)

    baroc_vort = (kwargs['u'] * bvgx + kwargs['v'] * bvgy) / np.hypot(kwargs['u'], kwargs['v'])
    return baroc_vort

def computeVH(**kwargs):
    vert_vort = computeVorticity(**kwargs)
    vert_hel = vert_vort * kwargs['w']

    return vert_hel

toRecArray = lambda **d: np.array(zip(*[ d[k].ravel() for k in sorted(d.keys()) ]), dtype=[ (k, float) for k in sorted(d.keys()) ]).reshape(d[d.keys()[0]].shape)
toDict = lambda r: dict( (f, r[f]) for f in r.dtype.fields.iterkeys() )

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
