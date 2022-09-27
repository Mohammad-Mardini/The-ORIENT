import os
from . import orient
data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
del os
orient.set_data_dir(data_dir)
from ._version import __version__
import numpy as np

try:
    import astropy.coordinates
    import astropy.units as u
    astropy.coordinates.galactocentric_frame_defaults.set('v4.0')
    has_astropy = True
except:
    has_astropy = False

Galaxy = orient.Galaxy

def icrs_transform(ra, dec, distance, pm_ra_cosdec, pm_dec, radial_velocity):
    icrs = astropy.coordinates.ICRS(
        ra              = ra              * u.deg,
        dec             = dec             * u.deg,
        distance        = distance        * u.kpc,
        pm_ra_cosdec    = pm_ra_cosdec    * u.mas/u.yr,
        pm_dec          = pm_dec          * u.mas/u.yr,
        radial_velocity = radial_velocity * u.km/u.s,
    )
    galactocentric = astropy.coordinates.Galactocentric()
    coords         = icrs.transform_to(galactocentric)
    result = np.empty(6)
    result[:3] = coords.data.xyz / u.kpc
    result[3:] = coords.velocity.d_xyz / (u.km/u.s)
    return result

class Integrate:
    kpcGyr = 1.0226911647958985
    def __init__(self,
                 galaxy:      Galaxy,
                 ic:          list,
                 t_min:       np.double,
                 t_max:       np.double,
                 stride_size: np.double,
                 max_size:    int = 100000,
                 kms:         bool = True,
                 icrs:        bool = False):
        """
        Construct an orbital trajectory object

        Parameters
        ----------
        galaxy : Galaxy
            A Galaxy object (ORIENT potential model)
        ic : list
            The initial conditions as a size 6 array [kpc for position and km/s for velocities].
        t_min : double
            Cosmic age at initial conditions [Gyr].
        t_max : double
            Cosmic age at end of integration.
        stride_size: double
            Record the orbit parameters at this time interval [Gyr].
        max_size: int, optional
            Maximum number of records to make during integration.
        kms: book, optional
            If true (default), velocity initial conditions are specified in km/s, otherwise kpc/Gyr.
        icrs: bool, optional
            If false (default), initial conditions specified in Cartesian coordinates, otherwise in the ICRS system (requires AstroPy).
        """
        _ic = np.array(ic, dtype=np.double)

        if icrs:
            if not has_astropy:
                raise RuntimeError('icrs==True requires astropy')
            if not kms:
                raise RuntimeError('icrs==True requires kms=True')
            _ic = icrs_transform(*_ic)

        if kms: _ic[3:] *= Integrate.kpcGyr
        self._data = orient.integrate(galaxy, _ic, t_min, t_max, stride_size, max_size)
        if kms: self._data[3:] /= Integrate.kpcGyr

        stride_count = int((t_max-t_min) / stride_size)
        if stride_count*stride_size > t_max-t_min: stride_count += 1
        self.t = np.linspace(t_min, t_max, stride_count)[:max_size]

    def __getitem__(self, key):
        return self._data[key]
    
    @property
    def x(self):  return self._data[:,0]

    @property
    def y(self):  return self._data[:,1]

    @property
    def z(self):  return self._data[:,2]

    @property
    def vx(self): return self._data[:,3]

    @property
    def vy(self): return self._data[:,4]

    @property
    def vz(self): return self._data[:,5]

    @property
    def R(self):
        return np.linalg.norm(self._data[:,:3], axis=1)

    @property
    def L(self):
        return np.cross(self._data[:,:3], self._data[:,3:])

    @property
    def Lmag(self):
        return np.linalg.norm(self.L, axis=1)
    
    @property
    def E(self):
        raise NotImplementedError

from typing import Optional
import scipy.interpolate
class Orbit:
    def __init__(self,
                 galaxy:      Galaxy,
                 ic:          list,
                 t_min:       np.double,
                 t_max:       np.double,
                 stride_size: Optional[np.double] = None,
                 max_size:    int = 100000,
                 kms:         bool = True,
                 icrs:        bool = False):
        t = None
        t_apsis = None
        if stride_size is None:
            stride_size = (t_max - t_min)/10000
            sign = np.sign(stride_size)
            stride_size = sign * 2**(np.floor(np.log(abs(stride_size))/np.log(2)))
        self.integrate_obj = Integrate(galaxy, ic, t_min, t_max, stride_size, max_size, kms, icrs)
        data = self.integrate_obj._data
        self._t = np.arange(t_min, t_max, stride_size)
        if t_min > t_max:
            self._t = self._t[::-1]
            data = data[::-1]
        self.data_interp = scipy.interpolate.interp1d(self._t, data, axis=0)

        r_vec = data[:,:3]
        x, y, z  = r_vec.T
        v_vec = data[:,3:]
        vx, vy, vz = v_vec.T
        
        r = np.linalg.norm(r_vec, axis=1)
        t = self._t

        sign_diff = np.diff(np.sign(np.diff(r)))
        peri_mask = sign_diff ==  2
        apo_mask  = sign_diff == -2
        self._t_peri = t[1:-1][peri_mask]
        self._r_peri = r[1:-1][peri_mask]
        self._t_apo  = t[1:-1][apo_mask]
        self._r_apo  = r[1:-1][apo_mask]

        zmax_mask = np.zeros(len(t), dtype=bool)
        params = np.array(galaxy.get_fit_params(t))
        galaxy_phi   = params[:,0]
        galaxy_theta = params[:,1]
        galaxy_direction = np.array([np.cos(galaxy_phi)*np.sin(galaxy_theta),
                                     np.sin(galaxy_phi)*np.sin(galaxy_theta),
                                     np.cos(galaxy_theta)]).T
        z_transformed = np.sum((r_vec * galaxy_direction[None,:])[0], axis=1)
        zmax_mask[1:-1] = np.diff(np.sign(np.diff(z_transformed)))!=0
        self._zmax = np.abs(z_transformed[zmax_mask])
        self._t_zmax = t[zmax_mask]

        t_apsis = np.empty((len(self._t_peri) + len(self._t_apo)))
        r_apsis = np.empty_like(t_apsis)
        if self._t_peri[0] < self._t_apo[0]:
            t_apsis[::2] = self._t_peri
            r_apsis[::2] = self._r_peri
            t_apsis[1::2] = self._t_apo
            r_apsis[1::2] = self._r_apo
        else:
            t_apsis[::2] = self._t_apo
            r_apsis[::2] = self._r_apo
            t_apsis[1::2] = self._t_peri
            r_apsis[1::2] = self._r_peri

        a = .5*(r_apsis[:-1] + r_apsis[1:])
        e = np.abs(r_apsis[:-1] - r_apsis[1:])/(2*a)
            
        self.e = scipy.interpolate.interp1d(.5*(t_apsis[:-1] + t_apsis[1:]), e, bounds_error=False)
        self.zmax = scipy.interpolate.interp1d(self._t_zmax, self._zmax, bounds_error=False)
        self.r_apo = scipy.interpolate.interp1d(self._t_apo, self._r_apo, bounds_error=False)

    @property
    def x(self):  return self.integrate_obj._data[:,0]
    
    @property
    def y(self):  return self.integrate_obj._data[:,1]
    
    @property
    def z(self):  return self.integrate_obj._data[:,2]
    
    @property
    def vx(self): return self.integrate_obj._data[:,3]
    
    @property
    def vy(self): return self.integrate_obj._data[:,4]
    
    @property
    def vz(self): return self.integrate_obj._data[:,5]