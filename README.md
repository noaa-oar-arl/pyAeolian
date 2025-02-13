# pyAeolian

A Python library to interface with [aeolian](https://en.wikipedia.org/wiki/Aeolian_processes)
dust parameterizations

### Installation

The easiest way to install is using `pip`.

`cd` into the `pyAeolian` directory where the `setup.py` file resides.

Then execute the command:

```
pip install .
```

### Usage

Simple example on how to use.  First import `pyaeolian`.

```python

import pyaeolian

```

Now all subroutines in the `aeolian.F90` file are available to you in Python.

```python
Help on package pyaeolian:

NAME
    pyaeolian

PACKAGE CONTENTS
    aeolian

FUNCTIONS
    draxler_hflux(ustar, threshold_velocity)

    fecan_dry_limit(clayfrac)

    fecan_moisture_correction(volumetric_soil_moisture, sandfrac, clayfrac)

    fengsha(rhoa, volumetric_soil_moisture, ssm, land, ustar, clayfrac, sandfrac, drag_partition, dry_threshold)

    fengsha_albedo(rhoa, volumetric_soil_moisture, ssm, land, ustar, clayfrac, sandfrac, drag_partition, dry_threshold)

    mackinnon_drag_partition(z0)

    mb95_drag_partition(z0)

    mb95_kvh(clay)

    modified_threshold_velocity(dry_threshold, moisture_correction, drag_partition)

    shao_1996_soil_moisture(volumetric_soil_moisture)

    shao_2004_soil_moisture(volumetric_soil_moisture)

    volumetric_to_gravimetric(volumetric_soil_moisture, sandfrac)

```


You can even wrap this with `xarray.apply_ufunc` to apply to a grid (without Python loops [ really it is just pure C loops ] , ie very fast) here `zz` is an `xarray.Dataset` with the clay sand silt threshold drag partition and sediment supply map, and `cc` is the CMAQ meteorology containing ustar and soil moisture.  Note that we just use constants for land and rhoa here.

```python
def xarray_fengsha(rhoa, volumetric_soil_moisture, ssm, land, ustar, clayfrac, sandfrac, drag_partition, dry_threshold):
    func = lambda a,b,c,d,e,f,g,h,i: pyfengsha.fengsha(a,b,c,d,e,f,g,h,i)
    return xr.apply_ufunc(func,rhoa, volumetric_soil_moisture, ssm, land, ustar, clayfrac, sandfrac, drag_partition, dry_threshold,vectorize=True)

result = xarray_fengsha(zz.SSM*0+1000,c.SOIM1.isel(time=40),zz.SSM,zz.SSM*0+1,c.USTAR.isel(time=40),zz.CLAY_FRAC,zz.SAND_FRAC,zz.DRAG_PART,zz.UTHRES)

```

### Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
