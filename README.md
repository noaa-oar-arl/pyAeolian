# pyAeolian
 A python library to interface with aelian dust paramterizations

 ### To Install

 The easiest way to install is using pip

 ```cd``` into the pyAeolian directory where the setup.py file exists

 Then execute the command:

 ```
 pip install .
 ```

 ### To use

 Simple example on how to use.  First import `pyaeolian`

 ```python

 import pyaeolian

 ```

 now all subroutines in the aeolian.F90 file is available to you in python

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


 you can even wrap this with xarray.apply_ufunc to apply to a grid (without python loops [ really it is just pure C loops ] , ie very fast) here zz is an xarray.Dataset with the clay sand silt threshold drag partition and sediment supply map and cc is the CMAQ meteorology containing ustar and soil moisture.  Note that we just use constants for land and rhoa here.


 ```python
 def xarray_fengsha(rhoa, volumetric_soil_moisture, ssm, land, ustar, clayfrac, sandfrac, drag_partition, dry_threshold):
     func = lambda a,b,c,d,e,f,g,h,i: pyfengsha.fengsha(a,b,c,d,e,f,g,h,i)
     return xr.apply_ufunc(func,rhoa, volumetric_soil_moisture, ssm, land, ustar, clayfrac, sandfrac, drag_partition, dry_threshold,vectorize=True)

 result = xarray_fengsha(zz.SSM*0+1000,c.SOIM1.isel(time=40),zz.SSM,zz.SSM*0+1,c.USTAR.isel(time=40),zz.CLAY_FRAC,zz.SAND_FRAC,zz.DRAG_PART,zz.UTHRES)

 ```
