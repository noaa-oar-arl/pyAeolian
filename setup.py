from numpy.distutils.core import Extension

ext = Extension(name="aeolian", sources=["pyaeolian/aeolian.pyf", "pyaeolian/aeolian.F90"])

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(
        name="pyAeolian",
        version="0.1",
        description="Wrapper around the NOAA ARL aeolian dust developmental parameterizations",
        author="Barry D. Baker",
        lisense="MIT",
        author_email="barry.baker@noaa.gov",
        packages=["pyaeolian"],
        ext_modules=[ext],
        install_requires=["numpy", "xarray"],
    )
