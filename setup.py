from numpy.distutils.core import Extension

ext2 = Extension(name="aeolian", sources=["pyaeolian/aeolian.pyf", "pyaeolian/aeolian.F90"])

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(
        name="pyAeolian",
        version="0.1",
        description="Wrapper around the NOAA ARL aeolian developmental parameterizations",
        author="Barry D. Baker",
        lisense="MIT",
        author_email="barry.baker@noaa.gov",
        source=["pyaeolian"],
        packages=["pyaeolian"],
        ext_modules=[ext2],
        install_requires=["numpy"],
    )
