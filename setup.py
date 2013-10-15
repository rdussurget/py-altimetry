from distutils.core import setup

setup(
    name='py-altimetry',
    version='0.3.0',
    author='R. Dussurget',
    author_email='renaud.dussurget@gmail.com',
    packages=['altimetry','altimetry.config','altimetry.data','altimetry.tools','altimetry.externals'],
    scripts=['bin/cdf_merge.py','altimetry/examples/spectral_tapering.py','altimetry/examples/j2_spectrum.py'],
    url='https://code.google.com/p/py-altiwaves/',
    license='LICENSE.TXT',
    description='PyAltimetry: Python tools for altimetry data and much more',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "NetCDF4"
    ],
)