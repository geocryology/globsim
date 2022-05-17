import re
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = re.search(
    '^__version__\\s*=\\s*"(.*)"',
    open('globsim/_version.py').read(),
    re.M
    ).group(1)

setuptools.setup(
    name="globsim",
    version=version,
    author="Stephan Gruber",
    author_email="stephan.gruber@carleton.ca",
    description="Using global reanalysis data for local permafrost simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/geocryology/globsim",
    packages=['globsim',
              'globsim.download',
              'globsim.interpolate',
              'globsim.scale'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Scientific/Engineering :: Atmospheric Science"
        ],
    entry_points={
        'console_scripts': [
            'globsim = globsim.globsim_cli:main',
        ]},

    install_requires=['lxml',
                      'numpy',
                      'pandas',
                      'netCDF4',
                      'scipy',
                      'pydap', 
                      'tomlkit',
                      'ecmwf-api-client',
                      'cdsapi',
                      'pysolar']
)
