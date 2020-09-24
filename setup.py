import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="globsim",
    version="0.0.1",
    author="Stephan Gruber",
    author_email="stephan.gruber@carleton.ca",
    description="Using global reanalysis data for local permafrost simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/geocryology/globsim",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Scientific/Engineering :: Atmospheric Science"
        
    ],
    install_requires=[
                    'numpy',
                    'tomlkit',
                    'netCDF4'
    ]
)