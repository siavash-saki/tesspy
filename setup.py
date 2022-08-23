import setuptools

with open("README.md", "r", encoding="utf8", errors="ignore") as fh:
    LONG_DESCRIPTION = fh.read()

# Packages that tesspy uses explicitly:
INSTALL_REQUIRES = [
    "geopandas >= 0.10.0",
    "scipy",
    "h3",  # conda h3-py
    "osmnx",
    "hdbscan",
    "mercantile",
    "scikit-learn",
]

setuptools.setup(
    name="tesspy",
    version="0.1.2",
    author="tesspy Developers",
    author_email="jonas.hamann@fb3.fra-uas.de",
    description="Tessellation of Urban Areas",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/siavash-saki/tesspy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: BSD License",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: GIS",
    ],
    python_requires=">=3.7",
    install_requires=INSTALL_REQUIRES,
)
