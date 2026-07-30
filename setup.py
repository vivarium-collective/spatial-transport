import re
from setuptools import setup


VERSION = '0.0.1'


with open("README.md", "r") as readme:
    description = readme.read()
    # Patch the relative links to absolute URLs that will work on PyPI.
    description2 = re.sub(
        r']\(([\w/.-]+\.png)\)',
        r'](https://github.com/vivarium-collective/spatial-transport/raw/main/\1)',
        description)
    long_description = re.sub(
        r']\(([\w/.-]+)\)',
        r'](https://github.com/vivarium-collective/spatial-transport/blob/main/\1)',
        description2)

setup(
    name="spatial-transport",
    version=VERSION,
    author="Tasnif Rahman",
    author_email="trahman@uchc.edu",
    description="Spatial transport processes for process-bigraph composites",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/vivarium-collective/spatial-transport",
    packages=[
        'spatial_transport',
        'spatial_transport.processes',
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.12",
    install_requires=[
        "cdFBA",
        # the schema types used here (dataclass `Node` types, `overwrite`)
        # landed in these versions
        "process_bigraph>=1.5.0",
        "bigraph-schema>=1.4.3",
        "numpy",
        "pandas",
        "seaborn",
        "matplotlib",
        "imageio",
        "pytest",
    ],
    extras_require={
        # only needed by spatial_transport.processes.tyssue_diffusion
        "tyssue": ["tyssue"],
    },
)
