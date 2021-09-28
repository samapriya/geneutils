from setuptools import setup
from setuptools import find_packages

setup(
    name="geneutils",
    version="0.0.5",
    packages=["geneutils"],
    data_files=[("", ["LICENSE"])],
    url="https://github.com/samapriya/geneutils",
    install_requires=["biopython>=1.77","pandas>=1.1.5","requests>=2.26.0","beautifulsoup4>=4.9.3"],
    license="MIT",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent"
    ],
    author="Samapriya Roy",
    author_email="samapriya.roy@gmail.com",
    description="CLI and utilities for Genetic analysis and database interface",
    entry_points={"console_scripts": ["geneutils=geneutils.geneutils:main",],},
)
