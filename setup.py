from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy as np
import sys

cpp_root = "mesh2pcf/Mesh2pcfCpp"
cpp_source = f"{cpp_root}/src/"
cpp_inc = f"{cpp_root}/include/"

clustering_stats_cpp_module = Extension(
    "mesh2pcf.Mesh2pcfCpp._Mesh2pcfCpp",
    sources=[f"{cpp_source}/Mesh2pcfCpp.cpp",
             f"{cpp_source}/Mesh2pcfCpp.i"],
    include_dirs = [cpp_inc, np.get_include()],
    
    swig_opts=["-threads",
            "-c++",
            f"-I{cpp_inc}",
            "-oh",
            cpp_inc+"Mesh2pcfCpp_wrap.h",
            "-outdir",
            cpp_root],
    extra_compile_args = ["-O3", "-fopenmp"],
    libraries=["gsl", "gslcblas", "gomp"]
)


# Carica il README in sicurezza
try:
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = ""

setup(
    name="mesh2pcf",
    version="0.1.2",
    author="Edoardo Maragliano",
    author_email="edoardo.maragliano@edu.unige.it",
    description="A simple package to compute the 2pcf on a gridded field.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EdoardoMaragliano/Mesh2PCF.git",
    packages=find_packages(),
    ext_modules=[clustering_stats_cpp_module],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],
    python_requires=">=3.6",
    install_requires=["numpy"],
)
