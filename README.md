# mesh2pcf

A simple package to measure the 2-point correlation function (2PCF) of a gridded density field, using a C++ backend wrapped with SWIG for Python.

---

## Installation

Make sure you have the required dependencies:

- Python 3.6+
- NumPy
- GSL and OpenMP libraries (required by the C++ code)
- SWIG (for building from source)
- A C++ compiler supporting OpenMP

### Installing from source

1. Clone the repository:
   ```bash
   git clone https://github.com/EdoardoMaragliano/Mesh2PCF.git
   cd mesh2pcf
   python setup.py build_ext --inplace
   python setup.py install

### Usage

See the notebooks and scripts in the [`examples`](examples) folder for basic usage and demonstrations of the package.

# Mesh2PCF
