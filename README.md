# mesh2pcf

`mesh2pcf` is a Python wrapper for a high-performance C++ library that computes the two-point correlation function (2PCF) multipoles from gridded density fields. It supports both isotropic and anisotropic measurements and uses OpenMP for fast execution.

---

## Features

- Computes monopole (ℓ=0), quadrupole (ℓ=2), and hexadecapole (ℓ=4) of the 2PCF
- Supports regular 3D gridded input (e.g., from density fields)
- Anisotropic correlation function with customizable μ-binning
- Multi-threaded execution using OpenMP
- C++ core with Python interface


---

## Installation

Make sure you have the required dependencies:

- Python 3.6+
- NumPy
- GSL and OpenMP libraries (required by the C++ code)
- SWIG (for building from source)
- A C++ compiler supporting OpenMP

---

### Installing from source

1. Clone the repository:
   ```bash
   git clone https://github.com/EdoardoMaragliano/Mesh2PCF.git
   cd mesh2pcf
   python setup.py build_ext --inplace
   python setup.py install
   ```

### Usage

See the notebooks and scripts in the [`examples`](examples) folder for basic usage and demonstrations of the package.

See the following for a quick example.

```python
import os

os.environ['OMP_NUM_THREADS'] = '24'

import numpy as np
from mesh2pcf.Mesh2pcfCpp.Mesh2pcfCpp import Mesh2PCF, Mesh2PCFAniso, DoubleVector, DoubleVectorVector, DoubleGrid

mesh = np.load('real_space_field.npy')
mesh = np.array(mesh, dtype=np.float64)
c_mesh = DoubleGrid(mesh)

test_2pcf_r_mu = Mesh2PCFAniso(c_mesh, 1000)
test_2pcf_r_mu.compute_mesh2PCF(
    max_distance = 200.0,
    bin_size = 5.0,
    num_mu_bins = 200,
)
r_array = test_2pcf_r_mu.get_r_array()
xi0_r_mu_128 = test_2pcf_r_mu.get_monopole()
xi2_r_mu_128 = test_2pcf_r_mu.get_quadrupole()
xi4_r_mu_128 = test_2pcf_r_mu.get_hexadecapole()

output = np.array([r_array, xi0_r_mu_128, xi2_r_mu_128, xi4_r_mu_128])

np.save('xi_r_mu_128.npy', output)
```

## Acknowledgments
I thank A. Veropalumbo for the OMP parallelization.