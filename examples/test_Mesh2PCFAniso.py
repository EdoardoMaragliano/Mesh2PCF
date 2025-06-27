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