#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>    // OpenMP for parallel computing, thanks to ChatGPT
#include "Mesh2pcfCpp.h"
#include <ctime>

using namespace std;

Mesh2PCF::Mesh2PCF(const vector<vector<vector<double>>>& density_field, const double box_size)
        : density_grid(density_field), box_size(box_size) {
        
        nmesh = density_field.size();
        cell_size = box_size / nmesh;

        cout << "Mesh2PCF object created" << endl;
        cout << "Density field size: " << nmesh << endl;
        cout << "Box size: " << box_size << endl;

        bool progress = true;
        cout << "Progress: " << progress << endl;
    }


void Mesh2PCF::set_progress(bool progress) {
    this->progress = progress;
}

vector<double> Mesh2PCF::get_xi() const {
    return xi;
}

vector<double> Mesh2PCF::get_r_array() const {
    return r_array;
}

vector<double> Mesh2PCF::get_counts() const {
    return counts;
}

vector<double> Mesh2PCF::get_prods() const {
    return prods;
}


void Mesh2PCF::compute_mesh2PCF(const double max_distance, const double bin_size) {
    this->max_distance = max_distance;
    this->bin_size = bin_size;
    num_bins = static_cast<int>(max_distance / bin_size) + 1;
    time_t start; time(&start); // start the time count

    // Initialize the counts and products
    counts.assign(num_bins, 0);
    prods.assign(num_bins, 0);

    precalculate_offsets();
    cout << "Offsets stored. Entering cycle on the 3D mesh..." << endl;

    int tid = 0;
    double fact_count = 100.0 / (nmesh * nmesh * nmesh);

    #pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
    {

        tid = omp_get_thread_num();
        vector<double> _counts(num_bins, 0);
        vector<double> _prod(num_bins, 0);

        // Compute the density products for each cell
        #pragma omp for schedule(dynamic)
        for (int index = 0; index < nmesh * nmesh * nmesh; ++index) {
            int i_c = index / (nmesh * nmesh); // Compute the i index
            int j_c = (index / nmesh) % nmesh; // Compute the j index
            int k_c = index % nmesh;            // Compute the k index

            vector<double> current_counts(num_bins, 0);
            vector<double> current_prod(num_bins, 0);

            tie(current_counts, current_prod) = compute_density_products(i_c, j_c, k_c);

            // Accumulate the counts and products
            for (int bin = 0; bin < num_bins; ++bin) {
                _counts[bin] += current_counts[bin];
                _prod[bin] += current_prod[bin];
            }

            // estimate the computational time and update the time count
            time_t end_temp; time(&end_temp); double diff_temp = difftime(end_temp, start);
            if (tid==0 && progress) { cout << "\r" << float(index)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }       
        }

        // Accumulate the counts and products
        #pragma omp critical
        {
            for (int bin = 0; bin < num_bins; ++bin) {
                counts[bin] += _counts[bin];
                prods[bin] += _prod[bin];
            }
        }
    }
    // Compute the correlation function
    xi.assign(num_bins, 0);
    for (int bin = 0; bin < num_bins; ++bin) {
        // Avoid division by zero
        // condition ? value_if_true : value_if_false
        xi[bin] = (counts[bin] != 0) ? (prods[bin] / counts[bin]) : 0;
    }

    // Create the separation array
    r_array.resize(num_bins);
    for (int bin = 0; bin < num_bins; ++bin) {
        r_array[bin] = (bin + 0.5) * bin_size; // Average distance for bin
    }
}

void Mesh2PCF::precalculate_offsets() {
    // Calculate the maximum index distance for a given max_distance
    int max_index_dist = static_cast<int>(ceil(max_distance / cell_size));

    for (int di = -max_index_dist; di <= max_index_dist; ++di) {
        for (int dj = -max_index_dist; dj <= max_index_dist; ++dj) {
            for (int dk = -max_index_dist; dk <= max_index_dist; ++dk) {
                // Calculate the distance
                double dist = cell_size * sqrt(di * di + dj * dj + dk * dk);

                if (dist <= max_distance && (di != 0 || dj != 0 || dk != 0)) {
                    // Store the offset and distance
                    offsets.emplace_back(di, dj, dk, dist); 
                }
            }
        }
    }
}

pair<vector<double>, vector<double>> Mesh2PCF::compute_density_products(int i_c, int j_c, int k_c) {
    vector<double> distance_bins(num_bins, 0);
    vector<double> count_bins(num_bins, 0);

    double density_central = density_grid[i_c][j_c][k_c];

    for (const auto& offset : offsets) {
        int di, dj, dk;
        double dist;
        tie(di, dj, dk, dist) = offset;

        // get the neighbor cell indices
        int ni = (i_c + di + nmesh) % nmesh; // Periodic boundary condition
        int nj = (j_c + dj + nmesh) % nmesh;
        int nk = (k_c + dk + nmesh) % nmesh;

        // Get the product of the densities
        double density_neighbor = density_grid[ni][nj][nk];
        double product = density_central * density_neighbor;

        int bin_index = static_cast<int>(dist / bin_size);
        if (bin_index < num_bins) {
            distance_bins[bin_index] += product;
            count_bins[bin_index] += 1;
        }
    }
    return {count_bins, distance_bins};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Anisotropic Mesh2PCF

Mesh2PCFAniso::Mesh2PCFAniso(const vector<vector<vector<double>>>& density_grid, double box_size)
        : density_grid(density_grid) {
        
        this->box_size = box_size;
        this->nmesh = density_grid.size();
        this->cell_size = box_size / nmesh;

        cout << "Mesh2PCFAniso object created" << endl;
        cout << "Density field size: " << nmesh << endl;
        cout << "Box size: " << box_size << endl;

        bool progress = true;
        cout << "Progress: " << progress << endl;
    }

void Mesh2PCFAniso::set_progress(bool progress) {
    this->progress = progress;
}

void Mesh2PCFAniso::compute_mesh2PCF(double max_distance, double bin_size, int num_mu_bins) {
    this->max_distance = max_distance;
    this->bin_size = bin_size;
    this->num_mu_bins = num_mu_bins;
    num_bins = static_cast<int>(max_distance / bin_size) + 1;

    // Create the separation arrays
    cout << "Creating separation arrays" << endl;
    r_array.resize(num_bins);
    mu_array.resize(num_mu_bins);
    for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
        r_array[r_bin] = (r_bin + 0.5) * bin_size;
    }
    for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
        mu_array[mu_bin] = -1.0 + (mu_bin + 0.5) * (2.0 / num_mu_bins);
    }

    time_t start; time(&start); // start the time count

    counts.assign(num_bins, vector<double>(num_mu_bins, 0));
    prods.assign(num_bins, vector<double>(num_mu_bins, 0));
    cout << "Counts and products initialized. Shape is " << num_bins << "x" << num_mu_bins << endl;
    precalculate_offsets();
    cout << "Offsets stored. Entering cycle on the 3D mesh..." << endl;

    int tid = 0;
    const double fact_count = 100.0 / (nmesh * nmesh * nmesh);
    #pragma omp parallel num_threads(omp_get_max_threads()) private(tid)
    {
        tid = omp_get_thread_num();

        //define an accumulator for each thread
        vector<vector<double>> _counts(num_bins, vector<double>(num_mu_bins, 0));
        vector<vector<double>> _prods(num_bins, vector<double>(num_mu_bins, 0));

        #pragma omp for schedule(dynamic)
        for (int index = 0; index < nmesh * nmesh * nmesh; ++index) {
            int i_c = index / (nmesh * nmesh);
            int j_c = (index / nmesh) % nmesh;
            int k_c = index % nmesh;

            // Compute the density products for each cell
            // The _counts and _prods arrays are then accumulated in the counts and prods arrays
            // Counts and prods are binned in r and mu bins
            vector<vector<double>> current_counts(num_bins, vector<double>(num_mu_bins, 0));
            vector<vector<double>> current_prods(num_bins, vector<double>(num_mu_bins, 0));
            tie(current_counts, current_prods) = compute_density_products(i_c, j_c, k_c);

            for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
                for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
                    _counts[r_bin][mu_bin] += current_counts[r_bin][mu_bin];
                    _prods[r_bin][mu_bin] += current_prods[r_bin][mu_bin];
                }
            }
            // estimate the computational time and update the time count
            time_t end_temp; time(&end_temp); double diff_temp = difftime(end_temp, start);
            if (tid==0 && progress) { cout << "\r" << float(index)*fact_count << "% completed (" << diff_temp << " seconds)\r"; cout.flush(); }
        }

        #pragma omp critical
        {
            for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
                for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
                    counts[r_bin][mu_bin] += _counts[r_bin][mu_bin];
                    prods[r_bin][mu_bin] += _prods[r_bin][mu_bin];
                }
            }
        }
        
    }

    // Compute the normalization factor, integrating counts over mu
    cout << "Computing normalization factor" << endl;
    xi_norm.assign(num_bins, 0);
    for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
        for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
            xi_norm[r_bin] += counts[r_bin][mu_bin];
        }
    }

    // compute the multipoles
    cout << "Computing multipoles" << endl;
    monopole = compute_monopole();
    quadrupole = compute_quadrupole();
    hexadecapole = compute_hexadecapole();

    cout << "Mesh2PCFAniso computation completed" << endl;
}

pair<vector<vector<double>>, vector<vector<double>>> Mesh2PCFAniso::compute_density_products(int i_c, int j_c, int k_c) {

    // Initialize the count and distance bins, with shape (num_bins, num_mu_bins)
    vector<vector<double>> count_bins(num_bins, vector<double>(num_mu_bins, 0));
    vector<vector<double>> distance_bins(num_bins, vector<double>(num_mu_bins, 0));

    // Get the central density
    double density_central = density_grid[i_c][j_c][k_c];

    for (const auto& offset : offsets) {
        // Get the offset values
        int di, dj, dk;
        double dist;
        tie(di, dj, dk, dist) = offset;

        // Get the neighbor cell indices, with periodic boundary conditions
        int ni = (i_c + di + nmesh) % nmesh;
        int nj = (j_c + dj + nmesh) % nmesh;
        int nk = (k_c + dk + nmesh) % nmesh;

        // Get the product of the densities
        double density_neighbor = density_grid[ni][nj][nk];
        double product = density_central * density_neighbor;

        // Compute the distance and mu bins
        int r_bin = static_cast<int>(dist / bin_size);
        if (r_bin >= num_bins) continue;

        if (r_bin < num_bins) {
            // Compute the mu value
            double mu = dk * cell_size / dist;
            // Compute the mu bin index
            int mu_bin = static_cast<int>((mu + 1.0) / 2.0 * num_mu_bins);
            // Ensure the mu bin is within bounds
            mu_bin = std::min(mu_bin, num_mu_bins - 1);
            if (mu_bin >= 0 && mu_bin < num_mu_bins) {
                // Accumulate the product and count
                distance_bins[r_bin][mu_bin] += product;
                count_bins[r_bin][mu_bin] += 1;
            }
        }        
    }
    return {count_bins, distance_bins};
}

void Mesh2PCFAniso::precalculate_offsets() {
    offsets.clear();
    for (int di = -nmesh / 2; di <= nmesh / 2; ++di) {
        for (int dj = -nmesh / 2; dj <= nmesh / 2; ++dj) {
            for (int dk = -nmesh / 2; dk <= nmesh / 2; ++dk) {
                double dist = sqrt(static_cast<double>(di * di + dj * dj + dk * dk)) * cell_size;

                if (dist > 0 && dist <= max_distance) {
                    offsets.emplace_back(di, dj, dk, dist);
                }
            }
        }
    }
}

vector<double> Mesh2PCFAniso::get_r_array() const {
    return r_array;
}

vector<double> Mesh2PCFAniso::get_mu_array() const {
    return mu_array;
}

vector<vector<double>> Mesh2PCFAniso::get_counts() const {
    return counts;
}

vector<vector<double>> Mesh2PCFAniso::get_prods() const {
    return prods;
}

// note about the multipoles: the normalization factor is computed as the integral of the counts over mu
// the multipoles are then computed as the sum of the products over mu, divided by the normalization factor
// this only works if RR does not depend on mu, which is the case here for Quiojote boxes, but not in general
// the general way is to compute products(r, mu) / counts(r, mu) and then integrate over mu

vector<double> Mesh2PCFAniso::compute_monopole() const {
    vector<double> monopole(num_bins, 0);
    for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
        for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
            monopole[r_bin] += prods[r_bin][mu_bin];
        }
        if (xi_norm[r_bin] != 0) {
            monopole[r_bin] /= xi_norm[r_bin];
        }
    }
    return monopole;
}

vector<double> Mesh2PCFAniso::compute_quadrupole() const {
    vector<double> quadrupole(num_bins, 0);
    for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
        for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
            quadrupole[r_bin] += prods[r_bin][mu_bin] * 0.5 * (3.0 * mu_array[mu_bin] * mu_array[mu_bin] - 1.0);
        }
        if (xi_norm[r_bin] != 0) {
            quadrupole[r_bin] /= xi_norm[r_bin];
        }
        quadrupole[r_bin] *= 5.0;   // 2l+1 normalization factor
    }
    return quadrupole;
}

vector<double> Mesh2PCFAniso::compute_hexadecapole() const {
    vector<double> hexadecapole(num_bins, 0);
    for (int r_bin = 0; r_bin < num_bins; ++r_bin) {
        for (int mu_bin = 0; mu_bin < num_mu_bins; ++mu_bin) {
            hexadecapole[r_bin] += prods[r_bin][mu_bin] * 0.125 * (35.0 * pow(mu_array[mu_bin], 4) - 30.0 * mu_array[mu_bin] * mu_array[mu_bin] + 3.0);
        }
        if (xi_norm[r_bin] != 0) {
            hexadecapole[r_bin] /= xi_norm[r_bin];
        }
        hexadecapole[r_bin] *= 9.0; // 2l+1     normalization factor
    }
    return hexadecapole;
}

vector<double> Mesh2PCFAniso::get_monopole() const {
    return monopole;
}

vector<double> Mesh2PCFAniso::get_quadrupole() const {
    return quadrupole;
}

vector<double> Mesh2PCFAniso::get_hexadecapole() const {
    return hexadecapole;
}