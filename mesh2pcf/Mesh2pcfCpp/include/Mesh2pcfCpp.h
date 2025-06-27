#ifndef MESH2PCFCPP_H
#define MESH2PCFCPP_H

#include <vector>
#include <tuple>

class Mesh2PCF {
public:
    // Constructor
    Mesh2PCF(const std::vector<std::vector<std::vector<double>>>& density_field, double box_size);

    // Method to compute the 2PCF
    void compute_mesh2PCF(double max_distance, double bin_size);
    std::vector<double> get_xi() const;
    std::vector<double> get_r_array() const;
    std::vector<double> get_counts() const;
    std::vector<double> get_prods() const;
    void set_progress(bool progress);

private:
    // Member variables
    std::vector<std::vector<std::vector<double>>> density_grid;
    double box_size;
    int nmesh;
    double cell_size;
    double max_distance;
    double bin_size;
    int num_bins;
    std::vector<double> counts;
    std::vector<double> prods;
    std::vector<double> r_array;
    std::vector<double> xi;
    std::vector<std::tuple<int, int, int, double>> offsets;
    bool progress;

    // Helper methods
    void precalculate_offsets();
    std::pair<std::vector<double>, std::vector<double>> compute_density_products(int i_c, int j_c, int k_c);
};

// Anisotropic Mesh2PCF
class Mesh2PCFAniso {
public:
    // Constructor
    Mesh2PCFAniso(const std::vector<std::vector<std::vector<double>>>& density_grid, double box_size);

    // Methods
    void compute_mesh2PCF(double max_distance, double bin_size, int num_mu_bins);
    std::vector<double> get_r_array() const;
    std::vector<double> get_mu_array() const; 
    std::vector<std::vector<double>> get_counts() const;
    std::vector<std::vector<double>> get_prods() const;
    std::vector<double> get_monopole() const;
    std::vector<double> get_quadrupole() const;
    std::vector<double> get_hexadecapole() const;
    void set_progress(bool progress);
        
private:
    // Member variables
    int nmesh;
    double box_size, max_distance, bin_size, cell_size;
    int num_bins, num_mu_bins;
    std::vector<std::vector<std::vector<double>>> density_grid;
    std::vector<std::vector<double>> counts, prods;
    std::vector<double> r_array, mu_array, xi_norm, monopole, quadrupole, hexadecapole;
    std::vector<std::tuple<int, int, int, double>> offsets;
    bool progress;

    // Helper methods
    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> compute_density_products(int i_c, int j_c, int k_c);
    void precalculate_offsets();
    std::vector<double> compute_monopole() const;
    std::vector<double> compute_quadrupole() const;
    std::vector<double> compute_hexadecapole() const;
};

#endif // MESH2PCFCPP_H
