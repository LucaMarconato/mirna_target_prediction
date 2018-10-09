#ifndef MATCHINGS_PREDICTOR_H
#define MATCHINGS_PREDICTOR_H

#include <utility>

#include "patient.hpp"

class Matchings_predictor {
    Patient & patient;
    std::unordered_map<Mirna_id, double> mirna_profile;
    std::unordered_map<Cluster *, double> cluster_profile;
    std::unordered_map<Cluster *, double> original_cluster_profile;
    // r_ic_values is described in the attached document
    std::unordered_map<std::pair<Mirna_id, Cluster *>, double> r_ic_values;
    std::unordered_map<std::pair<Mirna_id, Site *>, double> r_ijk_values;
    std::unordered_map<Cluster *, double> p_c_bound_values;
    std::unordered_map<Cluster *, double> sum_of_r_ijk_for_cluster;
    std::unordered_map<std::pair<Gene_id, Cluster *>, double> p_j_downregulated_given_c_bound_values;
    std::unordered_map<Gene_id, double> p_j_downregulated_values;

    void export_interaction_matrix();
    void compute_probabilities();
    void export_probabilities();
    void export_p_j_downregulated(int filename_suffix = 999999);
    inline void recusively_compute_p_j_downregulated(bool * b, int level, int max_level, double * sum, double p_j_downregulated_given_b, double p_b, double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened);
    inline double iteratively_compute_p_j_downregulated(double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened, int clusters_count);
public:
    Matchings_predictor(Patient & patient);
    void compute();
};

#endif // MATCHINGS_PREDICTOR_H
