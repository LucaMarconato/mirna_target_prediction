#ifndef MATCHINGS_PREDICTOR_H
#define MATCHINGS_PREDICTOR_H

#include <utility>

#include "patient.hpp"

class Perturbation_analyzer;
struct Binding_flattened;

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
    std::string simulation_id;
    std::list<Gene_id> genes_skipped_by_the_distance_based_predictor;

    std::string get_output_path();
    void export_interaction_matrix();
    void compute_probabilities();
    void export_probabilities();
    void export_p_j_downregulated(int filename_suffix = 999999);
    inline void recusively_compute_p_j_downregulated(bool * b, int level, int max_level, double * sum, double p_j_downregulated_given_b, double p_b, double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened);
    inline double iteratively_compute_p_j_downregulated(double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened, int clusters_count);
    void compute_distance_based_predictions();
    inline double distance_based_enchange(unsigned long distance);
    inline void recusively_compute_p_j_downregulated_distance_based(bool * b, int level, int max_level,
                                                                    double * sum, double p_j_downregulated_given_b, double p_b,
                                                                    Binding_flattened * bindings_flattened,
                                                                    double ** distance_based_enhance_matrix,
                                                                    unsigned long latest_bound_level,
                                                                    bool the_current_cluster_is_free_to_bind);
public:
    double lambda;
    
    Matchings_predictor(Patient & patient, std::string simulation_id = "original_data");
    Matchings_predictor(const Matchings_predictor & obj);
    friend void swap(Matchings_predictor & obj1, Matchings_predictor & obj2);
    Matchings_predictor & operator=(Matchings_predictor obj);
    void compute();
};

#endif // MATCHINGS_PREDICTOR_H
