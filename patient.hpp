#ifndef PATIENT_H
#define PATIENT_H

#include <string>

#include "mirna_expression_profile.hpp"
#include "gene_expression_profile.hpp"
#include "cluster_expression_profile.hpp"
#include "interaction_graph.hpp"

class Matchings_predictor;

class Patient {
    Patient();

public:
    std::string case_id;
    Mirna_expression_profile mirna_expression_profile;
    Gene_expression_profile gene_expression_profile;
    Cluster_expression_profile cluster_expression_profile;
    Interaction_graph interaction_graph;

    Patient(std::string case_id, bool export_data);
    Patient(const Patient & obj);
    friend void swap(Patient & obj1, Patient & obj2);
    Patient & operator=(Patient obj);
    void generate_cluster_expression_profiles();
    void export_expression_profiles(std::string patient_folder);

    // to access the default private constructor
    friend class Perturbation_analyzer;
};

#endif // PATIENT_H
