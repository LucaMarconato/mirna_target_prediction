#ifndef PATIENT_H
#define PATIENT_H

#include <string>

#include "mirna_expression_profile.hpp"
#include "gene_expression_profile.hpp"
#include "site_expression_profile.hpp"
#include "interaction_graph.hpp"

class Patient {
public:
    std::string case_id;
    Mirna_expression_profile normal_mirnas, tumor_mirnas;
    Gene_expression_profile normal_genes, tumor_genes;
    Site_expression_profile normal_sites, tumor_sites;
    // Site_expression_profile ...;
    Interaction_graph interaction_graph;

    Patient(std::string case_id);
    void export_expression_profiles(std::string patient_folder);
};

#endif // PATIENT_H
