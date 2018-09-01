#ifndef PATIENT_H
#define PATIENT_H

#include <string>

#include "mirna_expression_profile.hpp"
#include "gene_expression_profile.hpp"
#include "site_expression_profile.hpp"

class Patient {
public:
    std::string case_id;
    Mirna_expression_profile normal_mirna, tumor_mirna;
    Gene_expression_profile normal_mrna, tumor_mrna;
    // Site_expression_profile ...;

    Patient(std::string case_id);
};

#endif // PATIENT_H
