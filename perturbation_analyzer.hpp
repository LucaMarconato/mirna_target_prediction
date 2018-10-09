#ifndef PERTURBATION_ANALYZER_H
#define PERTURBATION_ANALYZER_H

#include "patient.hpp"

enum class Perturbation_type {
	Point_perturbation,
	Gaussian_pertubation
	// Bernoully_perturbation
};

class Perturbation_analyzer {
    Patient & patient;
    
public: 
    Perturbation_type mirna_perturbation_type;
    Perturbation_type gene_perturbation_type;
    int number_of_most_expressed_mirnas_to_perturb;
    int number_of_most_expressed_genes_to_perturb;
    double mirna_perturbation_amplifier;
    double gene_perturbation_amplifier;

    Perturbation_analyzer(Patient & patient);
    void configure(Perturbation_type mirna_perturbation_type,
                   Perturbation_type gene_perturbation_type,
                   int number_of_most_expressed_mirnas_to_perturb,
                   int number_of_most_expressed_genes_to_perturb,
                   double mirna_perturbation_amplifier,
                   double gene_perturbation_amplifier,
                   std::string output_name_suffix = "");
};

#endif // PERTURBATION_ANALYZER_H
