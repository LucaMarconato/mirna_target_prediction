#include "perturbation_analyzer.hpp"

Perturbation_analyzer::Perturbation_analyzer(Patient & patient) : patient(patient) {}

void Perturbation_analyzer::configure(Perturbation_type mirna_perturbation_type,
                                      Perturbation_type gene_perturbation_type,
                                      int number_of_most_expressed_mirnas_to_perturb,
                                      int number_of_most_expressed_genes_to_perturb,
                                      double mirna_perturbation_amplifier,
                                      double gene_perturbation_amplifier,
                                      std::string output_name_suffix)
{
    if(mirna_perturbation_type == Perturbation_type::Point_perturbation && gene_perturbation_type == Perturbation_type::Point_perturbation && number_of_most_expressed_mirnas_to_perturb * number_of_most_expressed_genes_to_perturb != 0) {
        std::cerr << "error: if you choose to perform a point perturbation to both mirnas and genes, then you can perturb either mirnas either genes at the same time.\n";
        exit(1);
    }
}
