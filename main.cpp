#include <iostream>

#include "global_parameters.hpp"
#include "mirna.hpp"
#include "gene.hpp"
#include "patient.hpp"
#include "matchings_predictor.hpp"
#include "perturbation_analyzer.hpp"

int main(int argc, char * argv[])
{
    Global_parameters::load_from_json();
    Mirna::initialize_mirna_dictionary();
    Gene::initialize_gene_dictionary();
    // Mirna::print_mirna_dictionary(10);
    // Gene::print_gene_dictionary(10);

    Patient patient("TCGA-CJ-4642", true);
    // Patient patient("artificial0", true);
    // Patient patient("artificial1", true);

    Matchings_predictor matching_predictor(patient);
    matching_predictor.compute();
    
    Matchings_predictor matching_predictor1(patient, "test");
    matching_predictor1.compute();

    Perturbation_analyzer perturbation_analyzer(patient);
    perturbation_analyzer.configure(Perturbation_type::Point_perturbation,
                                    Perturbation_type::Point_perturbation,
                                    3, 0, 1, 1);
    perturbation_analyzer.configure(Perturbation_type::Point_perturbation,
                                    Perturbation_type::Point_perturbation,
                                    0, 3, 1, 1);

    int experiments_count = 5;
    int mirnas_count = patient.tumor_mirnas.profile.size();
    int genes_count = patient.tumor_genes.profile.size();
    for(int i = 0; i < experiments_count; i++) {
        std::stringstream ss;
        ss << i;
        perturbation_analyzer.configure(Perturbation_type::Gaussian_pertubation,
                                        Perturbation_type::Gaussian_pertubation,
                                        mirnas_count, experiments_count, 1, 1, ss.str());
    }

    std::cout << "cleaning up\n";
    return 0;
}

