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

    // int mirnas_count = patient.tumor_mirnas.profile.size();
    // Perturbation_analyzer perturbation_analyzer(patient);

    // perturbation_analyzer.run(Perturbation_type::Gaussian_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Elements_from_nth_largest(mirnas_count)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           3, double());
    
    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           -1, double());
    
    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           -0.9, double());

    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           -0.5, double());

    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           0.5, double());

    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           1, double());
    
    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           10, double());
    
    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           100, double());
    
    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           0, double());
    // perturbation_analyzer.run(Perturbation_type::No_perturbation,
    //                           Perturbation_type::Point_perturbation,
    //                           Perturbation_target(Perturbation_target_type::Empty_target),
    //                           Perturbation_target(Perturbation_target_type::Nth_larget_element, 0),
    //                           NULL, 1);
    // perturbation_analyzer.run(Perturbation_type::No_perturbation,
    //                           Perturbation_type::Point_perturbation,
    //                           Perturbation_target(Perturbation_target_type::Empty_target),
    //                           Perturbation_target(Perturbation_target_type::Gene_id, 118),
    //                           NULL, 1);

    int experiments_count = 5;
    // int mirnas_count = patient.tumor_mirnas.profile.size();
    // int genes_count = patient.tumor_genes.profile.size();
    for(int i = 0; i < experiments_count; i++) {
        std::stringstream ss;
        ss << i;
        // perturbation_analyzer.run(Perturbation_type::Gaussian_pertubation,
        //                           Perturbation_type::Gaussian_pertubation,
        //                           Perturbation_target(Perturbation_target_type::Elements_from_nth_largest, mirnas_count),
        //                           Perturbation_target(Perturbation_target_type::Elements_from_nth_largest, genes_count),
        //                           1, 1, ss.str());
    }
    std::cout << "freeing pointers\n";
    patient.interaction_graph.free_pointers();
    std::cout << "other cleaning up\n";
    return 0;
}
