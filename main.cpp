#include <iostream>

#include <boost/filesystem.hpp>

#include "global_parameters.hpp"
#include "mirna.hpp"
#include "gene.hpp"
#include "patient.hpp"
#include "matchings_predictor.hpp"
#include "perturbation_analyzer.hpp"

int main(int argc, char * argv[])
{
    std::cout << "argv[0] = " << argv[0] << "\n";
    boost::filesystem::path full_path(boost::filesystem::current_path());
    std::cout << "Current path is: " << full_path << std::endl;
    
    Global_parameters::load_from_json();
    Mirna::initialize_mirna_dictionary();
    Gene::initialize_gene_dictionary();
    // Mirna::print_mirna_dictionary(10);
    // Gene::print_gene_dictionary(10);
    // Site::reduce_size_of_scored_interactions_file();

    Patient patient("TCGA-CJ-4642", true);
    // Patient patient("artificial0", true);    
    // Patient patient("artificial1", true);

    // Matchings_predictor matching_predictor(patient);
    // matching_predictor.compute();

    Perturbation_analyzer perturbation_analyzer(patient);

    perturbation_analyzer.run(Perturbation_type::No_perturbation,
                              Perturbation_type::No_perturbation,
                              Perturbation_target(Perturbation_target::Empty_target()),
                              Perturbation_target(Perturbation_target::Empty_target()),
                              Perturbation_extent(Perturbation_extent::No_perturbation()),
                              Perturbation_extent(Perturbation_extent::No_perturbation()));

    Global_parameters::consider_distance_for_predictions = true;
    
    perturbation_analyzer.run(Perturbation_type::No_perturbation,
                              Perturbation_type::No_perturbation,
                              Perturbation_target(Perturbation_target::Empty_target()),
                              Perturbation_target(Perturbation_target::Empty_target()),
                              Perturbation_extent(Perturbation_extent::No_perturbation()),
                              Perturbation_extent(Perturbation_extent::No_perturbation()),
                              "distance");

    // int mirnas_count = patient.tumor_mirnas.profile.size();
    // int first_n_mirnas_to_perturb = mirnas_count;
    // for(int i = 0; i < first_n_mirnas_to_perturb; i++) {
    //     perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                               Perturbation_type::No_perturbation,
    //                               Perturbation_target(Perturbation_target::Nth_largest_element(i)),
    //                               Perturbation_target(Perturbation_target::Empty_target()),
    //                               // Perturbation_extent(Perturbation_extent::Relative_perturbation(1)),
    //                               Perturbation_extent(Perturbation_extent::Absolute_perturbation(500000)),
    //                               Perturbation_extent(Perturbation_extent::No_perturbation()));
    // }

    // perturbation_analyzer.run(Perturbation_type::Gaussian_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Elements_from_nth_largest(mirnas_count)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           Perturbation_extent(Perturbation_extent::Relative_perturbation(3)),
    //                           Perturbation_extent(Perturbation_extent::No_perturbation()));
    
    // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    //                           Perturbation_type::No_perturbation,
    //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    //                           Perturbation_target(Perturbation_target::Empty_target()),
    //                           Perturbation_extent(Perturbation_extent::Relative_perturbation(-1)),
    //                           Perturbation_extent(Perturbation_extent::No_perturbation()));

    // int experiments_count = 5;
    // for(int i = 0; i < experiments_count; i++) {
    //     std::stringstream ss;
    //     ss << i;
    //     perturbation_analyzer.run(Perturbation_type::Gaussian_perturbation,
    //                               Perturbation_type::No_perturbation,
    //                               Perturbation_target(Perturbation_target::Elements_from_nth_largest(mirnas_count)),
    //                               Perturbation_target(Perturbation_target::Empty_target()),
    //                               Perturbation_extent(Perturbation_extent::Relative_perturbation(1)),
    //                               Perturbation_extent(Perturbation_extent::No_perturbation()));
    // }
    std::cout << "freeing pointers\n";
    patient.interaction_graph.free_pointers();
    std::cout << "other cleaning up\n";
    return 0;
}
