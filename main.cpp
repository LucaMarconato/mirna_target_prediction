#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>

#include "gene.hpp"
#include "global_parameters.hpp"
#include "matchings_predictor.hpp"
#include "mirna.hpp"
#include "patient.hpp"
#include "perturbation_analyzer.hpp"

int main(int argc, char* argv[])
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

    std::vector<std::string> patient_ids = {"artificial_ENCFF360IHM-hela"};
    // std::vector<std::string> patient_ids = {"artificial_ENCFF360IHM-hela", "artificial_ENCFF495ZXC-hela", "artificial_ENCFF612ZIR-hela", "artificial_ENCFF729EQX-hela", "artificial_ENCFF806EYY-hela", "artificial_ENCFF902KUU-hela"};
    std::vector<std::string> transfected_mirbase_ids = {"hsa-mir-122-5p", "hsa-mir-128-3p", "hsa-mir-132-3p", "hsa-mir-142-5p"};

    std::vector<Patient> patients;

    // patients.emplace_back(Patient("artificial_TCGA-CJ-4642", true));

    for (auto& patient_id : patient_ids) {
        for (auto& transfected_mirbase_id : transfected_mirbase_ids) {
            std::string transfected_patient_id = patient_id + "_" + transfected_mirbase_id + "_transfected";
            patients.emplace_back(Patient(transfected_patient_id, true));
            patients.emplace_back(Patient(patient_id, true));
        }
    }

    // for (auto& patient_id : patient_ids) {
    //     patients.emplace_back(Patient(patient_id, true));
    // }

    for (auto& patient : patients) {
        Matchings_predictor matching_predictor(patient);
        matching_predictor.compute();

        // Perturbation_analyzer perturbation_analyzer(patient);
        // int mirnas_count = patient.mirna_expression_profile.profile.size();
        // int first_n_mirnas_to_perturb = mirnas_count;
        // // int first_n_mirnas_to_perturb = 3;
        // for (int i = 0; i < first_n_mirnas_to_perturb; i++) {
        //     perturbation_analyzer.run(Perturbation_type::Point_perturbation,
        //                               Perturbation_type::No_perturbation,
        //                               Perturbation_target(Perturbation_target::Nth_largest_element(i)),
        //                               Perturbation_target(Perturbation_target::Empty_target()),
        //                               Perturbation_extent(Perturbation_extent::Absolute_perturbation(500000)),
        //                               Perturbation_extent(Perturbation_extent::No_perturbation()));
        // }
        patient.interaction_graph.free_pointers();
    }

    // Patient& patient = patients.at(0);
    // Perturbation_analyzer perturbation_analyzer(patient);

    // perturbation_analyzer.run(Perturbation_type::No_perturbation, Perturbation_type::No_perturbation, Perturbation_target(Perturbation_target::Empty_target()), Perturbation_target(Perturbation_target::Empty_target()),
    //                           Perturbation_extent(Perturbation_extent::No_perturbation()), Perturbation_extent(Perturbation_extent::No_perturbation()));

    // // Global_parameters::consider_distance_for_predictions = true;

    // // perturbation_analyzer.run(Perturbation_type::No_perturbation,
    // //                           Perturbation_type::No_perturbation,
    // //                           Perturbation_target(Perturbation_target::Empty_target()),
    // //                           Perturbation_target(Perturbation_target::Empty_target()),
    // //                           Perturbation_extent(Perturbation_extent::No_perturbation()),
    // //                           Perturbation_extent(Perturbation_extent::No_perturbation()),
    // //                           "distance");

    // int mirnas_count = patient.mirna_expression_profile.profile.size();
    // int first_n_mirnas_to_perturb = mirnas_count;
    // // int first_n_mirnas_to_perturb = 3;
    // for (int i = 0; i < first_n_mirnas_to_perturb; i++) {
    //     perturbation_analyzer.run(Perturbation_type::Point_perturbation, Perturbation_type::No_perturbation, Perturbation_target(Perturbation_target::Nth_largest_element(i)), Perturbation_target(Perturbation_target::Empty_target()),
    //                               // Perturbation_extent(Perturbation_extent::Relative_perturbation(-1)),
    //                               Perturbation_extent(Perturbation_extent::Absolute_perturbation(500000)), Perturbation_extent(Perturbation_extent::No_perturbation()));
    // }

    // // perturbation_analyzer.run(Perturbation_type::Gaussian_perturbation,
    // //                           Perturbation_type::No_perturbation,
    // //                           Perturbation_target(Perturbation_target::Elements_from_nth_largest(mirnas_count)),
    // //                           Perturbation_target(Perturbation_target::Empty_target()),
    // //                           Perturbation_extent(Perturbation_extent::Relative_perturbation(3)),
    // //                           Perturbation_extent(Perturbation_extent::No_perturbation()));

    // // perturbation_analyzer.run(Perturbation_type::Point_perturbation,
    // //                           Perturbation_type::No_perturbation,
    // //                           Perturbation_target(Perturbation_target::Nth_largest_element(0)),
    // //                           Perturbation_target(Perturbation_target::Empty_target()),
    // //                           Perturbation_extent(Perturbation_extent::Relative_perturbation(-1)),
    // //                           Perturbation_extent(Perturbation_extent::No_perturbation()));

    // // int max_perturbation_extent = 10;
    // // int replicates = 3;
    // // for(int i = 0; i < max_perturbation_extent; i++) {
    // //     for(int j = 0; j < replicates; j++) {
    // //         std::stringstream ss;
    // //         ss << j;
    // //         perturbation_analyzer.run(Perturbation_type::Gaussian_perturbation,
    // //                                   Perturbation_type::No_perturbation,
    // //                                   Perturbation_target(Perturbation_target::Elements_from_nth_largest(mirnas_count
    // //                                   - 1)),
    // //                                   Perturbation_target(Perturbation_target::Empty_target()),
    // //                                   Perturbation_extent(Perturbation_extent::Relative_perturbation(i+1)),
    // //                                   Perturbation_extent(Perturbation_extent::No_perturbation()),
    // //                                   ss.str());
    // //     }
    // // }
    std::cout << "freeing pointers\n";
    // patient.interaction_graph.free_pointers();
    std::cout << "other cleaning up\n";
    return 0;
}
