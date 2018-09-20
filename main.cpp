#include <iostream>

#include "global_parameters.hpp"
#include "mirna.hpp"
#include "gene.hpp"
#include "patient.hpp"
#include "matchings_predictor.hpp"

int main(int argc, char * argv[])
{
    Global_parameters::load_from_json();
    Mirna::initialize_mirna_dictionary();
    Gene::initialize_gene_dictionary();
    // Mirna::print_mirna_dictionary(10);
    // Gene::print_gene_dictionary(10);

    Patient patient("TCGA-CJ-4642");
    // Patient patient("artificial0");

    Matchings_predictor matching_predictor(patient);
    matching_predictor.compute();

    std::cout << "cleaning up\n";
    return 0;
}

