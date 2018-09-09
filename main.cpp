#include <iostream>

#include "global_parameters.hpp"
#include "mirna.hpp"
#include "gene.hpp"
#include "site.hpp"
#include "patient.hpp"

int main(int argc, char * argv[])
{
    Global_parameters::load_from_json();
    Mirna::initialize_mirna_dictionary();
    // Mirna::print_mirna_dictionary(10);
    Gene::initialize_gene_dictionary();
    // Gene::print_gene_dictionary(10);

    Patient patient("TCGA-CJ-4642");

    std::cout << "cleaning up\n";
    return 0;
}

