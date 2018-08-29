#include <iostream>
#include <fstream>
#include <unordered_map>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>

#include "mirna.hpp"
#include "gene.hpp"

int main(int argc, char * argv[])
{
    // 
    // std::unordered_map<Mirna, int> mirna_id_dictionary_left;
    Mirna::initialize_mirna_dictionary();
    Mirna::print_mirna_dictionary(10);

    Gene::initialize_gene_dictionary();
    Gene::print_gene_dictionary(100);

    return 0;
}

