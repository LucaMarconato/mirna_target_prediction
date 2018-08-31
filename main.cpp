#include <iostream>
#include <fstream>
#include <unordered_map>

#include <marconato/output_buffer/output_buffer.hpp>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>

#include "mirna.hpp"
#include "gene.hpp"
#include "site.hpp"
#include "interaction_graph.hpp"

int main(int argc, char * argv[])
{
    // 
    // std::unordered_map<Mirna, int> mirna_id_dictionary_left;
    Mirna::initialize_mirna_dictionary();
    Mirna::print_mirna_dictionary(10);

    Gene::initialize_gene_dictionary();
    Gene::print_gene_dictionary(10);

    // Site::reduce_size_of_scored_interactions_file();
    // Site::regenerate_original_file();

    Interaction_graph::build_interaction_graph();

    Site::delete_all_sites();
    return 0;
}

