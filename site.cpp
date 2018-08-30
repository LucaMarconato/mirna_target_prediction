#include "site.hpp"

#include <sstream>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

#include "gene.hpp"

std::ostream & operator<<(std::ostream & stream, const Seed_match_type & o)
{
    return stream << Mirna::format_mirna_matching(o.match_type);
}

//------------------------------

Site::Site() {}

void Site::reduce_size_of_scored_interactions_file()
{
    // ob.add_chunk() can only be called passing a number of bytes not bigger than the third argument of the constructor (1000 in this case)
    Output_buffer ob("./data/processed/scored_interactions_processed.tsv", 1024*1024*1024, 1000);
    std::string header = "mirna_id\tgene_id\tutr_start\tutr_end\tseed_match_type\tcontext_score\tweighted_context_score\tconserved\n";
    ob.add_chunk((void *)(header.c_str()), header.size());
    
    io::CSVReader<10, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/scored_interactions.tsv");
    in.read_header(io::ignore_extra_column, "mirna_family", "gene_id", "gene_symbol", "transcript_id", "utr_start", "utr_end", "seed_match_type", "context_score", "weighted_context_score", "conserved");
    std::string columns[10];
    while(in.read_row(columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7], columns[8], columns[9])) {
        Mirna mirna(columns[0]);
        int mirna_id = Mirna::mirna_id_dictionary.left.at(mirna);
        Gene gene(columns[1], columns[2], columns[3]);
        int gene_id = Gene::gene_id_dictionary.left.at(gene);
        Mirna::Mirna_matching mirna_matching;
        if(columns[6] == "7mer-a1") {
            mirna_matching = Mirna::Mirna_matching::canonical_7mer_A1;
        } else if(columns[6] == "7mer-m8") {
            mirna_matching = Mirna::Mirna_matching::canonical_7mer_m8;
        } else if(columns[6] == "8mer") {
            mirna_matching = Mirna::Mirna_matching::canonical_8mer;
        } else {
            std::cerr << "error: columns[6] = " << columns[6] << "\n";
            exit(1);
        }
        std::stringstream ss;
        ss << mirna_id << "\t" << gene_id << "\t" << columns[4] << "\t" << columns[5] << "\t" << mirna_matching << "\t" << columns[7] << "\t" << columns[8] << "\t" << columns[9] << "\n";
        std::string s = ss.str();
        ob.add_chunk((void *)(s.c_str()),s.size());
    }
}

void Site::regenerate_original_file()
{
    
}

std::ostream & operator<<(std::ostream & stream, const Site & o)
{
    return stream << "utr_start = " << o.utr_start << ", utr_end = " << o.utr_end << ", seed_match_type = " << o.seed_match_type << ", context_score = " << o.context_score << ", weighted_context_score = " << o.weighted_context_score << ", conserved = " << o.
conserved << "\n";
}
