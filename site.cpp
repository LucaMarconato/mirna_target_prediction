#include "site.hpp"

#include <sstream>
#include <fstream>
#include <algorithm>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

#include "seed_match_type.hpp"
#include "timer.hpp"

Site::Site() {}

Site::Site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end) :
    mirna_id(mirna_id),
    gene_id(gene_id),
    utr_start(utr_start),
    utr_end(utr_end)
{

}

void Site::reduce_size_of_scored_interactions_file()
{
    Timer::start();
    std::cout << "reducing the size of the scored interactions file\n";
    // ob.add_chunk() can only be called passing a number of bytes not bigger than the third argument of the constructor (1000 in this case)
    Output_buffer ob("./data/processed/scored_interactions_processed.tsv", 1024*1024*1024, 1000);
    ob.verbose = true;
    std::string header = "mirna_id\tgene_id\tutr_start\tutr_end\tseed_match_type\tcontext_score\tweighted_context_score\tconserved\n";
    ob.add_chunk(header);

    io::CSVReader<10, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/scored_interactions.tsv");
    in.read_header(io::ignore_extra_column, "mirna_family", "gene_id", "gene_symbol", "transcript_id", "utr_start", "utr_end", "seed_match_type", "context_score", "weighted_context_score", "conserved");
    std::string columns[10];
    while(in.read_row(columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7], columns[8], columns[9])) {
        std::string mirna_family = columns[0];
        std::transform(mirna_family.begin(), mirna_family.end(), mirna_family.begin(), ::tolower);
        Mirna mirna(mirna_family);
        Mirna_id mirna_id = Mirna::mirna_id_dictionary.left.at(mirna);
        Gene gene(columns[1], columns[2], columns[3]);
        Gene_id gene_id = Gene::gene_id_dictionary.left.at(gene);
        Seed_match_type::Mirna_site_matching mirna_matching;
        if(columns[6] == "7mer-a1") {
            mirna_matching = Seed_match_type::Mirna_site_matching::canonical_7mer_a1;
        } else if(columns[6] == "7mer-m8") {
            mirna_matching = Seed_match_type::Mirna_site_matching::canonical_7mer_m8;
        } else if(columns[6] == "8mer") {
            mirna_matching = Seed_match_type::Mirna_site_matching::canonical_8mer;
        } else {
            std::cerr << "error: columns[6] = " << columns[6] << "\n";
            exit(1);
        }
        bool conserved;
        if(columns[9] == "True") {
            conserved = true;
        } else if(columns[9] == "False") {
            conserved = false;
        } else {
            std::cerr << "error: columns[9] = " << columns[9] << "\n";
            exit(1);
        }
        std::stringstream ss;
        ss << mirna_id << "\t" << gene_id << "\t" << columns[4] << "\t" << columns[5] << "\t" << mirna_matching << "\t" << columns[7] << "\t" << columns[8] << "\t" << conserved << "\n";
        std::string s = ss.str();
        ob.add_chunk((void *)(s.c_str()),s.size());
    }
    ob.write_to_file();
    std::cout << "finished, \n";
    Timer::stop();
}

std::ostream & operator<<(std::ostream & stream, const Site & o)
{
    return stream << "utr_start = " << o.utr_start << ", utr_end = " << o.utr_end << "\n";
}
