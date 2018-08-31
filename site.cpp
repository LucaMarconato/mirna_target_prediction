#include "site.hpp"

#include <sstream>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/list.hpp>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

std::list<Site *> Site::all_sites;

Site::Site() {}

Site::Site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end, Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved) :
    mirna_id(mirna_id), 
    gene_id(gene_id), 
    utr_start(utr_start), 
    utr_end(utr_end), 
    seed_match_type(seed_match_type), 
    context_score(context_score), 
    weighted_context_score(weighted_context_score), 
    conserved(conserved)
{
    
}

template<class Archive>
void Site::serialize(Archive & ar, const unsigned int version)
{
    ar & this->mirna_id;
    ar & this->gene_id;
    ar & this->utr_start;
    ar & this->utr_end;
    ar & this->seed_match_type;
    ar & this->context_score;
    ar & this->weighted_context_score;
    ar & this->conserved;
}

void Site::reduce_size_of_scored_interactions_file()
{
    // ob.add_chunk() can only be called passing a number of bytes not bigger than the third argument of the constructor (1000 in this case)
    Output_buffer ob("./data/processed/scored_interactions_processed.tsv", 1024*1024*1024, 1000);
    std::string header = "mirna_id\tgene_id\tutr_start\tutr_end\tseed_match_type\tcontext_score\tweighted_context_score\tconserved\n";
    ob.add_chunk(header);
    
    io::CSVReader<10, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/scored_interactions.tsv");
    in.read_header(io::ignore_extra_column, "mirna_family", "gene_id", "gene_symbol", "transcript_id", "utr_start", "utr_end", "seed_match_type", "context_score", "weighted_context_score", "conserved");
    std::string columns[10];
    while(in.read_row(columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7], columns[8], columns[9])) {
        Mirna mirna(columns[0]);
        Mirna_id mirna_id = Mirna::mirna_id_dictionary.left.at(mirna);
        Gene gene(columns[1], columns[2], columns[3]);
        Gene_id gene_id = Gene::gene_id_dictionary.left.at(gene);
        Seed_match_type::Mirna_site_matching mirna_matching;
        if(columns[6] == "7mer-a1") {
            mirna_matching = Seed_match_type::Mirna_site_matching::canonical_7mer_A1;
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
}

// void Site::regenerate_original_file()
// {
//     Output_buffer ob("./data/processed/scored_interactions_test.tsv", 1024*1024*1024, 1000);
//     std::string header = "mirna_family\tgene_id\tgene_symbol\ttranscript_id\tutr_start\tutr_end\tseed_match_type\tcontext_score\tweighted_context_score\tconserved";
//     ob.add_chunk(header);
//     int mirna_count = Mirna::mirna_id_dictionary.size();
//     for(int i = 0; i < mirna_count; i++) {
//         std::string 
//     }
// }

Site * Site::new_site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end, Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved)
{
    Site * site = new Site(mirna_id, gene_id, utr_start, utr_end, seed_match_type, context_score, weighted_context_score, conserved);
    Site::all_sites.push_back(site);
    return site;
}

void Site::delete_all_sites()
{
    for(auto & e : Site::all_sites) {
        delete e;
    }
}

void Site::save_all_sites()
{
    std::ofstream out("all_sites.bin", std::ios::binary);
    boost::archive::binary_oarchive oa(out);
    oa << Site::all_sites;
    out.close();
}

void Site::load_all_sites()
{
    std::ifstream in("all_sites.bin", std::ios::binary);
    boost::archive::binary_iarchive ia(in);
    ia >> Site::all_sites;
    in.close();
}

std::ostream & operator<<(std::ostream & stream, const Site & o)
{
    return stream << "utr_start = " << o.utr_start << ", utr_end = " << o.utr_end << ", seed_match_type = " << o.seed_match_type << ", context_score = " << o.context_score << ", weighted_context_score = " << o.weighted_context_score << ", conserved = " << o.
conserved << "\n";
}
