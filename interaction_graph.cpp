#include "interaction_graph.hpp"

#include <cstdlib>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include <strasser/csv.h>
#include <sohail/serialize_tuple.hpp>

#include "seed_match_type.hpp"
#include "timer.hpp"

#define Ig Interaction_graph

std::unordered_map<boost::tuple<Mirna_id, Gene_id, unsigned int, unsigned int>, Site *> Ig::sites_by_location;

Mirna_site_arc::Mirna_site_arc(Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved) : seed_match_type(seed_match_type), context_score(context_score), weighted_context_score(weighted_context_score), conserved(conserved) {}

Mirna_site_arc::Mirna_site_arc(const Mirna_site_arc & obj)
{
    this->seed_match_type = obj.seed_match_type;
    this->context_score = obj.context_score;
    this->weighted_context_score = obj.weighted_context_score;
    this->conserved = obj.conserved;
}

void swap(Mirna_site_arc & obj1, Mirna_site_arc & obj2)
{
    std::swap(obj1.seed_match_type,obj2.seed_match_type);
    std::swap(obj1.context_score,obj2.context_score);
    std::swap(obj1.weighted_context_score,obj2.weighted_context_score);
    std::swap(obj1.conserved,obj2.conserved);    
}

Mirna_site_arc & Mirna_site_arc::operator=(Mirna_site_arc obj)
{
    swap(*this, obj);
    return *this;
}

// --------------------------------------------------

std::unordered_map<Gene_id, std::list<Site *>> Ig::gene_to_sites_arcs;
std::unordered_map<std::pair<Mirna_id, Site *>, Mirna_site_arc> Ig::mirna_site_arcs;
std::unordered_map<Mirna_id, std::list<Site *>> Ig::mirna_to_sites_arcs;
std::unordered_map<Site *, std::list<Mirna_id>> Ig::site_to_mirnas_arcs;
std::unordered_map<std::pair<Mirna_id, Gene_id>, std::list<Site *>> Ig::mirna_gene_arcs;
std::unordered_map<Mirna_id, std::list<Gene_id>> Ig::mirna_to_genes_arcs;
std::unordered_map<Gene_id, std::list<Mirna_id>> Ig::gene_to_mirnas_arcs;

void Ig::build_interaction_graph()
{
    if(!boost::filesystem::exists("all_sites.bin") || !boost::filesystem::exists("interaction_graph.bin")) {
        io::CSVReader<8, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/scored_interactions_processed.tsv");
        in.read_header(io::ignore_extra_column, "mirna_id",  "gene_id",  "utr_start",  "utr_end",  "seed_match_type",  "context_score",  "weighted_context_score", "conserved");
        std::string columns[8];
        while(in.read_row(columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7])) {
            Mirna_id mirna_id = std::strtoul(columns[0].c_str(), nullptr, 10);
            Gene_id gene_id = std::strtoul(columns[1].c_str(), nullptr, 10);
            unsigned int utr_start = std::strtoul(columns[2].c_str(), nullptr, 10);
            unsigned int utr_end = std::strtoul(columns[3].c_str(), nullptr, 10);
            Seed_match_type::Mirna_site_matching mirna_site_matching = static_cast<Seed_match_type::Mirna_site_matching>(atoi(columns[4].c_str()));
            Seed_match_type seed_match_type(mirna_site_matching);
            double context_score = std::strtod(columns[5].c_str(), nullptr);
            double weighted_context_score = std::strtod(columns[6].c_str(), nullptr);
            bool conserved = atoi(columns[7].c_str());
            Site * site = Ig::get_site(mirna_id, gene_id, utr_start, utr_end);

            // create gene-site arcs
            if(Ig::gene_to_sites_arcs.find(gene_id) == Ig::gene_to_sites_arcs.end()) {
                // this initialize the list as empty
                Ig::gene_to_sites_arcs[gene_id];
            }
            Ig::gene_to_sites_arcs[gene_id].push_back(site);       

            // create mirna-site arcs
            Mirna_site_arc mirna_site_arc(seed_match_type, context_score, weighted_context_score, conserved);
            Ig::mirna_site_arcs[std::make_pair(mirna_id, site)] = mirna_site_arc;
        
            if(Ig::mirna_to_sites_arcs.find(mirna_id) == Ig::mirna_to_sites_arcs.end()) {
                Ig::mirna_to_sites_arcs[mirna_id];
            }
            Ig::mirna_to_sites_arcs[mirna_id].push_back(site);
        
            if(Ig::site_to_mirnas_arcs.find(site) == Ig::site_to_mirnas_arcs.end()) {
                Ig::site_to_mirnas_arcs[site];
            }
            Ig::site_to_mirnas_arcs[site].push_back(mirna_id);
        
            // create mirna-gene arcs
            auto p = std::make_pair(mirna_id, gene_id);
            if(Ig::mirna_gene_arcs.find(p) == Ig::mirna_gene_arcs.end()) {
                Ig::mirna_gene_arcs[p];
            }
            Ig::mirna_gene_arcs[p].push_back(site);
        
            if(Ig::mirna_to_genes_arcs.find(mirna_id) == Ig::mirna_to_genes_arcs.end()) {
                Ig::mirna_to_genes_arcs[mirna_id];
            }
            Ig::mirna_to_genes_arcs[mirna_id].push_back(gene_id);
        
            if(Ig::gene_to_mirnas_arcs.find(gene_id) == Ig::gene_to_mirnas_arcs.end()) {
                Ig::gene_to_mirnas_arcs[gene_id];
            }
            Ig::gene_to_mirnas_arcs[gene_id].push_back(mirna_id);
        }
        // Site::save_all_sites();
        std::cout << "writing interaction_graph.bin\n";
        Timer::start();
        std::ofstream out("interaction_graph.bin", std::ios::binary);
        boost::archive::binary_oarchive oa(out);
        oa << Ig::sites_by_location;
        oa << Ig::gene_to_sites_arcs;
        oa << Ig::mirna_site_arcs;
        oa << Ig::mirna_to_sites_arcs;
        oa << Ig::site_to_mirnas_arcs;
        oa << Ig::mirna_gene_arcs;
        oa << Ig::mirna_to_genes_arcs;
        oa << Ig::gene_to_mirnas_arcs;
        out.close();
        std::cout << "written, ";
        Timer::stop();
    } else {
        // Site::load_all_sites();
        std::cout << "loading interaction_graph.bin\n";
        Timer::start();
        std::ifstream in("interaction_graph.bin", std::ios::binary);
        boost::archive::binary_iarchive ia(in);
        ia >> Ig::sites_by_location;
        ia >> Ig::gene_to_sites_arcs;        
        ia >> Ig::mirna_site_arcs;
        ia >> Ig::mirna_to_sites_arcs;
        ia >> Ig::site_to_mirnas_arcs;
        ia >> Ig::mirna_gene_arcs;
        ia >> Ig::mirna_to_genes_arcs;
        ia >> Ig::gene_to_mirnas_arcs;
        in.close();
        std::cout << "loaded, ";
        Timer::stop();
    }
}

Site * Ig::get_site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end)
{
    auto t = boost::make_tuple(mirna_id, gene_id, utr_start, utr_end);
    auto e = Ig::sites_by_location.find(t);
    if(e != Ig::sites_by_location.end()) {
        return e->second;
    } else {
        Site * site = new Site(mirna_id, gene_id, utr_start, utr_end);
        Ig::sites_by_location[boost::make_tuple(mirna_id, gene_id, utr_start, utr_end)] = site;
        return site;   
    }
}

void Ig::delete_all_sites()
{
    Timer::start();
    std::cout << "deleting all_sites\n";
    for(auto & e : Ig::sites_by_location) {
        delete e.second;
    }
    std::cout << "deleted, ";
    Timer::stop();
}
