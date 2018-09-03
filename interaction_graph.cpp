#include "interaction_graph.hpp"

#include <cstdlib>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

#include <strasser/csv.h>

#include "seed_match_type.hpp"
#include "timer.hpp"

#define Ig Interaction_graph

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

Ig::Interaction_graph() {}

void Ig::build_interaction_graph(std::set<Mirna_id> & mirnas, std::set<Gene_id> & genes)
{
    std::cout << "building the interaction graph\n";
    Timer::start();    
    io::CSVReader<8, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/scored_interactions_processed.tsv");
    in.read_header(io::ignore_extra_column, "mirna_id",  "gene_id",  "utr_start",  "utr_end",  "seed_match_type",  "context_score",  "weighted_context_score", "conserved");
    std::string columns[8];
    while(in.read_row(columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7])) {
        Mirna_id mirna_id = std::strtoul(columns[0].c_str(), nullptr, 10);
        Gene_id gene_id = std::strtoul(columns[1].c_str(), nullptr, 10);
        if(mirnas.find(mirna_id) == mirnas.end() || genes.find(gene_id) == genes.end()) {
            this->rows_skipped++;
            continue;
        }
        this->rows_processed++;            
        unsigned int utr_start = std::strtoul(columns[2].c_str(), nullptr, 10);
        unsigned int utr_end = std::strtoul(columns[3].c_str(), nullptr, 10);
        Seed_match_type::Mirna_site_matching mirna_site_matching = static_cast<Seed_match_type::Mirna_site_matching>(atoi(columns[4].c_str()));
        Seed_match_type seed_match_type(mirna_site_matching);
        double context_score = std::strtod(columns[5].c_str(), nullptr);
        double weighted_context_score = std::strtod(columns[6].c_str(), nullptr);
        bool conserved = atoi(columns[7].c_str());
        Site * site = this->get_site(mirna_id, gene_id, utr_start, utr_end);

        // create gene_site arcs
        if(this->gene_to_sites_arcs.find(gene_id) == this->gene_to_sites_arcs.end()) {
            // this initialize the list as empty
            this->gene_to_sites_arcs[gene_id];
        }
        this->gene_to_sites_arcs[gene_id].push_back(site);       

        // create mirna_site arcs
        Mirna_site_arc mirna_site_arc(seed_match_type, context_score, weighted_context_score, conserved);
        this->mirna_site_arcs[std::make_pair(mirna_id, site)] = mirna_site_arc;
        
        if(this->mirna_to_sites_arcs.find(mirna_id) == this->mirna_to_sites_arcs.end()) {
            this->mirna_to_sites_arcs[mirna_id];
        }
        this->mirna_to_sites_arcs[mirna_id].push_back(site);
        
        if(this->site_to_mirnas_arcs.find(site) == this->site_to_mirnas_arcs.end()) {
            this->site_to_mirnas_arcs[site];
        }
        this->site_to_mirnas_arcs[site].push_back(mirna_id);
        
        // create mirna_gene arcs
        auto p = std::make_pair(mirna_id, gene_id);
        if(this->mirna_gene_arcs.find(p) == this->mirna_gene_arcs.end()) {
            this->mirna_gene_arcs[p];
        }
        this->mirna_gene_arcs[p].push_back(site);
        
        if(this->mirna_to_genes_arcs.find(mirna_id) == this->mirna_to_genes_arcs.end()) {
            this->mirna_to_genes_arcs[mirna_id];
        }
        this->mirna_to_genes_arcs[mirna_id].push_back(gene_id);
        
        if(this->gene_to_mirnas_arcs.find(gene_id) == this->gene_to_mirnas_arcs.end()) {
            this->gene_to_mirnas_arcs[gene_id];
        }
        this->gene_to_mirnas_arcs[gene_id].push_back(mirna_id);
    }
    std::cout << "built, ";
    Timer::stop();
}

Site * Ig::get_site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end)
{
    auto t = boost::make_tuple(mirna_id, gene_id, utr_start, utr_end);
    auto e = this->sites_by_location.find(t);
    if(e != this->sites_by_location.end()) {
        return e->second;
    } else {
        Site * site = new Site(mirna_id, gene_id, utr_start, utr_end);
        this->sites_by_location[boost::make_tuple(mirna_id, gene_id, utr_start, utr_end)] = site;
        return site;   
    }
}

void Ig::print_statistics()
{
    std::cout << "rows_processed/rows_skipped = " << rows_processed << "/" << rows_skipped << " = " << ((double)rows_processed)/rows_skipped << "\n";
    std::cout << "sites_by_location.size() = " << sites_by_location.size() << "\n";
    std::cout << "arcs: gene_site\n";
    std::cout << "gene_to_sites_arcs.size() = " << gene_to_sites_arcs.size() << "\n";
    std::cout << "arcs: mirna_site\n";
    std::cout << "mirna_site_arcs.size() = " << mirna_site_arcs.size() << "\n";
    std::cout << "mirna_to_sites_arcs.size() = " << mirna_to_sites_arcs.size() << "\n";
    std::cout << "site_to_mirnas_arcs.size() = " << site_to_mirnas_arcs.size() << "\n";
    std::cout << "arcs: mirna_gene\n";
    std::cout << "mirna_gene_arcs.size() = " << mirna_gene_arcs.size() << "\n";
    std::cout << "mirna_to_genes_arcs.size() = " << mirna_to_genes_arcs.size() << "\n";
    std::cout << "gene_to_mirnas_arcs.size() = " << gene_to_mirnas_arcs.size() << "\n";
}

Ig::~Interaction_graph()
{
    Timer::start();
    std::cout << "deleting all_sites\n";
    for(auto & e : this->sites_by_location) {
        delete e.second;
    }
    std::cout << "deleted, ";
    Timer::stop();
}
