#ifndef INTERACTION_GRAPH_H
#define INTERACTION_GRAPH_H

#include <list>
#include <utility>
#include <unordered_map>
#include <set>

#include <boost/serialization/list.hpp>
#include <boost/serialization/unordered_map.hpp>

#include <sohail/serialize_tuple.hpp>

#include "mirna.hpp"
#include "gene.hpp"
#include "site.hpp"
#include "seed_match_type.hpp"

class Interaction_graph;

class Mirna_site_arc {    
public:
    // this is required by std::unordered_map, I should find a way to keep this private
    Mirna_site_arc() {};
    Seed_match_type seed_match_type;
    double context_score, weighted_context_score;
    bool conserved;

    Mirna_site_arc(Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved);
    Mirna_site_arc(const Mirna_site_arc & obj);
    Mirna_site_arc & operator=(Mirna_site_arc obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);

    friend void swap(Mirna_site_arc & obj1, Mirna_site_arc & obj2);
};

class Interaction_graph {
public:    
    std::unordered_map<boost::tuple<Mirna_id, Gene_id, unsigned int, unsigned int>, Site *> sites_by_location;
    // gene_site arcs
    std::unordered_map<Gene_id, std::list<Site *>> gene_to_sites_arcs;
    // mirna_site arcs
    std::unordered_map<std::pair<Mirna_id, Site *>, Mirna_site_arc> mirna_site_arcs;
    std::unordered_map<Mirna_id, std::list<Site *>> mirna_to_sites_arcs;
    std::unordered_map<Site *, std::list<Mirna_id>> site_to_mirnas_arcs;
    // mirna_gene arcs
    std::unordered_map<std::pair<Mirna_id, Gene_id>, std::list<Site *>> mirna_gene_arcs;
    std::unordered_map<Mirna_id, std::list<Gene_id>> mirna_to_genes_arcs;
    std::unordered_map<Gene_id, std::list<Mirna_id>> gene_to_mirnas_arcs;
    unsigned long long rows_processed = 0;
    unsigned long long rows_skipped = 0;

    Interaction_graph();
    void build_interaction_graph(std::set<Mirna_id> & mirnas, std::set<Gene_id> & genes);
    Site * get_site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end);
    void print_statistics();
    void export_interactions_data(std::string patient_folder);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    ~Interaction_graph();
};

#include "interaction_graph.tpp"

#endif // INTERACTION_GRAPH_H
