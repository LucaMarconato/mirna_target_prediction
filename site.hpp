#ifndef SITE_H
#define SITE_H

#include <iostream>
#include <string>
#include <list>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/functional/hash.hpp>

#include "seed_match_type.hpp"
#include "mirna.hpp"
#include "gene.hpp"

class Site {
private:
    Site();
public:
    Mirna_id mirna_id;
    Gene_id gene_id;
    unsigned int utr_start, utr_end;
    Seed_match_type seed_match_type;
    double context_score, weighted_context_score;
    bool conserved;
    
    static std::list<Site *> all_sites;
    
    Site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end, Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    /*
      this function rewrites the file data/processed/scored_interactions.tsv 
      into a new file in which the columns mirna_family, gene_id, gene_symbol, transcript_id and seed_match_type are substituted with integers
    */
    static void reduce_size_of_scored_interactions_file();
    /*
      this function is used to see if reduce_size_of_scored_interactions_file(),
      by checking if the generated data corresponds to the original one
    */
    // static void regenerate_original_file();
    static Site * new_site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end, Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved);
    static void delete_all_sites();
    static void save_all_sites();
    static void load_all_sites();

    friend std::ostream & operator<<(std::ostream & stream, const Site & o);
    template< class LeftType, class RightType, bool force_mutable >
    friend class boost::bimaps::relation::detail::relation_storage;
    // for calling the private default constructor
    friend class boost::serialization::access;
};

namespace std {
    template <>
    struct hash<std::pair<Mirna_id, Site *>>
    {
        size_t operator()(const std::pair<Mirna_id, Site *> & e) const noexcept
        {
            size_t seed = 0;
            boost::hash_combine(seed, e.first);
            boost::hash_combine(seed, e.second);
            return seed;
        }
    };
}

namespace std {
    template <>
    struct hash<std::pair<Mirna_id, Gene_id>>
    {
        size_t operator()(const std::pair<Mirna_id, Gene_id> & e) const noexcept
        {
            size_t seed = 0;
            boost::hash_combine(seed, e.first);
            boost::hash_combine(seed, e.second);
            return seed;
        }
    };
}

#endif // SITE_H

