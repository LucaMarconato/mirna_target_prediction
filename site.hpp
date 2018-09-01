#ifndef SITE_H
#define SITE_H

#include <iostream>
#include <string>
#include <unordered_map>
// #include <tuple>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "mirna.hpp"
#include "gene.hpp"

class Interaction_graph;

class Site {
private:
    Site();
    Site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end);
public:
    Mirna_id mirna_id;
    Gene_id gene_id;
    unsigned int utr_start, utr_end;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    /*
      this function rewrites the file data/processed/scored_interactions.tsv 
      into a new file in which the columns mirna_family, gene_id, gene_symbol, transcript_id and seed_match_type are substituted with integers
    */
    static void reduce_size_of_scored_interactions_file();

    friend std::ostream & operator<<(std::ostream & stream, const Site & o);
    template< class LeftType, class RightType, bool force_mutable >
    friend class boost::bimaps::relation::detail::relation_storage;
    // for calling the private default constructor
    friend class boost::serialization::access;
    // for calling the private non-default constructor
    friend class Interaction_graph;
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

namespace std {
    template <>
    struct hash<boost::tuple<Mirna_id, Gene_id, unsigned int, unsigned int>>
    {
        size_t operator()(const boost::tuple<Mirna_id, Gene_id, unsigned int, unsigned int> & e) const noexcept
        {
            size_t seed = 0;
            boost::hash_combine(seed, e.get<0>());
            boost::hash_combine(seed, e.get<1>());
            boost::hash_combine(seed, e.get<2>());
            boost::hash_combine(seed, e.get<3>());
            return seed;
        }
    };
}

#include "site.tpp"

#endif // SITE_H
