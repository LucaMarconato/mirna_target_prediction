#ifndef SITE_H
#define SITE_H

#include <iostream>
#include <string>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/functional/hash.hpp>

#include "mirna.hpp"

// TODO: move the enum inside Mirna, inside this class
class Seed_match_type {
public:
    Mirna::Mirna_matching match_type;
    friend std::ostream & operator<<(std::ostream & stream, const Seed_match_type & o);
};

class Site {
private:
    Site();
public:
    unsigned int utr_start, utr_end;
    Seed_match_type seed_match_type;
    double context_score, weighted_context_score;
    bool conserved;

    /*
      this function rewrites the file data/processed/scored_interactions.tsv 
      into a new file in which the columns mirna_family, gene_id, gene_symbol, transcript_id and seed_match_type are substituted with integers
    */
    static void reduce_size_of_scored_interactions_file();

    /*
      this function is used to see if reduce_size_of_scored_interactions_file(),
      by checking if the generated data corresponds to the original one
    */
    static void regenerate_original_file();
    
    friend std::ostream & operator<<(std::ostream & stream, const Site & o);
    template< class LeftType, class RightType, bool force_mutable >
    friend class boost::bimaps::relation::detail::relation_storage;   
};

#endif // SITE_H
