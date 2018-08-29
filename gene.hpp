#ifndef GENE_H
#define GENE_H

#include <string>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/functional/hash.hpp>

class Gene {
public:
    static boost::bimap<Gene, int> gene_id_dictionary;
    
    std::string gene_id;
    int gene_id_version;
    std::string gene_symbol;
    std::string transcript_id;
    int transcript_id_version;
    
    Gene(std::string gene_id, int gene_id_version, std::string gene_symbol, std::string transcript_id, int transcript_id_version);
    static void initialize_gene_dictionary();
    static void print_gene_dictionary(unsigned int max_rows = -1);
    friend bool operator<(Gene const & lhs, Gene const & rhs);
};

class Site {
    
};

namespace std {
    template<>
    struct hash<Gene>
    {
        size_t operator () (const Gene & e) const noexcept
        {
            size_t seed = 0;
            // for the moment we use only this key, in case we can add some of them
            boost::hash_combine(seed, e.gene_id);
            return seed;
        }
    };
}

#endif // GENE_H
