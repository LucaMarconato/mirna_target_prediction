#ifndef GENE_EXPRESSION_PROFILE_H
#define GENE_EXPRESSION_PROFILE_H

#include <unordered_map>
#include <string>

#include "expression_profile.hpp"
#include "gene.hpp"

class Gene_expression_profile : Expression_profile {
public:
    std::unordered_map<Gene_id, Expression> profile;
    unsigned long long recognized_distinct_genes = 0;
    unsigned long long recognized_reads = 0;
    unsigned long long not_recognized_distinct_genes = 0;
    unsigned long long not_recognized_reads = 0;
    unsigned long long total_distinct_genes;
    unsigned long long total_reads;
    // unsigned long long discarded_reads = 0;
    
    void load_from_gdc_file(std::string filename, std::string patient_folder);
    void print_statistics();
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#include "gene_expression_profile.tpp"

#endif // GENE_EXPRESSION_PROFILE_H
