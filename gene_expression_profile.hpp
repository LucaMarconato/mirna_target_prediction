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
    double recognized_reads = 0;
    unsigned long long not_recognized_distinct_genes = 0;
    double not_recognized_reads = 0;
    unsigned long long total_distinct_genes;
    double total_reads;
    unsigned long long filtered_out_distinct_genes = 0;
    double filtered_out_reads = 0;
    // unsigned long long discarded_reads = 0;

    Gene_expression_profile();
    Gene_expression_profile(const Gene_expression_profile & obj);
    friend void swap(Gene_expression_profile & obj1, Gene_expression_profile & obj2);
    Gene_expression_profile & operator=(Gene_expression_profile obj);
    void load_from_gdc_file(std::string tissue, std::string patient_folder);
    void print_statistics();
    void filter(double threshold_rpm);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};

#include "gene_expression_profile.tpp"

#endif // GENE_EXPRESSION_PROFILE_H
