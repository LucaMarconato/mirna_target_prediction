#ifndef MIRNA_EXPRESSION_PROFILE_H
#define MIRNA_EXPRESSION_PROFILE_H

#include <unordered_map>
#include <string>

#include "expression_profile.hpp"
#include "mirna.hpp"

class Mirna_expression_profile : Expression_profile {
public:
    std::unordered_map<Mirna_id, Expression> profile;
    unsigned long long recognized_distinct_mirnas = 0;
    unsigned long long recognized_reads = 0;
    unsigned long long not_recognized_distinct_mirnas = 0;
    unsigned long long not_recognized_reads = 0;
    unsigned long long total_distinct_mirnas;
    unsigned long long total_reads;
    unsigned long long filtered_out_distinct_mirnas = 0;
    unsigned long long filtered_out_reads = 0;

    // void load_from_gdc_file(std::string filename, std::string patient_folder);
    void load_from_gdc_file(std::string tissue, std::string patient_folder);
    void print_statistics();
    void filter(double threshold_rpm);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#include "mirna_expression_profile.tpp"

#endif // MIRNA_EXPRESSION_PROFILE_H
