#ifndef MIRNA_EXPRESSION_PROFILE_H
#define MIRNA_EXPRESSION_PROFILE_H

#include <string>
#include <unordered_map>

#include "expression_profile.hpp"
#include "mirna.hpp"

class Mirna_expression_profile : Expression_profile
{
public:
    std::unordered_map<Mirna_id, Expression> profile;
    unsigned long long distinct_mirnas = 0;
    double total_reads = 0;
    unsigned long long filtered_out_distinct_mirnas = 0;
    double filtered_out_reads = 0;

    Mirna_expression_profile();
    Mirna_expression_profile(const Mirna_expression_profile& obj);
    friend void swap(Mirna_expression_profile& obj1, Mirna_expression_profile& obj2);
    Mirna_expression_profile& operator=(Mirna_expression_profile obj);
    void load_from_file(std::string patient_folder);
    void print_statistics();
    void filter(double threshold_rpm);
    template <class Archive>
    void serialize(Archive& ar, const unsigned int);
};

#include "mirna_expression_profile.tpp"

#endif // MIRNA_EXPRESSION_PROFILE_H
