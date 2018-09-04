#ifndef SITE_EXPRESSION_PROFILE_H
#define SITE_EXPRESSION_PROFILE_H

#include <unordered_map>

#include "expression_profile.hpp"
#include "site.hpp"

class Gene_expression_profile;
class Interaction_graph;

class Site_expression_profile : Expression_profile {
public:
    std::unordered_map<Site *, Expression> profile;
    Site_expression_profile();
    Site_expression_profile(Gene_expression_profile & gep, Interaction_graph & ig);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#include "site_expression_profile.tpp"

#endif // SITE_EXPRESSION_PROFILE_H
