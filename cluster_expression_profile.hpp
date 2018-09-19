#ifndef CLUSTER_EXPRESSION_PROFILE_H
#define CLUSTER_EXPRESSION_PROFILE_H

#include <unordered_map>

#include "expression_profile.hpp"
#include "cluster.hpp"

class Gene_expression_profile;
class Interaction_graph;

class Cluster_expression_profile : Expression_profile {
public:
    std::unordered_map<Cluster *, Expression> profile;
    Cluster_expression_profile();
    Cluster_expression_profile(Gene_expression_profile & gep, Interaction_graph & ig);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#include "cluster_expression_profile.tpp"

#endif // CLUSTER_EXPRESSION_PROFILE_H
