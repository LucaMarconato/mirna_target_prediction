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
    Cluster_expression_profile(const Cluster_expression_profile & obj);
    friend void swap(Cluster_expression_profile & obj1, Cluster_expression_profile & obj2);
    Cluster_expression_profile & operator=(Cluster_expression_profile obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};

#include "cluster_expression_profile.tpp"

#endif // CLUSTER_EXPRESSION_PROFILE_H
