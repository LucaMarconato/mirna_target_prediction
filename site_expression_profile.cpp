#include "site_expression_profile.hpp"

#include <list>

#include "gene_expression_profile.hpp"
#include "interaction_graph.hpp"

#define Sep Site_expression_profile

Sep::Site_expression_profile() {}

Sep::Site_expression_profile(Gene_expression_profile & gep, Interaction_graph & ig)
{
    for(auto & e : gep.profile) {
        for(auto & l : ig.gene_to_sites_arcs[e.first]) {
            this->profile[l] = e.second;
        }
    }
}
