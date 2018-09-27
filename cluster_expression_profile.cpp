#include "cluster_expression_profile.hpp"

#include <list>

#include "gene_expression_profile.hpp"
#include "interaction_graph.hpp"

#define Cep Cluster_expression_profile

Cep::Cluster_expression_profile() {}

Cep::Cluster_expression_profile(Gene_expression_profile & gep, Interaction_graph & ig)
{
    unsigned long long total_reads = 0;
    for(auto & e : gep.profile) {
        if(ig.gene_to_clusters_arcs.find(e.first) != ig.gene_to_clusters_arcs.end()) {
            for(auto & l : ig.gene_to_clusters_arcs.at(e.first)) {
                // experiment, with no biological meaning
                // if(this->profile.find(l) == this->profile.end()) {
                //     this->profile[l] = Expression(Reads(0));
                // }
                // this->profile[l] = Reads(this->profile.at(l).to_reads() + e.second.to_reads());
                this->profile[l] = e.second;
                total_reads += e.second.to_reads();
            }
        }
    }

    for(auto & e : this->profile) {
        e.second.normalize_reads(total_reads);
        // note that it is improper to talk about RPM for clusters, so we set this variable
        e.second.valid_rpm = false;
    }
}
