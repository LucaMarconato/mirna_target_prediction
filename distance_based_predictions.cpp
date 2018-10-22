#include "matchings_predictor.hpp"

struct Site_flattened {
    unsigned int utr_start;
    // unsigned int utr_end; // not used, you may want to use it if employing it in the clusterization procedure in Interaction_graph
    double r_ijk;
    double p_j_downregulated_given_ijk;
};

void Matchings_predictor::compute_distance_based_predictions()
{
    // TODO
    // for gene to site arcs
    // determine the size of site_flattened
    // find the maximum size
    // create an array of site_flattened of length equal to the maximum size found
    // fill the array
    // recursively make the precictions using the array
    // store the prediction on the p_j_downregulated
}

// double compute_distance_based_cluster_contribution(std::unordered_map<Cluster *, double> )

// double Matchings_predictor::compute_p_j_downregulated_considering_distance(double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened, Cluster ** clusters_in_j_flattened, int clusters_count)
// {
    // std::unordered_map<Cluster *, double> cluster_contribution;

    // for(auto & e : this->patient.interaction_graph.gene_to_clusters_arcs) {
    //     Gene_id loop_gene_id = e.first;
    //     if(gene_id != loop_gene_id) {
    //         continue;
    //     }
    //     auto & clusters = e.second;
    //     for(Cluster * cluster : clusters) {
            
    //     }
    //     double sum = this->iteratively_compute_p_j_downregulated(p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened, clusters.size());
    //     this->p_j_downregulated_values[gene_id] = sum;
    // }
    // for(Cluster * cluster : this->p_c_bound_values) {
        
    // }
    // compute_distance_based_cluster_contribution(cluster_contribution, gene_id, p_j_downregulated_given_c_bound_values, p_c_bound_values;
    // return -1;
// }
