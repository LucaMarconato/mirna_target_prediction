#include "matchings_predictor.hpp"

struct Binding_flattened {
    unsigned long utr_start;
    // unsigned long utr_end; // not used, you may want to use it if employing it in the clusterization procedure in Interaction_graph
    double r_ijk;
    double p_j_downregulated_given_ijk;
};

typedef std::list<Binding_flattened> Distance_based_cluster;

void Matchings_predictor::compute_distance_based_predictions()
{
    unsigned long long total_possible_bindings_in_a_gene_max_value = 0;
    for(auto & e : this->patient.interaction_graph.gene_to_sites_arcs) {
        Gene_id gene_id = e.first;
        auto & sites = e.second;
        // j is tought to be the index of the current gene
        unsigned long total_possible_bindings_for_j = 0;
        for(Site * site : sites) {
            auto & mirnas_for_site = this->patient.interaction_graph.site_to_mirnas_arcs;
            total_possible_bindings_for_j += mirnas_for_site.size();
        }
        if(total_possible_bindings_in_a_gene_max_value <= total_possible_bindings_for_j) {
            total_possible_bindings_in_a_gene_max_value = total_possible_bindings_for_j;
        }        
    }

    Binding_flattened * bindings_in_j_flattened = new Binding_flattened [total_possible_bindings_in_a_gene_max_value];
    Binding_flattened * bindings_in_the_current_distance_based_cluster = new Binding_flattened [total_possible_bindings_in_a_gene_max_value];
    // debug purpose, to print the current call of the recursive procedure
    bool * b = new bool [total_possible_bindings_in_a_gene_max_value];
    const unsigned long threshold_for_exponential_algorithm = 20;
    double ** distance_based_enhance_matrix = new double * [threshold_for_exponential_algorithm];
    for(int k = 0; k < threshold_for_exponential_algorithm; k++) {
        distance_based_enhance_matrix[k] = new double [threshold_for_exponential_algorithm];
        for(int l = 0; l < threshold_for_exponential_algorithm; l++) {
            distance_based_enhance_matrix[k][l] = 0;
        }
    }
    /*
      Since the memory of the list is not freed upon calling clear(), placing these variables outside the loop reduces the number of malloc (unless the compiler does some optimization). 
      TODO: check if this affects the performance, if it does not, move this variable within the appropriate loops and not outside here
    */
    std::list<Distance_based_cluster> distance_based_clusters;
    Distance_based_cluster distance_based_cluster;
    for(auto & e : this->patient.interaction_graph.gene_to_sites_arcs) {
        Gene_id gene_id = e.first;
        auto & sites = e.second;

        // find all possible interactions with the gene and store them in a flattened form
        unsigned long i = 0;
        for(Site * site : sites) {
            auto & mirnas_for_site = this->patient.interaction_graph.site_to_mirnas_arcs;
            for(Mirna_id mirna_id : mirnas_for_site) {
                bindings_in_j[i].utr_start = site->utr_start;
                auto p = std::make_pair(mirna_id, site);
                bindings_in_j[i].r_ijk = this->r_ijk_values.at(p);
                Mirna_site_arc & mirna_site_arc = this->patient.interaction_graph.mirna_site_arcs.at(p);
                double context_score = mirna_site_arc.context_score;
                double probability_of_downregulation = 1 - std::pow(2, context_score);
                bindings_in_j[i].p_j_downregulated_given_ijk = probability_of_downregulation;
                i++;
            }            
        }

        // create all distance based clusters
        std::sort(bindings_in_the_current_distance_based_cluster.begin(),
                  bindings_in_the_current_distance_based_cluster.end(),
                  [](const Binding_flattened & a, const Binding_flattened & b)
                  {
                      return a.utr_start < b.utr_start;
                  });       
        // note that utr_end is not used there, you may want to use it
        distance_based_clusters.clear();
        distance_based_cluster.clear();
        unsigned long last_utr_start = 0;
        for(int k = 0; k < i; k++) {
            Binding_flattened * binding = bindings_in_j_flattened + k;
            if(binding->utr_start - last_utr_start > Global_parameters::threshold_for_overlapping_sites) {
                if(distance_based_cluster == 0) {
                    std::cerr << "error: distance_based_cluster = " << distance_based_cluster << "\n";
                    exit(1);
                }
                distance_based_clusters.push_back(distance_based_cluster);
                distance_based_cluster.clear();
            }
            distance_based_cluster.push_back(binding);
            last_utr_start = binding->last_utr_start;
        }
        if(distance_based_cluster.size() > 0) {
            distance_based_clusters.push_back(distance_based_cluster);
        }

        // if one of the distance based clusters is too long, then skip the computation since it would require to much time (exponential algorithm)
        bool skip_the_computation = false;
        for(auto & current_distance_based_cluster : distance_based_clusters) {
            if(current_distance_based_cluster.size() > threshold_for_exponential_algorithm) {
                this->genes_skipped_by_the_distance_based_predictor.push_back(gene_id);
                skip_the_computation = true;
            }
        }

        if(!skip_the_computation) {
            // compute the matrix describing how sites enhance depending on the distance
            /*
              Note that an apparently innocent modification as changing the conditions on Global_parameters::threshold_for_overlapping_sites from <= to < in interaction_graph.cpp and in this file would make the use of this matrix inconsisent, due to the start values for the variables k and l, that in that case would require to remove the "+ 1"
            */
            for(unsigned long k = Global_parameters::threshold_for_overlapping_sites + 1; k < threshold_for_exponential_algorithm; k++) {
                for(unsigned long l = Global_parameters::threshold_for_exponential_algorithm + 1; l < threshold_for_exponential_algorithm; l++) {
                    unsigned long utr_start0 = current_distance_based_cluster[k].utr_start;
                    unsigned long utr_start1 = current_distance_based_cluster[l].utr_start;
                    distance_based_enhance_matrix[k][l] = distance_based_enchange(std::abs(utr_start0, utr_start1));
                }
            }
            
            // for each distance based cluster, create a flattened version and process it to compute the probability of down-regualation
            double total_p_j_downregulated = 0;
            for(auto & current_distance_based_cluster : distance_based_clusters) {
                int k = 0;
                for(auto & binding_flattened : current_distance_based_cluster) {
                    bindings_in_the_current_distance_based_cluster[k] = binding_flattened;
                    k++;
                }
                this->recusively_compute_p_j_downregulated_distance_based(b, 0, current_distance_based_cluster.size(),
                                                                          &total_p_j_downregulated, 1, 1,
                                                                          bindings_in_the_current_distance_based_cluster,
                                                                          distance_based_enhance_matrix,
                                                                          -1 - Global_parameters::threshold_for_overlapping_sites);
            }
            this->p_j_downregulated_values[gene_id] = total_p_j_downregulated;   
        }
    }
    delete [] bindings_in_j;
    delete [] bindings_in_the_current_distance_based_cluster;
    delete [] b;
    for(int k = 0; k < threshold_for_exponential_algorithm; k++) {
        delete [] distance_based_enhance_matrix[k];
    }
    delete [] distance_based_enhance_matrix;
    // recursively make the precictions using the array
    // store the prediction on the p_j_downregulated
}

void Matchings_predictor::recusively_compute_p_j_downregulated_distance_based(bool * b, int level, int max_level, double * sum, double p_j_downregulated_given_b, double p_b, Binding_flattened * bindings_flattened, double ** distance_based_enhance_matrix, unsigned long latest_utr_start)
{
    if(level == max_level) {
        p_j_downregulated_given_b = 1 - p_j_downregulated_given_b;
        *sum += p_j_downregulated_given_b * p_b;
        // std::cout << "b:\n";
        // for(int i = 0; i < max_level; i++) {
        //     std::cout << b[i] << " ";
        // }
        // std::cout << p_j_downregulated_given_b * p_b << "\n\n";
    } else {
        // TODO: COMPLETE HERE!!!!!!!!
        // double p_j_downregulated_given_c_bound = p_j_downregulated_given_c_bound_values_flattened[level];
        // double p_c_bound = p_c_bound_values_flattened[level];

        // double new_p_j_downregulated_given_b;
        // double new_p_b;

        // // b is just used for debug purposes and can be removed
        // // b[level] = 0;
        // new_p_j_downregulated_given_b = p_j_downregulated_given_b;
        // new_p_b = p_b * (1 - p_c_bound);
        // recusively_compute_p_j_downregulated(b, level + 1, max_level, sum, new_p_j_downregulated_given_b, new_p_b, p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened);
        
        // // b[level] = 1;
        // new_p_j_downregulated_given_b = p_j_downregulated_given_b * (1 - p_j_downregulated_given_c_bound);
        // new_p_b = p_b * p_c_bound;
        // recusively_compute_p_j_downregulated(b, level + 1, max_level, sum, new_p_j_downregulated_given_b, new_p_b, p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened);
    }
}

double Matchings_predictor::distance_based_enchange(unsigned long distance)
{
    if(distance <= Global_parameters::threshold_for_overlapping_sites) {
        std::cerr << "error: distance = " << distance << "\n";
        exit(1);
    }
    // TODO: fill with real values, these are just very very rough approximations
    double additional_log_fold_change;
    if(distance < 15) {
        additional_log_fold_change = -1;
    } else if(distance < 30) {
        additional_log_fold_change = -2;
    } else if(distance < 50) {
        additional_log_fold_change = -1;
    } else {
        additional_log_fold_change = 0;
    }

    double additional_probability_of_downregulation = 1 - std::exp(2, additional_log_fold_change);
    return additional_probability_of_downregulation;
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
