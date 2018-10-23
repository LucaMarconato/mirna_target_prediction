#include "matchings_predictor.hpp"

#include <iomanip>

#include "global_parameters.hpp"

struct Binding_flattened {
    unsigned long utr_start;
    // unsigned long utr_end; // not used, you may want to use it if employing it in the clusterization procedure in Interaction_graph
    double p_ijk_bound;
    double p_j_downregulated_given_ijk_bound;
};

typedef std::vector<Binding_flattened> Distance_based_cluster;

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

    Binding_flattened * bindings_in_j = new Binding_flattened [total_possible_bindings_in_a_gene_max_value];
    Binding_flattened * bindings_in_the_current_distance_based_cluster = new Binding_flattened [total_possible_bindings_in_a_gene_max_value];
    // debug purpose, to print the current call of the recursive procedure
    short * b = new short [total_possible_bindings_in_a_gene_max_value];
    const unsigned long threshold_for_exponential_algorithm = 20;
    double ** distance_based_enhancement_matrix = new double * [threshold_for_exponential_algorithm];
    for(int k = 0; k < threshold_for_exponential_algorithm; k++) {
        distance_based_enhancement_matrix[k] = new double [threshold_for_exponential_algorithm];
        for(int l = 0; l < threshold_for_exponential_algorithm; l++) {
            distance_based_enhancement_matrix[k][l] = 0;
        }
    }
    /*
      Since the memory of the list is not freed upon calling clear(), placing these variables outside the loop reduces the number of malloc (unless the compiler does some optimization). 
      TODO: check if this affects the performance, if it does not, move this variable within the appropriate loops and not outside here
    */
    std::list<Distance_based_cluster> distance_based_clusters;
    Distance_based_cluster distance_based_cluster;
    // either performing the full computation, either skipping the computation because too many matchings are available
    unsigned long genes_processed = 0;
    for(auto & e : this->patient.interaction_graph.gene_to_clusters_arcs) {
        Gene_id gene_id = e.first;
        if(gene_id != 12484) {
            continue;            
        }        
        std::cout << "gene_id = " << gene_id << "\n";
        auto & clusters = e.second;

        // find all possible interactions with the gene and store them in a flattened form
        unsigned long i = 0;
        for(Cluster * cluster : clusters) {
            for(Site * site : cluster->sites) {
                auto & mirnas_for_site = this->patient.interaction_graph.site_to_mirnas_arcs.at(site);
                for(auto & mirna_id : mirnas_for_site) {
                    bindings_in_j[i].utr_start = site->utr_start;
                    auto p = std::make_pair(mirna_id, site);
                    bindings_in_j[i].p_ijk_bound = this->r_ijk_values.at(p)/this->original_cluster_profile.at(cluster);
                    Mirna_site_arc & mirna_site_arc = this->patient.interaction_graph.mirna_site_arcs.at(p);
                    double context_score = mirna_site_arc.context_score;
                    double probability_of_downregulation = 1 - std::pow(2, context_score);
                    bindings_in_j[i].p_j_downregulated_given_ijk_bound = probability_of_downregulation;
                    i++;
                }    
            }
        }

        // create all distance based clusters
        std::sort(bindings_in_j, bindings_in_j + i,
                  [](const Binding_flattened & a, const Binding_flattened & b)
                  {
                      return a.utr_start < b.utr_start;
                  });       
        // note that utr_end is not used there, you may want to use it
        distance_based_clusters.clear();
        distance_based_cluster.clear();
        unsigned long latest_utr_start = bindings_in_j[0].utr_start;
        for(int k = 0; k < i; k++) {
            Binding_flattened * binding = bindings_in_j + k;
            if(binding->utr_start - latest_utr_start > Global_parameters::threshold_for_overlapping_sites) {
                if(distance_based_cluster.size() == 0) {
                    std::cerr << "error: distance_based_cluster.size() = " << distance_based_cluster.size() << "\n";
                    exit(1);
                }
                distance_based_clusters.push_back(distance_based_cluster);
                distance_based_cluster.clear();
            }
            distance_based_cluster.push_back(*binding);
            latest_utr_start = binding->utr_start;
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
            // for each distance based cluster, create a flattened version, compute the distance based matrix, and finally compute the probability of down-regualation for that particular distance based cluster
            double ratio_j_not_downregulated = 1.0;
            for(auto & current_distance_based_cluster : distance_based_clusters) {
                double p_j_downregulated_given_current_distance_based_cluster = 0;
                // create the flattened version of the distance based cluster
                unsigned long k = 0;
                for(auto & binding_flattened : current_distance_based_cluster) {
                    bindings_in_the_current_distance_based_cluster[k] = binding_flattened;
                    k++;
                }
                
                // compute the upper triangular matrix describing how sites enhance depending on the distance
                std::cout << "current_distance_based_cluster.size() = " << current_distance_based_cluster.size() << "\n";
                for(unsigned long k = 0; k < current_distance_based_cluster.size(); k++) {
                    for(unsigned long l = k + 1; l < current_distance_based_cluster.size(); l++) {
                        unsigned long utr_start0 = current_distance_based_cluster[k].utr_start;
                        unsigned long utr_start1 = current_distance_based_cluster[l].utr_start;
                        distance_based_enhancement_matrix[k][l] = distance_based_enchange(std::abs((long long)utr_start0 - (long long)utr_start1));
                    }
                }

                std::cout << "starting the recursion, current_distance_based_cluster.size() = " << current_distance_based_cluster.size() << "\n";
                this->recusively_compute_p_j_downregulated_distance_based(b, 0, current_distance_based_cluster.size(),
                                                                          &p_j_downregulated_given_current_distance_based_cluster, 1, 1,
                                                                          bindings_in_the_current_distance_based_cluster,
                                                                          distance_based_enhancement_matrix,
                                                                          // WARNING: is the comparison of an unsigned long against -1 guaranteed by the standard? TODO: check
                                                                          -1,
                                                                          true);
                std::cout << "current_distance_based_cluster:\n";
                for(int k = 0; k < current_distance_based_cluster.size(); k++) {
                    std::cout << "k = 0\n";
                    std::cout << "utr_start = " << current_distance_based_cluster[k].utr_start << "\n";
                    std::cout << "p_ijk_bound = " << current_distance_based_cluster[k].p_ijk_bound << "\n";
                    std::cout << "p_j_downregulated_given_ijk_bound = " << current_distance_based_cluster[k].p_j_downregulated_given_ijk_bound << "\n";
                    std::cout << "\n";
                }
                std::cout << "distance_based_enhancement_matrix\n";
                for(int k = 0; k < current_distance_based_cluster.size(); k++) {
                    for(int l = 0; l < current_distance_based_cluster.size(); l++) {
                        std::cout << distance_based_enhancement_matrix[k][l] << " ";
                    }
                    std::cout << "\n";
                }
                std::cout << "\n\n";
                double newly_downregulated = ratio_j_not_downregulated * p_j_downregulated_given_current_distance_based_cluster;
                ratio_j_not_downregulated -= newly_downregulated;
            }
            double p_j_downregulated = 1 - ratio_j_not_downregulated;
            // test to check if the predictions made without considering the distance fit the one when the distance is considered but the distance contribution is set to zero
            if(std::abs(this->p_j_downregulated_values.at(gene_id) - p_j_downregulated) > Global_parameters::epsilon) {
                std::cout << "abs(this->p_j_downregulated_values.at(gene_id) - p_j_downregulated) = " << abs(this->p_j_downregulated_values.at(gene_id) - p_j_downregulated) << "\n";
                std::cout << "genes_processed = " << genes_processed << " (" << this->genes_skipped_by_the_distance_based_predictor.size() << " skipped)\n";
                std::cout << "distance_based_clusters.size() = " << distance_based_clusters.size() << "\n";
                exit(1);
            }
            this->p_j_downregulated_values[gene_id] = p_j_downregulated;
            genes_processed++;
            std::cout << "genes_processed = " << genes_processed << " (" << this->genes_skipped_by_the_distance_based_predictor.size() << " skipped)\n";
        }
    }

    delete [] bindings_in_j;
    delete [] bindings_in_the_current_distance_based_cluster;
    delete [] b;
    for(int k = 0; k < threshold_for_exponential_algorithm; k++) {
        delete [] distance_based_enhancement_matrix[k];
    }
    delete [] distance_based_enhancement_matrix;
    // recursively make the precictions using the array
    // store the prediction on the p_j_downregulated


    // TODO: TEST THE DISTANCE BASED PREDICTIONS IN THE CASE IN WHICH THE RESULT SHOULD COINCIDE TO THE USUAL PREDICITONS!!!!!
}

void Matchings_predictor::recusively_compute_p_j_downregulated_distance_based(short * b, int level, int max_level,
                                                                              double * sum, double p_j_downregulated_given_b, double p_b,
                                                                              Binding_flattened * bindings_flattened,
                                                                              double ** distance_based_enhancement_matrix,
                                                                              unsigned long latest_bound_level,
                                                                              bool the_current_cluster_is_free_to_bind)
{
    if(level == max_level) {
        p_j_downregulated_given_b = 1 - p_j_downregulated_given_b;
        *sum += p_j_downregulated_given_b * p_b;
        std::cout << "b: ";
        for(int i = 0; i < max_level; i++) {
            std::cout << std::setw(2) << b[i] << " ";
        }        
        std::cout << "| " << p_j_downregulated_given_b * p_b << "\n";
    } else {
        double p_j_downregulated_given_ijk_bound = bindings_flattened[level].p_j_downregulated_given_ijk_bound;
        double p_ijk_bound = bindings_flattened[level].p_ijk_bound;
        unsigned long utr_start = bindings_flattened[level].utr_start;
        
        bool overlapping;
        // TODO: BUG HERE!!!!!!!!!!! TO FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE VALUE OF DISTANE_BASED_ENHANCEMENT IS NOT CORRECT RIGHT NOW
        double distance_based_enhancement_factor = 0;
        if(latest_bound_level != -1) {
            unsigned long latest_utr_start = bindings_flattened[latest_bound_level].utr_start;
            overlapping = utr_start - latest_utr_start > Global_parameters::threshold_for_overlapping_sites;
            distance_based_enhancement_factor = 1 - distance_based_enhancement_matrix[latest_bound_level][level];
        } else {
            overlapping = false;
        }

        if(level > 0) {
            if(utr_start - bindings_flattened[level - 1].utr_start > Global_parameters::threshold_for_overlapping_sites) {
                the_current_cluster_is_free_to_bind = true;
            }
        }
        // TODO: BUG HERE!!!!!!!!!!!!!!!! WE NEED TO LOOP IN THE CLUSTER (CLUSTER, NOT DISTANCE BASED CLUSTER) TO FIND ALL THE INCOMPATIBLE EVENT AND SPAWN A NEW RECURSIVE CALL FOR EACH OF THEM (ALSO THE EVENT IN WHICH NONE OF THEM HAPPENS MUST BE CONSIDERED)
        if(!overlapping && the_current_cluster_is_free_to_bind) {
            double new_p_j_downregulated_given_b;
            double new_p_b;
            unsigned long new_bound_level;

            // b is just used for debug purposes and can be removed
            b[level] = 0;
            new_bound_level = latest_bound_level;
            new_p_j_downregulated_given_b = p_j_downregulated_given_b;
            new_p_b = p_b * (1 - p_ijk_bound);
            recusively_compute_p_j_downregulated_distance_based(b, level + 1, max_level,
                                                                sum, new_p_j_downregulated_given_b, new_p_b,
                                                                bindings_flattened,
                                                                distance_based_enhancement_matrix,
                                                                new_bound_level,
                                                                the_current_cluster_is_free_to_bind);
            b[level] = 1;
            new_bound_level = level;
            // not that the indexes of the matrix must be [latest_bound_level][bound_level] and not the opposite because the matrix is upper triangular
            new_p_j_downregulated_given_b = p_j_downregulated_given_b * (1 - p_j_downregulated_given_ijk_bound) * distance_based_enhancement_factor;
            new_p_b = p_b * p_ijk_bound;
            recusively_compute_p_j_downregulated_distance_based(b, level + 1, max_level,
                                                                sum, new_p_j_downregulated_given_b, new_p_b,
                                                                bindings_flattened,
                                                                distance_based_enhancement_matrix,
                                                                new_bound_level,
                                                                false);
        } else {
            b[level] = -1;
            recusively_compute_p_j_downregulated_distance_based(b, level + 1, max_level,
                                                                sum, p_j_downregulated_given_b, p_b,
                                                                bindings_flattened,
                                                                distance_based_enhancement_matrix,
                                                                latest_bound_level,
                                                                the_current_cluster_is_free_to_bind);
        }
    }
}

double Matchings_predictor::distance_based_enchange(unsigned long distance)
{    
    // TODO: fill with real values, these are just very very rough approximations
    double additional_log_fold_change;
    // if(distance < 15) {
    //     additional_log_fold_change = -1;
    // } else if(distance < 30) {
    //     additional_log_fold_change = -2;
    // } else if(distance < 50) {
    //     additional_log_fold_change = -1;
    // } else {
    //     additional_log_fold_change = 0;
    // }
    // to test if the distance based predictions give the same result of the usual method when distance does not influence predicitons
    additional_log_fold_change = 0;
    

    double additional_probability_of_downregulation = 1 - std::pow(2, additional_log_fold_change);
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
