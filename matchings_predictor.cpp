#include "matchings_predictor.hpp"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include <marconato/output_buffer/output_buffer.hpp>

#include "omp_disabler.h"
#ifndef OMP_TEMPORARILY_DISABLED
#include <omp.h>
#else
#endif

#include "global_parameters.hpp"
#include "timer.hpp"

Matchings_predictor::Matchings_predictor(Patient & patient) : patient(patient)
{
    // TODO: for the moment I am just considering tumor mirnas and tumor clusters. In the future I should consider the difference between normal and tumor expressions.
    // note that these are copies, we do not want to modify the original expressions
    for(auto & e : this->patient.tumor_mirnas.profile) {
        const Mirna_id & mirna_id = e.first;
        double relative_value = e.second.to_relative_expression();
        this->mirna_profile[mirna_id] = relative_value;
    }
    for(auto & e : this->patient.tumor_clusters.profile) {
        Cluster * cluster = e.first;
        double relative_value = e.second.to_relative_expression();
        this->cluster_profile[cluster] = relative_value;
    }

    /*
      The mirnas do not sum to one because of the filtering procedure, so we normalize the expression profile.
    */
    double sum = 0;
    for(auto & e : this->mirna_profile) {
        sum += e.second;
    }
    for(auto & e : this->mirna_profile) {
        e.second /= sum;
    }

    /*
      The mirna expression proile has been just normalized to 1.
      The clusters are already normalized to 1, by construtction.
      Anyway, we check this last condition explicitly.
    */      
    sum = 0;
    for(auto & e : this->cluster_profile) {
        sum += e.second;
    }
    if(std::abs(sum - 1.0) > Global_parameters::epsilon) {
        std::cerr << "error: std::abs(sum - 1.0) = " << std::abs(sum - 1.0) << "\n";
        exit(1);
    }

    this->original_cluster_profile = this->cluster_profile;
}

void Matchings_predictor::export_interaction_matrix()
{
    std::string filename = "./data/patients/" + this->patient.case_id + "/interaction_matrix.mat";
    std::ofstream out(filename);
    std::stringstream ss;
    unsigned long long cluster_count = this->cluster_profile.size();
    unsigned long long i = 0;
    for(auto & e : this->cluster_profile) {
        ss << "\"" << e.first << "\"" << ((i == cluster_count - 1) ? "\n" : "\t");
        i++;
    }
    ss << "\n";
    std::unordered_map<Mirna_id, unsigned long long> sum_of_rows;
    std::unordered_map<Cluster *, unsigned long long> sum_of_columns;
    for(auto & e : this->mirna_profile) {
        auto & mirna_id = e.first;
        sum_of_rows[mirna_id] = 0;
    }
    for(auto & e : this->cluster_profile) {
        Cluster * cluster = e.first;
        sum_of_columns[cluster] = 0;
    }
        
    for(auto & mirna : this->mirna_profile) {
        auto & mirna_id = mirna.first;
        ss << "\"" << mirna_id << "\"\t";
        unsigned long long i = 0;
        for(auto & e : this->cluster_profile) {
            Cluster * cluster = e.first;
            int matchings = 0;
            for(auto & site : cluster->sites) {
                auto p = std::make_pair(mirna_id, site);
                if(this->patient.interaction_graph.mirna_site_arcs.find(p) != this->patient.interaction_graph.mirna_site_arcs.end()) {
                    matchings++;
                }
            }

            sum_of_rows[mirna_id] = sum_of_rows.at(mirna_id) + matchings;
            sum_of_columns[cluster] = sum_of_columns.at(cluster) + matchings;
            ss << matchings << ((i == cluster_count - 1) ? "\n" : "\t");
            i++;
        }
    }
    out << ss.str();
    out.close();

    filename = "./data/patients/" + this->patient.case_id + "/sum_of_rows.tsv";
    out.open(filename);
    ss.str("mirna_id\tsum_of_row\n");    
    for(auto & e : this->mirna_profile) {
        auto & mirna_id = e.first;
        ss << mirna_id << "\t" << sum_of_rows.at(mirna_id) << "\n";
    }
    out << ss.str();
    out.close();

    filename = "./data/patients/" + this->patient.case_id + "/sum_of_columns.tsv";
    out.open(filename);
    ss.str("cluster_address\tsum_of_column\n");
    for(auto & e : this->cluster_profile) {
        Cluster * cluster = e.first;
        ss << cluster << "\t" << sum_of_columns.at(cluster) << "\n";
    }
    out << ss.str();
    out.close();
}

void Matchings_predictor::compute()
{
    std::cout << "for the moment, just a trivial explicit Euler scheme\n";
    unsigned long long max_steps = 1000;
    double mirna_lambda = 1;
    double cluster_lambda = 1;
    if(Global_parameters::lambda > 1) {
        mirna_lambda = 1.0/Global_parameters::lambda;
    } else {
        cluster_lambda = Global_parameters::lambda;
    }
    // h must be <= 1
    double h = 1;
    bool scaling = true;
    double cumulative_scaling = 1;
    bool logging = true;
    bool export_mirna_expression_profile = true;
    bool export_cluster_expression_profile = true;
    bool export_interaction_matrix = false;
    bool export_partial_predicted_downregulation = false;
    // just to be sure to avoid logical errors in the future
    if(!logging && (export_mirna_expression_profile || export_cluster_expression_profile || export_interaction_matrix || export_partial_predicted_downregulation)) {
        export_mirna_expression_profile = false;
        export_cluster_expression_profile = false;
        export_interaction_matrix = false;
        export_partial_predicted_downregulation = false;
        std::cout << "not logging (logging = false)\n";
    }
    std::stringstream ss0;
    std::ofstream out0;
    std::stringstream ss1_header;
    std::stringstream ss2_header;
    std::string mirna_log_filename_prefix;
    std::string cluster_log_filename_prefix;
    if(logging) {
        std::string filename0 = "./data/patients/" + this->patient.case_id + "/matchings_prediction.tsv";
        out0.open(filename0);
        ss0 << "cumulative_scaling\tavg_mirna_level\tavg_cluster_level\tmirna_total_exchange\tcluster_total_exchange\n";
        if(export_mirna_expression_profile) {
            std::string log_folder = "./data/patients/" + this->patient.case_id + "/mirna_expression_profiles";
            mirna_log_filename_prefix = log_folder + "/mirna_expression_profile_";
            std::string cmd = "mkdir -p " + log_folder;
            std::system(cmd.c_str());
            std::cout << "deleting old log from \"" + log_folder + "\"\n";
            cmd = "rm " + log_folder + "/mirna_expression_profile*";
            std::system(cmd.c_str());
            ss1_header << "mirna_id\trelative_expression\n";   
        }
        if(export_cluster_expression_profile) {
            std::string log_folder = "./data/patients/" + this->patient.case_id + "/cluster_expression_profiles";
            cluster_log_filename_prefix = log_folder + "/cluster_expression_profile_";
            std::string cmd = "mkdir -p " + log_folder;
            std::system(cmd.c_str());
            std::cout << "deleting old log from \"" + log_folder + "\"\n";
            cmd = "rm " + log_folder + "/cluster_expression_profile*";
            std::system(cmd.c_str());

            // you may want to export a file in which each Cluster * (i.e. what I call cluster_address below) is enriched with information to make you able to identify specific clusters
            ss2_header << "cluster_address\trelative_expression\n";
        }
        if(export_interaction_matrix) {
            this->export_interaction_matrix();
        }
        if(export_partial_predicted_downregulation) {
            std::string log_folder = "./data/patients/" + this->patient.case_id + "/predicted_downregulation";
            std::string cmd = "mkdir -p " + log_folder;
            std::system(cmd.c_str());
            std::cout << "deleting old log from \"" + log_folder + "\"\n";
            cmd = "rm " + log_folder + "/p_j_downregulated_values_*";
            std::system(cmd.c_str());
        }
    }    
    
    int latest_percentage = -1;
    // these four variables will be shared among the threads and set either only by the thread 0, either inside a #pragma omp critical
    double mirna_total_exchange;
    double cluster_total_exchange;
    std::unordered_map<Mirna_id, double> new_mirna_profile;
    std::unordered_map<Cluster *, double> new_cluster_profile;
    unsigned long long loop_size;
    int t;

    std::cout << "starting the simulation\n";
#ifndef OMP_TEMPORARILY_DISABLED
    unsigned long long my_num_threads = 4;
#pragma omp parallel num_threads(my_num_threads)
#endif
    {
#ifndef OMP_TEMPORARILY_DISABLED
        int rank = omp_get_thread_num();
#else
        unsigned long long my_num_threads = 1;
        int rank = 0;
#endif
        if(rank == 0) {
            std::cout << "my_num_threads = " << my_num_threads << "\n";
            t = 0;
        }
        for(t = 0; t < max_steps; ) {
            if(rank == 0) {
                int current_percentage = (int)(100*((double)(t)/max_steps));
                if(current_percentage > latest_percentage) {
                    latest_percentage = current_percentage;
                    std::cout << "t/max_steps: " << t << "/" << max_steps << " = " << current_percentage << "%\n";
                    if(logging) {
                        if(export_mirna_expression_profile) {
                            std::stringstream ss1;
                            ss1 << ss1_header.str();
                            for(auto & e : this->mirna_profile) {
                                double adjusted_value = e.second / cumulative_scaling;
                                ss1 << e.first << "\t" << adjusted_value << "\n";
                            }
                            std::stringstream mirna_log_filename;
                            mirna_log_filename << mirna_log_filename_prefix << std::setfill('0') << std::setw(6) << t << ".tsv";
                            std::ofstream out1(mirna_log_filename.str());
                            out1 << ss1.str();
                            out1.close();
                        }
                        if(export_cluster_expression_profile) {
                            std::stringstream ss2;
                            ss2 << ss2_header.str();
                            for(auto & e : this->cluster_profile) {
                                double adjusted_value = e.second / cumulative_scaling;
                                ss2 << e.first << "\t" << adjusted_value << "\n";
                            }
                            std::stringstream cluster_log_filename;
                            cluster_log_filename << cluster_log_filename_prefix << std::setfill('0') << std::setw(6) << t << ".tsv";
                            std::ofstream out2(cluster_log_filename.str());
                            out2 << ss2.str();
                            out2.close();
                        }
                        if(export_partial_predicted_downregulation) {
                            // if t == 0, then export_p_j_downregulated(t) will just export the initial expression profile
                            if(t > 0) {
                                this->compute_probabilities();
                            }
                            this->export_p_j_downregulated(t);
                        }
                    }
                }
                if(scaling) {
                    // perform the scaling
                    double mirna_sum = 0;
                    for(auto & e : this->mirna_profile) {
                        mirna_sum += e.second;
                    }
                    double cluster_sum = 0;
                    for(auto & e : this->cluster_profile) {
                        cluster_sum += e.second;
                    }
                    double sum = std::max(mirna_sum, cluster_sum);
                    for(auto & e : this->mirna_profile) {
                        e.second /= sum;
                    }
                    for(auto & e : this->cluster_profile) {
                        e.second /= sum;
                    }
                    if(sum > 1 + Global_parameters::epsilon) {
                        std::cerr << "error: sum = " << sum << "\n";
                        exit(1);
                    }
                    cumulative_scaling /= sum;
                }
                if(logging) {
                    ss0 << cumulative_scaling << "\t";
                    // compute the average level
                    double avg_mirna_level = 0;
                    for(auto & e : this->mirna_profile) {
                        avg_mirna_level += e.second;
                    }
                    avg_mirna_level /= this->mirna_profile.size();
                    double avg_cluster_level = 0;
                    for(auto & e : this->cluster_profile) {
                        double value = e.second;
                        avg_cluster_level += value;            
                    }
                    avg_cluster_level /= this->cluster_profile.size();
                    double adjusted_avg_mirna_level = avg_mirna_level / cumulative_scaling;
                    double adjusted_avg_cluster_level = avg_cluster_level / cumulative_scaling;
                    ss0 << adjusted_avg_mirna_level << "\t" << adjusted_avg_cluster_level << "\t";
                }
                // perform the exchange
                mirna_total_exchange = 0;
                cluster_total_exchange = 0;
                new_mirna_profile = this->mirna_profile;
                new_cluster_profile = this->cluster_profile;
                loop_size = this->patient.interaction_graph.gene_to_clusters_arcs.size();
            }
#ifndef OMP_TEMPORARILY_DISABLED
#pragma omp barrier
#endif
            auto rank_new_mirna_profile = std::unordered_map<Mirna_id, double>();
            auto rank_new_cluster_profile = std::unordered_map<Cluster *, double>();
            for(auto & e : this->mirna_profile) {
                rank_new_mirna_profile[e.first] = 0;
            }
            for(auto & e : this->cluster_profile) {
                rank_new_cluster_profile[e.first] = 0;
            }
            double rank_mirna_total_exchange = 0;
            double rank_cluster_total_exchange = 0;
            auto it = this->patient.interaction_graph.gene_to_clusters_arcs.begin();
            for(auto & e : rank_new_mirna_profile) {
                e.second = 0;
            }
            int i;            
            for(i = 0; i < loop_size && i < rank; i++) {
                it++;
            }
            
            for(; i < loop_size;) {
                auto & clusters = it->second;
                for(Cluster * cluster : clusters) {
                    for(Site * site : cluster->sites) {
                        for(const Mirna_id & mirna_id : this->patient.interaction_graph.site_to_mirnas_arcs.at(site)) {
                            double exchange = h * this->mirna_profile.at(mirna_id) * this->cluster_profile.at(cluster);
                            if(exchange < 0) {
                                static bool warning_already_presented = false;
                                if(abs(exchange) < Global_parameters::epsilon) {
                                    if(!warning_already_presented) {
                                        warning_already_presented = true;
                                        std::cerr << "warning: exchange = " << exchange << ", silencing further warnings like this\n";   
                                    }
                                } else {
                                    std::cout << "error: exchange = " << exchange << "\n";
                                    exit(1);
                                }
                            }
                            rank_new_mirna_profile.at(mirna_id) -= mirna_lambda * exchange;
                            rank_new_cluster_profile.at(cluster) -= cluster_lambda * exchange;
                            if(cumulative_scaling < 1 - Global_parameters::epsilon) {
                                std::cerr << "error: cumulative_scaling = " << cumulative_scaling << "\n";
                                exit(1);
                            }
                            rank_mirna_total_exchange += (mirna_lambda * exchange / cumulative_scaling);
                            rank_cluster_total_exchange += (cluster_lambda * exchange / cumulative_scaling);
                            // if(mirna_id == 397 && cluster->sites.front()->gene_id == 15242) {
                            //     std::cout << "rank_mirna_total_exchange = " << rank_mirna_total_exchange << ", rank_cluster_total_exchange = " << rank_cluster_total_exchange << "\n";
                            // }
                            auto p = std::make_pair(mirna_id, cluster);
                            if(this->r_ic_values.find(p) == this->r_ic_values.end()) {
                                this->r_ic_values[p] = 0;
                            }

                            // TODO: modify the formula for r_ic in the document
                            if(rank > 0) {
                                std::cerr << "error: I have not parallelized the part of the code accessing r_ic_value yet\n";
                                exit(1);
                            }
                            this->r_ic_values[p] += cluster_lambda * exchange / cumulative_scaling;
                        }
                    }
                }
                
                for(int j = 0; j < my_num_threads && i < loop_size; j++) {
                    it++;
                    i++;
                }
            }
#ifndef OMP_TEMPORARILY_DISABLED
#pragma omp critical
#endif
            {
                for(auto & e : rank_new_mirna_profile) {
                    new_mirna_profile.at(e.first) += e.second;
                }
                for(auto & e : rank_new_cluster_profile) {
                    new_cluster_profile.at(e.first) += e.second;
                }
                mirna_total_exchange += rank_mirna_total_exchange;
                cluster_total_exchange += rank_cluster_total_exchange;
            }
            if(rank == 0) {
                t++;
            }
#ifndef OMP_TEMPORARILY_DISABLED
#pragma omp barrier
#endif
            if(rank == 0) {
                this->mirna_profile = new_mirna_profile;
                this->cluster_profile = new_cluster_profile;
                if(logging) {
                    ss0 << mirna_total_exchange << "\t" << cluster_total_exchange << "\n";
                    out0 << ss0.str();
                    ss0.str("");
                }
            }
        } // end of step of the algorithm
    } // end of the parallelization block
    if(logging) {
        out0.close();   
    }
    this->compute_probabilities();
    this->export_probabilities();
}

void Matchings_predictor::compute_probabilities()
{
    // the function compute_probabilities is called many times if export_partial_predicted_downregulation is true, so we need to reset the probabilities before proceeding
    r_ijk_values.clear();
    p_c_bound_values.clear();
    sum_of_r_ijk_for_cluster.clear();
    p_j_downregulated_given_c_bound_values.clear();
    p_j_downregulated_values.clear();
    
    for(auto & e : this->r_ic_values) {
        // compute r_ijk_values for the sites in the current cluster and binding to the current miRNA
        auto & mirna_id = e.first.first;
        Cluster * cluster = e.first.second;
        double value = e.second;
        // sum_of_r_ijk_for_cluster is needed when computing the denominator of this->p_j_downregulated_given_c_bound_values
        if(this->sum_of_r_ijk_for_cluster.find(cluster) == this->sum_of_r_ijk_for_cluster.end()) {
            this->sum_of_r_ijk_for_cluster[cluster] = 0;
        }
        this->sum_of_r_ijk_for_cluster[cluster] = this->sum_of_r_ijk_for_cluster.at(cluster) + value;
        std::list<std::pair<Mirna_id, Site *>> just_added;
        for(auto & site : cluster->sites) {
            auto & m = this->patient.interaction_graph.mirna_site_arcs;
            auto p = std::make_pair(mirna_id, site);
            if(m.find(p) != m.end()) {
                this->r_ijk_values[p] = value;
                just_added.push_back(p);
            }
        }
        if(just_added.size() == 0) {
            std::cerr << "error: just_added.size() = " << just_added.size() << "\n";
            exit(1);
        }
        for(auto & p : just_added) {            
            this->r_ijk_values[p] = this->r_ijk_values.at(p) / just_added.size();            
        }

        // update p_c_bound_values, note that the values of a cluster are not definitive since we need to excecute the following code once per miRNA and the miRNA is given by the outer loop
        if(this->p_c_bound_values.find(cluster) == this->p_c_bound_values.end()) {
            this->p_c_bound_values[cluster] = 0;
        }        
        this->p_c_bound_values[cluster] = this->p_c_bound_values.at(cluster) + value / (this->original_cluster_profile.at(cluster) * cluster->sites.size());
    }

    // compute p_j_downregulated_given_c_bound_values for the current cluster
    for(auto & e : this->r_ic_values) {
        auto & mirna_id = e.first.first;
        Cluster * cluster = e.first.second;
        Gene_id gene_id = cluster->sites.front()->gene_id;
        std::unordered_map<Site *, std::list<std::pair<Mirna_id, Site *>>> just_added_by_site;
        for(auto & site : cluster->sites) {
            auto & m = this->patient.interaction_graph.mirna_site_arcs;
            auto p = std::make_pair(mirna_id, site);
            if(m.find(p) != m.end()) {
                just_added_by_site[site].push_back(p);
            }
        }
        if(just_added_by_site.size() == 0) {
            std::cerr << "error: just_added_by_site.size() = " << just_added_by_site.size() << "\n";
            exit(1);
        }        
        for(auto & e : just_added_by_site) {
            Site * site = e.first;
            for(auto & p : e.second) {
                auto key0 = std::make_pair(gene_id, cluster);
                if(this->p_j_downregulated_given_c_bound_values.find(key0) == this->p_j_downregulated_given_c_bound_values.end()) {
                    this->p_j_downregulated_given_c_bound_values[key0] = 0;
                }
                double r_ijk = this->r_ijk_values.at(p);
                Mirna_site_arc & mirna_site_arc = this->patient.interaction_graph.mirna_site_arcs.at(p);
                double context_score = mirna_site_arc.context_score;
                double probability_of_downregulation = 1 - std::pow(2, context_score);
                double to_add = r_ijk * probability_of_downregulation / sum_of_r_ijk_for_cluster.at(cluster);
                this->p_j_downregulated_given_c_bound_values[key0] = this->p_j_downregulated_given_c_bound_values.at(key0) + to_add;
            }
        }        
    }

    // compute p_j_downregulated_values
    // this code can be trivially parallelized, if needed
    // this variable is used as a debug purpose, it slows down the program and it can be removed
    // double cluster_limit = 20;
    // std::cout << "computing the probabilities of down-regulation for genes with at most " << cluster_limit << " clusters\n";
    // Timer::start();
    // int genes_debugged = 0;
    // int genes_not_debugged = 0;
    unsigned long max_number_of_clusters_per_gene = 0;
    for(auto & e : this->patient.interaction_graph.gene_to_clusters_arcs) {
        auto & clusters = e.second;
        max_number_of_clusters_per_gene = std::max(max_number_of_clusters_per_gene, clusters.size());
    }
    // b is just used for debug purposes and can be removed
    bool * b = new bool [max_number_of_clusters_per_gene];
    double * p_j_downregulated_given_c_bound_values_flattened = new double [max_number_of_clusters_per_gene];
    double * p_c_bound_values_flattened = new double [max_number_of_clusters_per_gene];    

    int latest_percentage = -1;
    for(auto & e : this->patient.interaction_graph.gene_to_clusters_arcs) {
        // int total_genes = this->patient.interaction_graph.gene_to_clusters_arcs.size();
        // int current_percentage = (int)(100*((double)(genes_debugged + genes_not_debugged)/total_genes));
        // if(current_percentage > latest_percentage) {
        //     latest_percentage = current_percentage;
            // std::cout << "genes_debugged/total_genes: " << genes_debugged + genes_not_debugged << "/" << total_genes << " = " << current_percentage << "%\n";
        // }
        Gene_id gene_id = e.first;
        auto & clusters = e.second;
        int i = 0;
        for(Cluster * cluster : clusters) {
            b[i] = false;
            p_j_downregulated_given_c_bound_values_flattened[i] = this->p_j_downregulated_given_c_bound_values.at(std::make_pair(gene_id, cluster));
            p_c_bound_values_flattened[i] = this->p_c_bound_values.at(cluster);
            i++;
        }
        double sum = this->iteratively_compute_p_j_downregulated(p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened, clusters.size());
        this->p_j_downregulated_values[gene_id] = sum;

        // // safety check for genes with at most "cluster_limit" clusters
        // if(clusters.size() <= cluster_limit) {
        //     double debug_sum = 0;
        //     this->recusively_compute_p_j_downregulated(b, 0, clusters.size(), &debug_sum, 1, 1, p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened);
        //     if(std::abs(sum - debug_sum) > Global_parameters::epsilon) {
        //         std::cerr << "error: sum = " << sum << ", debug_sum = " << debug_sum << ", std::abs(sum - debug_sum) = " << std::abs(sum - debug_sum) << "\n";
        //         std::cout << "describing the instance: clusters.size() = " << clusters.size() << "\n";
        //         std::cout << "p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened\n";
        //         for(int j = 0; j < clusters.size(); j++) {
        //             std::cout << p_j_downregulated_given_c_bound_values_flattened[j] << ", " << p_c_bound_values_flattened[j] << "\n";
        //         }
        //         exit(1);
        //     }            
        //     genes_debugged++;
        // } else {
        //     genes_not_debugged++;
        // }
    }
    delete [] b;
    delete [] p_j_downregulated_given_c_bound_values_flattened;
    delete [] p_c_bound_values_flattened;
    std::cout << "\n";

    // std::cout << "genes_debugged/total_genes: " << genes_debugged << "/" << (genes_debugged + genes_not_debugged) << " = " << ((double)genes_debugged)/(genes_debugged + genes_not_debugged) << "\n";
    // std::cout << "finished, ";
    // Timer::stop();
}

/*
  The algorithm described by this function is exponential in the number of clusters in a gene.
  The array b and the valus sum, p_j_downregulated_given_c_b and p_b are computed incrementally at every recursive call and hold the final value their name represents only when processed in the last level.
*/
void Matchings_predictor::recusively_compute_p_j_downregulated(bool * b, int level, int max_level, double * sum, double p_j_downregulated_given_b, double p_b, double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened)
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
        double p_j_downregulated_given_c_bound = p_j_downregulated_given_c_bound_values_flattened[level];
        double p_c_bound = p_c_bound_values_flattened[level];

        double new_p_j_downregulated_given_b;
        double new_p_b;

        // b is just used for debug purposes and can be removed
        b[level] = 0;
        new_p_j_downregulated_given_b = p_j_downregulated_given_b;
        new_p_b = p_b * (1 - p_c_bound);
        recusively_compute_p_j_downregulated(b, level + 1, max_level, sum, new_p_j_downregulated_given_b, new_p_b, p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened);
        
        b[level] = 1;
        new_p_j_downregulated_given_b = p_j_downregulated_given_b * (1 - p_j_downregulated_given_c_bound);
        new_p_b = p_b * p_c_bound;
        recusively_compute_p_j_downregulated(b, level + 1, max_level, sum, new_p_j_downregulated_given_b, new_p_b, p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened);
    }
}

/*
  The algoirthm described by this functions is polynomial in the number of clusters in a gene.
*/
double Matchings_predictor::iteratively_compute_p_j_downregulated(double * p_j_downregulated_given_c_bound_values_flattened, double * p_c_bound_values_flattened, int clusters_count)
{
    double ratio_j_not_downregulated = 1.0;
    for(int i = 0; i < clusters_count; i++) {
        double newly_downregualated = ratio_j_not_downregulated * p_c_bound_values_flattened[i] * p_j_downregulated_given_c_bound_values_flattened[i];
        ratio_j_not_downregulated -= newly_downregualated;
    }
    return 1 - ratio_j_not_downregulated;
}    

void Matchings_predictor::export_probabilities()
{
    std::ofstream out;
    std::stringstream ss;
    std::string filename;
    
    std::cout << "exporting r_ic_values\n";
    // Timer::start();
    ss.str("");
    ss << "mirna_id\tcluster_address\tr_ic\n";
    for(auto & e : this->r_ic_values) {
        auto & mirna_id = e.first.first;
        Cluster * cluster = e.first.second;
        double value = e.second;
        ss << mirna_id << "\t" << cluster << "\t" << value << "\n";
    }
    filename = "./data/patients/" + this->patient.case_id + "/r_ic_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();
    // std::cout << "finished, ";
    // Timer::stop();

    std::cout << "exporting r_ijk_values\n";
    // Timer::start();
    ss.str("");
    ss << "mirna_id\tsite_address\tr_ijk\n";
    for(auto & e : this->r_ijk_values) {
        auto & mirna_id = e.first.first;
        Site * site = e.first.second;
        double value = e.second;
        ss << mirna_id << "\t" << site << "\t" << value << "\n";
    }
    filename = "./data/patients/" + this->patient.case_id + "/r_ijk_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();
    // std::cout << "finished, ";
    // Timer::stop();

    std::cout << "exporting p_c_bound_values\n";
    // Timer::start();
    ss.str("");
    ss << "cluster_address\tp_c_bound\n";
    for(auto & e : this->p_c_bound_values) {
        Cluster * cluster = e.first;
        double value = e.second;
        ss << cluster << "\t" << value << "\n";
    }
    filename = "./data/patients/" + this->patient.case_id + "/p_c_bound_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();
    // std::cout << "finished, ";
    // Timer::stop();

    std::cout << "exporting p_j_downregulated_given_c_bound_values\n";
    // Timer::start();
    ss.str("");
    ss << "gene_id\tcluster_address\tp_j_downregulated_given_c_bound_values\n";
    for(auto & e : this->p_j_downregulated_given_c_bound_values) {
        auto & gene_id = e.first.first;
        Cluster * cluster = e.first.second;
        double value = e.second;
        ss << gene_id << "\t" << cluster << "\t" << value << "\n";
    }
    filename = "./data/patients/" + this->patient.case_id + "/p_j_downregulated_given_c_bound_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();
    // std::cout << "finished, ";
    // Timer::stop();

    this->export_p_j_downregulated();
}

void Matchings_predictor::export_p_j_downregulated(int filename_suffix)
{
    std::stringstream ss;
    std::cout << "exporting p_j_downregulated_values\n";
    // Timer::start();
    ss << "gene_id\tp_j_downregulated_values\n";
    if(filename_suffix > 0) {
        for(auto & e: this->p_j_downregulated_values) {
            auto & gene_id = e.first;
            double value = e.second;
            ss << gene_id << "\t" << value << "\n";
        }
    } else {
        for(auto & e: this->patient.tumor_genes.profile) {            
            auto & gene_id = e.first;
            ss << gene_id << "\t" << 0 << "\n";
        }
    }
    std::stringstream filename;
    filename << "./data/patients/" << this->patient.case_id << "/predicted_downregulation/p_j_downregulated_values_" << std::setfill('0') << std::setw(6) << filename_suffix << ".tsv";
    std::ofstream out(filename.str());
    out << ss.str();
    out.close();
    // Timer::stop();
}
