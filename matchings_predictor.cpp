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
    for(auto & mirna : this->mirna_profile) {
        auto & mirna_id = mirna.first;
        ss << "\"" << mirna_id << "\"\t";
        unsigned long long i = 0;
        for(auto & cluster : this->cluster_profile) {
            int matching = 0;
            for(auto & site : cluster.first->sites) {
                auto p = std::make_pair(mirna_id, site);
                if(this->patient.interaction_graph.mirna_site_arcs.find(p) != this->patient.interaction_graph.mirna_site_arcs.end()) {
                    matching++;
                }
            }
            ss << matching << ((i == cluster_count - 1) ? "\n" : "\t");
            i++;
        }
    }
    out << ss.str();
    out.close();
}

void Matchings_predictor::compute()
{
    std::cout << "for the moment, just a trivial explicit Euler scheme\n";
    unsigned long long max_steps = 100;
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
    bool logging = false;
    bool export_mirna_expression_profile = true;
    bool export_cluster_expression_profile = true;
    bool export_interaction_matrix = true;
    // just to be sure to avoid logical errors in the future
    if(!logging && (export_mirna_expression_profile || export_cluster_expression_profile || export_interaction_matrix)) {
        export_mirna_expression_profile = false;
        export_cluster_expression_profile = false;
        export_interaction_matrix = false;
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
    std::ofstream out("cluster_debugging");
    std::stringstream ss;
    
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
        ss << cluster << " " << cluster->sites.size() << " " << cluster->sites.front()->mirna_id << " " << cluster->sites.front()->gene_id << " " << cluster->sites.front()->utr_start << " " << this->p_c_bound_values.at(cluster) << " " << value << " " << this->original_cluster_profile.at(cluster) << " " << value / (this->original_cluster_profile.at(cluster) * cluster->sites.size()) << "\n";
        this->p_c_bound_values[cluster] = this->p_c_bound_values.at(cluster) + value / (this->original_cluster_profile.at(cluster) * cluster->sites.size());
    }
    out << ss.str();
    out.close();
    
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
                double probability = std::pow(2, 1 - context_score);
                double to_add = r_ijk * probability / sum_of_r_ijk_for_cluster.at(cluster);
                this->p_j_downregulated_given_c_bound_values[key0] = this->p_j_downregulated_given_c_bound_values.at(key0) + to_add;
            }
        }        
    }
}

void Matchings_predictor::export_probabilities()
{
    std::ofstream out;
    std::stringstream ss;
    std::string filename;
    
    std::cout << "exporting r_ic_values\n";
    Timer::start();
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
    std::cout << "finished, \n";
    Timer::stop();

    std::cout << "exporting r_ijk_values\n";
    Timer::start();
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
    std::cout << "finished, \n";
    Timer::stop();

    std::cout << "exporting p_c_bound_values\n";
    Timer::start();
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
    std::cout << "finished, \n";
    Timer::stop();

    std::cout << "exporting p_j_downregulated_given_c_bound_values\n";
    Timer::start();
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
    std::cout << "finished, \n";
    Timer::stop();
}
