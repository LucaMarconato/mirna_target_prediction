#include "matchings_predictor.hpp"

#include <sstream>
#include <fstream>
#include <iomanip>

#include "omp_disabler.h"
#ifndef OMP_TEMPORARILY_DISABLED
#include <omp.h>
#else
#endif

#include "global_parameters.hpp"

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
    
}

void Matchings_predictor::compute()
{
    std::cout << "for the moment, just a trivial explicit Euler scheme\n";
    // not used for the moment
    // double lambda = 1;
    unsigned long long max_steps = 100;
    // h must be <= 1
    double h = 1;
    bool scaling = true;
    double mirna_cumulative_scaling = 1;
    double cluster_cumulative_scaling = 1;
    bool logging = true;
    bool export_mirna_expression_profile = true;
    /*
      Use this only with very small interaction graphs.
      TODO: truncate logging if there are too many genes
    */
    // just to be sure to avoid logical errors in the future
    bool export_cluster_expression_profile = true;
    if(!logging && (export_mirna_expression_profile || export_cluster_expression_profile)) {
        export_mirna_expression_profile = false;
        export_cluster_expression_profile = false;
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
        ss0 << "mirna_cumulative_scaling\tcluster_cumulative_scaling\tavg_mirna_level\tavg_cluster_level\ttotal_exchange\n";
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
    }    
    
    int latest_percentage = -1;
    // these four variables will be shared among the threads and set either only by the thread 0, either inside a #pragma omp critical
    double total_exchange;
    std::unordered_map<Mirna_id, double> new_mirna_profile;
    std::unordered_map<Cluster *, double> new_cluster_profile;
    unsigned long long loop_size;
    int t;
    
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
                                double adjusted_value = e.second / mirna_cumulative_scaling;
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
                                double adjusted_value = e.second / cluster_cumulative_scaling;
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
                    auto scale = [](auto & profile, double & cumulative_scaling)
                                 {
                                     double sum = 0;
                                     for(auto & e : profile) {
                                         sum += e.second;
                                     }
                                     for(auto & e : profile) {
                                         e.second /= sum;
                                     }
                                     if(sum > 1 + Global_parameters::epsilon) {
                                         std::cerr << "error: sum = " << sum << "\n";
                                         exit(1);
                                     }
                                     cumulative_scaling /= sum;                
                                 };
                    
                    scale(this->mirna_profile, mirna_cumulative_scaling);
                    scale(this->cluster_profile, cluster_cumulative_scaling);
                }
                if(logging) {
                    ss0 << mirna_cumulative_scaling << "\t" << cluster_cumulative_scaling << "\t";
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
                    double adjusted_avg_mirna_level = avg_mirna_level / mirna_cumulative_scaling;
                    double adjusted_avg_cluster_level = avg_cluster_level / cluster_cumulative_scaling;
                    ss0 << adjusted_avg_mirna_level << "\t" << adjusted_avg_cluster_level << "\t";
                }
                // perform the exchange
                total_exchange = 0;
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
            double rank_total_exchange;
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
                                        std::cerr << "warning: exchange = " << exchange << ", silening further warnings like this\n";   
                                    }
                                } else {
                                    std::cout << "error: exchange = " << exchange << "\n";
                                    exit(1);   
                                }
                            }
                            double mirna_lambda = 1;
                            double cluster_lambda = 1;
                            if(Global_parameters::lambda < 1) {
                                mirna_lambda = Global_parameters::lambda;
                            } else {
                                cluster_lambda = 1.0/Global_parameters::lambda;
                            }
                            rank_new_mirna_profile.at(mirna_id) -= mirna_lambda * exchange;
                            rank_new_cluster_profile.at(cluster) -= cluster_lambda * exchange;
                            rank_total_exchange += exchange;
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
                total_exchange += rank_total_exchange;
            }
#ifndef OMP_TEMPORARILY_DISABLED
            if(rank == 0) {
                t++;
            }
#pragma omp barrier
#endif
            if(rank == 0) {
                this->mirna_profile = new_mirna_profile;
                this->cluster_profile = new_cluster_profile;
                if(abs(mirna_cumulative_scaling - cluster_cumulative_scaling) > Global_parameters::epsilon) {
                    std::cerr << "error: mirna_cumulative_scaling = " << mirna_cumulative_scaling << ", cluster_cumulative_scaling = " << cluster_cumulative_scaling << "\n";
                    exit(1);
                }
                if(logging) {
                    ss0 << total_exchange << "\n";
                    out0 << ss0.str();
                    ss0.str("");
                }
            }
        } // end of step of the algorithm
    } // end of the parallelization block
    if(logging) {
        out0.close();   
    }
}



// TODO: export the system of ODE
