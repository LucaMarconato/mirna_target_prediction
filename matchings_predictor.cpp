#include "matchings_predictor.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "omp_disabler.h"
#ifndef OMP_TEMPORARILY_DISABLED
#include <omp.h>
#else
#endif

#include <marconato/output_buffer/output_buffer.hpp>

#include "global_parameters.hpp"
#include "timer.hpp"

Matchings_predictor::Matchings_predictor(Patient& patient, std::string simulation_id) : patient(patient)
{
    double total_reads = 0;
    for (auto& e : this->patient.mirna_expression_profile.profile) {
        const Mirna_id& mirna_id = e.first;
        double relative_value = e.second.to_relative_expression();
        total_reads += e.second.to_reads();
        this->mirna_profile[mirna_id] = relative_value;
    }
    // grand_total also includes filtered reads, it will soon used to normalize reads when computing rpm.
    double grand_total = this->patient.mirna_expression_profile.profile.begin()->second.get_grand_total();
    for (auto& e : this->patient.cluster_expression_profile.profile) {
        Cluster* cluster = e.first;
        double relative_value = e.second.to_relative_expression();
        this->cluster_profile[cluster] = relative_value;
    }
    /*
      The mirnas do not sum to one because of the filtering procedure, so we
      normalize the expression profile.
    */
    double sum = 0;
    for (auto& e : this->mirna_profile) {
        sum += e.second;
    }
    this->lambda_adjustment_due_to_mirna_perturbation = 1;
    this->relative_mirna_to_reads_conversion_factor = total_reads;
    this->relative_mirna_to_rpm_conversion_factor = total_reads / grand_total * 1000000;
    // std::cout << "grand_total = " << grand_total << ", sum = " << sum << "\n";
    // std::cout << "this->relative_mirna_to_rpm_conversion_factor = " << this->relative_mirna_to_rpm_conversion_factor << "\n";
    // exit(1);
    for (auto& e : this->mirna_profile) {
        e.second /= sum;
    }

    /*
      The mirna expression profile has been just normalized to 1.
      The clusters are already normalized to 1, by construtction.
      Anyway, we check this last condition explicitly.
    */
    sum = 0;
    for (auto& e : this->cluster_profile) {
        sum += e.second;
    }
    if (std::abs(sum - 1.0) > Global_parameters::epsilon) {
        std::cerr << "error: std::abs(sum - 1.0) = " << std::abs(sum - 1.0) << "\n";
        exit(1);
    }

    this->original_cluster_profile = this->cluster_profile;

    this->simulation_id = simulation_id;
    std::string cmd = "mkdir -p ./data/patients/" + this->patient.case_id + "/matchings_predictor_output";
    std::system(cmd.c_str());
    cmd = "mkdir -p " + this->get_output_path();
    std::system(cmd.c_str());

    this->lambda = Global_parameters::lambda;
}

Matchings_predictor::Matchings_predictor(const Matchings_predictor& obj) : patient(obj.patient)
{
    this->mirna_profile = obj.mirna_profile;
    this->cluster_profile = obj.cluster_profile;
    this->original_cluster_profile = obj.original_cluster_profile;
    this->r_ic_values = obj.r_ic_values;
    this->r_ijk_values = obj.r_ijk_values;
    this->p_c_bound_values = obj.p_c_bound_values;
    this->sum_of_r_ijk_for_cluster = obj.sum_of_r_ijk_for_cluster;
    this->p_j_downregulated_given_c_bound_values = obj.p_j_downregulated_given_c_bound_values;
    this->p_j_downregulated_values = obj.p_j_downregulated_values;
    this->simulation_id = obj.simulation_id;
    this->lambda = obj.lambda;
    this->relative_mirna_to_reads_conversion_factor = obj.relative_mirna_to_reads_conversion_factor;
    this->relative_mirna_to_rpm_conversion_factor = obj.relative_mirna_to_rpm_conversion_factor;
    this->lambda_adjustment_due_to_mirna_perturbation = obj.lambda_adjustment_due_to_mirna_perturbation;
    // this->lambda_adjustment_due_to_gene_perturbation = obj.lambda_adjustment_due_to_gene_perturbation;
    this->genes_skipped_by_the_distance_based_predictor = obj.genes_skipped_by_the_distance_based_predictor;
}

void swap(Matchings_predictor& obj1, Matchings_predictor& obj2)
{
    std::swap(obj1.mirna_profile, obj2.mirna_profile);
    std::swap(obj1.cluster_profile, obj2.cluster_profile);
    std::swap(obj1.original_cluster_profile, obj2.original_cluster_profile);
    std::swap(obj1.r_ic_values, obj2.r_ic_values);
    std::swap(obj1.r_ijk_values, obj2.r_ijk_values);
    std::swap(obj1.p_c_bound_values, obj2.p_c_bound_values);
    std::swap(obj1.sum_of_r_ijk_for_cluster, obj2.sum_of_r_ijk_for_cluster);
    std::swap(obj1.p_j_downregulated_given_c_bound_values, obj2.p_j_downregulated_given_c_bound_values);
    std::swap(obj1.p_j_downregulated_values, obj2.p_j_downregulated_values);
    std::swap(obj1.simulation_id, obj2.simulation_id);
    std::swap(obj1.lambda, obj2.lambda);
    std::swap(obj1.relative_mirna_to_reads_conversion_factor, obj2.relative_mirna_to_reads_conversion_factor);
    std::swap(obj1.relative_mirna_to_rpm_conversion_factor, obj2.relative_mirna_to_rpm_conversion_factor);
    std::swap(obj1.lambda_adjustment_due_to_mirna_perturbation, obj2.lambda_adjustment_due_to_mirna_perturbation);
    // std::swap(obj1.lambda_adjustment_due_to_gene_perturbation, obj2.lambda_adjustment_due_to_gene_perturbation);
    std::swap(obj1.genes_skipped_by_the_distance_based_predictor, obj2.genes_skipped_by_the_distance_based_predictor);
}

Matchings_predictor& Matchings_predictor::operator=(Matchings_predictor obj)
{
    swap(*this, obj);
    return *this;
}

std::string Matchings_predictor::get_output_path()
{
    std::string path = "./data/patients/" + this->patient.case_id + "/matchings_predictor_output/" + simulation_id + "/";
    return path;
}

void Matchings_predictor::export_interaction_matrix()
{
    std::string filename = this->get_output_path() + "interaction_matrix.mat";
    std::ofstream out(filename);
    std::stringstream ss;
    unsigned long long cluster_count = this->cluster_profile.size();
    unsigned long long i = 0;
    for (auto& e : this->cluster_profile) {
        ss << "\"" << e.first << "\"" << ((i == cluster_count - 1) ? "\n" : "\t");
        i++;
    }
    ss << "\n";
    std::unordered_map<Mirna_id, unsigned long long> sum_of_rows;
    std::unordered_map<Cluster*, unsigned long long> sum_of_columns;
    for (auto& e : this->mirna_profile) {
        auto& mirna_id = e.first;
        sum_of_rows[mirna_id] = 0;
    }
    for (auto& e : this->cluster_profile) {
        Cluster* cluster = e.first;
        sum_of_columns[cluster] = 0;
    }

    for (auto& mirna : this->mirna_profile) {
        auto& mirna_id = mirna.first;
        ss << "\"" << mirna_id << "\"\t";
        unsigned long long i = 0;
        for (auto& e : this->cluster_profile) {
            Cluster* cluster = e.first;
            int matchings = 0;
            for (auto& site : cluster->sites) {
                auto p = std::make_pair(mirna_id, site);
                if (this->patient.interaction_graph.mirna_site_arcs.find(p) != this->patient.interaction_graph.mirna_site_arcs.end()) {
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

    filename = this->get_output_path() + "sum_of_rows.tsv";
    out.open(filename);
    ss.str("mirna_id\tsum_of_row\n");
    for (auto& e : this->mirna_profile) {
        auto& mirna_id = e.first;
        ss << mirna_id << "\t" << sum_of_rows.at(mirna_id) << "\n";
    }
    out << ss.str();
    out.close();

    filename = this->get_output_path() + "/sum_of_columns.tsv";
    out.open(filename);
    ss.str("cluster_address\tsum_of_column\n");
    for (auto& e : this->cluster_profile) {
        Cluster* cluster = e.first;
        ss << cluster << "\t" << sum_of_columns.at(cluster) << "\n";
    }
    out << ss.str();
    out.close();
}

void Matchings_predictor::compute()
{
    std::cout << "for the moment, just a trivial explicit Euler scheme\n";
    unsigned long long max_steps = 10;
    double mirna_lambda = 1;
    double cluster_lambda = 1;
    if (lambda > 1) {
        mirna_lambda = 1.0 / lambda;
    } else {
        cluster_lambda = lambda;
    }
    // h must be <= 1
    double h = 1;
    bool scaling = true;
    double cumulative_scaling = 1;
    bool logging = true;
    bool export_mirna_expression_profile = true;
    bool export_cluster_expression_profile = false;
    bool export_interaction_matrix = false;
    bool export_partial_predicted_downregulation = false;
    if (Global_parameters::test_parallelization) {
        std::cout << "generating only output for testing the parallelization\n";
        logging = true;
        export_mirna_expression_profile = true;
        export_cluster_expression_profile = false;
        export_interaction_matrix = false;
        export_partial_predicted_downregulation = false;
    }
    // just to be sure to avoid logical errors in the future
    if (!logging && (export_mirna_expression_profile || export_cluster_expression_profile || export_interaction_matrix || export_partial_predicted_downregulation)) {
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
    if (logging) {
        std::string filename0 = this->get_output_path() + "matchings_prediction.tsv";
        out0.open(filename0);
        ss0 << "cumulative_scaling\tavg_mirna_level\tavg_cluster_level\tmirna_"
               "total_exchange\tcluster_total_exchange\n";
        std::string log_folder;
        std::string cmd;
        log_folder = this->get_output_path() + "mirna_expression_profiles";
        mirna_log_filename_prefix = log_folder + "/mirna_expression_profile_";
        cmd = "mkdir -p " + log_folder;
        std::system(cmd.c_str());
        if (export_mirna_expression_profile) {
            std::cout << "deleting old log from \"" + log_folder + "\"\n";
            cmd = "rm " + log_folder + "/mirna_expression_profile*";
            std::system(cmd.c_str());
            std::cout << "deleted\n";
            ss1_header << "mirna_id\trelative_expression\treads\n";
        }

        log_folder = this->get_output_path() + "cluster_expression_profiles";
        cluster_log_filename_prefix = log_folder + "/cluster_expression_profile_";
        cmd = "mkdir -p " + log_folder;
        std::system(cmd.c_str());
        if (export_cluster_expression_profile) {
            std::cout << "deleting old log from \"" + log_folder + "\"\n";
            cmd = "rm " + log_folder + "/cluster_expression_profile*";
            std::system(cmd.c_str());
            std::cout << "deleted\n";
            // you may want to export a file in which each Cluster * (i.e. what I call
            // cluster_address below) is enriched with information to make you able to
            // identify specific clusters
            ss2_header << "cluster_address\trelative_expression\n";
        }

        if (export_interaction_matrix) {
            std::cout << "exporting the interaction matrix\n";
            this->export_interaction_matrix();
            std::cout << "done\n";
        }

        log_folder = this->get_output_path() + "predicted_downregulation";
        cmd = "mkdir -p " + log_folder;
        std::system(cmd.c_str());
        if (export_partial_predicted_downregulation) {
            std::cout << "deleting old log from \"" + log_folder + "\"\n";
            cmd = "rm " + log_folder + "/p_j_downregulated_values_*";
            std::system(cmd.c_str());
            std::cout << "deleted\n";
        }
    }

    int latest_percentage = -1;
    // these four variables will be shared among the threads and set either only
    // by the thread 0, either inside a #pragma omp critical
    double mirna_total_exchange;
    double cluster_total_exchange;
    std::unordered_map<Mirna_id, double> new_mirna_profile;
    std::unordered_map<Cluster*, double> new_cluster_profile;
    unsigned long long loop_size;
    unsigned int t;

    std::cout << "starting the simulation\n";
#ifndef OMP_TEMPORARILY_DISABLED
    unsigned long long my_num_threads = 4;
#pragma omp parallel num_threads(my_num_threads)
#endif
    {
#ifndef OMP_TEMPORARILY_DISABLED
        unsigned long long rank = omp_get_thread_num();
#else
        unsigned long long my_num_threads = 1;
        unsigned long long rank = 0;
#endif
        if (rank == 0) {
            std::cout << "my_num_threads = " << my_num_threads << "\n";
            t = 0;
        }
        for (t = 0; t < max_steps;) {
            if (rank == 0) {
                int current_percentage = (int)(100 * ((double)(t) / max_steps));
                if (current_percentage > latest_percentage) {
                    latest_percentage = current_percentage;
                    std::cout << "t/max_steps: " << t << "/" << max_steps << " = " << current_percentage << "%\n";
                    if (logging) {
                        if (export_mirna_expression_profile) {
                            std::stringstream ss1;
                            ss1 << ss1_header.str();
                            for (auto& e : this->mirna_profile) {
                                double adjusted_value = e.second / cumulative_scaling;
                                double reads = adjusted_value * this->relative_mirna_to_reads_conversion_factor;
                                ss1 << e.first << "\t" << adjusted_value << "\t" << reads << "\n";
                            }
                            std::stringstream mirna_log_filename;
                            mirna_log_filename << mirna_log_filename_prefix << std::setfill('0') << std::setw(6) << t << ".tsv";
                            std::ofstream out1(mirna_log_filename.str());
                            out1 << ss1.str();
                            out1.close();
                        }
                        if (export_cluster_expression_profile) {
                            std::stringstream ss2;
                            ss2 << ss2_header.str();
                            for (auto& e : this->cluster_profile) {
                                double adjusted_value = e.second / cumulative_scaling;
                                // TODO: export reads along with relative quantities, so to make different datasets comparable
                                ss2 << e.first << "\t" << adjusted_value << "\n";
                            }
                            std::stringstream cluster_log_filename;
                            cluster_log_filename << cluster_log_filename_prefix << std::setfill('0') << std::setw(6) << t << ".tsv";
                            std::ofstream out2(cluster_log_filename.str());
                            out2 << ss2.str();
                            out2.close();
                        }
                        if (export_partial_predicted_downregulation) {
                            // if t == 0, then export_p_j_downregulated(t) will just export
                            // the initial expression profile
                            if (t > 0) {
                                this->compute_probabilities();
                            }
                            this->export_p_j_downregulated(t);
                        }
                    }
                }
                if (scaling) {
                    // perform the scaling
                    double mirna_sum = 0;
                    for (auto& e : this->mirna_profile) {
                        mirna_sum += e.second;
                    }
                    double cluster_sum = 0;
                    for (auto& e : this->cluster_profile) {
                        cluster_sum += e.second;
                    }
                    double sum = std::max(mirna_sum, cluster_sum);
                    for (auto& e : this->mirna_profile) {
                        e.second /= sum;
                    }
                    for (auto& e : this->cluster_profile) {
                        e.second /= sum;
                    }
                    if (sum > 1 + Global_parameters::epsilon) {
                        std::cerr << "error: sum = " << sum << "\n";
                        exit(1);
                    }
                    cumulative_scaling /= sum;
                }
                if (logging) {
                    ss0 << cumulative_scaling << "\t";
                    // compute the average level
                    double avg_mirna_level = 0;
                    for (auto& e : this->mirna_profile) {
                        avg_mirna_level += e.second;
                    }
                    avg_mirna_level /= this->mirna_profile.size();
                    double avg_cluster_level = 0;
                    for (auto& e : this->cluster_profile) {
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
            // TODO: check if it is faster to copy the profile and have each rank
            // accessing its copy or to access the shared copy among the ranks auto
            // rank_mirna_profile = this->mirna_profile; auto rank_cluster_profile =
            // this->cluster_profile;
            auto rank_new_mirna_profile = std::unordered_map<Mirna_id, double>();
            auto rank_new_cluster_profile = std::unordered_map<Cluster*, double>();
            for (auto& e : this->mirna_profile) {
                rank_new_mirna_profile[e.first] = 0;
            }
            for (auto& e : this->cluster_profile) {
                rank_new_cluster_profile[e.first] = 0;
            }
            double rank_mirna_total_exchange = 0;
            double rank_cluster_total_exchange = 0;
            auto rank_r_ic_values = std::unordered_map<std::pair<Mirna_id, Cluster*>, double>();
            auto rank_r_ijk_values = std::unordered_map<std::pair<Mirna_id, Site*>, double>();

            auto it = this->patient.interaction_graph.gene_to_clusters_arcs.begin();

            unsigned long long int i;
            for (i = 0; i < loop_size && i < rank; i++) {
                it++;
            }
            for (; i < loop_size;) {
                auto& clusters = it->second;
                for (Cluster* cluster : clusters) {
                    for (Site* site : cluster->sites) {
                        for (const Mirna_id& mirna_id : this->patient.interaction_graph.site_to_mirnas_arcs.at(site)) {
                            // the computation of the following value is inefficient but I bet
                            // it does not impact the performance since there is no cluster
                            // with an high number of sites
                            int number_of_sites_for_the_mirna_in_the_cluster = 0;
                            for (Site* site : cluster->sites) {
                                if (site->mirna_id == mirna_id) {
                                    number_of_sites_for_the_mirna_in_the_cluster++;
                                }
                            }
                            double mirna_value = this->mirna_profile.at(mirna_id);
                            double cluster_value = this->cluster_profile.at(cluster);
                            double exchange = h * mirna_value * cluster_value / number_of_sites_for_the_mirna_in_the_cluster;

                            if (exchange < 0) {
                                static bool warning_already_presented = false;
                                if (abs(exchange) < Global_parameters::epsilon) {
                                    if (!warning_already_presented) {
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
                            if (cumulative_scaling < 1 - Global_parameters::epsilon) {
                                std::cerr << "error: cumulative_scaling = " << cumulative_scaling << "\n";
                                exit(1);
                            }
                            rank_mirna_total_exchange += (mirna_lambda * exchange / cumulative_scaling);
                            rank_cluster_total_exchange += (cluster_lambda * exchange / cumulative_scaling);

                            auto p = std::make_pair(mirna_id, cluster);
                            if (rank_r_ic_values.find(p) == rank_r_ic_values.end()) {
                                rank_r_ic_values[p] = 0;
                            }

                            rank_r_ic_values[p] += cluster_lambda * exchange / cumulative_scaling;

                            auto ijk = std::make_pair(mirna_id, site);
                            if (rank_r_ijk_values.find(ijk) == rank_r_ijk_values.end()) {
                                rank_r_ijk_values[ijk] = 0;
                            }
                            rank_r_ijk_values[ijk] += cluster_lambda * exchange / cumulative_scaling;
                        }
                    }
                }

                for (unsigned long long j = 0; j < my_num_threads && i < loop_size; j++) {
                    it++;
                    i++;
                }
            }
#ifndef OMP_TEMPORARILY_DISABLED
#pragma omp critical
#endif
            {
                for (auto& e : rank_new_mirna_profile) {
                    new_mirna_profile.at(e.first) += e.second;
                }
                for (auto& e : rank_new_cluster_profile) {
                    new_cluster_profile.at(e.first) += e.second;
                }
                mirna_total_exchange += rank_mirna_total_exchange;
                cluster_total_exchange += rank_cluster_total_exchange;
                for (auto& e : rank_r_ic_values) {
                    if (this->r_ic_values.find(e.first) == this->r_ic_values.end()) {
                        this->r_ic_values[e.first] = e.second;
                    } else {
                        this->r_ic_values[e.first] += e.second;
                    }
                }
                for (auto& e : rank_r_ijk_values) {
                    if (this->r_ijk_values.find(e.first) == this->r_ijk_values.end()) {
                        this->r_ijk_values[e.first] = e.second;
                    } else {
                        this->r_ijk_values[e.first] += e.second;
                    }
                }
            }
            if (rank == 0) {
                t++;
            }
#ifndef OMP_TEMPORARILY_DISABLED
#pragma omp barrier
#endif
            if (rank == 0) {
                this->mirna_profile = new_mirna_profile;
                this->cluster_profile = new_cluster_profile;
                if (logging) {
                    ss0 << mirna_total_exchange << "\t" << cluster_total_exchange << "\n";
                    out0 << ss0.str();
                    ss0.str("");
                }
            }
        } // end of step of the algorithm
    }     // end of the parallelization block
    if (logging) {
        out0.close();
    }
    this->compute_probabilities();
    if (!Global_parameters::test_parallelization) {
        this->export_probabilities();
    } else {
        this->export_p_j_downregulated();
    }
}

void Matchings_predictor::compute_probabilities()
{
    // the function compute_probabilities is called many times if
    // export_partial_predicted_downregulation is true, so we need to reset the
    // probabilities before proceeding r_ijk_values.clear();
    p_c_bound_values.clear();
    sum_of_r_ijk_for_cluster.clear();
    p_j_downregulated_given_c_bound_values.clear();
    p_j_downregulated_values.clear();

    for (auto& e : this->r_ic_values) {
        auto& mirna_id = e.first.first;
        Cluster* cluster = e.first.second;
        double value = e.second;
        // sum_of_r_ijk_for_cluster is needed when computing the denominator of
        // this->p_j_downregulated_given_c_bound_values
        if (this->sum_of_r_ijk_for_cluster.find(cluster) == this->sum_of_r_ijk_for_cluster.end()) {
            this->sum_of_r_ijk_for_cluster[cluster] = 0;
        }
        this->sum_of_r_ijk_for_cluster[cluster] = this->sum_of_r_ijk_for_cluster.at(cluster) + value;

        // update p_c_bound_values, note that the values of a cluster are not
        // definitive since we need to excecute the following code once per miRNA
        // and the miRNA is given by the outer loop
        if (this->p_c_bound_values.find(cluster) == this->p_c_bound_values.end()) {
            this->p_c_bound_values[cluster] = 0;
        }
        this->p_c_bound_values[cluster] = this->p_c_bound_values.at(cluster) + value / (this->original_cluster_profile.at(cluster));
    }

    // compute p_j_downregulated_given_c_bound_values for the current cluster
    for (auto& e : this->r_ic_values) {
        auto& mirna_id = e.first.first;
        Cluster* cluster = e.first.second;
        Gene_id gene_id = cluster->sites.front()->gene_id;
        std::unordered_map<Site*, std::list<std::pair<Mirna_id, Site*>>> just_added_by_site;
        for (auto& site : cluster->sites) {
            auto& m = this->patient.interaction_graph.mirna_site_arcs;
            auto p = std::make_pair(mirna_id, site);
            if (m.find(p) != m.end()) {
                just_added_by_site[site].push_back(p);
            }
        }
        if (just_added_by_site.size() == 0) {
            std::cerr << "error: just_added_by_site.size() = " << just_added_by_site.size() << "\n";
            exit(1);
        }
        for (auto& e : just_added_by_site) {
            // Site * site = e.first;
            for (auto& p : e.second) {
                auto key0 = std::make_pair(gene_id, cluster);
                if (this->p_j_downregulated_given_c_bound_values.find(key0) == this->p_j_downregulated_given_c_bound_values.end()) {
                    this->p_j_downregulated_given_c_bound_values[key0] = 0;
                }
                double r_ijk = this->r_ijk_values.at(p);
                Mirna_site_arc& mirna_site_arc = this->patient.interaction_graph.mirna_site_arcs.at(p);
                double context_score = mirna_site_arc.context_score;
                double probability_of_downregulation = 1 - std::pow(2, context_score);
                if (std::abs(sum_of_r_ijk_for_cluster.at(cluster)) < Global_parameters::epsilon) {
                    // if(std::abs(p_c_bound_values.at(cluster)) <
                    // Global_parameters::epsilon) {
                    static bool silence_warning = false;
                    // if(p.second->utr_start == 23 && gene_id == 7494) {
                    //     std::cerr << "error: we should never get here!\n";
                    //     exit(1);
                    // }
                    if (std::abs(sum_of_r_ijk_for_cluster.at(cluster)) < std::pow(2, -50)) {
                        if (!silence_warning) {
                            silence_warning = true;
                            std::cerr << "error: division by zero. If this arises in the "
                                         "context of a perturbation analysis in which some "
                                         "mirnas are knocked out then this warning is fine\n "
                                         "This warning will be silenced.";
                        }
                    } else {
                        std::cerr << "error: division by zero. Do you really need to multiply by "
                                     "r_ijk and divide by sum_of_r_ijk_for_clusters? Probably "
                                     "not! Review the code by printing those values and act "
                                     "accordingly. If the cluster has only one site in and the "
                                     "site bings to only one mirna, then those values will be "
                                     "identical! To reproduce this use some gaussian "
                                     "perturbations and you will eventually get this error\n";
                        // TODO: fix this. When considering gaussian perturbations we end
                        // here
                        // exit(1);
                    }
                } else {
                    double to_add = r_ijk * probability_of_downregulation / sum_of_r_ijk_for_cluster.at(cluster);
                    // double to_add = r_ijk * probability_of_downregulation /
                    // p_c_bound_values.at(cluster);
                    this->p_j_downregulated_given_c_bound_values[key0] = this->p_j_downregulated_given_c_bound_values.at(key0) + to_add;
                }
            }
        }
    }

    // compute p_j_downregulated_values
    // #define TEST_THE_ITERATIVE_ALGORITHM
#ifdef TEST_THE_ITERATIVE_ALGORITHM
    // begin debug code
    double cluster_limit = 20;
    std::cout << "computing the probabilities of down-regulation for genes with "
                 "at most "
              << cluster_limit << " clusters\n";
    Timer::start();
    int genes_debugged = 0;
    int genes_not_debugged = 0;
    // end debug code
#endif
    unsigned long max_number_of_clusters_per_gene = 0;
    for (auto& e : this->patient.interaction_graph.gene_to_clusters_arcs) {
        auto& clusters = e.second;
        max_number_of_clusters_per_gene = std::max(max_number_of_clusters_per_gene, clusters.size());
    }
#ifdef TEST_THE_ITERATIVE_ALGORITHM
    // begin debug code
    bool* b = new bool[max_number_of_clusters_per_gene];
    int latest_percentage = -1;
    // end debug code
#endif
    double* p_j_downregulated_given_c_bound_values_flattened = new double[max_number_of_clusters_per_gene];
    double* p_c_bound_values_flattened = new double[max_number_of_clusters_per_gene];

    for (auto& e : this->patient.interaction_graph.gene_to_clusters_arcs) {
#ifdef TEST_THE_ITERATIVE_ALGORITHM
        // begin debug code
        int total_genes = this->patient.interaction_graph.gene_to_clusters_arcs.size();
        int current_percentage = (int)(100 * ((double)(genes_debugged + genes_not_debugged) / total_genes));
        if (current_percentage > latest_percentage) {
            latest_percentage = current_percentage;
            std::cout << "genes_debugged/total_genes: " << genes_debugged + genes_not_debugged << "/" << total_genes << " = " << current_percentage << "%\n";
        }
        // end debug code
#endif
        Gene_id gene_id = e.first;
        auto& clusters = e.second;
        int i = 0;
        for (Cluster* cluster : clusters) {
#ifdef TEST_THE_ITERATIVE_ALGORITHM
            // begin debug code
            b[i] = false;
            // end debug code
#endif
            p_j_downregulated_given_c_bound_values_flattened[i] = this->p_j_downregulated_given_c_bound_values.at(std::make_pair(gene_id, cluster));
            p_c_bound_values_flattened[i] = this->p_c_bound_values.at(cluster);
            // std::cout << "utr_start values for the sites in the current
            // cluster:\n"; for(Site * site : cluster->sites) {
            //     std::cout << site->utr_start << " ";
            // }
            // std::cout << "\n";
            // std::cout << "p_j_downregulated_given_c_bound_values_flattened[" << i
            // << "], p_c_bound_values_flattened[" << i << "]\n"; std::cout <<
            // p_j_downregulated_given_c_bound_values_flattened[i] << " " <<
            // p_c_bound_values_flattened[i] << "\n\n";
            i++;
        }
        double sum = this->iteratively_compute_p_j_downregulated(p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened, clusters.size());
        this->p_j_downregulated_values[gene_id] = sum;

#ifdef TEST_THE_ITERATIVE_ALGORITHM
        // begin debug code
        // safety check for genes with at most "cluster_limit" clusters
        if (clusters.size() <= cluster_limit) {
            double debug_sum = 0;
            this->recusively_compute_p_j_downregulated(b, 0, clusters.size(), &debug_sum, 1, 1, p_j_downregulated_given_c_bound_values_flattened, p_c_bound_values_flattened);
            if (std::abs(sum - debug_sum) > Global_parameters::epsilon) {
                std::cerr << "error: sum = " << sum << ", debug_sum = " << debug_sum << ", std::abs(sum - debug_sum) = " << std::abs(sum - debug_sum) << "\n";
                std::cout << "describing the instance: clusters.size() = " << clusters.size() << "\n";
                std::cout << "p_j_downregulated_given_c_bound_values_flattened, "
                             "p_c_bound_values_flattened\n";
                for (int j = 0; j < clusters.size(); j++) {
                    std::cout << p_j_downregulated_given_c_bound_values_flattened[j] << ", " << p_c_bound_values_flattened[j] << "\n";
                }
                exit(1);
            }
            genes_debugged++;
        } else {
            genes_not_debugged++;
        }
#endif
    }

#ifdef TEST_THE_ITERATIVE_ALGORITHM
    std::cout << "genes_debugged/total_genes: " << genes_debugged << "/" << (genes_debugged + genes_not_debugged) << " = " << ((double)genes_debugged) / (genes_debugged + genes_not_debugged) << "\n";
    std::cout << "finished, ";
    Timer::stop();
    delete[] b;
    // end debug code
#endif
    delete[] p_j_downregulated_given_c_bound_values_flattened;
    delete[] p_c_bound_values_flattened;

    /*
      If considering distance then most of the values for p_j_downregualted will
      be overwritten, but not everyone, since computing sometimes the computation
      would require too much time (exponential algorithm in the worst case). For
      those values, the user will be notified that the value has not beed computed
      considering a distance-based prediction
    */
    this->genes_skipped_by_the_distance_based_predictor.clear();
    if (Global_parameters::consider_distance_for_predictions) {
        this->compute_distance_based_predictions();
    }
    std::cout << "\n";
}

/*
  The algorithm described by this function is exponential in the number of
  clusters in a gene. The array b and the valus sum, p_j_downregulated_given_c_b
  and p_b are computed incrementally at every recursive call and hold the final
  value their name represents only when processed in the last level.
*/
void Matchings_predictor::recusively_compute_p_j_downregulated(bool* b, int level, int max_level, double* sum, double p_j_downregulated_given_b, double p_b, double* p_j_downregulated_given_c_bound_values_flattened, double* p_c_bound_values_flattened)
{
    if (level == max_level) {
        p_j_downregulated_given_b = 1 - p_j_downregulated_given_b;
        *sum += p_j_downregulated_given_b * p_b;
        // std::cout << "b: ";
        // for(int i = 0; i < max_level; i++) {
        //     std::cout << b[i] << " ";
        // }
        // std::cout << "| " << p_j_downregulated_given_b * p_b << "\n";
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
  The algoirthm described by this functions is polynomial in the number of
  clusters in a gene.
*/
double Matchings_predictor::iteratively_compute_p_j_downregulated(double* p_j_downregulated_given_c_bound_values_flattened, double* p_c_bound_values_flattened, int clusters_count)
{
    double ratio_j_not_downregulated = 1.0;
    for (int i = 0; i < clusters_count; i++) {
        double newly_downregulated = ratio_j_not_downregulated * p_c_bound_values_flattened[i] * p_j_downregulated_given_c_bound_values_flattened[i];
        ratio_j_not_downregulated -= newly_downregulated;
    }
    return 1 - ratio_j_not_downregulated;
}

void Matchings_predictor::export_probabilities()
{
    std::ofstream out;
    std::stringstream ss;
    std::string filename;

    std::cout << "exporting r_ic_values\n";
    ss.str("");
    ss << "mirna_id\tcluster_address\tr_ic\n";
    for (auto& e : this->r_ic_values) {
        auto& mirna_id = e.first.first;
        Cluster* cluster = e.first.second;
        double value = e.second;
        ss << mirna_id << "\t" << cluster << "\t" << value << "\n";
    }
    filename = this->get_output_path() + "r_ic_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();

    std::cout << "exporting r_ijk_values\n";
    ss.str("");
    ss << "mirna_id\tsite_address\tr_ijk\n";
    for (auto& e : this->r_ijk_values) {
        auto& mirna_id = e.first.first;
        Site* site = e.first.second;
        double value = e.second;
        ss << mirna_id << "\t" << site << "\t" << value << "\n";
    }
    filename = this->get_output_path() + "r_ijk_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();

    std::cout << "exporting p_c_bound_values\n";
    ss.str("");
    ss << "cluster_address\tp_c_bound\n";
    for (auto& e : this->p_c_bound_values) {
        Cluster* cluster = e.first;
        double value = e.second;
        ss << cluster << "\t" << value << "\n";
    }
    filename = this->get_output_path() + "p_c_bound_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();

    std::cout << "exporting p_j_downregulated_given_c_bound_values\n";
    ss.str("");
    ss << "gene_id\tcluster_address\tp_j_downregulated_given_c_bound_values\n";
    for (auto& e : this->p_j_downregulated_given_c_bound_values) {
        auto& gene_id = e.first.first;
        Cluster* cluster = e.first.second;
        double value = e.second;
        ss << gene_id << "\t" << cluster << "\t" << value << "\n";
    }
    filename = this->get_output_path() + "p_j_downregulated_given_c_bound_values.tsv";
    out.open(filename);
    out << ss.str();
    out.close();

    this->export_p_j_downregulated();
}

void Matchings_predictor::export_p_j_downregulated(int filename_suffix)
{
    std::stringstream ss;
    std::cout << "exporting p_j_downregulated_values\n";

    ss << "gene_id\tp_j_downregulated_values\n";
    if (filename_suffix > 0) {
        for (auto& e : this->p_j_downregulated_values) {
            auto& gene_id = e.first;
            double value = e.second;
            ss << gene_id << "\t" << value << "\n";
        }
    } else {
        for (auto& e : this->patient.gene_expression_profile.profile) {
            auto& gene_id = e.first;
            ss << gene_id << "\t" << 0 << "\n";
        }
    }

    std::stringstream filename;
    filename << this->get_output_path() << "predicted_downregulation/p_j_downregulated_values_" << std::setfill('0') << std::setw(6) << filename_suffix << ".tsv";
    std::ofstream out(filename.str());
    out << ss.str();
    out.close();
}
