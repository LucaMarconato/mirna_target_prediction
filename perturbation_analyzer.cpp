#include "perturbation_analyzer.hpp"

#include <random>
#include <vector>

#include "timer.hpp"

// // --------------------Perturbation_target--------------------

std::string Perturbation_target::string_id()
{
    std::stringstream ss;
    if(this->specific_mirna.is_valid) {
        ss << "id" << this->specific_mirna.mirna_id;
    }
    if(this->specific_gene.is_valid) {
        ss << "id" << this->specific_gene.gene_id;
    }
    if(this->nth_largest_element.is_valid) {
        ss << this->nth_largest_element.n;
    }
    if(this->elements_from_nth_largest.is_valid) {
        ss << this->elements_from_nth_largest.n;
    }
    return ss.str();
}

// --------------------Perturbation_analyzer--------------------

Perturbation_analyzer::Perturbation_analyzer(Patient & patient) : patient(patient) {}

void Perturbation_analyzer::run(Perturbation_type mirna_perturbation_type,
                                Perturbation_type gene_perturbation_type,
                                Perturbation_target mirna_perturbation_target,
                                Perturbation_target gene_perturbation_target,
                                double mirna_perturbation_amplifier,
                                double gene_perturbation_amplifier,
                                std::string output_name_suffix)
{
    this->mirna_perturbation_type = mirna_perturbation_type;
    this->gene_perturbation_type = gene_perturbation_type;
    this->mirna_perturbation_target = mirna_perturbation_target;
    this->gene_perturbation_target = gene_perturbation_target;
    this->mirna_perturbation_amplifier = mirna_perturbation_amplifier;
    this->gene_perturbation_amplifier = gene_perturbation_amplifier;
    this->integrity_check();
    simulation_id.str("");
    simulation_id << (mirna_perturbation_type == Perturbation_type::Point_perturbation ? "p" : (mirna_perturbation_type == Perturbation_type::Gaussian_perturbation) ? "g" : "") << "_"
                  << (gene_perturbation_type == Perturbation_type::Point_perturbation ? "p" : (gene_perturbation_type == Perturbation_type::Gaussian_perturbation) ? "g" : "") << "_"
                  << mirna_perturbation_target.string_id() << "_"
                  << gene_perturbation_target.string_id() << "_"
                  << mirna_perturbation_amplifier << "_"
                  << gene_perturbation_amplifier << "_"
                  << output_name_suffix;
    this->perturb();
    this->matchings_predictor->compute();    
}

void Perturbation_analyzer::integrity_check()
{
    if(mirna_perturbation_type == Perturbation_type::Point_perturbation && gene_perturbation_type == Perturbation_type::Point_perturbation) {
        std::cerr << "error: can only use apply a point perturbation to only either one mirna, either one gene\n";
        exit(1);
    }
    
    bool error = false;
    if(mirna_perturbation_type == Perturbation_type::Point_perturbation &&
       !mirna_perturbation_target.specific_mirna.is_valid &&
       !mirna_perturbation_target.nth_largest_element.is_valid) {
        error = true;
    }
    if(gene_perturbation_type == Perturbation_type::Point_perturbation &&
       !gene_perturbation_target.specific_gene.is_valid &&
       !gene_perturbation_target.nth_largest_element.is_valid) {
        error = true;
    }
    if(mirna_perturbation_type == Perturbation_type::No_perturbation &&
       !mirna_perturbation_target.empty_target.is_valid) {
        error = true;
    }
    if(gene_perturbation_type == Perturbation_type::No_perturbation &&
       !gene_perturbation_target.empty_target.is_valid) {
        error = true;
    }
    if(mirna_perturbation_type == Perturbation_type::Gaussian_perturbation &&
       !mirna_perturbation_target.elements_from_nth_largest.is_valid) {
        error = true;
    }
    if(gene_perturbation_type == Perturbation_type::Gaussian_perturbation &&
       !gene_perturbation_target.elements_from_nth_largest.is_valid) {
        error = true;
    }    

    if(error) {
        std::cerr << "error: integrity check failed\n";
        exit(1);
    }
}

void Perturbation_analyzer::perturb()
{
    this->perturbed_patient = this->patient;
    std::map<double, std::list<Mirna_id>, std::greater<double>> mirna_profile_sorted;
    std::map<double, std::list<Gene_id>, std::greater<double>> gene_profile_sorted;
    for(auto & e : this->perturbed_patient.tumor_mirnas.profile) {
        const Mirna_id & mirna_id = e.first;
        double relative_value = e.second.to_relative_expression();
        mirna_profile_sorted[relative_value].push_back(mirna_id);
    }
    for(auto & e : this->perturbed_patient.tumor_genes.profile) {        
        const Gene_id & gene_id = e.first;
        double relative_value = e.second.to_relative_expression();
        gene_profile_sorted[relative_value].push_back(gene_id);
    }
    
    // perturb mirnas
    /*
      We compute total_reads_before_perturbation and to compare it with total_reads_after_perturbation because the grand_total value, which can be accessed for instance using this->patient.tumor_mirnas.profile.begin()->second.get_grand_total(), contains also information about the filtered reads.
      One could also look at this->patient.tumor_mirnas.filtered_out_reads.
      I prefer this approach because if one applies the filtering procedure more than once it is possible to include some subtle bugs
    */
    double total_reads_before_perturbation = 0;
    for(auto & e : this->perturbed_patient.tumor_mirnas.profile) {
        double reads = e.second.to_reads();
        total_reads_before_perturbation += reads;
    }
    
    if(this->mirna_perturbation_type == Perturbation_type::Point_perturbation) {        
        if(this->mirna_perturbation_target.specific_mirna.is_valid ||
           this->mirna_perturbation_target.nth_largest_element.is_valid) {
            Mirna_id mirna_id = -1;
            if(this->mirna_perturbation_target.specific_mirna.is_valid) {
                mirna_id = this->mirna_perturbation_target.specific_mirna.mirna_id;
            } else if(this->mirna_perturbation_target.nth_largest_element.is_valid) {
                unsigned int n = this->mirna_perturbation_target.nth_largest_element.n;
                if(n > this->perturbed_patient.tumor_mirnas.profile.size()) {
                    std::cerr << "error: n = " << n << ", this->perturbed_patient.tumor_mirnas.profile.size() = " << this->perturbed_patient.tumor_mirnas.profile.size() << "\n";
                    exit(1);
                }
                unsigned int i = 0;
                for(auto & e : mirna_profile_sorted) {
                    auto & list_of_ids = e.second;
                    if(i + list_of_ids.size() > n) {
                        mirna_id = list_of_ids.front();
                        break;
                    }
                    i += list_of_ids.size();
                }
            }
            std::cout << "mirna_id = " << mirna_id << "\n";
            double reads = this->perturbed_patient.tumor_mirnas.profile.at(mirna_id).to_reads();
            // when studying the gaussian perturbation you may to use (1 + sigma + exp(sigma/expression_level)) instead of the formula used
            reads *= (1 + this->mirna_perturbation_amplifier);
            this->perturbed_patient.tumor_mirnas.profile[mirna_id] = Reads(reads);
        } else {
            std::cerr << "error: exception in this->mirna_perturbation_type == Perturbation_type::Point_perturbation\n";
            exit(1);
        }
    } else if(this->mirna_perturbation_type == Perturbation_type::Gaussian_perturbation) {
        if(this->mirna_perturbation_target.elements_from_nth_largest.is_valid) {
            unsigned int n = this->mirna_perturbation_target.elements_from_nth_largest.n;
            if(n > this->perturbed_patient.tumor_mirnas.profile.size()) {
                std::cerr << "error: n = " << n << ", this->perturbed_patient.tumor_mirnas.profile.size() = " << this->perturbed_patient.tumor_mirnas.profile.size() << "\n";
                exit(1);
            }
                
            std::cout << "generating " << n << " standard normal values\n";
            Timer::start();
            std::default_random_engine generator;
            std::normal_distribution<double> distribution(0.0, 1.0);
            std::vector<double> samples;
            samples.reserve(n);
            for(int i = 0; i < n; i++) {
                samples[i] = distribution(generator);
            }
            std::cout << "done, \n";
            Timer::stop();

            unsigned int i = 0;
            std::vector<Mirna_id> to_process;
            to_process.reserve(n);
            for(auto & e : mirna_profile_sorted) {
                auto & list_of_ids = e.second;
                if(i < n) {
                    for(auto & mirna_id : list_of_ids) {
                        if(to_process.size() < n) {
                            to_process.push_back(mirna_id);
                        }
                    }
                }
                i += list_of_ids.size();
            }
            if(to_process.size() != n) {
                std::cout << "to_process.size() = " << to_process.size() << ", n = " << n << "\n";
                exit(1);
            }
            for(auto & mirna_id : to_process) {
                double reads = this->perturbed_patient.tumor_mirnas.profile.at(mirna_id).to_reads();
                // when studying the gaussian perturbation you may to use (1 + sigma + exp(sigma/expression_level)) instead of the formula used
                reads *= (1 + this->mirna_perturbation_amplifier);
                this->perturbed_patient.tumor_mirnas.profile[mirna_id] = Reads(reads);   
            }
        } else {
            std::cerr << "error: exception in this->mirna_perturbation_type == Perturbation_type::Gaussian_perturbation\n";
            exit(1);
        }
    }

    double total_reads_after_perturbation = 0;
    for(auto & e : this->perturbed_patient.tumor_mirnas.profile) {
        double reads = e.second.to_reads();
        total_reads_after_perturbation += reads;
    }
    double total_reads = this->patient.tumor_mirnas.profile.begin()->second.get_grand_total() - (total_reads_before_perturbation - total_reads_after_perturbation);
    for(auto & e : this->perturbed_patient.tumor_mirnas.profile) {
        e.second.normalize_reads(total_reads);
    }
    /*
      Note that this is computed on the values of the reads that comprise the filtered reads.
      You could be interested in leaving the lambda as the one it was before.
    */
    double lambda_adjustment_due_to_mirnas = total_reads / this->patient.tumor_mirnas.profile.begin()->second.get_grand_total();
    std::cout << "lambda_adjustment_due_to_mirnas = " << lambda_adjustment_due_to_mirnas << "\n";

    // perturb genes
    std::cout << "TODO: do the same for genes, better to insert all the code in a lambda function and call first for mirnas and then for genes\n";
    
    this->matchings_predictor = std::unique_ptr<Matchings_predictor>(new Matchings_predictor(this->perturbed_patient, simulation_id.str()));
    // it is very important to adjust the value of lambda stored inside the Matchings_predictor object
    this->matchings_predictor->lambda *= lambda_adjustment_due_to_mirnas;
}
