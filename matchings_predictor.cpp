#include "matchings_predictor.hpp"

#include "global_parameters.hpp"

Matchings_predictor::Matchings_predictor(Patient & patient) : patient(patient)
{
    // TODO: for the moment I am just considering tumor mirnas and tumor clusters. In the future I should consider the difference between normal and tumor expressions.
    // note that these are copies, we do not want to modify the original expressions
    this->mirnas = patient.tumor_mirnas;
    this->clusters = patient.tumor_clusters;
}

void Matchings_predictor::compute()
{
    std::cout << "for the moment, just a trivial explicit Euler scheme\n";
    double lambda = 1;
    unsigned long long max_steps = 10000;

    /*
      The mirnas do not sum to one because of the filtering procedure, so we normalize the expression profile.
      After this normalization it will be improper to talk about RPM.
    */
    unsigned long long total_considered_reads = 0;
    for(auto & e : this->mirnas.profile) {
        total_considered_reads += e.second.to_reads();
        e.second.valid_rpm = false;        
    }
    for(auto & e : this->mirnas.profile) {
        e.second.normalize_reads(total_considered_reads);
    }

    /*
      The mirna expression proile has been just normalized to 1.
      The clusters are already normalized to 1, by construtction.
      Anyway, we check these condition explicitly
    */      
    double sum0 = 0, sum1 = 0;
    for(auto & e : this->mirnas.profile) {
        sum0 += e.second.to_relative_expression();
    }
    for(auto & e : this->clusters.profile) {
        sum1 += e.second.to_relative_expression();
    }
    if(std::abs(sum0 - 1.0) > Global_parameters::epsilon) {
        std::cerr << "error: std::abs(sum0 - 1.0) = " << std::abs(sum0 - 1.0) << "\n";
        exit(1);
    }
    if(std::abs(sum1 - 1.0) > Global_parameters::epsilon) {
        std::cerr << "error: std::abs(sum1 - 1.0) = " << std::abs(sum1 - 1.0) << "\n";
        exit(1);
    }
    
    for(unsigned long long t = 0; t < max_steps; t++) {
        
    }
}
