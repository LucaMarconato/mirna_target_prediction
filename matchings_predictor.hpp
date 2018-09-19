#ifndef MATCHINGS_PREDICTOR_H
#define MATCHINGS_PREDICTOR_H

#include "patient.hpp"

class Matchings_predictor {
    Patient & patient;
    std::unordered_map<Mirna_id, double> mirna_profile;
    std::unordered_map<Cluster *, double> cluster_profile;
public:
    Matchings_predictor(Patient & patient);
    void compute();
    
};

#endif // MATCHINGS_PREDICTOR_H
