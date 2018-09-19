#ifndef MATCHINGS_PREDICTOR_H
#define MATCHINGS_PREDICTOR_H

#include "patient.hpp"

class Matchings_predictor {
    Patient & patient;
    Mirna_expression_profile mirnas;
    Cluster_expression_profile clusters;    
public:
    Matchings_predictor(Patient & patient);
    void compute();
    
};

#endif // MATCHINGS_PREDICTOR_H
