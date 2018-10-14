#include "global_parameters.hpp"

#include <string>
#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

double Global_parameters::mirna_threshold_rpm;
double Global_parameters::gene_threshold_rpm;
// TODO: this has to be chosen wisely looking at the relative expression profiles of mirnas, genes, sites and clusters
double Global_parameters::epsilon = 0.00000001;
double Global_parameters::lambda;
bool Global_parameters::test_parallelization;

void Global_parameters::load_from_json()
{
    std::ifstream in("global_parameters.json");
    std::stringstream buffer;
    buffer << in.rdbuf();
    auto j = json::parse(buffer);
    in.close();
    
    Global_parameters::mirna_threshold_rpm = j["mirna_threshold_rpm"].get<double>();
    Global_parameters::gene_threshold_rpm = j["gene_threshold_rpm"].get<double>();
    Global_parameters::lambda = j["lambda"].get<double>();
    Global_parameters::test_parallelization = j["test_parallelization"].get<bool>();
}
