#include "global_parameters.hpp"

#include <string>
#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

double Global_parameters::mirna_threshold_rpm;
double Global_parameters::gene_threshold_rpm;

void Global_parameters::load_from_json()
{
    std::ifstream in("global_parameters.json");
    std::stringstream buffer;
    buffer << in.rdbuf();
    auto j = json::parse(buffer);
    in.close();
    
    Global_parameters::mirna_threshold_rpm = j["mirna_threshold_rpm"].get<double>();
    Global_parameters::gene_threshold_rpm = j["gene_threshold_rpm"].get<double>();
}
