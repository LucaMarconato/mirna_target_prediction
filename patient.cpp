#include "patient.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <set>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "timer.hpp"

Patient::Patient(std::string case_id) : case_id(case_id)
{
    std::string patient_folder = "./data/patients/" + case_id + "/";
    if(!boost::filesystem::exists(patient_folder)) {
        std::cerr << "error: \"" << patient_folder << "\" not found\n";
        exit(1);
    } else {
        std::string patient_file = patient_folder + case_id + ".bin";
        if(!boost::filesystem::exists(patient_file)) {
            std::ifstream in(patient_folder + "info.json");
            std::stringstream buffer;
            buffer << in.rdbuf();
            // std::cout << "buffer.str() = " << buffer.str() << "\n";
            auto j = json::parse(buffer.str());
            in.close();
            
            std::string file_prefix = patient_folder;
            std::string normal_mirnas_file = j["normal_mirnas"]["uuid"].get<std::string>() + "/" + j["normal_mirnas"]["file"].get<std::string>();
            std::string normal_genes_file = j["normal_genes"]["uuid"].get<std::string>() + "/" + j["normal_genes"]["file"].get<std::string>();
            std::string tumor_mirnas_file = j["tumor_mirnas"]["uuid"].get<std::string>() + "/" + j["tumor_mirnas"]["file"].get<std::string>();
            std::string tumor_genes_file = j["tumor_genes"]["uuid"].get<std::string>() + "/" + j["tumor_genes"]["file"].get<std::string>();
            // std::cout << "normal_mirnas_file = " << normal_mirnas_file << ", normal_genes_file = " << normal_genes_file << ", tumor_mirnas_file = " << tumor_mirnas_file << ", tumor_genes_file = " << tumor_genes_file << "\n";

            this->normal_mirnas.load_from_gdc_file(normal_mirnas_file, patient_folder);
            this->normal_genes.load_from_gdc_file(normal_genes_file, patient_folder);
            this->tumor_mirnas.load_from_gdc_file(tumor_mirnas_file, patient_folder);
            this->tumor_genes.load_from_gdc_file(tumor_genes_file, patient_folder);

            std::set<Mirna_id> mirnas;
            for(auto & e : this->normal_mirnas.profile) {
                mirnas.insert(e.first);
            }
            for(auto & e : this->tumor_mirnas.profile) {
                mirnas.insert(e.first);
            }
            std::set<Gene_id> genes;
            for(auto & e : this->normal_genes.profile) {
                genes.insert(e.first);
            }
            for(auto & e : this->tumor_genes.profile) {
                genes.insert(e.first);
            }
            this->interaction_graph.build_interaction_graph(mirnas, genes);
            
            std::cout << "writing " << patient_file << "\n";
            Timer::start();
            std::ofstream out(patient_file, std::ios::binary);
            boost::archive::binary_oarchive oa(out);
            oa << this->case_id;
            oa << this->normal_mirnas << this->normal_genes << this->tumor_mirnas << this->tumor_genes;
            oa << this->interaction_graph;
            out.close();
            std::cout << "written, ";
            Timer::stop();
        } else {
            std::cout << "loading " << patient_file << "\n";
            Timer::start();
            std::ifstream in(patient_file, std::ios::binary);
            boost::archive::binary_iarchive ia(in);
            in >> this->case_id;
            ia >> this->normal_mirnas >> this->normal_genes >> this->tumor_mirnas >> this->tumor_genes;            
            ia >> this->interaction_graph;
            in.close();
            std::cout << "loaded, ";
            Timer::stop();
        }
        this->normal_mirnas.print_statistics();
        this->normal_genes.print_statistics();
        this->tumor_mirnas.print_statistics();
        this->tumor_genes.print_statistics();
        this->interaction_graph.print_statistics();
    }
}
