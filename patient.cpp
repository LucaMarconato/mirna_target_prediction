#include "patient.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

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
            std::string normal_mirna_file = j["normal_mirna"]["uuid"].get<std::string>() + "/" + j["normal_mirna"]["file"].get<std::string>();
            std::string normal_mrna_file = j["normal_mrna"]["uuid"].get<std::string>() + "/" + j["normal_mrna"]["file"].get<std::string>();
            std::string tumor_mirna_file = j["tumor_mirna"]["uuid"].get<std::string>() + "/" + j["tumor_mirna"]["file"].get<std::string>();
            std::string tumor_mrna_file = j["tumor_mrna"]["uuid"].get<std::string>() + "/" + j["tumor_mrna"]["file"].get<std::string>();
            // std::cout << "normal_mirna_file = " << normal_mirna_file << ", normal_mrna_file = " << normal_mrna_file << ", tumor_mirna_file = " << tumor_mirna_file << ", tumor_mrna_file = " << tumor_mrna_file << "\n";

            this->normal_mirna.load_from_gdc_file(normal_mirna_file, patient_folder);
            this->normal_mrna.load_from_gdc_file(normal_mrna_file, patient_folder);
            this->tumor_mirna.load_from_gdc_file(tumor_mirna_file, patient_folder);
            this->tumor_mrna.load_from_gdc_file(tumor_mrna_file, patient_folder);

            // TODO: build interaction graph, which must not be static!!!!!!!!!!!
            
            std::cout << "writing " << patient_file << "\n";
            Timer::start();
            std::ofstream out(patient_file, std::ios::binary);
            boost::archive::binary_oarchive oa(out);
            oa << this->normal_mirna << this->normal_mrna << this->tumor_mirna << this->tumor_mrna;
            // TODO: save the interaction graph
            out.close();
            std::cout << "written, ";
            Timer::stop();
        } else {
            std::cout << "loading " << patient_file << "\n";
            Timer::start();
            std::ifstream in(patient_file, std::ios::binary);
            boost::archive::binary_iarchive ia(in);
            ia >> this->normal_mirna >> this->normal_mrna >> this->tumor_mirna >> this->tumor_mrna;
            // TODO: load the interaction graph
            in.close();
            std::cout << "loaded, ";
            Timer::stop();
        }
        this->normal_mirna.print_statistics();
        this->normal_mrna.print_statistics();
        this->tumor_mirna.print_statistics();
        this->tumor_mrna.print_statistics();
    }
}
