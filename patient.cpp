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

Patient::Patient() {}

Patient::Patient(const Patient & obj)
{
    this->case_id = obj.case_id;
    this->normal_mirnas = obj.normal_mirnas;
    this->tumor_mirnas = obj.tumor_mirnas;
    this->normal_genes = obj.normal_genes;
    this->tumor_genes = obj.tumor_genes;
    this->normal_clusters = obj.normal_clusters;
    this->tumor_clusters = obj.tumor_clusters;
    this->interaction_graph = obj.interaction_graph;
}

void swap(Patient & obj1, Patient & obj2)
{
    std::swap(obj1.case_id, obj2.case_id);
    std::swap(obj1.normal_mirnas, obj2.normal_mirnas);
    std::swap(obj1.tumor_mirnas, obj2.tumor_mirnas);
    std::swap(obj1.normal_genes, obj2.normal_genes);
    std::swap(obj1.tumor_genes, obj2.tumor_genes);
    std::swap(obj1.normal_clusters, obj2.normal_clusters);
    std::swap(obj1.tumor_clusters, obj2.tumor_clusters);
    std::swap(obj1.interaction_graph, obj2.interaction_graph);
}

Patient & Patient::operator=(Patient obj)
{
    swap(*this, obj);
    return *this;
}

Patient::Patient(std::string case_id, bool export_data) : case_id(case_id)
{
    std::string patient_folder = "./data/patients/" + case_id + "/";
    if(!boost::filesystem::exists(patient_folder)) {
        std::cerr << "error: \"" << patient_folder << "\" not found\n";
        exit(1);
    } else {
        std::string patient_file = patient_folder + case_id + ".bin";
        // TODO: bug with boost serialization
        bool dont_use_boost_serialization_here = !boost::filesystem::exists(patient_file) || true;
        if(dont_use_boost_serialization_here) {
            std::ifstream in(patient_folder + "info.json");
            std::stringstream buffer;
            buffer << in.rdbuf();
            // std::cout << "buffer.str() = " << buffer.str() << "\n";
            auto j = json::parse(buffer.str());
            in.close();

            std::string file_prefix = patient_folder;
            // std::string normal_mirnas_file = j["normal_mirnas"]["uuid"].get<std::string>() + "/" + j["normal_mirnas"]["file"].get<std::string>();
            std::string normal_genes_file = j["normal_genes"]["uuid"].get<std::string>() + "/" + j["normal_genes"]["file"].get<std::string>();
            // std::string tumor_mirnas_file = j["tumor_mirnas"]["uuid"].get<std::string>() + "/" + j["tumor_mirnas"]["file"].get<std::string>();
            std::string tumor_genes_file = j["tumor_genes"]["uuid"].get<std::string>() + "/" + j["tumor_genes"]["file"].get<std::string>();
            // std::cout << "normal_mirnas_file = " << normal_mirnas_file << ", normal_genes_file = " << normal_genes_file << ", tumor_mirnas_file = " << tumor_mirnas_file << ", tumor_genes_file = " << tumor_genes_file << "\n";

            std::cout << "\n***loading expression profiles***\n";
            // this->normal_mirnas.load_from_gdc_file(normal_mirnas_file, patient_folder);
            this->normal_mirnas.load_from_gdc_file("normal", patient_folder);
            std::cout << "\n";
            this->normal_genes.load_from_gdc_file(normal_genes_file, patient_folder);
            std::cout << "\n";
            // this->tumor_mirnas.load_from_gdc_file(tumor_mirnas_file, patient_folder);
            this->tumor_mirnas.load_from_gdc_file("tumor", patient_folder);
            std::cout << "\n";
            this->tumor_genes.load_from_gdc_file(tumor_genes_file, patient_folder);
            std::cout << "\n";

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

            if(genes.find(19322) != genes.end()) {
                std::cout << "all right\n";
            }

            std::cout << "\n***building interaction graph***\n";
            this->interaction_graph.build_interaction_graph(mirnas, genes);
            this->generate_cluster_expression_profiles();

            if(!dont_use_boost_serialization_here) {
                std::cout << "writing " << patient_file << "\n";
                Timer::start();
                std::ofstream out(patient_file, std::ios::binary);
                boost::archive::binary_oarchive oa(out);
                oa << this->normal_mirnas << this->normal_genes << this->tumor_mirnas << this->tumor_genes << this->normal_clusters << this->tumor_clusters;
                oa << this->interaction_graph;
                out.close();
                std::cout << "written, ";
                Timer::stop();
            }
        } else {
            std::cout << "loading " << patient_file << "\n";
            Timer::start();
            std::ifstream in(patient_file, std::ios::binary);
            boost::archive::binary_iarchive ia(in);
            ia >> this->normal_mirnas >> this->normal_genes >> this->tumor_mirnas >> this->tumor_genes >> this->normal_clusters >> this->tumor_clusters;
            ia >> this->interaction_graph;
            in.close();
            std::cout << "loaded, ";
            Timer::stop();
        }
        std::cout << "\n***expression profiles statistics***\n";
        this->normal_mirnas.print_statistics();
        std::cout << "\n";
        this->normal_genes.print_statistics();
        std::cout << "\n";
        this->tumor_mirnas.print_statistics();
        std::cout << "\n";
        this->tumor_genes.print_statistics();
        std::cout << "\n";
        this->interaction_graph.print_statistics();
        std::cout << "\n";
        if(export_data) {
            this->interaction_graph.export_interactions_data(patient_folder);
            std::cout << "exporting expression profiles\n";
            this->export_expression_profiles(patient_folder);
            std::cout << "finished\n";
        }
    }
}

void Patient::generate_cluster_expression_profiles()
{
    this->normal_clusters = Cluster_expression_profile(this->normal_genes, this->interaction_graph);
    this->tumor_clusters = Cluster_expression_profile(this->tumor_genes, this->interaction_graph);
}

void Patient::export_expression_profiles(std::string patient_folder)
{
    // export expression profiles for the mirnas and the genes involved in the interaction graph
    auto export_mirna_expression_profile = [](Patient * patient, std::string filename, Mirna_expression_profile & ep) {
		std::stringstream ss;
        ss << "mirna_id\trpm\n";
        for(auto & e : ep.profile) {
            Mirna_id mirna_id = e.first;
            if(patient->interaction_graph.mirna_to_genes_arcs.find(mirna_id) != patient->interaction_graph.mirna_to_genes_arcs.end()) {
                ss << mirna_id << "\t" << e.second.to_rpm() << "\n";
            }
        }
        std::ofstream out(filename);
        out << ss.str();
        out.close();
	};
    auto export_gene_expression_profile = [](Patient * patient, std::string filename, Gene_expression_profile & ep) {
		std::stringstream ss;
        ss << "gene_id\trpm\n";
        for(auto & e : ep.profile) {
            Gene_id gene_id = e.first;
            if(patient->interaction_graph.gene_to_mirnas_arcs.find(gene_id) != patient->interaction_graph.gene_to_mirnas_arcs.end()) {
                ss << gene_id << "\t" << e.second.to_rpm() << "\n";
            }
        }
        std::ofstream out(filename);
        out << ss.str();
        out.close();
	};

    export_mirna_expression_profile(this, patient_folder + "normal_mirna_expression_profile.tsv", normal_mirnas);
    export_mirna_expression_profile(this, patient_folder + "tumor_mirna_expression_profile.tsv", tumor_mirnas);
    export_gene_expression_profile(this, patient_folder + "normal_gene_expression_profile.tsv", normal_genes);
    export_gene_expression_profile(this, patient_folder + "tumor_gene_expression_profile.tsv", tumor_genes);
}
