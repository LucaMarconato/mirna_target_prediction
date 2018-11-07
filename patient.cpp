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
    this->mirna_expression_profile = obj.mirna_expression_profile;
    this->gene_expression_profile = obj.gene_expression_profile;
    this->cluster_expression_profile = obj.cluster_expression_profile;
    this->interaction_graph = obj.interaction_graph;
}

void swap(Patient & obj1, Patient & obj2)
{
    std::swap(obj1.case_id, obj2.case_id);
    std::swap(obj1.mirna_expression_profile, obj2.mirna_expression_profile);
    std::swap(obj1.gene_expression_profile, obj2.gene_expression_profile);
    std::swap(obj1.cluster_expression_profile, obj2.cluster_expression_profile);
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
            std::cout << "\n***loading expression profiles***\n";
            this->mirna_expression_profile.load_from_file(patient_folder);
            std::cout << "\n";
            this->gene_expression_profile.load_from_file(patient_folder);
            std::cout << "\n";

            std::set<Mirna_id> mirnas;
            for(auto & e : this->mirna_expression_profile.profile) {
                mirnas.insert(e.first);
            }
            std::set<Gene_id> genes;
            for(auto & e : this->gene_expression_profile.profile) {
                genes.insert(e.first);
            }

            std::cout << "\n***building interaction graph***\n";
            this->interaction_graph.build_interaction_graph(mirnas, genes);
            this->generate_cluster_expression_profiles();

            if(!dont_use_boost_serialization_here) {
                std::cout << "writing " << patient_file << "\n";
                Timer::start();
                std::ofstream out(patient_file, std::ios::binary);
                boost::archive::binary_oarchive oa(out);
                oa << this->mirna_expression_profile << this->gene_expression_profile << this->cluster_expression_profile;
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
            ia >> this->mirna_expression_profile >> this->gene_expression_profile >> this->cluster_expression_profile;
            ia >> this->interaction_graph;
            in.close();
            std::cout << "loaded, ";
            Timer::stop();
        }
        std::cout << "\n***expression profiles statistics***\n";
        this->mirna_expression_profile.print_statistics();
        std::cout << "\n";
        this->gene_expression_profile.print_statistics();
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
    this->cluster_expression_profile = Cluster_expression_profile(this->gene_expression_profile, this->interaction_graph);
}

void Patient::export_expression_profiles(std::string patient_folder)
{
    // export expression profiles for the mirnas and the genes involved in the interaction graph
    // I have used two lambda functions because in a previous version of the code there were two expression profiles for mirnas and two for genes
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

    export_mirna_expression_profile(this, patient_folder + "mirna_expression_profile.tsv", this->mirna_expression_profile);
    export_gene_expression_profile(this, patient_folder + "gene_expression_profile.tsv", this->gene_expression_profile);
}
