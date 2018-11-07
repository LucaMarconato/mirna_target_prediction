#include "gene_expression_profile.hpp"

#include <cstdlib>
#include <sstream>

#include <boost/filesystem.hpp>

#include <marconato/output_buffer/output_buffer.hpp>
#include <strasser/csv.h>

#include "gene.hpp"
#include "global_parameters.hpp"

#define Gep Gene_expression_profile

Gep::Gene_expression_profile()
{
}

Gep::Gep(const Gep& obj)
{
    this->profile = obj.profile;
    this->distinct_genes = obj.distinct_genes;
    this->total_reads = obj.total_reads;
    this->filtered_out_distinct_genes = obj.filtered_out_distinct_genes;
    this->filtered_out_reads = obj.filtered_out_reads;
}

void swap(Gep& obj1, Gep& obj2)
{
    std::swap(obj1.profile, obj2.profile);
    std::swap(obj1.distinct_genes, obj2.distinct_genes);
    std::swap(obj1.total_reads, obj2.total_reads);
    std::swap(obj1.filtered_out_distinct_genes, obj2.filtered_out_distinct_genes);
    std::swap(obj1.filtered_out_reads, obj2.filtered_out_reads);
}

Gep& Gep::operator=(Gep obj)
{
    swap(*this, obj);
    return *this;
}

void Gep::load_from_file(std::string patient_folder)
{
    std::string filename = patient_folder + "gene_expression_profile.tsv";
    if (!boost::filesystem::exists(filename)) {
        std::cout << "\"" << filename << "\" does not exists, aborting\n";
        exit(1);
    } else {
        std::cout << "parsing \"" << filename << "\"\n";
        io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(filename);
        in.read_header(io::ignore_extra_column, "gene_id", "reads");

        std::string gene_id_and_version, reads;
        while (in.read_row(gene_id_and_version, reads)) {
            Gene gene(gene_id_and_version, "", "");
            double reads_ull = std::strtod(reads.c_str(), nullptr);
            auto e = Gene::gene_id_dictionary.left.find(gene);
            if (e == Gene::gene_id_dictionary.left.end()) {
                std::cerr << "error: gene not recognized, gene_id_and_version = " << gene_id_and_version << "\n";
                exit(1);
            } else {
                this->distinct_genes++;
                this->total_reads += reads_ull;
                Gene_id gene_id = e->second;
                Reads reads(reads_ull);
                Expression expression(reads);
                if (this->profile.find(gene_id) != this->profile.end()) {
                    std::cerr << "error: gene_id " << gene_id << " already present in this->profile\n";
                    exit(1);
                }
                this->profile[gene_id] = expression;
            }
        }
        for (auto& e : this->profile) {
            e.second.normalize_reads(total_reads);
        }
        // TODO: the filter threshold shuold be decided by considering many patients
        this->filter(Global_parameters::gene_threshold_rpm);
        this->initialized = true;
    }
}

void Gep::print_statistics()
{
    if (this->initialized) {
        std::cout << "distinct_genes = " << distinct_genes << ", total_genes = " << distinct_genes << "\n";
        std::cout << "(distinct_genes - filtered_out_distinct_genes)/distinct_genes: " << (distinct_genes - filtered_out_distinct_genes) << "/" << distinct_genes << " = "
                  << ((double)(distinct_genes - filtered_out_distinct_genes)) / distinct_genes << "\n";
        std::cout << "(total_reads - filtered_out_reads)/total_reads: " << (total_reads - filtered_out_reads) << "/" << total_reads << " = " << ((double)(total_reads - filtered_out_reads)) / total_reads << "\n";
    }
}

void Gep::filter(double threshold_rpm)
{
    unsigned long long newly_filtered_out_distinct_genes = 0;
    double newly_filtered_out_reads = 0;
    for (std::unordered_map<Gene_id, Expression>::iterator it = this->profile.begin(); it != this->profile.end();) {
        double rpm = it->second.to_rpm();
        double reads = it->second.to_reads();
        if (rpm < threshold_rpm) {
            newly_filtered_out_distinct_genes++;
            newly_filtered_out_reads += reads;
            it = this->profile.erase(it++);
        } else {
            ++it;
        }
    }
    // if we are performing the filtering process only once then the following lines print redundant information
    // std::cout << "newly_filtered_out_distinct_genes/recognized_distinct_genes: " << newly_filtered_out_distinct_genes << "/" << recognized_distinct_genes << " = " << ((double)newly_filtered_out_distinct_genes)/recognized_distinct_genes
    // << "\n"; std::cout << "newly_filtered_out_reads/recognized_reads: " << newly_filtered_out_reads << "/" << recognized_reads << " = " << ((double)newly_filtered_out_reads)/recognized_reads << "\n";

    this->filtered_out_distinct_genes += newly_filtered_out_distinct_genes;
    this->filtered_out_reads += newly_filtered_out_reads;
    std::cout << "(distinct_genes - filtered_out_distinct_genes)/distinct_genes: " << distinct_genes - filtered_out_distinct_genes << "/" << distinct_genes << " = "
              << ((double)distinct_genes - filtered_out_distinct_genes) / distinct_genes << "\n";
    std::cout << "(total_reads - filtered_out_reads)/recognized_reads: " << total_reads - filtered_out_reads << "/" << total_reads << " = " << ((double)total_reads - filtered_out_reads) / total_reads << "\n";
}
