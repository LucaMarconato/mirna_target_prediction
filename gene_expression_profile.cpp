#include "gene_expression_profile.hpp"

#include <cstdlib>
#include <sstream>

#include <boost/filesystem.hpp>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

#include "global_parameters.hpp"
#include "gene.hpp"

#define Gep Gene_expression_profile

Gep::Gene_expression_profile() {}

Gep::Gep(const Gep & obj)
{
    this->profile = obj.profile;
    this->recognized_distinct_genes = obj.recognized_distinct_genes;
    this->recognized_reads = obj.recognized_reads;
    this->not_recognized_distinct_genes = obj.not_recognized_distinct_genes;
    this->not_recognized_reads = obj.not_recognized_reads;
    this->total_distinct_genes = obj.total_distinct_genes;
    this->total_reads = obj.total_reads;
    this->filtered_out_distinct_genes = obj.filtered_out_distinct_genes;
    this->filtered_out_reads = obj.filtered_out_reads;
}

void swap(Gep & obj1, Gep & obj2)
{
    std::swap(obj1.profile, obj2.profile);
    std::swap(obj1.recognized_distinct_genes, obj2.recognized_distinct_genes);
    std::swap(obj1.recognized_reads, obj2.recognized_reads);
    std::swap(obj1.not_recognized_distinct_genes, obj2.not_recognized_distinct_genes);
    std::swap(obj1.not_recognized_reads, obj2.not_recognized_reads);
    std::swap(obj1.total_distinct_genes, obj2.total_distinct_genes);
    std::swap(obj1.total_reads, obj2.total_reads);
    std::swap(obj1.filtered_out_distinct_genes, obj2.filtered_out_distinct_genes);
    std::swap(obj1.filtered_out_reads, obj2.filtered_out_reads);
}

Gep & Gep::operator=(Gep obj)
{
    swap(*this, obj);
    return *this;
}

void Gep::load_from_gdc_file(std::string filename, std::string patient_folder)
{
    if(!boost::filesystem::exists(patient_folder + filename)) {
        std::cout << "\"" << patient_folder + filename << "\" does not exists, skipping it\n";
    } else {
        std::cout << "parsing \"" << patient_folder + filename << "\"\n";
        io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(patient_folder + filename);
        
        Output_buffer ob(patient_folder + "gene_not_recognized.tsv", 10000, 1000);
        std::string s = "gene_id\tgene_id_version\treads\n";

        std::string gene_id_and_version, reads;
        while(in.read_row(gene_id_and_version, reads)) {
            if(gene_id_and_version == "__no_feature" || gene_id_and_version == "__ambiguous" || gene_id_and_version == "__too_low_aQual" || gene_id_and_version == "__not_aligned" || gene_id_and_version == "__alignment_not_unique") {
                // WARNING! the information about these reads could be relevant, so mind to consider them
                // this sum is an overestimation, since some reads could belong to many of the above classes
                // discarded_reads += reads_ull;
                // TODO: the ratios of the reads considered should be comparable across all samples
                continue;
            }
            Gene gene(gene_id_and_version, "", "");
            double reads_ull = std::strtod(reads.c_str(), nullptr);
            auto e = Gene::gene_id_dictionary.left.find(gene);
            if(e == Gene::gene_id_dictionary.left.end()) {
                this->not_recognized_distinct_genes++;
                this->not_recognized_reads += reads_ull;
                std::stringstream ss;
                ss << gene.gene_id << "\t" << gene.gene_id_version << "\t" << reads << "\n";
                std::string line = ss.str();
                ob.add_chunk(line);
            } else {
                this->recognized_distinct_genes++;
                this->recognized_reads += reads_ull;
                Gene_id gene_id = e->second;
                Reads reads(reads_ull);
                Expression expression(reads);
                if(this->profile.find(gene_id) != this->profile.end()) {
                    std::cerr << "error: gene_id " << gene_id << " already present in this->profile\n";
                    exit(1);
                }
                this->profile[gene_id] = expression;
            }
        }
        this->total_distinct_genes = this->recognized_distinct_genes + this->not_recognized_distinct_genes;
        this->total_reads = this->recognized_reads + this->not_recognized_reads;
        for(auto & e : this->profile) {
            e.second.normalize_reads(total_reads);
        }
        // TODO: the filter threshold shuold be decided by considering many patients
        this->filter(Global_parameters::gene_threshold_rpm);
        this->initialized = true;
    }
}

void Gep::print_statistics()
{
    if(this->initialized) {
        std::cout << "recognized_distinct_genes/total_distinct_genes: " << recognized_distinct_genes << "/" << total_distinct_genes << " = " << ((double)recognized_distinct_genes)/total_distinct_genes << "\n";
        std::cout << "recognized_reads/total_reads: " << recognized_reads << "/" << total_reads << " = " << ((double)recognized_reads)/total_reads << "\n";
        // std::cout << "filtered_out_distinct_genes/total_distinct_genes: " << filtered_out_distinct_genes << "/" << total_distinct_genes << " = " << ((double)filtered_out_distinct_genes)/total_distinct_genes << "\n";
        // std::cout << "filtered_out_reads/total_reads: " << filtered_out_reads << "/" << total_reads << " = " << ((double)filtered_out_reads)/total_reads << "\n";
        std::cout << "remaining = recognized - filtered_out\n";
        std::cout << "remaining/total_distinct_genes: " << (recognized_distinct_genes - filtered_out_distinct_genes) << "/" << total_distinct_genes << " = " << ((double)(recognized_distinct_genes - filtered_out_distinct_genes))/total_distinct_genes << "\n";
        std::cout << "remaining/total_reads: " << (recognized_reads - filtered_out_reads) << "/" << total_reads << " = " << ((double)(recognized_reads - filtered_out_reads))/total_reads << "\n";
        // std::cout << "discarded_reads/total_reads: " << discarded_reads << "/" << total_reads << " = " << ((double)discarded_reads)/total_reads << "\n";
    }
}

void Gep::filter(double threshold_rpm)
{
    unsigned long long newly_filtered_out_distinct_genes = 0;
    double newly_filtered_out_reads = 0;
    for(std::unordered_map<Gene_id, Expression>::iterator it = this->profile.begin(); it != this->profile.end();) {
//        if(it->first == 19322) {
//            int a = 0;
//        }
        double rpm = it->second.to_rpm();
        double reads = it->second.to_reads();
        if(rpm < threshold_rpm) {
            newly_filtered_out_distinct_genes++;
            newly_filtered_out_reads += reads;
            it = this->profile.erase(it++);
        } else {
            ++it;
        }
    }
    // if we are performing the filtering process only once then the following lines print redundant information
    // std::cout << "newly_filtered_out_distinct_genes/recognized_distinct_genes: " << newly_filtered_out_distinct_genes << "/" << recognized_distinct_genes << " = " << ((double)newly_filtered_out_distinct_genes)/recognized_distinct_genes << "\n";
    // std::cout << "newly_filtered_out_reads/recognized_reads: " << newly_filtered_out_reads << "/" << recognized_reads << " = " << ((double)newly_filtered_out_reads)/recognized_reads << "\n";
    
    this->filtered_out_distinct_genes += newly_filtered_out_distinct_genes;
    this->filtered_out_reads += newly_filtered_out_reads;
    std::cout << "remaining = recognized - filtered_out\n";
    std::cout << "remaining/recognized_distinct_genes: " << recognized_distinct_genes - filtered_out_distinct_genes << "/" << recognized_distinct_genes << " = " << ((double)recognized_distinct_genes - filtered_out_distinct_genes)/recognized_distinct_genes << "\n";
    std::cout << "remaining/recognized_reads: " << recognized_reads - filtered_out_reads << "/" << recognized_reads << " = " << ((double)recognized_reads - filtered_out_reads)/recognized_reads << "\n";
}
