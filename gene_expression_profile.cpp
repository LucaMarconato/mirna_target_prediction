#include "gene_expression_profile.hpp"

#include <cstdlib>
#include <sstream>

#include <boost/filesystem.hpp>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

#include "gene.hpp"

#define Gep Gene_expression_profile

void Gep::load_from_gdc_file(std::string filename, std::string patient_folder)
{
    if(!boost::filesystem::exists(patient_folder + filename)) {
        std::cout << "\"" << patient_folder + filename << "\" does not exists, skipping it\n";
    } else {
        std::cout << "parsing \"" << patient_folder + filename << "\"\n";
        io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(patient_folder + filename);
        
        Output_buffer ob(patient_folder + "mrna_not_recognized.tsv", 10000, 1000);
        std::string s = "gene_id\tgene_id_version\treads\n";

        std::string gene_id_and_version, reads;
        while(in.read_row(gene_id_and_version, reads)) {
            unsigned long long reads_ull = std::strtoull(reads.c_str(), nullptr, 10);
            if(gene_id_and_version == "__no_feature" || gene_id_and_version == "__ambiguos" || gene_id_and_version == "__too_low_aQual" || gene_id_and_version == "__not_aligned" || gene_id_and_version == "__alignment_not_unique" || gene_id_and_version == "__ambiguous") {
                // WARNING! the information about these reads could be relevant, so mind to consider them
                // this sum is an overestimation, since some reads could belong to many of the above classes
                // discarded_reads += reads_ull;
                continue;
            }
            Gene gene(gene_id_and_version, "", "");

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

        for(auto & e : this->profile) {
            e.second.normalize_reads(total_reads);
        }
        this->initialized = true;
    }
}

void Gep::print_statistics()
{
    if(this->initialized) {
        total_distinct_genes = this->recognized_distinct_genes + this->not_recognized_distinct_genes;
        total_reads = this->recognized_reads + this->not_recognized_reads;
        std::cout << "recognized_distinct_genes/total_distinct_genes: " << recognized_distinct_genes << "/" << total_distinct_genes << " = " << ((double)recognized_distinct_genes)/total_distinct_genes << "\n";
        std::cout << "recognized_reads/total_reads: " << recognized_reads << "/" << total_reads << " = " << ((double)recognized_reads)/total_reads << "\n";
        // std::cout << "discarded_reads/total_reads: " << discarded_reads << "/" << total_reads << " = " << ((double)discarded_reads)/total_reads << "\n";        
    }
}
