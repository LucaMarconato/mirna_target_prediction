#include "mirna_expression_profile.hpp"

#include <cstdlib>

#include <boost/filesystem.hpp>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

#include "mirna.hpp"

#define Mep Mirna_expression_profile

void Mep::load_from_gdc_file(std::string filename, std::string patient_folder)
{
    if(!boost::filesystem::exists(patient_folder + filename)) {
        std::cout << "\"" << patient_folder + filename << "\" does not exists, skipping it\n";
    } else {
        std::cout << "parsing \"" << patient_folder + filename << "\"\n";
        io::CSVReader<4, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(patient_folder + filename);
        in.read_header(io::ignore_extra_column, "miRNA_ID", "read_count", "reads_per_million_miRNA_mapped", "cross-mapped");

        Output_buffer ob(patient_folder + "mirna_not_recognized.tsv", 10000, 1000);
        std::string s = "mirna_family\treads\n";

        std::string mirna_family, reads, unused0, unused1;
        while(in.read_row(mirna_family, reads, unused0, unused1)) {
            Mirna mirna(mirna_family);
            unsigned long long reads_ull = std::strtoull(reads.c_str(), nullptr, 10);
            auto e = Mirna::mirna_id_dictionary.left.find(mirna);
            if(e == Mirna::mirna_id_dictionary.left.end()) {
                this->not_recognized_distinct_mirnas++;
                this->not_recognized_reads += reads_ull;
                std::string line = mirna_family + "\t" + reads + "\n";
                ob.add_chunk(line);
            } else {
                this->recognized_distinct_mirnas++;
                this->recognized_reads += reads_ull;
                Mirna_id mirna_id = e->second;
                Reads reads(reads_ull);
                Expression expression(reads);
                if(this->profile.find(mirna_id) != this->profile.end()) {
                    std::cerr << "error: mirna_id " << mirna_id << " already present in this->profile\n";
                    exit(1);
                }
                this->profile[mirna_id] = expression;
            }
        }
        this->total_distinct_mirnas = this->recognized_distinct_mirnas + this->not_recognized_distinct_mirnas;
        this->total_reads = this->recognized_reads + this->not_recognized_reads;
        for(auto & e : this->profile) {
            e.second.normalize_reads(total_reads);
        }
        this->filter(1.0);
        this->initialized = true;
    }
}

void Mep::print_statistics()
{
    if(this->initialized) {
        std::cout << "recognized_distinct_mirnas/total_distinct_mirnas: " << recognized_distinct_mirnas << "/" << total_distinct_mirnas << " = " << ((double)recognized_distinct_mirnas)/total_distinct_mirnas << "\n";
        std::cout << "recognized_reads/total_reads: " << recognized_reads << "/" << total_reads << " = " << ((double)recognized_reads)/total_reads << "\n";
        std::cout << "filtered_out_distinct_mirnas/total_distinct_mirnas: " << filtered_out_distinct_mirnas << "/" << total_distinct_mirnas << " = " << ((double)filtered_out_distinct_mirnas)/total_distinct_mirnas << "\n";
        std::cout << "filtered_out_reads/total_reads: " << filtered_out_reads << "/" << total_reads << " = " << ((double)filtered_out_reads)/total_reads << "\n";
        std::cout << "remaining = recognized - filtered_out\n";
        std::cout << "remaining/total_distinct_mirnas: " << (recognized_distinct_mirnas - filtered_out_distinct_mirnas) << "/" << total_distinct_mirnas << " = " << ((double)(recognized_distinct_mirnas - filtered_out_distinct_mirnas))/total_distinct_mirnas << "\n";
        std::cout << "remaining/total_reads: " << (recognized_reads - filtered_out_reads) << "/" << total_reads << " = " << ((double)(recognized_reads - filtered_out_reads))/total_reads << "\n";
    }
}

void Mep::filter(double threshold_rpm)
{
    unsigned long long newly_filtered_out_distinct_mirnas = 0;
    unsigned long long newly_filtered_out_reads = 0;
    for(std::unordered_map<Mirna_id, Expression>::iterator it = this->profile.begin(); it != this->profile.end();) {
        double rpm = it->second.to_rpm();
        unsigned long long reads = it->second.to_reads();
        if(rpm < threshold_rpm) {
            newly_filtered_out_distinct_mirnas++;
            newly_filtered_out_reads += reads;
            it = this->profile.erase(it++);
        } else {
            ++it;
        }
    }
    std::cout << "newly_filtered_out_distinc_mirnas/recognized_distinct_mirnas: " << newly_filtered_out_distinct_mirnas << "/" << recognized_distinct_mirnas << " = " << ((double)newly_filtered_out_distinct_mirnas)/recognized_distinct_mirnas << "\n";
    std::cout << "newly_filtered_out_reads/recognized_reads: " << newly_filtered_out_reads << "/" << recognized_reads << " = " << ((double)newly_filtered_out_reads)/recognized_reads << "\n";
    this->filtered_out_distinct_mirnas += newly_filtered_out_distinct_mirnas;
    this->filtered_out_reads += newly_filtered_out_reads;
}
