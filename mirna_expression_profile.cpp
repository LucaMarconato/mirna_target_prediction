#include "mirna_expression_profile.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <strasser/csv.h>

#include "global_parameters.hpp"
#include "mirna.hpp"

#define Mep Mirna_expression_profile

Mep::Mirna_expression_profile()
{
}

Mep::Mep(const Mep& obj)
{
    this->profile = obj.profile;
    this->distinct_mirnas = obj.distinct_mirnas;
    this->total_reads = obj.total_reads;
    this->filtered_out_distinct_mirnas = obj.filtered_out_distinct_mirnas;
    this->filtered_out_reads = obj.filtered_out_reads;
}

void swap(Mep& obj1, Mep& obj2)
{
    std::swap(obj1.profile, obj2.profile);
    std::swap(obj1.distinct_mirnas, obj2.distinct_mirnas);
    std::swap(obj1.total_reads, obj2.total_reads);
    std::swap(obj1.filtered_out_distinct_mirnas, obj2.filtered_out_distinct_mirnas);
    std::swap(obj1.filtered_out_reads, obj2.filtered_out_reads);
}

Mep& Mep::operator=(Mep obj)
{
    swap(*this, obj);
    return *this;
}

void Mep::load_from_file(std::string patient_folder)
{
    std::string filename = patient_folder + "mirna_expression_profile.tsv";
    if (!boost::filesystem::exists(filename)) {
        std::cout << "\"" << filename << "\" does not exists, aborting\n";
        exit(1);
    } else {
        std::cout << "parsing \"" << filename << "\"\n";
        io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(filename);
        in.read_header(io::ignore_extra_column, "mirbase_id", "reads");

        std::string mirbase_id, reads;
        while (in.read_row(mirbase_id, reads)) {
            Mirna mirna(mirbase_id);
            double reads_ull = std::strtod(reads.c_str(), nullptr);
            auto e = Mirna::mirna_id_dictionary.left.find(mirna);
            if (e == Mirna::mirna_id_dictionary.left.end()) {
                std::cerr << "error: mirna not recognized, mirbase_id = " << mirbase_id << "\n";
                exit(1);
            } else {
                this->distinct_mirnas++;
                this->total_reads += reads_ull;
                Mirna_id mirna_id = e->second;
                Reads reads(reads_ull);
                Expression expression(reads);
                if (this->profile.find(mirna_id) != this->profile.end()) {
                    std::cerr << "error: mirna_id " << mirna_id << ", mirbase_id = " << mirbase_id << " already present in this->profile\n";
                    exit(1);
                }
                this->profile[mirna_id] = expression;
            }
        }
        for (auto& e : this->profile) {
            e.second.normalize_reads(total_reads);
        }
        // TODO: the filter threshold should be decided by considering many patients
        this->filter(Global_parameters::mirna_threshold_rpm);
        this->initialized = true;
    }
}

void Mep::print_statistics()
{
    if (this->initialized) {
        std::cout << "distinct_mirnas = " << distinct_mirnas << ", total_reads = " << total_reads << "\n";
        std::cout << "(distinct_mirnas - filtered_out_distinct_mirnas)/distinct_mirnas: " << (distinct_mirnas - filtered_out_distinct_mirnas) << "/" << distinct_mirnas << " = "
                  << ((double)(distinct_mirnas - filtered_out_distinct_mirnas)) / distinct_mirnas << "\n";
        std::cout << "(total_reads - filtered_out_reads)/total_reads: " << (total_reads - filtered_out_reads) << "/" << total_reads << " = " << ((double)(total_reads - filtered_out_reads)) / total_reads << "\n";
    }
}

void Mep::filter(double threshold_rpm)
{
    unsigned long long newly_filtered_out_distinct_mirnas = 0;
    double newly_filtered_out_reads = 0;
    for (std::unordered_map<Mirna_id, Expression>::iterator it = this->profile.begin(); it != this->profile.end();) {
        double rpm = it->second.to_rpm();
        double reads = it->second.to_reads();
        if (rpm < threshold_rpm) {
            newly_filtered_out_distinct_mirnas++;
            newly_filtered_out_reads += reads;
            it = this->profile.erase(it++);
        } else {
            ++it;
        }
    }
    // if we are performing the filtering process only once then the following lines print redundant information
    // std::cout << "newly_filtered_out_distinct_mirnas/recognized_distinct_mirnas: " << newly_filtered_out_distinct_mirnas << "/" << recognized_distinct_mirnas << " = " << ((double)newly_filtered_out_distinct_mirnas)/recognized_distinct_mirnas << "\n";
    // std::cout << "newly_filtered_out_reads/recognized_reads: " << newly_filtered_out_reads << "/" << recognized_reads << " = " << ((double)newly_filtered_out_reads)/recognized_reads << "\n";

    this->filtered_out_distinct_mirnas += newly_filtered_out_distinct_mirnas;
    this->filtered_out_reads += newly_filtered_out_reads;
    std::cout << "(distinct_mirnas - filtered_out_distinct_mirnas)/distinct_mirnas: " << distinct_mirnas - filtered_out_distinct_mirnas << "/" << distinct_mirnas << " = "
              << ((double)distinct_mirnas - filtered_out_distinct_mirnas) / distinct_mirnas << "\n";
    std::cout << "(total_reads - filtered_out_reads)/recognized_reads: " << total_reads - filtered_out_reads << "/" << total_reads << " = " << ((double)total_reads - filtered_out_reads) / total_reads << "\n";
}
