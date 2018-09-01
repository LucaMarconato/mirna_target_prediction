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
        std::string s = "mirna_family\tRPM\n";

        std::string mirna_family, unused0, rpm, unused1;
        while(in.read_row(mirna_family, unused0, rpm, unused1)) {
            Mirna mirna(mirna_family);
            double rpm_double = std::strtod(rpm.c_str(), nullptr);
            auto e = Mirna::mirna_id_dictionary.left.find(mirna);
            if(e == Mirna::mirna_id_dictionary.left.end()) {
                this->not_recognized_distinct_mirnas++;
                this->not_recognized_rpm += rpm_double;
                std::string line = mirna_family + "\t" + rpm + "\n";
                ob.add_chunk(line);
            } else {
                this->recognized_distinct_mirnas++;
                this->recognized_rpm += rpm_double;
                Mirna_id mirna_id = e->second;
                Rpm rpm(rpm_double);
                Expression expression(rpm);
                if(this->profile.find(mirna_id) != this->profile.end()) {
                    std::cerr << "error: mirna_id " << mirna_id << " already present in this->profile\n";
                    exit(1);
                }
                this->profile[mirna_id] = expression;
            }
        }
        this->initialized = true;
    }
}

void Mep::print_statistics()
{
    if(this->initialized) {
        this->total_distinct_mirnas = this->recognized_distinct_mirnas + this->not_recognized_distinct_mirnas;
        this->total_rpm = this->recognized_rpm + this->not_recognized_rpm;
        std::cout << "recognized_distinct_mirnas/total_distinct_mirnas: " << recognized_distinct_mirnas << "/" << total_distinct_mirnas << " = " << ((double)recognized_distinct_mirnas)/total_distinct_mirnas << "\n";
        std::cout << "recognized_rpm/total_rpm: " << recognized_rpm << "/" << total_rpm << " = " << ((double)recognized_rpm)/total_rpm << "\n";        
    }
}
