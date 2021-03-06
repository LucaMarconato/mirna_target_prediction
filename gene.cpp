#include "gene.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

#include <strasser/csv.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>

#include "timer.hpp"

boost::bimap<Gene, Gene_id> Gene::gene_id_dictionary;
std::unordered_set<std::string> Gene::ensembl_ids_used_in_targetscan;

Gene::Gene() {}

Gene::Gene(std::string gene_id, int gene_id_version, std::string gene_symbol, std::string transcript_id, int transcript_id_version)
    : gene_id(gene_id), gene_id_version(gene_id_version), gene_symbol(gene_symbol), transcript_id(transcript_id), transcript_id_version(transcript_id_version)
{
}

Gene::Gene(std::string gene_id_and_version, std::string gene_symbol, std::string transcript_id_and_version) : gene_symbol(gene_symbol)
{
    auto split_into_id_and_version = [](std::string& to_split, std::string& id, int& version) {
        size_t pos = to_split.find(".");
        if (pos != std::string::npos) {
            id = to_split.substr(0, pos);
            version = atoi(to_split.substr(pos + 1, to_split.size() - 1 - id.size()).c_str());
        } else {
            id = to_split;
            version = -1;
        }
    };
    split_into_id_and_version(gene_id_and_version, this->gene_id, this->gene_id_version);
    split_into_id_and_version(transcript_id_and_version, this->transcript_id, this->transcript_id_version);
    Gene gene(gene_id, gene_id_version, gene_symbol, transcript_id, transcript_id_version);
}

void Gene::initialize_gene_dictionary()
{
    if (!boost::filesystem::exists("gene_id_dictionary.bin")) {
        Gene_id i = 0;
        io::CSVReader<3, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/sites_with_scored_interactions.tsv");
        in.read_header(io::ignore_extra_column, "gene_id", "gene_symbol", "transcript_id");
        std::string column0, column1, column2;
        while (in.read_row(column0, column1, column2)) {
            Gene gene(column0, column1, column2);
            if (Gene::gene_id_dictionary.left.find(gene) == Gene::gene_id_dictionary.left.end()) {
                Gene::gene_id_dictionary.insert(boost::bimap<Gene, Gene_id>::value_type(gene, i));
                i++;
            }
        }
        std::cout << "writing gene_id_dictionary.bin\n";
        Timer::start();
        std::ofstream out("gene_id_dictionary.bin");
        boost::archive::binary_oarchive oa(out);
        oa << Gene::gene_id_dictionary;
        out.close();
        std::cout << "written, ";
        Timer::stop();

        std::stringstream ss;
        ss << "gene_id\tgene_id_cpp\n";
        for (auto& e : Gene::gene_id_dictionary.left) {
            ss << e.first.gene_id << "\t" << e.second << "\n";
        }
        out.open("data/processed/gene_id_dictionary.tsv");
        out << ss.str();
        out.close();
    } else {
        std::cout << "loading gene_id_dictionary.bin\n";
        Timer::start();
        std::ifstream in("gene_id_dictionary.bin");
        boost::archive::binary_iarchive ia(in);
        ia >> Gene::gene_id_dictionary;
        in.close();
        std::cout << "loaded, ";
        Timer::stop();
    }

    io::CSVReader<1, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/targetscan_genes.tsv");
    in.read_header(io::ignore_extra_column, "ensembl_id");
    std::string ensembl_id;
    while (in.read_row(ensembl_id)) {
        if (Gene::ensembl_ids_used_in_targetscan.find(ensembl_id) != Gene::ensembl_ids_used_in_targetscan.end()) {
            if(ensembl_id == "") {
                continue;
            }
            std::cerr << "error: ensembl_id = " << ensembl_id << " already present\n";
            exit(1);
        }
        Gene::ensembl_ids_used_in_targetscan.insert(ensembl_id);
    }
}

void Gene::print_gene_dictionary(unsigned int max_rows)
{
    if (max_rows != (unsigned int)-1) {
        std::cout << "printing at most " << max_rows << " rows\n";
    }
    unsigned int j = 0;
    for (auto& e : Gene::gene_id_dictionary.left) {
        if (j++ < max_rows) {
            Gene& gene = const_cast<Gene&>(e.first);
            Gene_id& i = const_cast<Gene_id&>(e.second);
            std::cout << "gene.gene_id = " << gene.gene_id << ", gene.gene_id_version = " << gene.gene_id_version << ", gene.gene_symbol = " << gene.gene_symbol << ", gene.transcript_id = " << gene.transcript_id
                      << ", gene.transcript_id_version = " << gene.transcript_id_version << ", i = " << i << "\n";
        }
    }
}

std::ostream& operator<<(std::ostream& stream, const Gene& o)
{
    return stream << "gene_id = " << o.gene_id << ", gene_id_version = " << o.gene_id_version << ", gene_symbol = " << o.gene_symbol << ", transcript_id = " << o.transcript_id << ", transcript_id_version = " << o.transcript_id_version << "\n";
}

bool operator<(Gene const& lhs, Gene const& rhs) { return lhs.gene_id < rhs.gene_id; }
