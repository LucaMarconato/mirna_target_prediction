#include "gene.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#include <strasser/csv.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

boost::bimap<Gene, Gene_id> Gene::gene_id_dictionary;

Gene::Gene() {}

Gene::Gene(std::string gene_id, int gene_id_version, std::string gene_symbol, std::string transcript_id, int transcript_id_version) :
    gene_id(gene_id), gene_id_version(gene_id_version), gene_symbol(gene_symbol), transcript_id(transcript_id), transcript_id_version(transcript_id_version) {}

Gene::Gene(std::string gene_id_and_version, std::string gene_symbol, std::string transcript_id_and_version) : gene_symbol(gene_symbol)
{
    auto split_into_id_and_version = [](std::string & to_split, std::string & id, int & version)
    {
        int pos = to_split.find(".");
        if(pos != std::string::npos) {
            id = to_split.substr(0, pos);
            version = atoi(to_split.substr(pos + 1, 1).c_str());                                             
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
    if(!boost::filesystem::exists("gene_id_dictionary.bin")) {
        Gene_id i = 0;
    	io::CSVReader<3, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/sites_with_scored_interactions.tsv");
    	in.read_header(io::ignore_extra_column, "gene_id", "gene_symbol", "transcript_id");
    	std::string column0, column1, column2;
    	while(in.read_row(column0, column1, column2)) {
            Gene gene(column0, column1, column2);
        	Gene::gene_id_dictionary.insert( boost::bimap<Gene, Gene_id>::value_type(gene, i) );
        	i++;
    	}
        std::ofstream out("gene_id_dictionary.bin");
        boost::archive::binary_oarchive oa(out);
        oa << Gene::gene_id_dictionary;
        out.close();
    } else {
        std::ifstream in("gene_id_dictionary.bin");
        boost::archive::binary_iarchive ia(in);
        ia >> Gene::gene_id_dictionary;
        in.close();
    }
}

void Gene::print_gene_dictionary(unsigned int max_rows)
{
    if(max_rows != -1) {
        std::cout << "printing at most " << max_rows << " rows\n";
    }
    int j = 0;
    for(auto & e : Gene::gene_id_dictionary.left) {
        if(j++ < max_rows) {
            Gene & gene = const_cast<Gene &>(e.first);
            Gene_id & i = const_cast<Gene_id &>(e.second);
            std::cout << "gene.gene_id = " << gene.gene_id << ", gene.gene_id_version = " << gene.gene_id_version << ", gene.gene_symbol = " << gene.gene_symbol << ", gene.transcript_id = " << gene.transcript_id << ", gene.transcript_id_version = " << gene.transcript_id_version << ", i = " << i << "\n";
        }
    }
}

template<class Archive>
void Gene::serialize(Archive & ar, const unsigned int version)
{
    ar & this->gene_id & this->gene_id_version & this->gene_symbol & this->transcript_id & this->transcript_id_version;
}

std::ostream & operator<<(std::ostream & stream, const Gene & o)
{
    return stream << "gene_id = " << o.gene_id << ", gene_id_version = " << o.gene_id_version << ", gene_symbol = " << o.gene_symbol << ", transcript_id = " << o.transcript_id << ", transcript_id_version = " << o.transcript_id_version << "\n";
}

bool operator<(Gene const & lhs, Gene const & rhs)
{
    return lhs.gene_id < rhs.gene_id;
}
