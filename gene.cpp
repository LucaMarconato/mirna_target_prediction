#include "gene.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

boost::bimap<Gene, int> Gene::gene_id_dictionary;

Gene::Gene(std::string gene_id, int gene_id_version, std::string gene_symbol, std::string transcript_id, int transcript_id_version)
{
	this->gene_id = gene_id;
    this->gene_id_version = gene_id_version;
    this->gene_symbol = gene_symbol;
    this->transcript_id = transcript_id;
    this->transcript_id_version = transcript_id_version;
}

void Gene::initialize_gene_dictionary()
{
    std::ifstream in("./data/processed/sites_with_scored_interactions.tsv");
    if(!in.is_open()) {
        std::cerr << "error: unable to open file\n";
        exit(1);
    }
    std::string line;
    int i = 0;
    while(!in.eof()) {
        getline(in, line);
        // this should happen at most in the end of the file
        if(line == "") {
            continue;
        }
        std::string gene_id;
        int gene_id_version;
        std::string gene_symbol;
        std::string transcript_id;
        int transcript_id_version;
        int field = 0;
        int pos = 0;
        int end_of_last_field = 0;
        bool end_of_field = false;
        while(field < 5) {
            if(end_of_field) {
                std::string temp = line.substr(end_of_last_field, pos-end_of_last_field);
                switch(field) {
                case 0:
                    gene_id = temp;
                    break;
                case 1:
                    gene_id_version = atoi(temp.c_str());
                    break;
                case 2:
                    gene_symbol = temp;
                    break;
                case 3:
                    transcript_id = temp;
                    break;
                case 4:
                    transcript_id_version = atoi(temp.c_str());
                    break;
                default:
                    std::cerr << "error: field = " << field << "\n";
                    exit(1);
                    break;
                }
                pos++;
                end_of_last_field = pos;
                end_of_field = false;
                field++;
            } else {
                if(line[pos] == '\t') {
                    end_of_field = true;
                } else {
                    pos++; 
                }
            }
        }
        // std::stringstream ss(line);
        // ss >> gene_id >> gene_id_version >> gene_symbol >> transcript_id >> transcript_id_version;
        Gene gene(gene_id, gene_id_version, gene_symbol, transcript_id, transcript_id_version);
        Gene::gene_id_dictionary.insert( boost::bimap<Gene, int>::value_type(gene, i) );
        i++;
    }
    in.close();
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
            int & i = const_cast<int &>(e.second);
            std::cout << "gene.gene_id = " << gene.gene_id << ", gene.gene_id_version = " << gene.gene_id_version << ", gene.gene_symbol = " << gene.gene_symbol << ", gene.transcript_id = " << gene.transcript_id << ", gene.transcript_id_version = " << gene.transcript_id_version << ", i = " << i << "\n";
        }
    }
}

bool operator<(Gene const & lhs, Gene const & rhs)
{
    return lhs.gene_id < rhs.gene_id;
}
