#include "mirna.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

#include "timer.hpp"

boost::bimap<Mirna, Mirna_id> Mirna::mirna_id_dictionary;

Mirna::Mirna() {}

Mirna::Mirna(std::string mirna_family) : mirna_family(mirna_family) {}

void Mirna::initialize_mirna_dictionary()
{
    if(!boost::filesystem::exists("mirna_id_dictionary.bin")) {
        std::ifstream in("./data/processed/mirnas_with_scored_interactions.tsv");
        if(!in.is_open()) {
            std::cerr << "error: unable to open file\n";
            exit(1);
        }
        std::string line;
        Mirna_id i = 0;
        // get the header
        getline(in, line);
        if(line != "mirna_family") {
            std::cerr << "error: line = " << line << "\n";
            exit(1);
        }
        while(!in.eof()) {
            getline(in, line);
            // this addresses the case in which the file ends with a newline
            if(line != "") {
                std::transform(line.begin(), line.end(), line.begin(), ::tolower);
                Mirna mirna(line);
                if(Mirna::mirna_id_dictionary.left.find(mirna) == Mirna::mirna_id_dictionary.left.end()) {
                    Mirna::mirna_id_dictionary.insert( boost::bimap<Mirna, Mirna_id>::value_type(mirna, i) );
                    i++;
                }
            }
        }
        in.close();        
        std::cout << "writing mirna_id_dictionary.bin\n";
        Timer::start();
        std::ofstream out("mirna_id_dictionary.bin", std::ios::binary);
        boost::archive::binary_oarchive oa(out);
        oa << Mirna::mirna_id_dictionary;
        out.close();
        std::cout << "written, ";
        Timer::stop();
        
        std::stringstream ss;
        ss << "mirna_family\tmirna_id_cpp\n";
        for(auto & e : Mirna::mirna_id_dictionary.left) {
            ss << e.first.mirna_family << "\t" << e.second << "\n";
        }
        out.open("data/processed/mirna_id_dictionary.tsv");
        out << ss.str();
        out.close();
    } else {
        std::cout << "loading mirna_id_dictionary.bin\n";
        Timer::start();
        std::ifstream in("mirna_id_dictionary.bin", std::ios::binary);
        boost::archive::binary_iarchive ia(in);        
        ia >> Mirna::mirna_id_dictionary;
        in.close();
        std::cout << "loaded, ";
        Timer::stop();
    }
}

void Mirna::print_mirna_dictionary(unsigned int max_rows)
{
    if(max_rows != -1) {
        std::cout << "printing at most " << max_rows << " rows\n";
    }
    int j = 0;
    for(auto & e : Mirna::mirna_id_dictionary.left) {
        if(j++ < max_rows) {
            Mirna & mirna = const_cast<Mirna &>(e.first);
            Mirna_id & i = const_cast<Mirna_id &>(e.second);
            std::cout << "mirna.mirna_family = " << mirna.mirna_family << ", i = " << i << "\n";            
        }
    }
}

bool operator<(Mirna const & lhs, Mirna const & rhs)
{
    return lhs.mirna_family < rhs.mirna_family;
}

std::ostream & operator<<(std::ostream & stream, const Mirna & o)
{
    return stream << "mirna_family = " << o.mirna_family << "\n";
}

// Mirna::Mirna_matching Mirna::matches_with_string(char * site)
// {
//     int len = strlen(site);
//     if(len < 8) {
//         std::cerr << "error: len = " << len << "\n";
//         exit(1);        
//     }
//     char * mirna = const_cast<char *>(this->sequence.c_str());
//     bool A_at_pos_1 = site[0] == 'A';
//     int mismatch = rna_rna_mismatch(mirna+1, site+1, 7, true);
//     if(A_at_pos_1 && mismatch == -1) {
//         return canonical_8mer;
//     }
//     if(mismatch == -1) {
//         return canonical_7mer_m8;
//     }
//     if(A_at_pos_1 && mismatch == 6) {
//         return canonical_7mer_A1;
//     }
//     if(mismatch == 6) {
//         return canonical_6mer;
//     }
//     if(rna_rna_mismatch(mirna+2, site+2, 6, true) == -1) {
//         return canonical_offset_6mer;
//     }       
//     return no_matching;
// }
