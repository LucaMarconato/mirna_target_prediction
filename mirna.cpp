#include "mirna.hpp"

#include <iostream>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

boost::bimap<Mirna, int> Mirna::mirna_id_dictionary;

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
        int i = 0;
        while(!in.eof()) {
            getline(in, line);
            Mirna mirna(line);
            Mirna::mirna_id_dictionary.insert( boost::bimap<Mirna, int>::value_type(mirna, i) );
            i++;
        }
        in.close();
        std::ofstream out("mirna_id_dictionary.bin", std::ios::binary);
        boost::archive::binary_oarchive oa(out);
        oa << Mirna::mirna_id_dictionary;
        out.close();
    } else {
        std::ifstream in("mirna_id_dictionary.bin", std::ios::binary);
        boost::archive::binary_iarchive ia(in);        
        ia >> Mirna::mirna_id_dictionary;
        in.close();
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
            int & i = const_cast<int &>(e.second);
            std::cout << "mirna.mirna_family = " << mirna.mirna_family << ", i = " << i << "\n";            
        }
    }
}

std::string Mirna::format_mirna_matching(Mirna::Mirna_matching mirna_matching)
{    
    switch(mirna_matching) {
    case canonical_8mer:
        return "canonical_8mer";
    case canonical_7mer_m8:
        return "canonical_7mer_m8";
    case canonical_7mer_A1:
        return "canonical_7mer_A1";
    case canonical_6mer:
        return "canonical_6mer";
    case canonical_offset_6mer:
        return "canonical_offset_6mer";
    case non_canonical:
        return "non_canonical";
    case no_matching:
        return "no_matching";
    default:
        std::cerr << "error: mirna_matching = " << mirna_matching << "\n";
        exit(1);
    }
}

template<class Archive>
void Mirna::serialize(Archive & ar, const unsigned int)
{
    ar & this->mirna_family;
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
