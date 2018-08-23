#include "mirna.hpp"

#include <iostream>
#include <fstream>
#include <cstdio>

#include "matching.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

bool Mirna::human_mirnas_initialized = false;

void Mirna::initialize_human_mirnas()
{
    std::ifstream in("human_mirnas.json");
    std::filebuf* pbuf = in.rdbuf();
    std::size_t size = pbuf->pubseekoff (0, in.end, in.in);
    pbuf->pubseekpos (0,in.in);
    char * buffer=new char[size+1];
    pbuf->sgetn(buffer,size);
    buffer[size] = '\0';
    in.close();
    auto j = json::parse(buffer);
    delete [] buffer;

    int mirna_added = 0;
    for(auto && j_mirna : j["list"]) {
        std::string mirna_id = j_mirna["id"].get<std::string>();
        std::string mirna_sequence = j_mirna["s"].get<std::string>();
        Mirna::human_mirnas[mirna_added].id = mirna_id;
        Mirna::human_mirnas[mirna_added].sequence = mirna_sequence;
        Mirna::mirna_id_to_int[mirna_id] = mirna_added;
        mirna_added++;
    }
    // TODO, initialize also the mirna_id_to_int variable
}

Mirna::Mirna()
{
    if(!Mirna::human_mirnas_initialized) {
        Mirna::human_mirnas_initialized = true;
        Mirna::initialize_human_mirnas();
    }    
}

Mirna Mirna::human_mirnas [3000];

std::unordered_map<std::string, int> Mirna::mirna_id_to_int;

Mirna::Mirna_matching Mirna::matches_with_string(char * site)
{
    int len = strlen(site);
    if(len < 8) {
        std::cerr << "error: len = " << len << "\n";
        exit(1);        
    }
    char * mirna = const_cast<char *>(this->sequence.c_str());
    bool A_at_pos_1 = site[0] == 'A';
    int mismatch = rna_rna_mismatch(mirna+1, site+1, 7, true);
    if(A_at_pos_1 && mismatch == -1) {
        return canonical_8mer;
    }
    if(mismatch == -1) {
        return canonical_7mer_m8;
    }
    if(A_at_pos_1 && mismatch == 6) {
        return canonical_7mer_A1;
    }
    if(mismatch == 6) {
        return canonical_6mer;
    }
    if(rna_rna_mismatch(mirna+2, site+2, 6, true) == -1) {
        return canonical_offset_6mer;
    }       
    return no_matching;
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
    case no_matching:
        return "no_matching";
    default:
        std::cerr << "error: mirna_matching = " << mirna_matching << "\n";
        exit(1);
    }
}
