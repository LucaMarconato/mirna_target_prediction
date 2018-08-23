#ifndef MIRNA_H
#define MIRNA_H

#include <string>
#include <unordered_map>

class Mirna {
private:
    static bool human_mirnas_initialized;
    static void initialize_human_mirnas();
public:
    enum Mirna_matching {
		canonical_8mer,
		canonical_7mer_m8,
		canonical_7mer_A1,
		canonical_6mer,
		canonical_offset_6mer,
        no_matching
    };
    
    std::string sequence;
    std::string id;    
    static Mirna human_mirnas [3000];
    static std::unordered_map<std::string, int> mirna_id_to_int;

    Mirna();
    /*
      Only the first characters of site are used.
      The string must be at least 7 nt long, otherwise an error will be raise.
    */
    Mirna_matching matches_with_string(char * site);
    static std::string format_mirna_matching(Mirna::Mirna_matching mirna_matching);
};

#endif // MIRNA_H
