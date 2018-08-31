#ifndef SEED_MATCH_TYPE_H
#define SEED_MATCH_TYPE_H

#include <iostream>
#include <string>

// TODO: I should change the names in this class to make them more intuitive

class Seed_match_type {
public:
    enum Mirna_site_matching {
		canonical_8mer,
		canonical_7mer_m8,
		canonical_7mer_A1,
		canonical_6mer,
		canonical_offset_6mer,
		non_canonical,
		no_matching,
        undefined
    };
    
    Seed_match_type::Mirna_site_matching match_type;
    // int match_type;

    Seed_match_type();
    Seed_match_type(Seed_match_type::Mirna_site_matching mirna_site_matching);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    static std::string format_mirna_matching(Seed_match_type::Mirna_site_matching mirna_matching);
    friend std::ostream & operator<<(std::ostream & stream, const Seed_match_type & o);
};

#endif // SEED_MATCH_TYPE_H
