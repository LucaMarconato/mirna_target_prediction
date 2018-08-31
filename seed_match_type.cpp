#include "seed_match_type.hpp"

// #include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

Seed_match_type::Seed_match_type()
{
    // TODO: this is needed by Site, which must call the default constructor: maybe I should make this constructor private and make Site a friend class
    this->match_type = Seed_match_type::Mirna_site_matching::undefined;
}

Seed_match_type::Seed_match_type(Seed_match_type::Mirna_site_matching mirna_site_matching)
{
    this->match_type = mirna_site_matching;
}

std::string Seed_match_type::format_mirna_matching(Seed_match_type::Mirna_site_matching mirna_matching)
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

std::ostream & operator<<(std::ostream & stream, const Seed_match_type & o)
{
    return stream << Seed_match_type::format_mirna_matching(o.match_type);
}
