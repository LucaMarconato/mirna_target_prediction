#ifndef MIRNA_H
#define MIRNA_H

#include <iostream>
#include <string>

// #include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/functional/hash.hpp>

typedef unsigned int Mirna_id;

class Mirna {    
private:
    Mirna();
public:
    static boost::bimap<Mirna, Mirna_id> mirna_id_dictionary;
    std::string mirna_family;

    Mirna(std::string mirna_family);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    static void initialize_mirna_dictionary();
    static void print_mirna_dictionary(unsigned int max_rows = -1);
    
    friend bool operator<(Mirna const & lhs, Mirna const & rhs);
    friend std::ostream & operator<<(std::ostream & stream, const Mirna & o);
    /* // for accessing the private default constructor */
    template< class LeftType, class RightType, bool force_mutable >
    friend class boost::bimaps::relation::detail::relation_storage;

    // std::string sequence;
    // std::string id;
    // static Mirna human_mirnas [3000];
    // static std::unordered_map<std::string, int> mirna_id_to_int;

    // Mirna();
    // /\*
    //   Only the first characters of site are used.
    //   The string must be at least 7 nt long, otherwise an error will be raise.
    // *\/
    // Mirna_matching matches_with_string(char * site);
    // static std::string format_mirna_matching(Mirna::Mirna_matching mirna_matching);
};

namespace std {
    template<>
    struct hash<Mirna>
    {
        size_t operator()(const Mirna & e) const noexcept
        {
            return boost::hash<std::string>()(e.mirna_family);
        }
    };
}

#endif // MIRNA_H
