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
    std::string mirbase_id;

    Mirna(std::string mirbase_id);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
    static void initialize_mirna_dictionary();
    static void print_mirna_dictionary(unsigned int max_rows = -1);

    friend bool operator<(Mirna const & lhs, Mirna const & rhs);
    friend std::ostream & operator<<(std::ostream & stream, const Mirna & o);
    /* // for accessing the private default constructor */
    template< class LeftType, class RightType, bool force_mutable >
    friend class boost::bimaps::relation::detail::relation_storage;
};

namespace std {
    template<>
    struct hash<Mirna>
    {
        size_t operator()(const Mirna & e) const noexcept
        {
            return boost::hash<std::string>()(e.mirbase_id);
        }
    };
}

#include "mirna.tpp"

#endif // MIRNA_H
