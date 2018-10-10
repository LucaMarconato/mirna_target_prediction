#ifndef CLUSTER_H
#define CLUSTER_H

#include <list>

#include <boost/functional/hash.hpp>
#include <boost/serialization/list.hpp>

#include "site.hpp"

class Cluster {
public:
    std::list<Site *> sites;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};

namespace std {
    template <>
    struct hash<std::pair<Mirna_id, Cluster *>>
    {
        size_t operator()(const std::pair<Mirna_id, Cluster *> & e) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, boost::hash_value(e.first));
            boost::hash_combine(seed, boost::hash_value(e.second));
            return seed;
        }
    };
}

// namespace std {
//     template <>
//     struct hash<std::pair<, Cluster *>>
//     {
//         size_t operator()(const std::pair<Mirna_id, Cluster *> & e) const noexcept
//         {
//             std::size_t seed = 0;
//             boost::hash_combine(seed, boost::hash_value(e.first));
//             boost::hash_combine(seed, boost::hash_value(e.second));
//             return seed;
//         }
//     };
// }

#include "cluster.tpp"

#endif // CLUSTER_H
