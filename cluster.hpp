#ifndef CLUSTER_H
#define CLUSTER_H

#include <list>

#include <boost/serialization/list.hpp>

#include "site.hpp"

class Cluster {
public:
    std::list<Site *> sites;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#include "cluster.tpp"

#endif // CLUSTER_H
