#ifndef CLUSTER_H
#define CLUSTER_H

#include <list>

#include "site.hpp"

class Cluster {
public:
    std::list<Site *> sites;
};

#endif // CLUSTER_H
