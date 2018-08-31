#ifndef MIRNA_EXPRESSION_PROFILE_H
#define MIRNA_EXPRESSION_PROFILE_H

#include <unordered_map>

#include "expression_profile.hpp"
#include "mirna.hpp"

class Mirna_expression_profile : Expression_profile {
    std::unordered_map<Mirna_id, Expression> profile;
};

#endif // MIRNA_EXPRESSION_PROFILE_H
