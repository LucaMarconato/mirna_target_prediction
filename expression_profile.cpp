#include "expression_profile.hpp"

#include <cmath>

Expression::Expression() {}

Expression::Expression(Relative_expression relative_expression) : relative_expression(relative_expression) {}

Expression::Expression(RPM rpm)
{
    this->relative_expression = rpm/1000000.0;   
}

RPM Expression::to_RPM()
{
    return std::round(this->relative_expression*1000000.0);
}

Relative_expression Expression::to_relative_expression()
{
    return this->relative_expression;
}

// --------------------------------------------------
