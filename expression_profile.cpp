#include "expression_profile.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

// --------------------------------------------------

Relative_expression::Relative_expression() {}

Relative_expression::Relative_expression(const Relative_expression & obj) {
    this->value = obj.value;
    this->undefined = obj.undefined;
    this->normalized = obj.normalized;
}

Relative_expression::Relative_expression(double value) : value(value) {
    undefined = false;
}

void Relative_expression::swap(Relative_expression & obj)
{
    std::swap(this->value, obj.value);
    std::swap(this->undefined, obj.undefined);
    std::swap(this->normalized, obj.normalized);
}

void Relative_expression::normalize(double normalization_factor)
{
    this->value /= normalization_factor;
    this->normalized = true;
}

// --------------------------------------------------

Reads::Reads() {}

Reads::Reads(const Reads & obj) {
    this->value = obj.value;
    this->undefined = obj.undefined;
}

Reads::Reads(unsigned long long value) : value(value) {
    undefined = false;
}

void Reads::swap(Reads & obj)
{
    std::swap(this->value, obj.value);
    std::swap(this->undefined, obj.undefined);
}

// --------------------------------------------------

Rpm::Rpm() {}

Rpm::Rpm(const Rpm & obj) {
    this->value = obj.value;
    this->undefined = obj.undefined;
}

Rpm::Rpm(double value) : value(value) {
    undefined = false;
}

void Rpm::swap(Rpm & obj)
{
    std::swap(this->value, obj.value);
    std::swap(this->undefined, obj.undefined);
}

// --------------------------------------------------

Expression::Expression() {}

Expression::Expression(Relative_expression relative_expression) : relative_expression(relative_expression) {}

void swap(Expression & obj1, Expression & obj2)
{
    obj1.relative_expression.swap(obj2.relative_expression);
    // I would prefer the following syntax, which requires some extra code
    // std::swap(obj1.relative_expression, obj2.relative_expression);
    std::swap(obj1.to_normalize, obj2.to_normalize);
}

Expression & Expression::operator=(Expression obj)
{
    swap(*this, obj);
    return *this;
}

Expression::Expression(Reads reads)
{
    this->relative_expression.value = (double)reads.value;
    this->to_normalize = true;
}

Expression::Expression(Rpm rpm)
{
    this->relative_expression.value = rpm.value/1000000.0;   
}

Rpm Expression::to_rpm()
{
    if(this->to_normalize) {
        std::cerr << "error: to_rpm(), this->to_normalize = " << this->to_normalize << "\n";
        exit(1);
    }
    return std::round(this->relative_expression.value*1000000.0);
}

Relative_expression Expression::to_relative_expression()
{
    if(this->to_normalize) {
        std::cerr << "error: to_relative_expression(), this->to_normalize = " << this->to_normalize << "\n";
        exit(1);
    }
    return this->relative_expression;
}

void Expression::normalize_reads(double total_reads)
{
    if(!this->to_normalize) {
        std::cerr << "error: normalize_reads(" << total_reads << "), this->to_normalize = " << this->to_normalize << "\n";
        exit(1);
    }
    this->relative_expression.value /= (double)total_reads;
}

double Expression::get_normalization_factor()
{
    if(!this->to_normalize) {
        std::cerr << "error: normalize_reads(" << total_reads << "), this->to_normalize = " << this->to_normalize << "\n";
        exit(1);
    }
    return this->normalization_factor;
}
