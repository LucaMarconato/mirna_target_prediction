#include "expression_profile.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

// --------------------------------------------------

Relative_expression::Relative_expression() {}

Relative_expression::Relative_expression(const Relative_expression & obj) {
    this->value = obj.value;
    this->undefined = obj.undefined;
}

Relative_expression::Relative_expression(double value) : value(value) {
    undefined = false;
}

void Relative_expression::swap(Relative_expression & obj)
{
    std::swap(this->value, obj.value);
    std::swap(this->undefined, obj.undefined);
}

// --------------------------------------------------

Reads::Reads() {}

Reads::Reads(const Reads & obj) {
    this->value = obj.value;
    this->undefined = obj.undefined;
    this->grand_total = obj.grand_total;
}

Reads::Reads(double value) : value(value) {
    undefined = false;
}

void Reads::swap(Reads & obj)
{
    std::swap(this->value, obj.value);
    std::swap(this->undefined, obj.undefined);
    std::swap(this->grand_total, obj.grand_total);
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

Expression::Expression(Relative_expression relative_expression) {
    this->relative_expression = Relative_expression(relative_expression);
    this->reads = Reads();
    this->rpm = Rpm();
    this->valid_rpm = false;
}

Expression::Expression(Reads reads)
{
    this->relative_expression = Relative_expression();
    this->reads = Reads(reads);
    this->rpm = Rpm();
    this->valid_rpm = true;
}

Expression::Expression(Rpm rpm)
{
    this->relative_expression = Relative_expression();
    this->reads = Reads();
    this->rpm = Rpm(rpm);
    this->valid_rpm = true;
}

void swap(Expression & obj1, Expression & obj2)
{
    // I prefer the syntax with std::swap, I will make the following conform to it
    obj1.relative_expression.swap(obj2.relative_expression);
    obj1.reads.swap(obj2.reads);
    obj1.rpm.swap(obj2.rpm);
    std::swap(obj1.valid_rpm, obj2.valid_rpm);
}

Expression & Expression::operator=(Expression obj)
{
    swap(*this, obj);
    return *this;
}

double Expression::to_relative_expression()
{
    integrity_check();
    if(!this->relative_expression.undefined) {
        return this->relative_expression.value;
    } else if(!this->reads.undefined) {
        if(this->reads.grand_total == 0) {
            std::cerr << "this->reads is not normalized\n";
            exit(1);
        } else {
            return ((double)this->reads.value)/this->reads.grand_total;
        }
    } else if(!this->rpm.undefined) {
        return ((double)this->rpm.value)/1000000;
    } else {
        // this case is actually caught in integrity_check()
        std::cerr << "Relative_expression not initialized\n";
        exit(1);
    }
}

double Expression::to_reads()
{
    integrity_check();
    if(!this->relative_expression.undefined) {
        std::cerr << "error: unable to infer the number of reads from the relative expression\n";
        exit(1);
    } else if(!this->reads.undefined) {
        return this->reads.value;
    } else if(!this->rpm.undefined) {
        std::cerr << "error: unable to infer the number of reads from the rpm\n";
        exit(1);
    } else {
        // this case is actually caught in integrity_check()
        std::cerr << "Relative_expression not initialized\n";
        exit(1);
    }
}

double Expression::to_rpm()
{
    integrity_check();
    if(!this->valid_rpm) {
        std::cerr << "error: requesting an rpm value which would be artificial\n";
        exit(1);
    }
    if(!this->relative_expression.undefined) {
        return this->relative_expression.value/1000000;            
    } else if(!this->reads.undefined) {
        if(this->reads.grand_total == 0) {
            std::cerr << "this->reads is not normalized\n";
            exit(1);
        } else {
            return ((double)this->reads.value)/this->reads.grand_total*1000000;
        }
    } else if(!this->rpm.undefined) {
        return this->rpm.value;
    } else {
        // this case is actually caught in integrity_check()
        std::cerr << "Relative_expression not initialized\n";
        exit(1);
    }
}

void Expression::integrity_check()
{
    int true_count = !this->relative_expression.undefined + !this->reads.undefined + !this->rpm.undefined;   
    if(true_count != 1) {
        std::cerr << "error: this->relative_expression.undefined = " << this->relative_expression.undefined << ", this->reads.undefined = " << this->reads.undefined << ", this->rpm.undefined = " << this->rpm.undefined << "\n";
        exit(1);
    }    
}

void Expression::normalize_reads(double grand_total)
{
    integrity_check();
    if(!this->reads.undefined) {
        this->reads.grand_total = grand_total;
    } else {
        std::cerr << "error: normalized_reads, reads not initialized; this->relative_expression.undefined = " << this->relative_expression.undefined << ", this->rpm.undefined = " << this->rpm.undefined << "\n";
        exit(1);
    }
}

double Expression::get_grand_total()
{
    integrity_check();
    if(!this->reads.undefined) {
        return this->reads.grand_total;
    } else {
        std::cerr << "error: normalized_reads, reads not initialized; this->relative_expression.undefined = " << this->relative_expression.undefined << ", this->rpm.undefined = " << this->rpm.undefined << "\n";
        exit(1);
    }
}
