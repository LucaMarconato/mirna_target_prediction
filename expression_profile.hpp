#ifndef EXPRESSION_PROFILE_H
#define EXPRESSION_PROFILE_H

struct Relative_expression {
    double value;
    bool undefined = true;

    Relative_expression();
    Relative_expression(const Relative_expression & obj);
    Relative_expression(double value);
    void swap(Relative_expression & obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};
struct Reads {
    double value;
    bool undefined = true;
    double grand_total = 0;
    
    Reads();
    Reads(const Reads & obj);
    Reads(double value);
    void swap(Reads & obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};
// reads per million
struct Rpm {
    double value;
    bool undefined = true;
    
    Rpm();
    Rpm(const Rpm & obj);
    Rpm(double value);
    void swap(Rpm & obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};

class Expression {
private:
    Relative_expression relative_expression;
    Reads reads;
    Rpm rpm;
public:
    bool valid_rpm = false;
    
    Expression();
    Expression(Relative_expression relative_expression);
    Expression(Reads reads);
    Expression(Rpm rpm);
    double to_relative_expression();
    double to_reads();
    double to_rpm();
    void integrity_check();
    Expression & operator=(Expression obj);
    bool is_valid();
    void clear();
    void normalize_reads(double grand_total);
    double get_grand_total();
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
    
    friend void swap(Expression & obj1, Expression & obj2);
};

class Expression_profile {
public:
    bool initialized = false;    
};

#include "expression_profile.tpp"

#endif // EXPRESSION_PROFILE_H
