#ifndef EXPRESSION_PROFILE_H
#define EXPRESSION_PROFILE_H

struct Relative_expression {
    double value;
    
    Relative_expression() {}
    Relative_expression(const Relative_expression & obj) {
        this->value = obj.value;
    }
    Relative_expression(double value) : value(value) {}
    void swap(Relative_expression & obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};
struct Reads {
    unsigned long long value;
    
    Reads() {}
    Reads(Reads & obj) {
        this->value = obj.value;
    }
    Reads(unsigned long long value) : value(value) {}
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};
// reads per million
struct Rpm {
    double value;
    
    Rpm() {}
    Rpm(Reads & obj) {
        this->value = obj.value;
    }
    Rpm(double value) : value(value) {}
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

class Expression {
private:
    Relative_expression relative_expression;
    bool to_normalize = false;    
public:    
    Expression();
    Expression(Relative_expression relative_expression);
    Expression(Reads reads);
    Expression(Rpm rpm);
    Rpm to_rpm();
    Relative_expression to_relative_expression();
    void normalize_reads(double total_reads);
    Expression & operator=(Expression obj);
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    
    friend void swap(Expression & obj1, Expression & obj2);
};

class Expression_profile {
public:
    bool initialized = false;
};

#include "expression_profile.tpp"

#endif // EXPRESSION_PROFILE_H
