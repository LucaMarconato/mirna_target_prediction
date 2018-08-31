#ifndef EXPRESSION_PROFILE_H
#define EXPRESSION_PROFILE_H

typedef double Relative_expression;
// reads per million
typedef unsigned int RPM;

class Expression {
private:
    Relative_expression relative_expression;    
public:    
    Expression();
    Expression(Relative_expression relative_expression);
    Expression(RPM rpm);
    RPM to_RPM();
    Relative_expression to_relative_expression();
};

class Expression_profile {
    
};

#endif // EXPRESSION_PROFILE_H
