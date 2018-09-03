template<class Archive>
void Relative_expression::serialize(Archive & ar, const unsigned int version)
{
    ar & this->value;
    ar & this->undefined;
}

template<class Archive>
void Reads::serialize(Archive & ar, const unsigned int version)
{
    ar & this->value;
    ar & this->undefined;
    ar & this->grand_total;
}

template<class Archive>
void Rpm::serialize(Archive & ar, const unsigned int version)
{
    ar & this->value;
    ar & this->undefined;
}

template<class Archive>
void Expression::serialize(Archive & ar, const unsigned int version)
{
    ar & this->relative_expression;
    ar & this->reads;
    ar & this->rpm;
}