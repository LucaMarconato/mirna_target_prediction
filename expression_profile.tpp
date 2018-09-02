template<class Archive>
void Relative_expression::serialize(Archive & ar, const unsigned int version)
{
    ar & this->value;
    ar & this->undefined;
    ar & this->normalized;
}

template<class Archive>
void Reads::serialize(Archive & ar, const unsigned int version)
{
    ar & this->value;
    ar & this->undefined;
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
    ar & this->to_normalize;
    ar & this->normalization_factor;
}