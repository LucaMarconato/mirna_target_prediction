template<class Archive>
void Relative_expression::serialize(Archive & ar, const unsigned int version)
{
    ar & value;
}

template<class Archive>
void Reads::serialize(Archive & ar, const unsigned int version)
{
    ar & value;
}

template<class Archive>
void Rpm::serialize(Archive & ar, const unsigned int version)
{
    ar & value;
}

template<class Archive>
void Expression::serialize(Archive & ar, const unsigned int version)
{
    ar & this->relative_expression;
    ar & to_normalize;
}