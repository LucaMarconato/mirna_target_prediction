template<class Archive>
void Seed_match_type::serialize(Archive & ar, const unsigned int version)
{
    ar & this->match_type;
}