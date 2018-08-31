template<class Archive>
void Mirna::serialize(Archive & ar, const unsigned int version)
{
    ar & this->mirna_family;
}