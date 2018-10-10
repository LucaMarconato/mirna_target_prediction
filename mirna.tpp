template<class Archive>
void Mirna::serialize(Archive & ar, const unsigned int)
{
    ar & this->mirna_family;
}