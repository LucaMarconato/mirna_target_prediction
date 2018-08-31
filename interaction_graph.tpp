template<class Archive>
void Mirna_site_arc::serialize(Archive & ar, const unsigned int)
{
    ar & this->context_score & this->weighted_context_score & this->conserved;
}