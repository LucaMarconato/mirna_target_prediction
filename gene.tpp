template<class Archive>
void Gene::serialize(Archive & ar, const unsigned int version)
{
    ar & this->gene_id & this->gene_id_version & this->gene_symbol & this->transcript_id & this->transcript_id_version;
}