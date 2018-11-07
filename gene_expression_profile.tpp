template<class Archive>
void Gene_expression_profile::serialize(Archive & ar, const unsigned int)
{
    ar & this->initialized;
    ar & this->profile;
    ar & this->distinct_genes;
    ar & this->total_reads;
    ar & this->filtered_out_distinct_genes;
    ar & this->filtered_out_reads;
}