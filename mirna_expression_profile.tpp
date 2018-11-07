template<class Archive>
void Mirna_expression_profile::serialize(Archive & ar, const unsigned int)
{
    ar & this->initialized;
    ar & this->profile;
    ar & this->distinct_mirnas;
    ar & this->total_reads;
    ar & this->filtered_out_distinct_mirnas;
    ar & this->filtered_out_reads;
}