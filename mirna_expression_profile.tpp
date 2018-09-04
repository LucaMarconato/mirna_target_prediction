template<class Archive>
void Mirna_expression_profile::serialize(Archive & ar, const unsigned int version)
{
    ar & this->initialized;
    ar & this->profile;
    ar & this->recognized_distinct_mirnas;
    ar & this->recognized_reads;
    ar & this->not_recognized_distinct_mirnas;
    ar & this->not_recognized_reads;
    ar & this->total_distinct_mirnas;
    ar & this->total_reads;
    ar & this->filtered_out_distinct_mirnas;
    ar & this->filtered_out_reads;
}