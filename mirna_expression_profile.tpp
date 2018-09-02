template<class Archive>
void Mirna_expression_profile::serialize(Archive & ar, const unsigned int version)
{
    ar & this->profile;
    ar & this->recognized_distinct_mirnas;
    ar & this->recognized_rpm;
    ar & this->not_recognized_distinct_mirnas;
    ar & this->not_recognized_rpm;
    ar & this->total_distinct_mirnas;
    ar & this->total_rpm;
    ar & this->initialized;
    ar & this->filtered_out_distinct_mirnas;
    ar & this->filtered_out_rpm;
}