template<class Archive>
void Gene_expression_profile::serialize(Archive & ar, const unsigned int version)
{
    ar & this->initialized;
    ar & this->profile;
    ar & this->recognized_distinct_genes;
    ar & this->recognized_reads;
    ar & this->not_recognized_distinct_genes;
    ar & this->not_recognized_reads;
    ar & this->total_distinct_genes;
    ar & this->total_reads;
    ar & this->filtered_out_distinct_genes;
    ar & this->filtered_out_reads;
}