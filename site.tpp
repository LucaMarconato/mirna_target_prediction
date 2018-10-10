template<class Archive>
void Site::serialize(Archive & ar, const unsigned int)
{
    ar & this->mirna_id;
    ar & this->gene_id;
    ar & this->utr_start;
    ar & this->utr_end;
}