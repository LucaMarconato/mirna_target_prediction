template<class Archive>
void Site_expression_profile::serialize(Archive & ar, const unsigned int version)
{
    ar & profile;
}