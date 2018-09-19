template<class Archive>
void Cluster_expression_profile::serialize(Archive & ar, const unsigned int version)
{
    ar & profile;
}