template<class Archive>
void Cluster::serialize(Archive & ar, const unsigned int version)
{
    ar & sites;    
}