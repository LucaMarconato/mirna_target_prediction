template<class Archive>
void Cluster::serialize(Archive & ar, const unsigned int)
{
    ar & sites;    
}