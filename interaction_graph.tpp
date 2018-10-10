template<class Archive>
void Mirna_site_arc::serialize(Archive & ar, const unsigned int)
{
    ar & this->context_score & this->weighted_context_score & this->conserved;
}

template<class Archive>
void Interaction_graph::serialize(Archive & ar, const unsigned int)
{
    // TODO: if you want to use this you have to update this list
    ar & this->sites_by_location;
    ar & this->gene_to_sites_arcs;
    ar & this->mirna_site_arcs;
    ar & this->mirna_to_sites_arcs;
    ar & this->site_to_mirnas_arcs;
    ar & this->mirna_gene_arcs;
    ar & this->mirna_to_genes_arcs;
    ar & this->gene_to_mirnas_arcs;
    ar & this->rows_processed;
    ar & this->rows_skipped;
}