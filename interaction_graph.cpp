#include "interaction_graph.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <set>
#include <queue>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

#include <strasser/csv.h>
#include <marconato/output_buffer/output_buffer.hpp>

#include "seed_match_type.hpp"
#include "timer.hpp"
#include "global_parameters.hpp"

#define Ig Interaction_graph

Mirna_site_arc::Mirna_site_arc(Seed_match_type seed_match_type, double context_score, double weighted_context_score, bool conserved) : seed_match_type(seed_match_type), context_score(context_score), weighted_context_score(weighted_context_score), conserved(conserved) {}

Mirna_site_arc::Mirna_site_arc(const Mirna_site_arc & obj)
{
    this->seed_match_type = obj.seed_match_type;
    this->context_score = obj.context_score;
    this->weighted_context_score = obj.weighted_context_score;
    this->conserved = obj.conserved;
}

void swap(Mirna_site_arc & obj1, Mirna_site_arc & obj2)
{
    std::swap(obj1.seed_match_type,obj2.seed_match_type);
    std::swap(obj1.context_score,obj2.context_score);
    std::swap(obj1.weighted_context_score,obj2.weighted_context_score);
    std::swap(obj1.conserved,obj2.conserved);    
}

Mirna_site_arc & Mirna_site_arc::operator=(Mirna_site_arc obj)
{
    swap(*this, obj);
    return *this;
}

// --------------------------------------------------

Ig::Interaction_graph() {}

Ig::Ig(const Ig & obj)
{
    this->sites_by_location = obj.sites_by_location;
    this->gene_to_sites_arcs = obj.gene_to_sites_arcs;
    this->gene_to_clusters_arcs = obj.gene_to_clusters_arcs;
    this->mirna_site_arcs = obj.mirna_site_arcs;
    this->mirna_to_sites_arcs = obj.mirna_to_sites_arcs;
    this->site_to_mirnas_arcs = obj.site_to_mirnas_arcs;
    this->mirna_gene_arcs = obj.mirna_gene_arcs;
    this->mirna_to_genes_arcs = obj.mirna_to_genes_arcs;
    this->gene_to_mirnas_arcs = obj.gene_to_mirnas_arcs;
    this->site_to_overlapping_sites = obj.site_to_overlapping_sites;
    this->rows_processed = obj.rows_processed;
    this->rows_skipped = obj.rows_skipped;
}

void swap(Ig & obj1, Ig & obj2)
{
    std::swap(obj1.sites_by_location, obj2.sites_by_location);
    std::swap(obj1.gene_to_sites_arcs, obj2.gene_to_sites_arcs);
    std::swap(obj1.gene_to_clusters_arcs, obj2.gene_to_clusters_arcs);
    std::swap(obj1.mirna_site_arcs, obj2.mirna_site_arcs);
    std::swap(obj1.mirna_to_sites_arcs, obj2.mirna_to_sites_arcs);
    std::swap(obj1.site_to_mirnas_arcs, obj2.site_to_mirnas_arcs);
    std::swap(obj1.mirna_gene_arcs, obj2.mirna_gene_arcs);
    std::swap(obj1.mirna_to_genes_arcs, obj2.mirna_to_genes_arcs);
    std::swap(obj1.gene_to_mirnas_arcs, obj2.gene_to_mirnas_arcs);
    std::swap(obj1.site_to_overlapping_sites, obj2.site_to_overlapping_sites);
    std::swap(obj1.rows_processed, obj2.rows_processed);
    std::swap(obj1.rows_skipped, obj2.rows_skipped);
}

Ig & Ig::operator=(Ig obj)
{
    swap(*this, obj);
    return *this;
}

void Ig::build_interaction_graph(std::set<Mirna_id> & mirnas, std::set<Gene_id> & genes)
{
    std::cout << "building the interaction graph\n";
    Timer::start();    
    io::CSVReader<8, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in("./data/processed/scored_interactions_processed.tsv");
    in.read_header(io::ignore_extra_column, "mirna_id",  "gene_id",  "utr_start",  "utr_end",  "seed_match_type",  "context_score",  "weighted_context_score", "conserved");
    std::string columns[8];
    while(in.read_row(columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], columns[6], columns[7])) {
        Mirna_id mirna_id = std::strtoul(columns[0].c_str(), nullptr, 10);
        Gene_id gene_id = std::strtoul(columns[1].c_str(), nullptr, 10);
        if(mirnas.find(mirna_id) == mirnas.end() || genes.find(gene_id) == genes.end()) {
            this->rows_skipped++;
            continue;
        }
        this->rows_processed++;            
        unsigned int utr_start = std::strtoul(columns[2].c_str(), nullptr, 10);
        unsigned int utr_end = std::strtoul(columns[3].c_str(), nullptr, 10);
        Seed_match_type::Mirna_site_matching mirna_site_matching = static_cast<Seed_match_type::Mirna_site_matching>(atoi(columns[4].c_str()));
        Seed_match_type seed_match_type(mirna_site_matching);
        double context_score = std::strtod(columns[5].c_str(), nullptr);
        double weighted_context_score = std::strtod(columns[6].c_str(), nullptr);
        bool conserved = atoi(columns[7].c_str());
        Site * site = this->get_site(mirna_id, gene_id, utr_start, utr_end);

        // create gene-site arcs
        this->gene_to_sites_arcs[gene_id].push_back(site);        

        // create mirna-site arcs
        Mirna_site_arc mirna_site_arc(seed_match_type, context_score, weighted_context_score, conserved);
        auto p0 = std::make_pair(mirna_id, site);
        if(this->mirna_site_arcs.find(p0) != this->mirna_site_arcs.end()) {
            std::cout << "warning: considering only one binding between a mirna and one sites\n";
        }
        this->mirna_site_arcs[p0] = mirna_site_arc;        
        this->mirna_to_sites_arcs[mirna_id].push_back(site);        
        this->site_to_mirnas_arcs[site].push_back(mirna_id);
        
        // create mirna-gene arcs
        auto p1 = std::make_pair(mirna_id, gene_id);
        this->mirna_gene_arcs[p1].push_back(site);        
        this->mirna_to_genes_arcs[mirna_id].push_back(gene_id);        
        this->gene_to_mirnas_arcs[gene_id].push_back(mirna_id);
    }
    // create site-site arcs
    for(auto & e : this->gene_to_sites_arcs) {
        // auto & gene_id = e.first;
        auto & sites = e.second;
        for(auto & site0 : sites) {
            for(auto & site1 : sites) {
                if(site0 != site1) {
                    // WARNING: you may be interested in using utr_start and utr_end instead of the following condition
                    if(std::abs((((long long)site0->utr_start) - ((long long)site1->utr_start))) <= Global_parameters::threshold_for_overlapping_sites) {
                        this->site_to_overlapping_sites[site0].push_back(site1);
                    }                    
                }
            }
        }
    }

    // extract clusters
    for(auto & e : this->gene_to_sites_arcs) {
        auto & gene_id = e.first;
        std::set<Site *> processed;
        for(auto & site : e.second) {
            if(processed.find(site) == processed.end()) {
                Cluster * c = new Cluster();
                std::queue<Site *> to_explore;
                to_explore.push(site);                
                while(to_explore.size() > 0) {
                    // std::cout << "processed.size() = " << processed.size() << ", c.sites.size() = " << c.sites.size() << ", to_explore.size() = " << to_explore.size() << "\n";
                    auto & e = to_explore.front();
                    to_explore.pop();
                    if(processed.find(e) != processed.end()) {
                        continue;
                    }
                    processed.insert(e);
                    c->sites.push_back(e);
                    if(this->site_to_overlapping_sites.find(e) != this->site_to_overlapping_sites.end()) {
                        for(auto & overlapping : this->site_to_overlapping_sites.at(e)) {
                            if(processed.find(overlapping) == processed.end()) {
                                to_explore.push(overlapping);            
                            }
                        }   
                    }
                }
                this->gene_to_clusters_arcs[gene_id].push_back(c);
            }
        }
    }
    
    /*
      WARNING: the model has a limitation (which can be addressed if needed): it underestimate the effect of a miRNA in the case in which it has more than one site available withing the same cluster. This scenario happens in approximately 0.5% of the cases, has it can be deduced running the following command 
      cut -f 6 -d " " | uniq | sort | uniq | wc
      to the output of the code commented below
      // for(auto & e : this->gene_to_clusters_arcs) {
      //     for(Cluster * cluster : e.second) {
      //         std::unordered_map<Mirna_id, int> count;
      //         for(Site * site : cluster->sites) {
      //             for(auto & mirna : this->site_to_mirnas_arcs.at(site)) {
      //                 if(count.find(mirna) == count.end()) {
      //                     count[mirna] = 0;
      //                 }
      //                 count[mirna] = count.at(mirna) + 1;
      //             }
      //         }
      //         for(auto e : count) {
      //             if(e.second > 1) {
      //                 std::cout << "e.first = " << e.first << ", cluster = " << cluster << "\n";
      //             }
      //         }
      //     }
      // }
     */
    
    std::cout << "built, ";
    Timer::stop();
}

Site * Ig::get_site(Mirna_id mirna_id, Gene_id gene_id, unsigned int utr_start, unsigned int utr_end)
{
    auto t = boost::make_tuple(mirna_id, gene_id, utr_start, utr_end);
    auto e = this->sites_by_location.find(t);
    if(e != this->sites_by_location.end()) {
        return e->second;
    } else {
        Site * site = new Site(mirna_id, gene_id, utr_start, utr_end);
        this->sites_by_location[boost::make_tuple(mirna_id, gene_id, utr_start, utr_end)] = site;
        return site;   
    }
}

void Ig::print_statistics()
{
    std::cout << "rows_processed/rows_skipped = " << rows_processed << "/" << rows_skipped << " = " << ((double)rows_processed)/rows_skipped << "\n";
    std::cout << "sites_by_location.size() = " << sites_by_location.size() << "\n";
    std::cout << "arcs: gene-site\n";
    std::cout << "gene_to_sites_arcs.size() = " << gene_to_sites_arcs.size() << "\n";
    std::cout << "arcs: mirna-site\n";
    std::cout << "mirna_site_arcs.size() = " << mirna_site_arcs.size() << "\n";
    std::cout << "mirna_to_sites_arcs.size() = " << mirna_to_sites_arcs.size() << "\n";
    std::cout << "site_to_mirnas_arcs.size() = " << site_to_mirnas_arcs.size() << "\n";
    std::cout << "arcs: mirna-gene\n";
    std::cout << "mirna_gene_arcs.size() = " << mirna_gene_arcs.size() << "\n";
    std::cout << "mirna_to_genes_arcs.size() = " << mirna_to_genes_arcs.size() << "\n";
    std::cout << "gene_to_mirnas_arcs.size() = " << gene_to_mirnas_arcs.size() << "\n";
    std::cout << "site_to_overlapping_sites.size() = " << site_to_overlapping_sites.size() << "\n";

    if(this->gene_to_sites_arcs.size() != this->gene_to_mirnas_arcs.size()) {
        std::cerr << "error: gene_to_sites_arcs.size() = " << gene_to_sites_arcs.size() << ", gene_to_mirnas_arcs.size() = " << gene_to_mirnas_arcs.size() << "\n";
        exit(1);
    }
    if(this->mirna_to_sites_arcs.size() != this->mirna_to_genes_arcs.size()) {
        std::cerr << "error: mirna_to_sites_arcs.size() = " << mirna_to_sites_arcs.size() << ", mirna_to_genes_arcs.size() = " << mirna_to_genes_arcs.size() << "\n";
        exit(1);
    }
    // note that if the condition is false (so that the equality holds), we are not considering sites in which the same mirna can bind in different ways
    if(this->sites_by_location.size() != this->mirna_site_arcs.size()) {
        std::cerr << "error: sites_by_location.size() = " << sites_by_location.size() << ", mirna_site_arcs.size() = " << mirna_site_arcs.size() << "\n";
        exit(1);        
    }    
    unsigned long long total_clusters = 0;
    for(auto & e : this->gene_to_clusters_arcs) {
        total_clusters += e.second.size();
    }
    std::cout << "total_clusters/total_sites: " << total_clusters << "/" << mirna_site_arcs.size() << " = " << ((double)total_clusters)/mirna_site_arcs.size() << "\n";
    
}

void Ig::export_interactions_data(std::string patient_folder)
{
    // matrix with one row for each mirna, one column for each gene and with the cells showing the number of sites for the mirna-gene combination
    std::string filename = patient_folder + "mirna_gene_adjacency_matrix.mat";
    unsigned long long mirna_count = this->mirna_to_genes_arcs.size();
    unsigned long long gene_count = this->gene_to_mirnas_arcs.size();
    // Timer::start();
    std::cout << "writing an " << mirna_count << "x" << gene_count << " sparse matix to " << filename << "\n";
    unsigned long long ** m = new unsigned long long * [mirna_count];
    for(unsigned long long i = 0; i < mirna_count; i++) {
        m[i] = new unsigned long long [gene_count];
        for(unsigned long long j = 0; j < gene_count; j++) {
            m[i][j] = 0;
        }
    }
    boost::bimap<unsigned long long, Mirna_id> i_to_mirna_id;
    boost::bimap<unsigned long long, Gene_id> j_to_gene_id;
    for(auto & e : this->mirna_to_genes_arcs) {
        Mirna_id mirna_id = e.first;
        if(i_to_mirna_id.right.find(mirna_id) == i_to_mirna_id.right.end()) {
            unsigned long long new_i = i_to_mirna_id.size();
            i_to_mirna_id.insert(boost::bimap<unsigned long long, Mirna_id>::value_type(new_i, mirna_id));
        }
    }
    for(auto & e : this->gene_to_mirnas_arcs) {
        Gene_id gene_id = e.first;
        if(j_to_gene_id.right.find(gene_id) == j_to_gene_id.right.end()) {
            unsigned long long new_j = j_to_gene_id.size();
            j_to_gene_id.insert(boost::bimap<unsigned long long, Gene_id>::value_type(new_j, gene_id));
        }
    }
    if(i_to_mirna_id.size() != mirna_count) {
        std::cerr << "error: i_to_mirna_id.size() = " << i_to_mirna_id.size() << ", mirna_count = " << mirna_count << "\n";
        exit(1);
    }
    if(j_to_gene_id.size() != gene_count) {
        std::cerr << "error: j_to_gene_id.size() = " << j_to_gene_id.size() << ", gene_count = " << gene_count << "\n";
        exit(1);
    }
    
    for(auto & e : this->mirna_to_genes_arcs) {
        Mirna_id mirna_id = e.first;
        unsigned long long i = i_to_mirna_id.right.at(mirna_id);
        for(auto & f: e.second) {
            Gene_id gene_id = f;
            unsigned long long j = j_to_gene_id.right.at(gene_id);
            m[i][j] = this->mirna_gene_arcs[std::make_pair(mirna_id, gene_id)].size();
        }
    }
    std::stringstream ss;
    for(unsigned long long j = 0; j < gene_count; j++) {
        // the \" are used by the R function read.table() to distinguish between entries and the column and row names
        ss << "\"" << j_to_gene_id.left.at(j) << "\"" << ((j == gene_count - 1) ? "\n" : "\t");
    }
    for(unsigned long long i = 0; i < mirna_count; i++) {
        ss << "\"" << i_to_mirna_id.left.at(i) << "\"\t";
        for(unsigned long long j = 0; j < gene_count; j++) {
            ss << m[i][j] << ((j == gene_count - 1) ? "\n" : "\t");
        }
    }
    std::ofstream out(filename);
    out << ss.str();
    out.close();

    for(unsigned long long i = 0; i < mirna_count; i++) {
        delete [] m[i];
    }
    delete [] m;
    // std::cout << "written, ";
    // Timer::stop();
    
    // information about overlapping sites
    // Timer::start();
    std::cout << "writing information about overlapping sites\n";
    Output_buffer ob0(patient_folder + "overlapping_sites.tsv", 100000, 1000);
    std::string s = "gene_id\tutr_start\toverlapping_sites_count\n";
    ob0.add_chunk(s);
    for(auto & e : this->gene_to_sites_arcs) {
        auto & gene_id = e.first;
        auto & sites = e.second;
        for(auto & site : sites) {
            std::stringstream ss;
            ss << gene_id << "\t" << site->utr_start << "\t" << this->site_to_overlapping_sites[site].size() << "\n";
            s = ss.str();
            ob0.add_chunk(s);
        }
    }
    // std::cout << "written, ";
    // Timer::stop();

    // information about clusters
    // Timer::start();
    std::cout << "writing information about clusters\n";
    Output_buffer ob1(patient_folder + "clusters.tsv", 100000, 1000);
    s = "gene_id\tcluster_size\n";
    ob1.add_chunk(s);
    for(auto & e : this->gene_to_clusters_arcs) {
        auto & gene_id = e.first;
        auto & clusters = e.second;
        for(auto & cluster : clusters) {
            std::stringstream ss;
            ss << gene_id << "\t" << cluster->sites.size() << "\n";
            s = ss.str();
            ob1.add_chunk(s);
        }
    }
    // std::cout << "written, ";
    // Timer::stop();
}

void Interaction_graph::free_pointers()
{
    // Timer::start();
    std::cout << "deleting all sites\n";
    for(auto & e : this->sites_by_location) {
        delete e.second;
    }
    // std::cout << "deleted, ";
    // Timer::stop();
    // Timer::start();
    std::cout << "deleting all clusters\n";
    for(auto & e : this->gene_to_clusters_arcs) {
        for(auto & cluster : e.second) {
            delete cluster;
        }
    }
    // std::cout << "deleted, ";
    // Timer::stop();    
}

Interaction_graph::~Interaction_graph()
{
    // this destructor is called by the copy-constructor triggered by Perturbation_analyzer when deleting the copy of the patient from the variable perturbed_patient when a new patient is being copied
    // so pointers must be freed manually or smart pointers, like std::shared_ptr, must be used
//    this->free_pointers();
}
