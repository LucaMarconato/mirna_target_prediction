#ifndef PERTURBATION_ANALYZER_H
#define PERTURBATION_ANALYZER_H

#include <memory>

#include "matchings_predictor.hpp"

enum class Perturbation_type {
	No_perturbation,
	Point_perturbation,
	Gaussian_perturbation
	// Bernoully_perturbation
};

class Perturbation_target {
public:
    class Specific_mirna {
    public:
        Specific_mirna() {}
        Specific_mirna(Mirna_id mirna_id) : mirna_id(mirna_id) {}
        Mirna_id mirna_id;
        bool is_valid = false;        
    };
    class Specific_gene {
    public:
        Specific_gene() {}
        Specific_gene(Gene_id gene_id) : gene_id(gene_id) {}
        Gene_id gene_id;
        bool is_valid = false;
    };
    class Nth_largest_element {
    public:
        Nth_largest_element() {}
        Nth_largest_element(int n) : n(n) {}
        int n;
        bool is_valid = false;
    };
    class Elements_from_nth_largest {
    public:
        Elements_from_nth_largest() {}
        Elements_from_nth_largest(int n) : n(n) {}
        int n;
        bool is_valid = false;        
    };
    class Empty_target {
    public:
        Empty_target() {}
        bool is_valid = false;
    };

    Specific_mirna specific_mirna;
    Specific_gene specific_gene;
    Nth_largest_element nth_largest_element;
    Elements_from_nth_largest elements_from_nth_largest;
    Empty_target empty_target;

    Perturbation_target() {}
    Perturbation_target(Specific_mirna specific_mirna) : specific_mirna(specific_mirna.mirna_id) {this->specific_mirna.is_valid = true;}
    Perturbation_target(Specific_gene specific_gene) : specific_gene(specific_gene.gene_id) {this->specific_gene.is_valid = true;}
    Perturbation_target(Nth_largest_element nth_largest_element) : nth_largest_element(nth_largest_element.n) {this->nth_largest_element.is_valid = true;}
    Perturbation_target(Elements_from_nth_largest elements_from_nth_largest) : elements_from_nth_largest(elements_from_nth_largest.n) {this->elements_from_nth_largest.is_valid = true;}
    Perturbation_target(Empty_target empty_target) {this->empty_target.is_valid = empty_target.is_valid = true; /* just to avoid an unused parameter warning */}
    std::string string_id();
};

class Perturbation_extent {
public:
    /*
      If the perturbation type is set to a point perturbation, then the variable "value" of the "relative_perturbation" specify the relative change to apply to the target.
      examples:
      mirna_perturbation_amplifier = 0 means no perturbation
      mirna_perturbation_amplifier = -1 means to delete all copies of the target mirna
      mirna_perturbation_amplifier = 1 means to double the copies of the target mirna

      If the pertrubation type is set to a gaussian perturbation, then "value" specify how to amplify the gaussian perturbation.
      examples:
      mirna_perturbation_amplifier = 0 means no perturbation
      mirna_perturbation_amplifier = -1, 1 means that the perturbation is given by a standard normal
    */
    class Relative_perturbation {
    public:
        Relative_perturbation() {}
        Relative_perturbation(double value) : value(value) {}
        double value;
        bool is_valid = false;
    };
    class Absolute_perturbation {
    public:
        Absolute_perturbation() {}
        Absolute_perturbation(double value) : value(value) {}
        double value;
        bool is_valid = false;
    };
    class No_perturbation {
    public:
        No_perturbation() {}
        bool is_valid = false;
    };

    Relative_perturbation relative_perturbation;
    Absolute_perturbation absolute_perturbation;
    No_perturbation no_perturbation;

    Perturbation_extent() {}
    Perturbation_extent(Relative_perturbation relative_perturbation) : relative_perturbation(relative_perturbation.value) {this->relative_perturbation.is_valid = true;}
    Perturbation_extent(Absolute_perturbation absolute_perturbation) : absolute_perturbation(absolute_perturbation.value) {this->absolute_perturbation.is_valid = true;}
    Perturbation_extent(No_perturbation no_perturbation) {this->no_perturbation.is_valid = no_perturbation.is_valid = true; /* just to avoid an unused parameter warning */}
    std::string string_id();
};

class Perturbation_analyzer {
    Patient & patient;
    Patient perturbed_patient;
    std::unique_ptr<Matchings_predictor> matchings_predictor;
    std::stringstream simulation_id;

    void integrity_check();
    void perturb();
    
public: 
    Perturbation_type mirna_perturbation_type;
    Perturbation_type gene_perturbation_type;
    Perturbation_target mirna_perturbation_target;
    Perturbation_target gene_perturbation_target;
    Perturbation_extent mirna_perturbation_extent;
    Perturbation_extent gene_perturbation_extent;

    Perturbation_analyzer(Patient & patient);
    void run(Perturbation_type mirna_perturbation_type,
             Perturbation_type gene_perturbation_type,
             Perturbation_target mirna_perturbation_target,
             Perturbation_target gene_perturbation_target,
             Perturbation_extent mirna_perturbation_extent,
             Perturbation_extent gene_perturbation_extent,
             std::string output_name_suffix = "");
};

#endif // PERTURBATION_ANALYZER_H
