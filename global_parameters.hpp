#ifndef GLOBAL_PARAMETERS_H
#define GLOBAL_PARAMETERS_H

class Global_parameters {
public:
    static double mirna_threshold_rpm;
    static double gene_threshold_rpm;    
    static double epsilon;
    /*
      This parameter quantifies the ratio between the average absolute number of miRNAs total copies that one expects to be involved in the regulation of a cell and the average absolute number of mRNAs total copies that one finds in the cells.

      So, a lambda of 0.5 means that there are on average twice as many mRNAs than miRNAs; a lambda of 2 means that tehre are on average twice as many miRNAs than mRNAs.
      If lambda is very low, then most of the miRNAs will bind to some mRNAs; if lambda is very high, then most of the miRNAs will remain unbound because there will be no more mRNA to bind to.
      
      Please do not confuse the number of target sites with the number of mRNAs copies, lambda refers to the number of mRNA copies
     */
    static double lambda;
    static bool test_parallelization;

    static void load_from_json();
};

#endif // GLOBAL_PARAMETERS_H
