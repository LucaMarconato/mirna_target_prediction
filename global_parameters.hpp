#ifndef GLOBAL_PARAMETERS_H
#define GLOBAL_PARAMETERS_H

class Global_parameters {
public:
    static double mirna_threshold_rpm;
    static double gene_threshold_rpm;
    static double epsilon;

    static void load_from_json();
};

#endif // GLOBAL_PARAMETERS_H
