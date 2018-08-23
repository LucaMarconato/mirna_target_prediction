#include <iostream>
#include <fstream>

#include "mirna.hpp"

int main(int argc, char * argv[])
{
    char * gene [4] = {"AAAUUGUGUG", "AAAAAAAAAAA", "AUAUAUUAUA", "UUUUUUUUU"};
    for(int i = 0; i < 4; i++) {
        std::cout << Mirna::format_mirna_matching(Mirna::human_mirnas[0].matches_with_string(gene[i])) << "\n";
    }
    
    return 0;
}
