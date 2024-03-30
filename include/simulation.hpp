#pragma once
#include "asc_file.hpp" 

namespace Flowtastic {

    class Simulation{
        public:
        Simulation(const AscFile & asc_file): asc_file(asc_file){};
        
        AscFile asc_file;

        private:

        double avg_lobe_thickness = 0;
        double lobe_area = 0;
        double max_semiaxis = 0;
        int max_cells = 0;
        double thickness_min = 0; 
    };

}