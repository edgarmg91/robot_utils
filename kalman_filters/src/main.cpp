//Basic libraries
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>

//Import Kalman FIlter
#include "kalman_filter.hpp"

//Main loop
int main(int argc, char **argv)
{
 
    //Kalman filters params
    int dim = 1;
    int n_filters = 5;
    int k_rate = 15.0;
    float R_gain = 0.1;
    float Q_gain = 0.1; 

    //Create kalman struct
    kalman_struct* signal_filter = new kalman_struct(n_filters, dim, k_rate, Q_gain, R_gain);

    //Test signal
    int iters = 1000;
    for(int i = 0; i < iters; i++)
    {
        //Update filters
        for(int j = 0; j < signal_filter->filters.size(); j++)
        {
            //Create false signal
            float signal = 3.1416*sin(i/100);

            std::vector<float> out = signal_filter->update_filter(j, signal);

            std::cout << signal << ", " << out[0] << std::endl;
        }
    }

 
 return 0;
}
