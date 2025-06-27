/*
POM Kalman Filter implementation
Authors: Edgar Macias Garcia (edgar.macias.garcia@intel.com)
Human Robot Collaboration (HRC), Intel Labs
*/

//Std libs
#include <stdio.h>
#include <iostream>
#include <vector>

//Eigen
#include <Eigen/Dense>
#include <cmath>

class kalman_vf 
{
    public:

        //Init parameters
        float k_rate;
        float Q_gain, R_gain; 
        
        //Status
        bool status, first_det;
        int counter;

   
    private:

        //Init Model parameters
        Eigen::MatrixXf A, At;
        Eigen::MatrixXf H, Ht;

        //Init Filter parameters
        Eigen::MatrixXf P;
        Eigen::MatrixXf K;
        Eigen::MatrixXf Q;
        Eigen::MatrixXf R;

        //Init Signals
        Eigen::VectorXf X;
        Eigen::VectorXf F;

        //Time step
        float T;

        //Filter size
        int dim;

        //Output
        Eigen::VectorXf signal_filtered;
        Eigen::VectorXf signal_history;

    public:

        //Constructor
        kalman_vf(int size, float time);
    
        //Configuration functions
        void set_parameters(std::vector<float> R_deff, std::vector<float> Q_deff, std::vector<float> position, float T_int);
    
        //Processing functions
        Eigen::VectorXf upgrade(std::vector<float> senial, bool arrived);
        void set_position(std::vector<float> position);
        void reset();
    
        //Status function
        std::vector<float> get_position();
    	std::vector<float> get_velocity();

};

//Kalman struct
class kalman_struct
{
    public:

        //Kalman parameters
        int dim;
        float k_rate;
        float Q_gain, R_gain; 
        
        //Struct params
        int n_filters; 

        //Filters
        std::vector<kalman_vf*> filters;

        //Constructor
        kalman_struct(int n_filters, int dim, float k_rate, float Q_gain, float R_gain);

        //Update filters
        std::vector<float> update_filter(int id, float value);

};
