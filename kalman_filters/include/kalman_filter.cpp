/*
Constant Velocity Kalman Filter implementation
Authors: Edgar Macias Garcia (edgar.macias.garcia@intel.com)
Human Robot Collaboraiton (HRC), Intel Labs
*/

//Headers
#include "kalman_filter.hpp"

//Constructor
kalman_vf::kalman_vf(int size, float time)
{

    //Init Model parameters
    A = Eigen::MatrixXf::Identity(2 * size, 2 * size);
    H = Eigen::MatrixXf::Zero(2 * size, 2 * size);

    //Init Filter parameters
    P = Eigen::MatrixXf::Identity(2 * size, 2 * size);
    K = Eigen::MatrixXf::Identity(2 * size, 2 * size);
    Q = Eigen::MatrixXf::Identity(2 * size, 2 * size);
    R = Eigen::MatrixXf::Identity(2 * size, 2 * size);

    //Init Signals
    X = Eigen::VectorXf::Zero(2 * size);
    F = Eigen::VectorXf::Zero(2 * size);

    //Time step
    T = time;

    //Kalman size
    dim = size;

    //Setup time step
    for (int i = 0; i < dim; i++)
    {
        A(2 * i, 2 * i + 1) = time;
    }

    //Change status
    status = true;
    
    //Default values
    k_rate = 30;
    Q_gain = 0.01;
    R_gain = 0.01;
    first_det = false;
    counter = 0;

    //Output
    signal_filtered = Eigen::VectorXf::Zero(2 * size);
    signal_history = Eigen::VectorXf::Zero(size);
}

//Reset
void kalman_vf::reset()
{
    //Init Model parameters
    A = Eigen::MatrixXf::Identity(2 * dim, 2 * dim);
    H = Eigen::MatrixXf::Zero(2 * dim, 2 * dim);

    //Init Filter parameters
    P = Eigen::MatrixXf::Identity(2 * dim, 2 * dim);
    K = Eigen::MatrixXf::Identity(2 * dim, 2 * dim);
    Q = Eigen::MatrixXf::Identity(2 * dim, 2 * dim);
    R = Eigen::MatrixXf::Identity(2 * dim, 2 * dim);

    //Init Signals
    X = Eigen::VectorXf::Zero(2 * dim);
    F = Eigen::VectorXf::Zero(2 * dim);

    //Setup time step
    for (int i = 0; i < dim; i++)
    {
        A(2 * i, 2 * i + 1) = T;
    }

    //Change status
    status = true;

    //Output
    signal_filtered = Eigen::VectorXf::Zero(2 * dim);
    signal_history = Eigen::VectorXf::Zero(dim);
}

Eigen::VectorXf kalman_vf::upgrade(std::vector<float> senial, bool arrived)
{
    //Auxiliary vectors
    Eigen::MatrixXf I = Eigen::MatrixXf::Identity(2 * dim, 2 * dim);

    //std::cout << "Upgrading Kalman..." << std::endl;

    //Main loop
    for (int i = 0; i < dim; i++)
    {
        //Get previous filtered signal
        F(2 * i) = signal_filtered(2 * i);
        F(2 * i + 1) = signal_filtered(2 * i + 1);

        //std::cout << "F: " << F << std::endl;

        //Get goal signal
        if(arrived)
        {
            //Build states measurement vector
            X(2*i) = senial[i];

            //Save signal history
            signal_history[i] = senial[i];
        }
        else 
        {
            //Build states measurement vector
            X(2*i) = signal_history[i];
        }
    }

    //std::cout << "F_old: " << F << std::endl;
    //std::cout << "P_inc: " << X << std::endl;
    //std::cout << "P1: " << At << std::endl;

    P = A * P * At + Q;
    //std::cout << "P2: " << P << std::endl;

    //Compute the Kalman gain
    K = P * Ht * (H * P * Ht + R).inverse();

    //std::cout << "K: " << K << std::endl;

    //Control to kalman states
    F = F + K * (X - H * F);

    //std::cout << "F: " << F << std::endl;
    //std::cout << "F_new: " << F << std::endl;

    //Update the error covariance
    P = P - K * H * P;

    //std::cout << "P: " << P << std::endl;

    /*
    //Update the estimate via Zk
    if (!arrived)
    {
        //std::cout << "F_old: " << F << std::endl;

        //std::cout << A << std::endl;
        //std::cout << F << std::endl;

        //Time Update
        F = A * F;

        //std::cout << "F2: " << F << std::endl;

        //std::cout << "F: " << F << std::endl;
        //std::cout << "F_new: " << F << std::endl;
    }
    */

    //Output signal
    //std::cout << "---------" << std::endl;
    for (int i = 0; i < 2 * dim; i++)
    {
        signal_filtered(i) = F(i);
    }

    return signal_filtered;
}

void kalman_vf::set_parameters(std::vector<float> R_deff, std::vector<float> Q_deff, std::vector<float> position, float T_int)
{
    //Set time interval
    T = T_int;

    //Save init gains
    k_rate = T_int;
    Q_gain = Q_deff[0]; 
    R_gain = R_deff[0]; 

    //std::cout << "Setting parameters" << std::endl;

    //Main loop
    for (int i = 0; i < dim; i++)
    {
        //Set transformation matrix
        A(2 * i, 2 * i + 1) = T;

        //Set coefficients
        R(2*i, 2*i) = R_deff[i];
        R(2*i + 1, 2*i + 1) = 2.0*R_deff[i];

        //H Matrix
        H(2*i, 2 * i) = 1.0;
        H(2*i + 1, 2*i + 1) = 1.0;

        //Time-dependet arrays
        Q(2 * i, 2 * i) = 0.25 * pow(T, 4) * Q_deff[i];
        Q(2 * i, 2 * i + 1) = 0.5 * pow(T, 3) * Q_deff[i];
        Q(2 * i + 1, 2 * i) = 0.5 * pow(T, 3) * Q_deff[i];
        Q(2 * i + 1, 2 * i + 1) = pow(T, 2) * Q_deff[i];

        //Set position
        signal_filtered(2 * i) = position[i];
    }

    //Optimize calculus
    At = A.transpose();
    Ht = H.transpose();

    //std::cout << "A: " << A << std::endl;
    //std::cout << "At: " << At << std::endl;

    //std::cout << "H: " << H << std::endl; 

    //std::cout << "R:" << R << std::endl;

    //std::cout << "Q:" << Q << std::endl;

}

void kalman_vf::set_position(std::vector<float> position)
{
    //Main loop
    for (int i = 0; i < dim; i++)
    {
        //Set position
        signal_filtered(2 * i) = position[i];
    }
}

std::vector<float> kalman_vf::get_position()
{
    std::vector<float> position;

    for (int i = 0; i < dim; i++)
    {
        position.push_back(signal_filtered(2 * i));
    }

    return position;
}

std::vector<float> kalman_vf::get_velocity()
{
    std::vector<float> vel;

    for (int i = 0; i < dim; i++)
    {
        vel.push_back(signal_filtered(2 * i + 1));
    }

    return vel;
}

//Kalman filter struct constructor
kalman_struct::kalman_struct(int n_filters, int dim, float k_rate, float Q_gain, float R_gain)
{
    //Init Kalman filters
    for(int i = 0; i < n_filters; i++)
    {
        //Create Kalman filter
        this->filters.push_back(new kalman_vf(dim, 1.0/k_rate));
        
        //Configure filter
        this->filters[i]->set_parameters(std::vector<float>{R_gain}, std::vector<float>{Q_gain}, std::vector<float>{0.0}, 1.0 / k_rate);
    }
}

//Update filter
std::vector<float> kalman_struct::update_filter(int id, float value)
{
    //Update specific filter
    this->filters[id]->upgrade(std::vector<float>{value}, true);

    //Get filtered signal
    std::vector<float> out = this->filters[id]->get_position();

    return out;

}