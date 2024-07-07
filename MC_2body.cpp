//Simulation of Charged Pion decay
//Author: Avinandan Mondal, Dept. of Physics, IIT-M

#include<iostream>
#include<cmath>
#include<random>
#include<vector>
#include<fstream>
#include<algorithm>

using namespace std;

//Decay is A -> B + C

double p_CM(double mA, double mB, double mC){
    return(sqrt((pow(mA,2)-pow((mB+mC),2)) * (pow(mA,2) - pow((mB-mC),2)))/(2*mA));
}
double E(double m, double p){
    return(sqrt(pow(m,2) + pow(p,2)));
}

double* Lab(double EA_Lab, double mA, double mB, double mC, double x){ //x is cos(theta_CM)
    double* arr = new double[11];
    double p = p_CM(mA, mB, mC);
    double EB_CM = E(mB,p);
    double EC_CM = E(mC,p);
    double y = sqrt(1-pow(x,2));
    double pB_z_CM = p*x;
    double pB_x_CM = p*y;
    double pC_z_CM = -pB_z_CM;
    double pC_x_CM = -pB_x_CM;
    double gamma = EA_Lab/mA;
    double beta = sqrt(1-1/pow(gamma,2));
    double pB_z_Lab = gamma*pB_z_CM + gamma*beta*EB_CM;
    double pB_x_Lab = pB_x_CM;
    double EB_Lab = gamma*EB_CM + gamma*beta*pB_z_CM;
    double pB_Lab = sqrt(pow(pB_z_Lab,2) + pow(pB_x_Lab,2));
    double theta_B = acos(pB_z_Lab/pB_Lab)*180/M_PI; //Angles in degrees
    double pC_z_Lab = gamma*pC_z_CM + gamma*beta*EC_CM;
    double pC_x_Lab = pC_x_CM;
    double EC_Lab = gamma*EC_CM + gamma*beta*pC_z_CM;
    double pC_Lab = sqrt(pow(pC_z_Lab,2) + pow(pC_x_Lab,2));
    double theta_C = acos(pC_z_Lab/pC_Lab)*180/M_PI; //Angles in degrees
    arr[0] = pB_z_Lab;
    arr[1] = pB_x_Lab;
    arr[2] = pB_Lab;
    arr[3] = EB_Lab;
    arr[4] = theta_B;
    arr[5] = pC_z_Lab;
    arr[6] = pC_x_Lab;
    arr[7] = pC_Lab;
    arr[8] = EC_Lab;
    arr[9] = theta_C;
    arr[10] = theta_B + theta_C;
    return arr;
}

const double m_pion_mean = 139.57018;
const double m_pion_sd = 0.00035;
const double m_muon = 105.65837;
const double m_positron = 0.51100;
const double PDW = 0.999876; //Partial decay width of muon decay channel of pi+ 

double* pion_generator(){
    double* arr = new double[3];
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> distA(m_pion_mean, m_pion_sd);
    double m_pion = distA(gen); //Mass of generated pion
    double mom_decay_constant = 60.0; //Custom value, depends on beam parameters
    exponential_distribution<double> distB(1.0 / mom_decay_constant);
    double mom_pion = distB(gen);    
    double E_pi = E(m_pion, mom_pion);
    arr[0] = E_pi;
    arr[1] = m_pion;
    arr[2] = mom_pion;
    return arr;
}

double* decay_generator(){
    double* arr = new double[14];
    
    random_device ran;
    mt19937 gen(ran());
    uniform_real_distribution<double> distC(0.0, 1.0);
    double r = distC(gen);
    int chan_det;
    if(r <= PDW){
        chan_det = 0;  //Muon decay channel
    }
    else{
        chan_det = 1; //Positron decay channel
    }
    
    double* pion = pion_generator();
    double E_pi = pion[0];
    double m_pi = pion[1];
    double mom_pi = pion[2];
    delete[] pion;
    
    random_device rd;
    mt19937 gen1(rd());
    uniform_real_distribution<double> distD(-1.0, 1.0);
    double x = distD(gen1);    //Cosine of decay angle in CM frame
    
    double* lab_data;
    if(chan_det == 0){
        lab_data = Lab(E_pi, m_pi, m_muon, 0, x);  //Muon neutrino assumed massless 
    }
    else{
        lab_data = Lab(E_pi, m_pi, m_positron, 0, x);    //Electron neutrino assumed massless
    }

    arr[0] = chan_det;
    arr[1] = E_pi;
    arr[2] = mom_pi;
    for(int i=3; i<=13; i++){
        arr[i] = lab_data[i-3];
    }
    delete[] lab_data;
    
    return arr;
}

void charged_pion_decay_simulator(int N){ //N is number of pions simulated
    vector<double> pion_mom_muon_channel;
    vector<double> pion_mom_positron_channel;
    vector<double> muon_mom;
    vector<double> positron_mom;
    vector<double> muon_neutrino_mom;
    vector<double> electron_neutrino_mom;
    vector<double> muon_angle;
    vector<double> muon_neutrino_angle;
    vector<double> positron_angle;
    vector<double> electron_neutrino_angle;
    vector<double> muon_decay_opening_angle;
    vector<double> positron_decay_opening_angle;
    vector<double> muon_channel_avg_opening_angle;
    vector<double> positron_channel_avg_opening_angle;
    for(int j=0; j<N; j++){
        double* master_data = decay_generator();
        if(master_data[0] == 0){
            pion_mom_muon_channel.push_back(master_data[2]);
            muon_mom.push_back(master_data[5]);
            muon_neutrino_mom.push_back(master_data[10]);
            muon_angle.push_back(master_data[7]);
            muon_neutrino_angle.push_back(master_data[12]);
            muon_decay_opening_angle.push_back(master_data[13]);
        }
        else{
            pion_mom_positron_channel.push_back(master_data[2]);
            positron_mom.push_back(master_data[5]);
            electron_neutrino_mom.push_back(master_data[10]);
            positron_angle.push_back(master_data[7]);
            electron_neutrino_angle.push_back(master_data[12]);
            positron_decay_opening_angle.push_back(master_data[13]);      
        }
        delete[] master_data;
    }
    
    auto max_it_muon = max_element(pion_mom_muon_channel.begin(), pion_mom_muon_channel.end());
    auto max_it_positron = max_element(pion_mom_positron_channel.begin(), pion_mom_positron_channel.end());
    double max_pion_mom_muon = *max_it_muon;
    double max_pion_mom_positron = *max_it_positron;
    
    double n = 0;
    while(n<max_pion_mom_muon){
        double counter = 0;
        double sum_holder = 0;
        for(int i=0; i<pion_mom_muon_channel.size(); i++){
            if(n <= pion_mom_muon_channel[i] && pion_mom_muon_channel[i] < n+10){
                sum_holder += muon_decay_opening_angle[i];
                counter += 1;           
            }   
        }
        if(counter!=0){
            muon_channel_avg_opening_angle.push_back(sum_holder/counter);
        }
        n += 10;       
    }    

    double m = 0;
    while(m<max_pion_mom_positron){
        double counter1 = 0;
        double sum_holder_1 = 0;
        for(int i=0; i<pion_mom_positron_channel.size(); i++){
            if(m <= pion_mom_positron_channel[i] && pion_mom_positron_channel[i] < m+10){
                sum_holder_1 += positron_decay_opening_angle[i];
                counter1 += 1;           
            }   
        }
        if(counter1!=0){
            positron_channel_avg_opening_angle.push_back(sum_holder_1/counter1);
        }
        m += 10;       
    } 
    
    //Data storing in .csv files
    ofstream myfile("muon_channel_vals_pion_momentum.csv");
    for(auto elem : pion_mom_muon_channel){
        myfile << elem << "\n" ;
    }
    myfile.close();

    ofstream myfile1("muon_momentum_vals.csv");
    for(auto elem : muon_mom){
        myfile1 << elem << "\n" ;
    }
    myfile1.close();

    ofstream myfile11("muon_neutrino_momentum_vals.csv");
    for(auto elem : muon_neutrino_mom){
        myfile11 << elem << "\n" ;
    }
    myfile11.close();

    ofstream myfile2("muon_angle_vals.csv");
    for(auto elem : muon_angle){
        myfile2 << elem << "\n" ;
    }
    myfile2.close();

    ofstream myfile3("muon_neutrino_angle_vals.csv");
    for(auto elem : muon_neutrino_angle){
        myfile3 << elem << "\n" ;
    }
    myfile3.close();

    ofstream myfile4("muon_channel_opening_angle_vals.csv");
    for(auto elem : muon_decay_opening_angle){
        myfile4 << elem << "\n" ;
    }
    myfile4.close();

    ofstream myfile0("positron_channel_vals_pion_momentum.csv");
    for(auto elem : pion_mom_positron_channel){
        myfile0 << elem << "\n" ;
    }
    myfile0.close();

    ofstream myfile5("positron_momentum_vals.csv");
    for(auto elem : positron_mom){
        myfile5 << elem << "\n" ;
    }
    myfile5.close();

    ofstream myfile51("electron_neutrino_momentum_vals.csv");
    for(auto elem : electron_neutrino_mom){
        myfile51 << elem << "\n" ;
    }
    myfile51.close();

    ofstream myfile6("positron_angle_vals.csv");
    for(auto elem : positron_angle){
        myfile6 << elem << "\n" ;
    }
    myfile6.close();

    ofstream myfile7("electron_neutrino_angle_vals.csv");
    for(auto elem : electron_neutrino_angle){
        myfile7 << elem << "\n" ;
    }
    myfile7.close();

    ofstream myfile8("positron_channel_opening_angle_vals.csv");
    for(auto elem : positron_decay_opening_angle){
        myfile8 << elem << "\n" ;
    }
    myfile8.close();

    ofstream myfile9("muon_channel_avg_opening_angle_vals.csv");
    for(auto elem : muon_channel_avg_opening_angle){
        myfile9 << elem << "\n" ;
    }
    myfile9.close();

    ofstream myfile10("positron_channel_avg_opening_angle_vals.csv");
    for(auto elem : positron_channel_avg_opening_angle){
        myfile10 << elem << "\n" ;
    }
    myfile10.close();

}



int main(){
    cout << "Simulation in Progress-------------------------------" << endl;
    charged_pion_decay_simulator(10000); // !! SIMULATE AT LEAST 10000 PIONS TO PREVENT SEGMENTATION FAULT !! 
    cout << "Data written in corresponding .csv files." << endl;
    cout << "End of Simulation------------------------------------" << endl;
    return 0;
}






