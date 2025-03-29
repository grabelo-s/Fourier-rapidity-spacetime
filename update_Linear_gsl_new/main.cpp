#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <complex>

#include "fft.hpp"
#include <gsl/gsl_integration.h>



int main() {

	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	
	
	double params;
	double real_result;
	double imag_result;
	double real_abserr;  
	//size_t real_neval;  
	double imag_abserr;  
	//size_t imag_neval;
	
	gsl_function F_real;
	gsl_function F_imag;  
	F_real.function=&real_integrand;  
	F_imag.function=&imag_integrand; 
	F_real.params=&params;
	F_imag.params=&params;


    const double tau_m = 1.0;

    std::string input_filename = "xv_ud_sea_ud.dat";

    std::vector<Data> data = read_data(input_filename);
    std::vector<Data> data_resampled;

    const int new_size = data.size(); 
    data_resampled = resample_data(data, new_size);
    
    data_resampled.erase(data_resampled.end()-1, data_resampled.end());
    

    double eta_s_initial = -2.0;
    double eta_s_final = 2.0;
    double abs_F_eta_s;
    const double delta_eta_s = (eta_s_final - eta_s_initial) / 200;

    if (data.empty()) {
        std::cerr << "Error: No valid data read!" << std::endl;
        return 1;
    }
    
    /*

    std::ofstream file(".dat");
    if (!file) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    file.close();

	*/
	std::ofstream fourier_abs_file("Abs_Fourier_udud.dat");
	
    for (double eta_s = eta_s_initial; eta_s < eta_s_final; eta_s += delta_eta_s) {
        //abs_F_eta_s = std::sqrt(std::norm(gauss_quadrature(data_resampled, tau_m, eta_s)));
     
  		params = eta_s;
       gsl_integration_qag(&F_real, eta_s_initial, eta_s_final, 1e-1, 1e-2, 1000,6,w, &real_result, &real_abserr);
       gsl_integration_qag(&F_imag, eta_s_initial, eta_s_final, 1e-1, 1e-2, 1000,6,w, &imag_result, &imag_abserr);
      //std::cout << eta_s << " " << real_result << " " << imag_result << " " << sqrt(real_result*real_result + imag_result*imag_result) << "\n";
      
      double abs_F = sqrt(real_result * real_result + imag_result * imag_result);
      
       fourier_abs_file << eta_s << " " << abs_F << "\n";
       std::cout << eta_s << " " << abs_F << "\n";
       real_result = 0;

    
   }
   fourier_abs_file.close();
   


    std::ofstream transformed_file("y_" + input_filename);
    if (!transformed_file) {
        std::cerr << "Error opening transformed file!" << std::endl;
        return 1;
    }

    for (const auto& [y_val, f_val] : data) {
        transformed_file << y_val << " " << f_val << "\n";
    }

    transformed_file.close();

	gsl_integration_workspace_free(w); 
    return 0;
}

