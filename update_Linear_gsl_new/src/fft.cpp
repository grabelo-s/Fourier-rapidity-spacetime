#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <complex>
#include <algorithm>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include "fft.hpp"

std::vector<Data> resample_data(const std::vector<Data>& data, int new_size) {

    std::vector<Data> resampled_data;
    
    size_t size = data.size();

    double data_y[size];
    double data_fy[size];

    for(int i=0; i < size; i++) {
        data_y[i] = data[i].y; 
        data_fy[i] = data[i].fy; 
    }

    gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, size);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline_init(spline, data_y, data_fy, size);

    double y_min = data.front().y; 
    double y_max = data.back().y; 

    double step  = (y_max - y_min) / (new_size - 1);

    double new_y = y_min;

    for (int i = 0; i < new_size - 1; i++) {
        new_y += step;
        double new_fy = (new_y < y_min || new_y > y_max) ? 0 : gsl_spline_eval(spline, new_y, acc);
        resampled_data.emplace_back(Data(new_y, new_fy));

        //std::cout << new_y << "\t" << new_fy << "\t" <<  data[i].y << "\t" << data[i].fy << std::endl;
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    
    std::ofstream interpfile("interp_quarks_udud.dat");
    
    for (auto [y,fy]: resampled_data){
         interpfile << y << " " << fy << "\n";
     
    } 
	
    interpfile.close();
    
    return resampled_data;
}

double interpolation(double yp) {

    double fyp = 0.0;
    double p  = 1.0;
    double y0 = 0.0;
    double fy0 = 0.0;
    double y1 = 0.0;
    double fy1 = 0.0;

	std::vector<Data> data = read_data("xv_ud_sea_ud.dat");
    std::vector<Data> data_resampled;
    
    std::vector<Data> data2;

    const int new_size = data.size(); 
    data_resampled = resample_data(data, new_size);
    
    data_resampled.erase(data_resampled.end()-1, data_resampled.end());

    for (auto [y, fy] : data_resampled) {
        if (yp > y) {
            y0 = y;
            fy0 = fy;
        } else if (yp == y) {
            return fy;
        } else if (yp < y) {
            y1 = y;
            fy1 = fy;
            break;
        }
    }
   
    data2.push_back({y0, fy0});
    data2.push_back({y1, fy1});

    for (auto [yi, fyi] : data2) {
        p = 1.0;
        for (auto [yj, fyj] : data2) {
            if (yi != yj) {
                p *= (yp - yj) / (yi - yj);
            }
        }
        fyp += p * fyi;
    }
	return fyp;
}   
double real_integrand(double yp, void *params) {

	double eta_s=*((double *) params);
	
	//std::cout << eta_s << "\r" << std::flush;
	double tau_m=1.0;


 	return interpolation(yp) * std::cos(std::sinh(eta_s) * std::sinh(yp))*std::cosh(yp);
    
}
double imag_integrand(double y, void *params) {

	double eta_s=*((double *) params);
	double tau_m=1.0;

 	return interpolation(y)*std::sin(std::sinh(eta_s) * std::sinh(y))*std::cosh(y);
    
}
 
/*
std::complex<double> gauss_quadrature(const std::vector<Data>& data, double tau_m, double eta_s){ 
    int n = data.size();
    std::complex<double> integral = 0.0;
    const int Nq = 1;
    const double xg[] = {-0.9739, -0.8650, -0.6794, -0.4333, -0.1488, 0.1488, 0.4333, 0.6794, 0.8650, 0.9739};
    const double wg[] = {0.0666, 0.1494, 0.2190, 0.2692, 0.2955, 0.2955, 0.2692, 0.2190, 0.1494, 0.0666};

    for (int j = 0; j < n - 1; ++j) {

        double a = data[j].y; 
        double b = data[j + 1].y;
        double mid = 0.5 * (a + b);
        double half_length = 0.5 * (b - a);

        for (int k = 0; k < Nq; ++k) {
            double y_prime = mid + half_length * xg[k];
            double f_interp = (data[j].fy + data[j + 1].fy) / 2.0;
            integral += wg[k] * f_interp * custom_exponent(tau_m, eta_s, y_prime) * half_length * std::cosh(y_prime);
        }
    }
    return integral;
}
*/

std::complex<double> custom_exponent(double tau_m, double eta_s, double y_prime) {
    return std::exp(-std::complex<double>(0, 1) * tau_m * std::sinh(eta_s) * std::sinh(y_prime));
    // return std::exp(-std::complex<double>(0, 1) * eta_s * y_prime);
}

std::vector<Data> read_data(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Data> data;
    double y, fy;

    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;
    }

    while (file >> y >> fy) {
        data.emplace_back(Data(-std::log(y), fy));
        //data.emplace_back(Data(y, fy));
    }

    file.close();

    std::sort(data.begin(), data.end());
    data.erase(std::unique(data.begin(), data.end()), data.end());

    return data;
}
