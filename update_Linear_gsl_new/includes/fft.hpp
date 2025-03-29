#include <complex>
#include <vector>
#include <cmath>

struct Data {
    double y;
    double fy;
    Data();
    Data(double m_y, double m_fy) : y(m_y), fy(m_fy) {};

    bool operator==(const Data& data) {
        return (y == data.y);
    }

    bool operator<(const Data& data) {
        return (y < data.y);
    }
};

static std::vector<Data> data_resampled;


double interpolation(double yp);

double real_integrand(double y, void *params);
double imag_integrand(double y, void *params);
std::complex<double> gauss_quadrature(const std::vector<Data>& data, double tau_m, double eta_s);
std::complex<double> custom_exponent(double tau_m, double eta_s, double y_prime);
std::vector<Data> read_data(const std::string& filename);

std::vector<Data> resample_data(const std::vector<Data>& data, int new_size);
















