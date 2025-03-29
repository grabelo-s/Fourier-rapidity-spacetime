#include "interpolation.hpp"

double linerp(const std::vector<Point> &data, double xp)
{
    std::vector<Point> data2;

    double yp = 0.0;
    double p  = 1.0;
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;

    for (auto [x, y] : data) {
        if (xp > x) {
            x0 = x;
            y0 = y;
        } else if (xp == x) {
            return y;
        } else if (xp < x) {
            x1 = x;
            y1 = y;
            break;
        }
    }
   
    data2.push_back({x0, y0});
    data2.push_back({x1, y1});

    for (auto [xi, yi] : data2) {
        p = 1.0;
        for (auto [xj, yj] : data2) {
            if (xi != xj) {
                p *= (xp - xj) / (xi - xj);
            }
        }
        yp += p * yi;
    }
    return yp;
}
