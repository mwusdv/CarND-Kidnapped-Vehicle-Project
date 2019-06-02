#include <math.h>
#include <string>

#include "helper_functions.h"
#include "utils.h"

using std::vector;
using std::string;

double multiv_gauss(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}



wheel_sampler::wheel_sampler(const std::vector<double>& sampling_probs) 
 : _gen(0)
 , _uniform_dist(0.0, 1.0) {
    int N = sampling_probs.size();

    // build the wheel for sampling
    _wheel.resize(N+1);
    _wheel[0] = 0;
    for (int i = 0; i < N; ++i) {
        _wheel[i+1] = _wheel[i] + sampling_probs[i];
    }
    _wheel[N] = 1.01;  // left close, right open
}

wheel_sampler::~wheel_sampler() {
    // empty
}

// sample once, return index
int wheel_sampler::sample() {
    double random_value = _uniform_dist(_gen);

    // search range
    int left = 0, right = _wheel.size() - 2;
    while(left <= right) {
        int mid = (left + right) / 2;
        if (_wheel[mid] <= random_value && random_value < _wheel[mid+1]) {
            return mid;
        }
        else if (random_value < _wheel[mid]) {
            right = mid - 1;
        }
        else {
            left = mid + 1;
        }
    }

    // something is wrong if we reach here
    throw string("Cannot make a proper sampling");
}
