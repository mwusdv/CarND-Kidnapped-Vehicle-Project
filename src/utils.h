#ifndef UTILS_H_   
#define UTILS_H_

#include <vector>
#include <random>

// multivaraiage gaussian
double multiv_gauss(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y);

// randomly sample one accodring to the given probabilities
class wheel_sampler {
 public:
    wheel_sampler(const std::vector<double>& sampling_probs);
    ~wheel_sampler();

    // sample once, return index
    int sample();

 protected:
    std::vector<double> _wheel;
    std::default_random_engine _gen;
    std::uniform_real_distribution<double> _uniform_dist;
};

#endif