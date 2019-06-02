/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <limits>

#include "helper_functions.h"
#include "map.h"
#include "utils.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles

  // normal distributions for x, y and theta
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  // initialize particles
  particles.resize(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0;
  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; ++i) {
    double theta = particles[i].theta;
    particles[i].x += (velocity/yaw_rate*(sin(theta + yaw_rate*delta_t) - sin(theta)) + dist_x(gen));
    particles[i].y += (velocity/yaw_rate*(cos(theta) - cos(theta + yaw_rate*delta_t)) + dist_y(gen));
    particles[i].theta += (yaw_rate * delta_t + dist_theta(gen));
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (size_t i = 0; i < observations.size(); ++i) {
    double xo = observations[i].x;
    double yo = observations[i].y;
    double minDist = std::numeric_limits<double>::max();
    size_t nearest_id = -1;

    for (size_t j = 0; j < predicted.size(); ++j) {
      double xp = predicted[j].x;
      double yp = predicted[j].y;

      double d = dist(xo, yo, xp, yp);
      if (d < minDist) {
        minDist = d;
        nearest_id = j;
      }

      observations[i].id = nearest_id;
    }
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  for (size_t p = 0; p < particles.size(); ++p) {
    double xp = particles[p].x;
    double yp = particles[p].y;
    double theta = particles[p].theta;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double weight = 1.0;

    // scan each observation
    for (size_t o = 0; o < observations.size(); ++o) {
      double xo = observations[o].x;
      double yo = observations[o].y;

      // coordinate transformation
      double xt = cos_theta*xo - sin_theta*yo + xp;
      double yt = sin_theta*xo + cos_theta*yo + yp;

      // find nearest landmark
      const vector<Map::single_landmark_s>& landmark_list = map_landmarks.landmark_list;
      double min_dist = std::numeric_limits<double>::max();
      size_t nearest_id = -1;
      for (size_t l = 0; l < landmark_list.size(); ++l) {
        double d = dist(xt, yt, landmark_list[l].x_f, landmark_list[l].y_f);
        if (d < min_dist) {
          min_dist = d;
          nearest_id = l;
        }
      }

      // Gaussian
      double prob = multiv_gauss(std_landmark[0], std_landmark[1], xt, yt, 
                                 landmark_list[nearest_id].x_f, landmark_list[nearest_id].y_f);
      weight *= prob;
    }

    // update the weight of current particile
    particles[p].weight = weight;
  }
 
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // compute sampling probabilities
  vector<double> sampling_probs(particles.size());
  int sum = 0;
  for (size_t p = 0; p < particles.size(); ++p) {
    sampling_probs[p] = particles[p].weight;
    sum += sampling_probs[p];
  }
  for (size_t p = 0; p < sampling_probs[p]; ++p) {
    sampling_probs[p] /= sum;
  }

  // resampling
  wheel_sampler sampler(sampling_probs);
  vector<Particle> new_particles(particles.size());
  for (size_t i = 0; i < new_particles.size(); ++i) {
    int idx = sampler.sample();
    new_particles[i] = particles[idx];
  }
  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}