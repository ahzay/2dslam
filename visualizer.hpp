#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP
#include <Eigen/Core>
class Visualizer {
public:
  Visualizer();
  void add_ellipse(const double p[6], const char *label);
  void add_points(Eigen::Matrix<double, Eigen::Dynamic, 2> pts,
                  const char *label);
  void show();
};

#endif // VISUALIZER_HPP
