#include "visualizer.hpp"
#include "plotty/matplotlibcpp.hpp"
namespace plt = plotty;
Visualizer::Visualizer() {
  plt::figure();
  plt::axis("equal");
}
void Visualizer::add_ellipse(const double p[6], const char *label) {
  Eigen::ArrayXd t = Eigen::ArrayXd::LinSpaced(100, 0, 2 * M_PI);
  Eigen::ArrayXd X = t.cos().abs().pow(p[5]) * p[3] * t.cos().sign();
  Eigen::ArrayXd Y = t.sin().abs().pow(p[5]) * p[4] * t.sin().sign();
  Eigen::MatrixXd rotm(2, 2);
  rotm << cos(p[2]), -sin(p[2]), sin(p[2]), cos(p[2]);
  Eigen::MatrixXd d(X.size(), 2);
  d.col(0) = X.matrix();
  d.col(1) = Y.matrix();
  Eigen::MatrixXd rotd = (rotm * d.transpose()).transpose();
  rotd.col(0) = rotd.col(0).array() + p[0];
  rotd.col(1) = rotd.col(1).array() + p[1];
  plt::plot(rotd.col(0), rotd.col(1), label);
  // plt::axis("equal");
}
void Visualizer::add_points(Eigen::Matrix<double, Eigen::Dynamic, 2> pts,
                            const char *label) {
  plt::plot(pts(Eigen::all, 0), pts(Eigen::all, 1), label);
}
void Visualizer::show() { plt::show(); }
