#ifndef LS_HPP
#define LS_HPP
#include "preproc.hpp"
#include <ceres/autodiff_cost_function.h>
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>

double cov(VectorXd x, VectorXd y) {
  auto xm = x.array() - x.mean();
  auto ym = y.array() - y.mean();
  auto m = xm.array() * ym.array();
  return m.sum() / (m.size() - 1);
}

VectorXd init(VectorXd x, VectorXd y, Vector2d loc) {
  double *p = new double[9];
  MatrixXd m({{cov(x, x), cov(x, y)}, {cov(y, x), cov(y, y)}});
  JacobiSVD<MatrixXd> svd(m, ComputeFullU | ComputeFullV);
  auto V = svd.matrixV();
  auto D = svd.singularValues();
  double cond = D(0) / D(D.size() - 1);
  double th0 = atan2(V(1, 0), V(0, 0));
  // rotation of ptcloud
  MatrixXd d(x.size(), 2);
  d.col(0) = x;
  d.col(1) = y;
  MatrixXd rotm(2, 2);
  rotm << cos(-th0), -sin(-th0), sin(-th0), cos(-th0);
  MatrixXd rotd = (rotm * d.transpose()).transpose();
  double x0 = d.col(0).mean();
  double y0 = d.col(1).mean();
  double a0 = (rotd.col(0).maxCoeff() - rotd.col(0).minCoeff()) / 2;
  double b0 = (rotd.col(1).maxCoeff() - rotd.col(1).minCoeff()) / 2;
  double angle = atan2(y0 - loc(1), x0 - loc(0));
  p[0] = x0 + sqrt(a0 + b0) * cos(angle);
  p[1] = y0 + sqrt(a0 + b0) * sin(angle);
  p[6] = x0 + sqrt(a0 * a0 + b0 * b0) * cos(angle); // alt x0 and y0
  p[7] = y0 + sqrt(a0 * a0 + b0 * b0) * sin(angle);
  // using old a0 b0
  a0 = sqrt(2 * D(0)); // a0
  b0 = sqrt(2 * D(1)); // b0
  //
  if (b0 < 1)
    b0 = 1;
  if (a0 < 1)
    a0 = 1;
  p[2] = th0;
  p[3] = a0;
  p[4] = b0;
  p[5] = 1;
  p[8] = th0 + M_PI / 2; // alt th0
  VectorXd p_i = Map<Vector<double, 9>>(p);
  return p_i;
}

struct MyLossFunction {
  MyLossFunction(double x, double y) : _x(x), _y(y) {}

  template <typename T> bool operator()(const T *const p, T *residual) const {
    T f1 = ((_x - p[0]) * cos(p[2]) + (_y - p[1]) * sin(p[2])) / p[3];
    T f2 = ((_x - p[0]) * sin(p[2]) - (_y - p[1]) * cos(p[2])) / p[4];
    residual[0] =
        (p[3] * p[4]) *
        (pow(pow(pow(f1, 2.0), (1.0 / p[5])) + pow(pow(f2, 2.0), (1.0 / p[5])),
             p[5]) -
         1.0);
    return true;
  }

private:
  double _x;
  double _y;
};

Vector<double, 6> LS(VectorXd xdata, VectorXd ydata, Vector2d loc) {
  auto p0 = init(xdata, ydata, loc);
  double *p1 = new double[]{p0[0], p0[1], p0[2], p0[3], p0[4], 1.98};
  double *p2 = new double[]{p0[6], p0[7], p0[2], p0[3], p0[4], 1};
  double *p3 = new double[]{p0[6], p0[7], p0[2], p0[3], p0[4], 1.98};
  // theta variation
  double *p4 = new double[]{p0[0], p0[1], p0[8], p0[3], p0[4], p0[5]};
  double *p5 = new double[]{p1[0], p1[1], p1[8], p1[3], p1[4], p1[5]};
  double *p6 = new double[]{p2[0], p2[1], p2[8], p2[3], p2[4], p2[5]};
  double *p7 = new double[]{p3[0], p3[1], p3[8], p3[3], p3[4], p3[5]};
  double *pa[] = {p0.data(), p1, p2, p3, p4, p5, p6, p7};
  ceres::Problem *problem = new ceres::Problem[8]; // destructor segfaults
  for (unsigned i = 0; i < xdata.size(); i++) {
    ceres::CostFunction *cost_function =
        new ceres::AutoDiffCostFunction<MyLossFunction, 1, 6>(
            new MyLossFunction(xdata(i), ydata(i)));
    for (unsigned i = 0; i < 8; i++)
      problem[i].AddResidualBlock(cost_function, nullptr, pa[i]);
  }
  // bounds
  for (unsigned i = 0; i < 8; i++) {
    problem[i].SetParameterLowerBound(pa[i], 3, 1);    // a
    problem[i].SetParameterLowerBound(pa[i], 4, 1);    // b
    problem[i].SetParameterLowerBound(pa[i], 5, 0.1);  // eps <- very important
    problem[i].SetParameterUpperBound(pa[i], 3, 30);   // a
    problem[i].SetParameterUpperBound(pa[i], 4, 30);   // b
    problem[i].SetParameterUpperBound(pa[i], 5, 1.99); // eps
    // avoids auto-occlusion
    problem[i].SetParameterUpperBound(pa[i], 0, max(p0[6], p0[0])); // x
    problem[i].SetParameterLowerBound(pa[i], 0, min(p0[6], p0[0])); // x
    problem[i].SetParameterUpperBound(pa[i], 1, max(p0[7], p0[1])); // y
    problem[i].SetParameterLowerBound(pa[i], 1, min(p0[7], p0[1])); // y
  }
  ceres::Solver::Options options;
  options.num_threads = 24;
  options.minimizer_progress_to_stdout = false;
  options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY; //<- slower
  options.preconditioner_type = ceres::CLUSTER_JACOBI;       //<- also slower
  options.initial_trust_region_radius = 1e8;                 // important
  options.max_num_iterations = 1000;
  ceres::Solver::Summary summary[8];
  for (unsigned i = 0; i < 8; i++)
    ceres::Solve(options, &problem[i], &summary[i]);
  unsigned opt = 0;
  double min_cost = 10000;
  for (unsigned i = 0; i < 8; i++)
    if (summary[i].final_cost < min_cost) {
      min_cost = summary[i].final_cost;
      opt = i;
    }
  Vector<double, 6> ans = Map<Vector<double, 6>>(pa[opt]);
  return ans;
}

#endif // LS_HPP
