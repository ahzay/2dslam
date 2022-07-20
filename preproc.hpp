#ifndef PREPROC_HPP
#define PREPROC_HPP
#include "simulate.hpp"
#include <Eigen/SVD>
using namespace Eigen;
struct Scan {
  Matrix<double, Dynamic, 2> pts, mes;
  Vector<double, 2> loc;
};
vector<Scan> data_to_scans(vector<tuple<vecext<double>, vecext<double>,
                                        vecext<double>, vecext<double>>>
                               data,
                           vector<vecext<double>> ls) {
  vector<Scan> scans;
  for (int i = 0; i < data.size(); i++) {
    MatrixXd pts(get<0>(data[i]).size(), 2);
    MatrixXd mes(get<0>(data[i]).size(), 2);
    pts.col(0) = Map<VectorXd>(get<0>(data[i]).data(), get<0>(data[i]).size());
    pts.col(1) = Map<VectorXd>(get<1>(data[i]).data(), get<1>(data[i]).size());
    mes.col(0) = Map<VectorXd>(get<2>(data[i]).data(), get<2>(data[i]).size());
    mes.col(1) = Map<VectorXd>(get<3>(data[i]).data(), get<3>(data[i]).size());
    Vector<double, 2> loc = Map<Vector<double, 2>>(ls[i].data());
    scans.push_back({pts, mes, loc});
  }
  return scans;
}
#endif // PREPROC_HPP
