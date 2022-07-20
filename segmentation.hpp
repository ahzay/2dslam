#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP
#include "ekf.hpp"

struct Aggregate {
  Matrix<double, Dynamic, 2> pts, mes;
  Vector<double, 6> p;
  int idx = -1; // idx in ps global state vector
};

vector<Aggregate> segment_scan(Scan scan, double tol_flx, double tol_str,
                               int npts) {
  vector<Aggregate> ans;

  int i = 0; // measure counter
  int lasti; // i of start of last segment
}

#endif // SEGMENTATION_HPP
