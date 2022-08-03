#include "segmentation.hpp"
#include "visualizer.hpp"
#include <iostream>
using namespace std;

int main() {

  Visualizer v;

  cout << "Define environment and path then simulate" << endl;
  vector<vecext<double>> ps = {//{-5, 10, M_PI / 2, 2.5, 2.5, 0.5},
                               //{8, 8, M_PI / 3, 3, 3, 1.5},
                               //{4, 0, 0, 0.2, 0.2, 1.5},
                               {0, 0, M_PI / 4, 2, 3, 1.7}};
  vector<vecext<double>> ls;
  for (auto &p : ps)
    v.add_ellipse(p.data(), "r-");

  for (auto a : linspace<double>(0, M_PI / 2, 10)) {
    ls.push_back({6 * cos(a), 6 * sin(a)});
  }
  // for (auto a : linspace<double>(1, 6, 6)) {
  //   ls.push_back({0, 6 + a});
  // }
  // for (auto a : linspace<double>(0, M_PI * -2, 40)) {
  //   ls.push_back({6 * cos(a), 6 * sin(a)});
  // }

  // ls.resize(2); // testing first n scans
  for (auto &l : ls) {
    Matrix<double, 1, 2> m;
    m << l[0], l[1];
    v.add_points(m, "r.");
  }
  auto data = full_sim(ps, ls, 0.025, 0.02);

  cout << "Preprocessing (to structs)" << endl;

  cout << "test" << endl;

  vector<Scan> scans = data_to_scans(data, ls);

  // vector<Scan> scans = import_scans("/home/user/Documents/2dslam/testdata/");
  cout << "Iterate over scans" << endl;
  // Filter variables
  DiagonalMatrix<double, 6> Q1(pow(0.1, 2), pow(0.1, 2), pow(0.05, 2),
                               pow(0.1, 2), pow(0.1, 2), pow(0.05, 2));
  MatrixXd W1(2, 2);
  W1 << 1e-2, 0, 0, 0;
  DiagonalMatrix<double, 6> Q2(pow(0.1, 2), pow(0.1, 2), pow(0.1, 2),
                               pow(0.1, 2), pow(0.1, 2), pow(1, 2));
  MatrixXd W2(2, 2);
  W2 << 1, 0, 0, 0;
  // global detected ps variable
  vector<Vector<double, 6>> P;
  for (auto &scan : scans) {
    v.add_points(scans.back().pts, "r.");
    cout << "    Segmenting scan" << endl;
    // perfect for multiple objects
    // auto aggregates = segment_scan(scan, P, Q, W, 1300, 3, 20, 700);
    auto aggregates = segment_scan(scan, P, Q1, W1, 1300, 10, 20, 10000);
    // iterate over aggregates for augmentation
    for (auto &a : aggregates) {
      if (a.idx > -1) {
        cout << "        Augmenting scan: ";
        augment_scan(a, P, Q2, W2);
      } else {
        P.push_back(a.p);
      }
      // v.add_points(a.pts, "g.");
    }
  }
  for (auto p : P)
    v.add_ellipse(p.data(), "");
  v.show();
  return 0;
}
