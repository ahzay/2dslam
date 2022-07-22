#include "segmentation.hpp"
#include "visualizer.hpp"
#include <iostream>
using namespace std;

int main() {
  cout << "Define environment and path then simulate" << endl;
  vector<vecext<double>> ps = {{8, 8, M_PI / 3, 3, 3, 1.3},
                               //{-5, 10, 0, 2, 2, 1.0},
                               {0, 0, 0, 2, 2, 0.4}};
  vector<vecext<double>> ls;
  Visualizer v;
  for (auto &p : ps)
    v.add_ellipse(p.data(), "r-");

  for (auto a : linspace<double>(0, M_PI / 2, 10)) {
    Matrix<double, 1, 2> m;
    m << 6 * cos(a), 6 * sin(a);
    v.add_points(m, "r.");
    ls.push_back({6 * cos(a), 6 * sin(a)});
  }
  for (auto a : linspace<double>(1, 6, 6)) {
    Matrix<double, 1, 2> m;
    m << 0, 6 + a;
    v.add_points(m, "r.");
    ls.push_back({0, 6 + a});
  }
  auto data = full_sim(ps, ls, 0.05, 0.02);

  cout << "Preprocessing (to structs)" << endl;
  vector<Scan> scans = data_to_scans(data, ls);

  cout << "Iterate over scans" << endl;
  // Filter variables
  DiagonalMatrix<double, 6> Q(pow(0.2, 2), pow(0.2, 2), pow(0.1, 2),
                              pow(0.2, 2), pow(0.2, 2), pow(0.1, 2));
  MatrixXd W(2, 2);
  W << 1e-2, 0, 0, 0;
  // global detected ps variable
  vector<Vector<double, 6>> P;
  for (auto &scan : scans) {
    v.add_points(scan.pts, "r.");
    cout << "    Segmenting scan" << endl;
    auto aggregates = segment_scan(scan, P, Q, W, 1300, 1, 20, 500);
    // iterate over aggregates for augmentation
    for (auto &a : aggregates) {
      if (a.idx > -1) {
        cout << "        Augmenting scan" << endl;
        augment_scan(a, P, Q, W);
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
