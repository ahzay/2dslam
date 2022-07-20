#include "segmentation.hpp"
#include <iostream>
using namespace std;

int main() {
  cout << "Define environment and path then simulate:" << endl;
  vector<vecext<double>> ps = {{0, 0, 0, 2, 2, 0.4},
                               {8, 8, M_PI / 3, 3, 3, 1.3}};
  vector<vecext<double>> ls;
  for (auto a : linspace<double>(0, M_PI, 5))
    ls.push_back({6 * cos(a), 6 * sin(a)});
  auto data = full_sim(ps, ls, 0.05, 0.02);
  cout << "Preprocessing (to structs):" << endl;
  vector<Scan> scans = data_to_scans(data, ls);
  cout << "Iterate over scans" << endl;
  for (auto scan : scans) {
    cout << "   Segmenting scan" << endl;
  }
  return 0;
}
