#ifndef PREPROC_HPP
#define PREPROC_HPP
#include "simulate.hpp"
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <fstream>
#include <iostream>
namespace fsm = boost::filesystem;
using namespace Eigen;
struct Pose {
  Vector<float, 3> p;
  Quaternionf o;
};
struct Meas {
  float t, p, d; // theta(horiz), phi(vert), dist
};

struct Scan {
  Matrix<double, Dynamic, 2> pts, mes;
  Vector<double, 2> loc;
  double ori;
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
    scans.push_back({pts, mes, loc, 0.0});
  }
  return scans;
}

vector<Scan> import_scans(string folder) {
  vector<Scan> scans;
  vector<string> files;
  float fbuf;
  for (auto &entry : boost::make_iterator_range(
           fsm::directory_iterator(fsm::path(folder)), {}))
    files.push_back(entry.path().string());
  // iterate files <files.size() / 2
  for (unsigned i = 1; i < 2; i++) {
    Pose p;
    Meas m;
    std::ifstream ifs;
    string filename;
    // reading single scan, define postfix (.dat not included)
    char num[6];
    snprintf(num, 6, "%05d", i);
    // READING POSE FIRST
    filename = folder + "pose_" + string(num) + ".dat";
    ifs.open(filename, ios::in | ios::binary);
    ifs.read((char *)&fbuf, sizeof(fbuf)); // px
    p.p(0) = fbuf;
    ifs.read((char *)&fbuf, sizeof(fbuf)); // py
    p.p(1) = fbuf;
    ifs.read((char *)&fbuf, sizeof(fbuf)); // pz
    p.p(2) = fbuf;
    ifs.read((char *)&fbuf, sizeof(fbuf)); // ox
    p.o.x() = fbuf;
    ifs.read((char *)&fbuf, sizeof(fbuf)); // oy
    p.o.y() = fbuf;
    ifs.read((char *)&fbuf, sizeof(fbuf)); // oz
    p.o.z() = fbuf;
    ifs.read((char *)&fbuf, sizeof(fbuf)); // ow
    p.o.w() = fbuf;
    ifs.close();
    // READING SCAN DATA
    filename = folder + "data_" + string(num) + ".dat";
    cout << filename << endl;
    ifs.open(filename, ios::in | ios::binary);
    size_t sz; // read filesize written at the beginning (in bytes)
    ifs.read((char *)&sz, sizeof(sz));
    vector<Meas> ms;
    ms.resize(sz / 12);
    for (unsigned j = 0; j < sz / 12; j++) { // iterate (12bytes/mes)
      ifs.read((char *)&fbuf, sizeof(fbuf)); // theta
      m.t = fbuf;                            // normalize
      ifs.read((char *)&fbuf, sizeof(fbuf)); // phi
      m.p = fbuf;                            // normalize
      cout << m.p << ", ";
      ifs.read((char *)&fbuf, sizeof(fbuf)); // dist
      m.d = fbuf;
      ms.push_back(m);
    }
    ifs.close();
    cout << endl << "number of pts: " << ms.size() << endl;
    // NOW TRANSLATING INTO OUR SCAN STRUCTURE
    Scan s;
    vector<double> d, an, x, y;
    MatrixXd pts(ms.size() / 10, 2);
    MatrixXd mes(ms.size() / 10, 2);
    s.loc(0) = p.p(0); // location
    s.loc(1) = p.p(1);
    s.ori = p.o.toRotationMatrix().eulerAngles(2, 1, 0)(0); // yaw
    for (unsigned j = 0; j < ms.size(); j += 10) {
      // eliminate vertical angle
      d.push_back(ms[j].d * cos(ms[j].p));
      an.push_back(ms[j].t);
      // coords
      x.push_back(d.back() * cos(an.back()));
      y.push_back(d.back() * sin(an.back()));
    }
    // map vectors
    pts.col(0) = Map<VectorXd>(x.data(), x.size());
    pts.col(1) = Map<VectorXd>(y.data(), y.size());
    mes.col(0) = Map<VectorXd>(d.data(), d.size());
    mes.col(1) = Map<VectorXd>(an.data(), an.size());
    s.pts = pts;
    cout << mes << endl;
    s.mes = mes;
    scans.push_back(s);
  }
  // exit(0);
  return scans;
}
#endif // PREPROC_HPP
