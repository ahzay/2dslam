#ifndef DIFF_HPP
#define DIFF_HPP
#include "ls.hpp"
#include <Eigen/LU>
#include <autodiff/reverse/var.hpp>
namespace ad = autodiff;

template <typename T>
T ft2(T xc, T yc, T th, T a, T b, T e, T xp, T yp, T d, T an, T mud, T muan) {
  T x = xp + (d + mud) * cos(an + muan);
  T y = yp + (d + mud) * sin(an + muan);
  T f1 = ((x - xc) * cos(th) + (y - yc) * sin(th)) / a;
  T f2 = ((x - xc) * sin(th) - (y - yc) * cos(th)) / b;
  T buf = pow(pow(f1, 2.), 1. / e) + pow(pow(f2, 2.), 1. / e);
  return pow(buf, e) - 1.;
}

MatrixXd dftd(Vector<double, 6> p, Vector2d pos, MatrixX2d m) {
  ad::var xc(p(0)), yc(p(1)), th(p(2)), a(p(3)), b(p(4)), e(p(5)), xp(pos(0)),
      yp(pos(1)), d(0.), an(0.), mud(0.), muan(0.);
  ad::var f = ft2(xc, yc, th, a, b, e, xp, yp, d, an, mud, muan);

  MatrixXd ans(m.rows(), 12);
  for (unsigned i = 0; i < m.rows(); i++) {
    d.update(m.row(i)(0));
    an.update(m.row(i)(1));
    f.update();
    auto [dxc, dyc, dth, da, db, de, dxp, dyp, dd, dan, dmud, dmuan] =
        ad::derivatives(f,
                        ad::wrt(xc, yc, th, a, b, e, xp, yp, d, an, mud, muan));
    ans(i, 0) = dxc;
    ans(i, 1) = dyc;
    ans(i, 2) = dth;
    ans(i, 3) = da;
    ans(i, 4) = db;
    ans(i, 5) = de;
    ans(i, 6) = dxp;
    ans(i, 7) = dyp;
    ans(i, 8) = dd;
    ans(i, 9) = dan;
    ans(i, 10) = dmud;
    ans(i, 11) = dmuan;
  }
  return ans;
}

#endif // DIFF_HPP
