#ifndef DOP_HPP
#define DOP_HPP
#include "diff.hpp"
// wrapper for use in EKF
template <typename T>
T fs(const Vector<double, 6> &p, const Vector2d &pos, const Vector2d &m) {
  // return ft<T, T, T, T, T, T, T, T, T, T, T, T>(xc, yc, th, a, b, e, xp, yp,
  // d,
  //                                              an, 0, 0);
  return ft2<T>(p(0), p(1), p(2), p(3), p(4), p(5), pos(0), pos(1), m(0), m(1),
                0, 0);
}
// for Jmes->Emes
MatrixX2d dftdmu(VectorXd p, Vector2d pos, MatrixX2d m) {
  MatrixX2d ans(m.rows(), 2);
  for (unsigned i = 0; i < m.rows(); i++) {
    ans(i, 0) = dftdn(p, pos, m.row(i), 10); // dftd(mud)
    ans(i, 1) = dftdn(p, pos, m.row(i), 11); // dftd(muan)
  }
  return ans;
}
// for H
MatrixXd dftdp(VectorXd p, Vector2d pos, MatrixX2d m) {
  MatrixXd ans(m.rows(), 6);
  for (unsigned i = 0; i < m.rows(); i++) {
    ans(i, 0) = dftdn(p, pos, m.row(i), 0); // dftd(xc)
    ans(i, 1) = dftdn(p, pos, m.row(i), 1); // dftd(yc)
    ans(i, 2) = dftdn(p, pos, m.row(i), 2); // dftd(phi)
    ans(i, 3) = dftdn(p, pos, m.row(i), 3); // dftd(a)
    ans(i, 4) = dftdn(p, pos, m.row(i), 4); // dftd(b)
    ans(i, 5) = dftdn(p, pos, m.row(i), 5); // dftd(eps)
  }
  return ans;
}

Matrix<double, 6, 6> DOP(VectorXd p, Vector2d pos, MatrixX2d m,
                         double sigma_d) {
  Matrix<double, Dynamic, 12> df = dftd(p, pos, m);
  MatrixX2d Jmes = df(all, {10, 11});
  MatrixXd Emes = Jmes * sigma_d * sigma_d * Jmes.transpose();
  Matrix<double, Dynamic, 6> H = df(all, {0, 1, 2, 3, 4, 5});
  Matrix<double, 6, 6> En =
      (H.transpose() * Emes.completeOrthogonalDecomposition().pseudoInverse() *
       H)
          .completeOrthogonalDecomposition()
          .pseudoInverse();
  return En;
}

#endif // DOP_HPP
