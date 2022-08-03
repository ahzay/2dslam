#ifndef EKF_HPP
#define EKF_HPP
#include "dop.hpp"
#include <time.h>
class Ekf {
public:
  // E is Sigma
  Ekf(const VectorXd &p, const MatrixXd &E, const MatrixXd &Q,
      const MatrixXd &W)
      : _p(p), _E(E), _Q(Q), _W(W) {}

  // i for debugging purposes
  int update(const Vector2d &pos, const Vector2d &m, const double tol,
             const int i) {
    // system propagation
    _E = _E + _Q;
    _r = -fs<double>(_p, pos, m);
    MatrixX2d _m(1, 2);
    _m << m(0), m(1);
    Matrix<double, Dynamic, 12> df = dftd(_p, pos, _m);
    H = df(all, {0, 1, 2, 3, 4, 5});
    J = df(all, {10, 11});
    S = H * _E * H.transpose() + J * _W * J.transpose();
    K = _E * H.transpose() / S(0);
    // m_dist
    m_dist = _r * _r / S(0);
    cout << i << " mdist: " << m_dist << endl;
    if (m_dist > tol)
      return 1;
    // update
    _p = _p + K * _r;
    _E = (I - K * H) * _E;
    return 0; // no bp or aberrant measure
  }

  // private:
  Vector<double, 6> _p;
  Matrix<double, 6, 6> _E;
  Matrix<double, 2, 2> _W;
  double _r;
  const Matrix<double, 6, 6> _Q;
  Matrix<double, 6, 6> I = Matrix<double, 6, 6>::Identity(6, 6);
  // uninit params
  Matrix<double, 1, 6> H = Matrix<double, 1, 6>::Zero(1, 6);
  Matrix<double, 1, 2> J = Matrix<double, 1, 2>::Zero(1, 2);
  Matrix<double, 1, 1> S = Matrix<double, 1, 1>::Zero(1, 1);
  Matrix<double, 6, 1> K = Matrix<double, 6, 1>::Zero(6, 1);
  double m_dist;
};

class Iekf {
public:
  // E is Sigma
  Iekf(const VectorXd &p, const MatrixXd &E, const MatrixXd &Q,
       const MatrixXd &W)
      : _p(p), _E(E), _Q(Q), _W(W) {}

  int update(const Vector2d &pos, const Matrix<double, 2, 1> m0) {
    // system propagation
    _E = _E + _Q;
    mk = m0;
    pk = _p;
    Matrix<double, Dynamic, 12> df;
    for (int j = 0; j < 10; j++) {
      // jacobians
      df = dftd(pk, pos, mk.transpose());
      H = df(all, {0, 1, 2, 3, 4, 5});
      J = df(all, {10, 11}); // dfdmu ? dfdmes?
      //
      r = fs<double>(pk, pos, mk);
      // cout << "r: " << r << endl;
      dr = (J * (m0 - mk) + H * (_p - pk))(0);
      // cout << "dr: " << dr << endl;
      S = (H * _E * H.transpose() + J * _W * J.transpose())(0);
      if ((_E * H.transpose() / S * (r + dr)).array().isNaN().sum())
        break;
      K = _E * H.transpose() / S;
      // cout << "K: " << K.transpose() << endl;
      pk = _p - K * (r + dr);
      mk = m0 - (_W * J.transpose() / S) * (r + dr);
      // cout << "mk: " << mk.transpose() << endl;
      // cout << "pk: " << pk.transpose() << endl;
      // cout << "sum: " << _E.sum() << endl << endl;
      if (abs(r) < 1)
        break;
      // usleep(200e3);
    }
    // << endl
    //     <<
    //     "break!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    //       "!!!!!!!!!!!"
    //    << endl;
    // testing
    if (pk(5) < 0.1)
      pk(5) = 0.1;
    if (pk(5) > 1.99)
      pk(5) = 1.99;
    _p = pk;
    // jacs
    df = dftd(_p, pos, mk.transpose());
    H = df(all, {0, 1, 2, 3, 4, 5});
    J = df(all, {10, 11});
    //
    L = I - K * H;
    _E = L * _E * L.transpose() + K * J * _W * J.transpose() * K.transpose();
    return 0;
  }

  // private:
  // Matrix<double, Dynamic, 12> df;
  Matrix<double, 2, 1> mk;
  Vector<double, 6> _p, pk;
  Matrix<double, 6, 6> _E, L;
  Matrix<double, 2, 2> _W;
  double r, dr, S;
  const Matrix<double, 6, 6> _Q;
  Matrix<double, 6, 6> I = Matrix<double, 6, 6>::Identity(6, 6);
  // uninit params
  Matrix<double, 1, 6> H = Matrix<double, 1, 6>::Zero(1, 6);
  Matrix<double, 1, 2> J = Matrix<double, 1, 2>::Zero(1, 2);
  Matrix<double, 6, 1> K = Matrix<double, 6, 1>::Zero(6, 1);
  double m_dist;
};
#endif // EKF_HPP
