#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP
#include "ekf.hpp"

struct Aggregate {
  Matrix<double, Dynamic, 2> pts, mes;
  Matrix<double, 6, 6> E;
  Vector<double, 6> p;
  Vector<double, 2> loc;
  int idx = -1; // idx in ps global state vector
};

bool isContinuous(Matrix<double, Dynamic, 2> mes, double d_tol, double an_tol) {
  for (int i = 1; i < mes.rows(); i++) {
    double var_d = abs(mes(i, 0) - mes(i - 1, 0));
    double var_an = abs(remainder(mes(i, 1) - mes(i - 1, 1), M_PI / 2));
    if (var_d > d_tol || var_an > an_tol) {
      cout << "not continuous!" << var_d << " " << var_an << endl;
      return 0;
    }
  }
  return 1;
}

int associate_index(Scan scan, int b, int e, vector<Vector<double, 6>> &P,
                    double thres) {
  Matrix<double, 2, 2> W;
  W << 0.02, 0, 0, 0;
  double min_m = 100000;
  int min_i = -1;
  for (unsigned i = 0; i < P.size(); i++) {
    auto p = P[i];
    VectorXd r = fsn(p, scan.loc, scan.mes(seq(b, e), all));
    auto E = DOP(p, scan.loc, scan.mes(seq(b, e), all), 0.02);
    Matrix<double, Dynamic, 12> df =
        dftd(p, scan.loc, scan.mes(seq(b, e), all));
    Matrix<double, Dynamic, 6> H = df(all, {0, 1, 2, 3, 4, 5});
    Matrix<double, Dynamic, 2> J = df(all, {10, 11});
    Matrix<double, Dynamic, Dynamic> S =
        H * E * H.transpose() + J * W * J.transpose();
    ArrayXd m_dist = r.array() * r.array() / S.diagonal().array();
    double sum = m_dist.matrix().sum() / m_dist.size();
    cout << "sum: " << sum << endl;
    if (isnan(sum)) {
      cout << "isnan - S: " << S.diagonal().transpose();
    }
    if (abs(sum) < min_m) {
      min_m = abs(sum);
      min_i = i;
    }
  }
  if (min_m <= thres)
    return min_i;
  else
    return -1;
}

vector<Aggregate> segment_scan(Scan scan, vector<Vector<double, 6>> &P,
                               DiagonalMatrix<double, 6> Q,
                               Matrix<double, 2, 2> W, double tol_flx,
                               double tol_str, int npts, double tol_aug) {
  vector<Aggregate> ans;
  int i = 0, LScnt = 0; // measure counter
  int lasti, idx;       // i of start of last segment
  // temps
  Vector<double, 6> p_i = Vector<double, 6>::Zero();
  Matrix<double, 6, 6> En = Matrix<double, 6, 6>::Zero();
  // state machine
init:
  cout << "init! " << i << endl;
  lasti = i;
  if (i < scan.mes.rows() - npts) { // STATE MACHINE
    // basic continuous check
    if (!isContinuous(scan.mes(seq(i, i + npts), all), 0.5, 0.1)) {
      i += 1;
      goto init;
    }
    // association with P global vector
    idx = associate_index(scan, i, i + npts, P, tol_aug);
    cout << "associated with: " << idx << endl;
    if (idx > -1)
      p_i = P[idx];
    else
      p_i = LS(scan.pts(seq(i, i + npts), 0), scan.pts(seq(i, i + npts), 1),
               scan.loc);
    En = DOP(p_i, scan.loc, scan.mes, 0.02); // init En with DOP
    Ekf ekf(p_i, En, Q, W);
    goto flexible;
  flexible:
    cout << "flexible! " << i << endl;
    //  aberrant measure
    if (ekf.update(scan.loc, scan.mes.row(i), tol_flx, i)) {
      i += 1; // ITERATE
      goto init;
    } else {
      i += npts;
      goto strict;
    }
  strict: // break point maybe
    cout << "strict! " << i << endl;
    if (ekf.update(scan.loc, scan.mes.row(i), tol_str, i))
      goto LS;
    else {
      if (i >= scan.mes.rows() - 1) {
        LScnt++;
        cout << "LS! " << LScnt << endl;
        ekf._p = LS(scan.pts(seq(lasti, i), 0), scan.pts(seq(lasti, i), 1),
                    //         scan.loc, ekf._p.data());
                    scan.loc);
        ekf._E = DOP(ekf._p, scan.loc, scan.mes(seq(lasti, i), all), 0.02);
        goto close;
      }
      i += 1;
      goto strict;
    }

  LS:
    LScnt++;
    cout << "LS! " << LScnt << " i:" << i << endl;
    ekf._p = LS(scan.pts(seq(lasti, i - 2), 0), scan.pts(seq(lasti, i - 2), 1),
                //           scan.loc, ekf._p.data());
                scan.loc);
    ekf._E = DOP(ekf._p, scan.loc, scan.mes(seq(lasti, i - 2), all), 0.02);
    if (ekf.update(scan.loc, scan.mes.row(i), tol_str, i) ||
        i >= scan.mes.rows() - 1)
      goto close;
    else {
      i += 1;
      goto strict;
    }
  close:
    cout << "close! ls count: " << LScnt << " i:" << i << endl;
    ans.push_back(Aggregate());
    ans.back().p = ekf._p;
    ans.back().E = ekf._E;
    ans.back().loc = scan.loc;
    ans.back().pts = scan.pts(seq(lasti, i - 1), all);
    ans.back().mes = scan.mes(seq(lasti, i - 1), all);
    ans.back().idx = idx;
    cout << ekf._p.transpose() << endl;
    goto init;
  }
  return ans;
}

void augment_scan(Aggregate a, vector<Vector<double, 6>> &P,
                  DiagonalMatrix<double, 6> Q, Matrix<double, 2, 2> W) {
  Iekf iekf(P[a.idx], a.E, Q, W);
  for (unsigned i = 0; i < a.mes.rows(); i++) {
    iekf.update(a.loc, a.mes(i, all).transpose());
  }
  P[a.idx] = iekf._p;
  cout << iekf._p.transpose() << endl;
}

#endif // SEGMENTATION_HPP
