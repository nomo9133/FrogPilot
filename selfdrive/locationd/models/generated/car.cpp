#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1506063965276717887) {
   out_1506063965276717887[0] = delta_x[0] + nom_x[0];
   out_1506063965276717887[1] = delta_x[1] + nom_x[1];
   out_1506063965276717887[2] = delta_x[2] + nom_x[2];
   out_1506063965276717887[3] = delta_x[3] + nom_x[3];
   out_1506063965276717887[4] = delta_x[4] + nom_x[4];
   out_1506063965276717887[5] = delta_x[5] + nom_x[5];
   out_1506063965276717887[6] = delta_x[6] + nom_x[6];
   out_1506063965276717887[7] = delta_x[7] + nom_x[7];
   out_1506063965276717887[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3364971332114261188) {
   out_3364971332114261188[0] = -nom_x[0] + true_x[0];
   out_3364971332114261188[1] = -nom_x[1] + true_x[1];
   out_3364971332114261188[2] = -nom_x[2] + true_x[2];
   out_3364971332114261188[3] = -nom_x[3] + true_x[3];
   out_3364971332114261188[4] = -nom_x[4] + true_x[4];
   out_3364971332114261188[5] = -nom_x[5] + true_x[5];
   out_3364971332114261188[6] = -nom_x[6] + true_x[6];
   out_3364971332114261188[7] = -nom_x[7] + true_x[7];
   out_3364971332114261188[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4094953735274879016) {
   out_4094953735274879016[0] = 1.0;
   out_4094953735274879016[1] = 0;
   out_4094953735274879016[2] = 0;
   out_4094953735274879016[3] = 0;
   out_4094953735274879016[4] = 0;
   out_4094953735274879016[5] = 0;
   out_4094953735274879016[6] = 0;
   out_4094953735274879016[7] = 0;
   out_4094953735274879016[8] = 0;
   out_4094953735274879016[9] = 0;
   out_4094953735274879016[10] = 1.0;
   out_4094953735274879016[11] = 0;
   out_4094953735274879016[12] = 0;
   out_4094953735274879016[13] = 0;
   out_4094953735274879016[14] = 0;
   out_4094953735274879016[15] = 0;
   out_4094953735274879016[16] = 0;
   out_4094953735274879016[17] = 0;
   out_4094953735274879016[18] = 0;
   out_4094953735274879016[19] = 0;
   out_4094953735274879016[20] = 1.0;
   out_4094953735274879016[21] = 0;
   out_4094953735274879016[22] = 0;
   out_4094953735274879016[23] = 0;
   out_4094953735274879016[24] = 0;
   out_4094953735274879016[25] = 0;
   out_4094953735274879016[26] = 0;
   out_4094953735274879016[27] = 0;
   out_4094953735274879016[28] = 0;
   out_4094953735274879016[29] = 0;
   out_4094953735274879016[30] = 1.0;
   out_4094953735274879016[31] = 0;
   out_4094953735274879016[32] = 0;
   out_4094953735274879016[33] = 0;
   out_4094953735274879016[34] = 0;
   out_4094953735274879016[35] = 0;
   out_4094953735274879016[36] = 0;
   out_4094953735274879016[37] = 0;
   out_4094953735274879016[38] = 0;
   out_4094953735274879016[39] = 0;
   out_4094953735274879016[40] = 1.0;
   out_4094953735274879016[41] = 0;
   out_4094953735274879016[42] = 0;
   out_4094953735274879016[43] = 0;
   out_4094953735274879016[44] = 0;
   out_4094953735274879016[45] = 0;
   out_4094953735274879016[46] = 0;
   out_4094953735274879016[47] = 0;
   out_4094953735274879016[48] = 0;
   out_4094953735274879016[49] = 0;
   out_4094953735274879016[50] = 1.0;
   out_4094953735274879016[51] = 0;
   out_4094953735274879016[52] = 0;
   out_4094953735274879016[53] = 0;
   out_4094953735274879016[54] = 0;
   out_4094953735274879016[55] = 0;
   out_4094953735274879016[56] = 0;
   out_4094953735274879016[57] = 0;
   out_4094953735274879016[58] = 0;
   out_4094953735274879016[59] = 0;
   out_4094953735274879016[60] = 1.0;
   out_4094953735274879016[61] = 0;
   out_4094953735274879016[62] = 0;
   out_4094953735274879016[63] = 0;
   out_4094953735274879016[64] = 0;
   out_4094953735274879016[65] = 0;
   out_4094953735274879016[66] = 0;
   out_4094953735274879016[67] = 0;
   out_4094953735274879016[68] = 0;
   out_4094953735274879016[69] = 0;
   out_4094953735274879016[70] = 1.0;
   out_4094953735274879016[71] = 0;
   out_4094953735274879016[72] = 0;
   out_4094953735274879016[73] = 0;
   out_4094953735274879016[74] = 0;
   out_4094953735274879016[75] = 0;
   out_4094953735274879016[76] = 0;
   out_4094953735274879016[77] = 0;
   out_4094953735274879016[78] = 0;
   out_4094953735274879016[79] = 0;
   out_4094953735274879016[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_1631448948235285413) {
   out_1631448948235285413[0] = state[0];
   out_1631448948235285413[1] = state[1];
   out_1631448948235285413[2] = state[2];
   out_1631448948235285413[3] = state[3];
   out_1631448948235285413[4] = state[4];
   out_1631448948235285413[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1631448948235285413[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1631448948235285413[7] = state[7];
   out_1631448948235285413[8] = state[8];
}
void F_fun(double *state, double dt, double *out_518627117871004420) {
   out_518627117871004420[0] = 1;
   out_518627117871004420[1] = 0;
   out_518627117871004420[2] = 0;
   out_518627117871004420[3] = 0;
   out_518627117871004420[4] = 0;
   out_518627117871004420[5] = 0;
   out_518627117871004420[6] = 0;
   out_518627117871004420[7] = 0;
   out_518627117871004420[8] = 0;
   out_518627117871004420[9] = 0;
   out_518627117871004420[10] = 1;
   out_518627117871004420[11] = 0;
   out_518627117871004420[12] = 0;
   out_518627117871004420[13] = 0;
   out_518627117871004420[14] = 0;
   out_518627117871004420[15] = 0;
   out_518627117871004420[16] = 0;
   out_518627117871004420[17] = 0;
   out_518627117871004420[18] = 0;
   out_518627117871004420[19] = 0;
   out_518627117871004420[20] = 1;
   out_518627117871004420[21] = 0;
   out_518627117871004420[22] = 0;
   out_518627117871004420[23] = 0;
   out_518627117871004420[24] = 0;
   out_518627117871004420[25] = 0;
   out_518627117871004420[26] = 0;
   out_518627117871004420[27] = 0;
   out_518627117871004420[28] = 0;
   out_518627117871004420[29] = 0;
   out_518627117871004420[30] = 1;
   out_518627117871004420[31] = 0;
   out_518627117871004420[32] = 0;
   out_518627117871004420[33] = 0;
   out_518627117871004420[34] = 0;
   out_518627117871004420[35] = 0;
   out_518627117871004420[36] = 0;
   out_518627117871004420[37] = 0;
   out_518627117871004420[38] = 0;
   out_518627117871004420[39] = 0;
   out_518627117871004420[40] = 1;
   out_518627117871004420[41] = 0;
   out_518627117871004420[42] = 0;
   out_518627117871004420[43] = 0;
   out_518627117871004420[44] = 0;
   out_518627117871004420[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_518627117871004420[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_518627117871004420[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_518627117871004420[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_518627117871004420[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_518627117871004420[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_518627117871004420[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_518627117871004420[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_518627117871004420[53] = -9.8000000000000007*dt;
   out_518627117871004420[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_518627117871004420[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_518627117871004420[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_518627117871004420[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_518627117871004420[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_518627117871004420[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_518627117871004420[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_518627117871004420[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_518627117871004420[62] = 0;
   out_518627117871004420[63] = 0;
   out_518627117871004420[64] = 0;
   out_518627117871004420[65] = 0;
   out_518627117871004420[66] = 0;
   out_518627117871004420[67] = 0;
   out_518627117871004420[68] = 0;
   out_518627117871004420[69] = 0;
   out_518627117871004420[70] = 1;
   out_518627117871004420[71] = 0;
   out_518627117871004420[72] = 0;
   out_518627117871004420[73] = 0;
   out_518627117871004420[74] = 0;
   out_518627117871004420[75] = 0;
   out_518627117871004420[76] = 0;
   out_518627117871004420[77] = 0;
   out_518627117871004420[78] = 0;
   out_518627117871004420[79] = 0;
   out_518627117871004420[80] = 1;
}
void h_25(double *state, double *unused, double *out_621159288242216091) {
   out_621159288242216091[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8707689189734903757) {
   out_8707689189734903757[0] = 0;
   out_8707689189734903757[1] = 0;
   out_8707689189734903757[2] = 0;
   out_8707689189734903757[3] = 0;
   out_8707689189734903757[4] = 0;
   out_8707689189734903757[5] = 0;
   out_8707689189734903757[6] = 1;
   out_8707689189734903757[7] = 0;
   out_8707689189734903757[8] = 0;
}
void h_24(double *state, double *unused, double *out_3522597052295566447) {
   out_3522597052295566447[0] = state[4];
   out_3522597052295566447[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3434194906821891163) {
   out_3434194906821891163[0] = 0;
   out_3434194906821891163[1] = 0;
   out_3434194906821891163[2] = 0;
   out_3434194906821891163[3] = 0;
   out_3434194906821891163[4] = 1;
   out_3434194906821891163[5] = 0;
   out_3434194906821891163[6] = 0;
   out_3434194906821891163[7] = 0;
   out_3434194906821891163[8] = 0;
   out_3434194906821891163[9] = 0;
   out_3434194906821891163[10] = 0;
   out_3434194906821891163[11] = 0;
   out_3434194906821891163[12] = 0;
   out_3434194906821891163[13] = 0;
   out_3434194906821891163[14] = 1;
   out_3434194906821891163[15] = 0;
   out_3434194906821891163[16] = 0;
   out_3434194906821891163[17] = 0;
}
void h_30(double *state, double *unused, double *out_7260113362511833094) {
   out_7260113362511833094[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6189356231227655130) {
   out_6189356231227655130[0] = 0;
   out_6189356231227655130[1] = 0;
   out_6189356231227655130[2] = 0;
   out_6189356231227655130[3] = 0;
   out_6189356231227655130[4] = 1;
   out_6189356231227655130[5] = 0;
   out_6189356231227655130[6] = 0;
   out_6189356231227655130[7] = 0;
   out_6189356231227655130[8] = 0;
}
void h_26(double *state, double *unused, double *out_5988695097038667037) {
   out_5988695097038667037[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5997551565100591635) {
   out_5997551565100591635[0] = 0;
   out_5997551565100591635[1] = 0;
   out_5997551565100591635[2] = 0;
   out_5997551565100591635[3] = 0;
   out_5997551565100591635[4] = 0;
   out_5997551565100591635[5] = 0;
   out_5997551565100591635[6] = 0;
   out_5997551565100591635[7] = 1;
   out_5997551565100591635[8] = 0;
}
void h_27(double *state, double *unused, double *out_6604083127548357844) {
   out_6604083127548357844[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3965762160043711913) {
   out_3965762160043711913[0] = 0;
   out_3965762160043711913[1] = 0;
   out_3965762160043711913[2] = 0;
   out_3965762160043711913[3] = 1;
   out_3965762160043711913[4] = 0;
   out_3965762160043711913[5] = 0;
   out_3965762160043711913[6] = 0;
   out_3965762160043711913[7] = 0;
   out_3965762160043711913[8] = 0;
}
void h_29(double *state, double *unused, double *out_4079395689340521331) {
   out_4079395689340521331[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5679124886913262946) {
   out_5679124886913262946[0] = 0;
   out_5679124886913262946[1] = 1;
   out_5679124886913262946[2] = 0;
   out_5679124886913262946[3] = 0;
   out_5679124886913262946[4] = 0;
   out_5679124886913262946[5] = 0;
   out_5679124886913262946[6] = 0;
   out_5679124886913262946[7] = 0;
   out_5679124886913262946[8] = 0;
}
void h_28(double *state, double *unused, double *out_8805874894716861087) {
   out_8805874894716861087[0] = state[0];
}
void H_28(double *state, double *unused, double *out_7685220169726758096) {
   out_7685220169726758096[0] = 1;
   out_7685220169726758096[1] = 0;
   out_7685220169726758096[2] = 0;
   out_7685220169726758096[3] = 0;
   out_7685220169726758096[4] = 0;
   out_7685220169726758096[5] = 0;
   out_7685220169726758096[6] = 0;
   out_7685220169726758096[7] = 0;
   out_7685220169726758096[8] = 0;
}
void h_31(double *state, double *unused, double *out_2979304202071846475) {
   out_2979304202071846475[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5371343462867240159) {
   out_5371343462867240159[0] = 0;
   out_5371343462867240159[1] = 0;
   out_5371343462867240159[2] = 0;
   out_5371343462867240159[3] = 0;
   out_5371343462867240159[4] = 0;
   out_5371343462867240159[5] = 0;
   out_5371343462867240159[6] = 0;
   out_5371343462867240159[7] = 0;
   out_5371343462867240159[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1506063965276717887) {
  err_fun(nom_x, delta_x, out_1506063965276717887);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3364971332114261188) {
  inv_err_fun(nom_x, true_x, out_3364971332114261188);
}
void car_H_mod_fun(double *state, double *out_4094953735274879016) {
  H_mod_fun(state, out_4094953735274879016);
}
void car_f_fun(double *state, double dt, double *out_1631448948235285413) {
  f_fun(state,  dt, out_1631448948235285413);
}
void car_F_fun(double *state, double dt, double *out_518627117871004420) {
  F_fun(state,  dt, out_518627117871004420);
}
void car_h_25(double *state, double *unused, double *out_621159288242216091) {
  h_25(state, unused, out_621159288242216091);
}
void car_H_25(double *state, double *unused, double *out_8707689189734903757) {
  H_25(state, unused, out_8707689189734903757);
}
void car_h_24(double *state, double *unused, double *out_3522597052295566447) {
  h_24(state, unused, out_3522597052295566447);
}
void car_H_24(double *state, double *unused, double *out_3434194906821891163) {
  H_24(state, unused, out_3434194906821891163);
}
void car_h_30(double *state, double *unused, double *out_7260113362511833094) {
  h_30(state, unused, out_7260113362511833094);
}
void car_H_30(double *state, double *unused, double *out_6189356231227655130) {
  H_30(state, unused, out_6189356231227655130);
}
void car_h_26(double *state, double *unused, double *out_5988695097038667037) {
  h_26(state, unused, out_5988695097038667037);
}
void car_H_26(double *state, double *unused, double *out_5997551565100591635) {
  H_26(state, unused, out_5997551565100591635);
}
void car_h_27(double *state, double *unused, double *out_6604083127548357844) {
  h_27(state, unused, out_6604083127548357844);
}
void car_H_27(double *state, double *unused, double *out_3965762160043711913) {
  H_27(state, unused, out_3965762160043711913);
}
void car_h_29(double *state, double *unused, double *out_4079395689340521331) {
  h_29(state, unused, out_4079395689340521331);
}
void car_H_29(double *state, double *unused, double *out_5679124886913262946) {
  H_29(state, unused, out_5679124886913262946);
}
void car_h_28(double *state, double *unused, double *out_8805874894716861087) {
  h_28(state, unused, out_8805874894716861087);
}
void car_H_28(double *state, double *unused, double *out_7685220169726758096) {
  H_28(state, unused, out_7685220169726758096);
}
void car_h_31(double *state, double *unused, double *out_2979304202071846475) {
  h_31(state, unused, out_2979304202071846475);
}
void car_H_31(double *state, double *unused, double *out_5371343462867240159) {
  H_31(state, unused, out_5371343462867240159);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
