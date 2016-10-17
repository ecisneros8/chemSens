#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cppad/cppad.hpp>

using namespace std;

using std::vector;
using CppAD::AD;
using CppAD::Value;
using CppAD::Var2Par;

extern "C" {
  void nspecies_(int& ns);
  void nreactions_(int& nr);
  void nelements_(int& ne);
  void getmolecularweights_(double* mw);
  void getmassfractions_(double* x, double* y);
  void getenthalpymass_(double* T, double* y, double* h);
  void getdensity_(double* T, double* p, double* y, double* rho);
  void gettemperature_(double* h, double* y, double* Told, double* Tnew);
  void gettemperaturesource_(double* p, double* h, double* T, double* y, double* fT);
  void getsource_(double* p, double* h, double* Told, double* y, double* fy);
  void getjacobian_(double* p, double* h, double* Told, double* y, double* dfy);
  void readreactionrates_(char* filename);
  void correctname(char* fname);
}

namespace mech {

  /***********************************************************************************/
  /* some constants                                                                  */
  /***********************************************************************************/
  double GasConstant = 8314.4621;
  double OneAtm      = 1.01325e5;
  double OneThird    = 1.0/3.0;

  /***********************************************************************************/
  /* mech definitions                                                                */
  /***********************************************************************************/
  int kk = 9;
  int ii = 24;
  int mm = 3;

  struct Kfentry {
    double p1, p2, p3;
  };
  vector<Kfentry> kfwd_coeffs;

  vector<double> getMolecularWeights() {

    vector<double> mw(kk,0.0);
    mw[0] = 2.015880e+00;
    mw[1] = 1.007940e+00;
    mw[2] = 3.199880e+01;
    mw[3] = 1.599940e+01;
    mw[4] = 1.700734e+01;
    mw[5] = 3.300674e+01;
    mw[6] = 3.401468e+01;
    mw[7] = 1.801528e+01;
    mw[8] = 2.801348e+01;

    return(mw);
    
  };

  /***********************************************************************************/
  /* thermodynamics                                                                  */
  /***********************************************************************************/
  template <class Type>
  void getSpecificHeats_R(Type& T, vector<Type>& cp0_R) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      cp0_R[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 + 2.015720940000000e-08 * tt2 + -7.376117610000001e-12 * tt3;
    } else {
      cp0_R[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 + -1.795663940000000e-10 * tt2 + 2.002553760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 + 2.300816320000000e-18 * tt2 + -9.277323320000001e-22 * tt3;
    } else {
      cp0_R[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 + -4.735152350000000e-18 * tt2 + 4.981973570000000e-22 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[2] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 + -9.681295090000001e-09 * tt2 + 3.243728370000000e-12 * tt3;
    } else {
      cp0_R[2] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 + 2.094705550000000e-10 * tt2 + -2.167177940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[3] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 + -6.128066240000000e-09 * tt2 + 2.112659710000000e-12 * tt3;
    } else {
      cp0_R[3] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 + -1.001777990000000e-11 * tt2 + 1.228336910000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[4] = 4.125305610000000e+00 + -3.225449390000000e-03 * tt0 + 6.527646910000000e-06 * tt1 + -5.798536430000000e-09 * tt2 + 2.062373790000000e-12 * tt3;
    } else {
      cp0_R[4] = 2.864728860000000e+00 + 1.056504480000000e-03 * tt0 + -2.590827580000000e-07 * tt1 + 3.052186740000000e-11 * tt2 + -1.331958760000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[5] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 + -2.427638940000000e-08 * tt2 + 9.292251240000000e-12 * tt3;
    } else {
      cp0_R[5] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 + 1.142463700000000e-10 * tt2 + -1.079085350000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[6] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 + -2.157708130000000e-08 * tt2 + 8.624543630000000e-12 * tt3;
    } else {
      cp0_R[6] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 + 3.711859860000000e-10 * tt2 + -2.879083050000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[7] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 + -5.487970620000000e-09 * tt2 + 1.771978170000000e-12 * tt3;
    } else {
      cp0_R[7] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 + -9.704198700000000e-11 * tt2 + 1.682009920000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      cp0_R[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };


  };

  template <class Type>
  void getEnthalpies_RT(Type& T, vector<Type>& h0_RT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;
    Type tt4 = 1.0 / T;

    if(tt0 < 1.0000e+03) {
      h0_RT[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 * 0.50 + -1.947815100000000e-05 * tt1 * OneThird + 2.015720940000000e-08 * tt2 * 0.25 + -7.376117610000001e-12 * tt3 * 0.20 + -9.179351730000000e+02 * tt4;
    } else {
      h0_RT[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 * 0.50 + 4.994567780000000e-07 * tt1 * OneThird + -1.795663940000000e-10 * tt2 * 0.25 + 2.002553760000000e-14 * tt3 * 0.20 + -9.501589220000000e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 * 0.50 + -1.995919640000000e-15 * tt1 * OneThird + 2.300816320000000e-18 * tt2 * 0.25 + -9.277323320000001e-22 * tt3 * 0.20 + 2.547365990000000e+04 * tt4;
    } else {
      h0_RT[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 * 0.50 + 1.615619480000000e-14 * tt1 * OneThird + -4.735152350000000e-18 * tt2 * 0.25 + 4.981973570000000e-22 * tt3 * 0.20 + 2.547365990000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[2] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 * 0.50 + 9.847302010000000e-06 * tt1 * OneThird + -9.681295090000001e-09 * tt2 * 0.25 + 3.243728370000000e-12 * tt3 * 0.20 + -1.063943560000000e+03 * tt4;
    } else {
      h0_RT[2] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 * 0.50 + -7.579666690000000e-07 * tt1 * OneThird + 2.094705550000000e-10 * tt2 * 0.25 + -2.167177940000000e-14 * tt3 * 0.20 + -1.088457720000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[3] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 * 0.50 + 6.643063960000000e-06 * tt1 * OneThird + -6.128066240000000e-09 * tt2 * 0.25 + 2.112659710000000e-12 * tt3 * 0.20 + 2.912225920000000e+04 * tt4;
    } else {
      h0_RT[3] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 * 0.50 + 4.194845890000000e-08 * tt1 * OneThird + -1.001777990000000e-11 * tt2 * 0.25 + 1.228336910000000e-15 * tt3 * 0.20 + 2.921757910000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[4] = 4.125305610000000e+00 + -3.225449390000000e-03 * tt0 * 0.50 + 6.527646910000000e-06 * tt1 * OneThird + -5.798536430000000e-09 * tt2 * 0.25 + 2.062373790000000e-12 * tt3 * 0.20 + 3.381538120000000e+03 * tt4;
    } else {
      h0_RT[4] = 2.864728860000000e+00 + 1.056504480000000e-03 * tt0 * 0.50 + -2.590827580000000e-07 * tt1 * OneThird + 3.052186740000000e-11 * tt2 * 0.25 + -1.331958760000000e-15 * tt3 * 0.20 + 3.718857740000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[5] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 * 0.50 + 2.115828910000000e-05 * tt1 * OneThird + -2.427638940000000e-08 * tt2 * 0.25 + 9.292251240000000e-12 * tt3 * 0.20 + 2.948080400000000e+02 * tt4;
    } else {
      h0_RT[5] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 * 0.50 + -6.336581500000000e-07 * tt1 * OneThird + 1.142463700000000e-10 * tt2 * 0.25 + -1.079085350000000e-14 * tt3 * 0.20 + 1.118567130000000e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[6] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 * 0.50 + 1.673357010000000e-05 * tt1 * OneThird + -2.157708130000000e-08 * tt2 * 0.25 + 8.624543630000000e-12 * tt3 * 0.20 + -1.770258210000000e+04 * tt4;
    } else {
      h0_RT[6] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 * 0.50 + -1.901392250000000e-06 * tt1 * OneThird + 3.711859860000000e-10 * tt2 * 0.25 + -2.879083050000000e-14 * tt3 * 0.20 + -1.786178770000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[7] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 * 0.50 + 6.520402110000000e-06 * tt1 * OneThird + -5.487970620000000e-09 * tt2 * 0.25 + 1.771978170000000e-12 * tt3 * 0.20 + -3.029372670000000e+04 * tt4;
    } else {
      h0_RT[7] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 * 0.50 + -1.640725180000000e-07 * tt1 * OneThird + -9.704198700000000e-11 * tt2 * 0.25 + 1.682009920000000e-14 * tt3 * 0.20 + -3.000429710000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 * 0.50 + -3.964481129109908e-06 * tt1 * OneThird + 5.642920880408571e-09 * tt2 * 0.25 + -2.445407041148433e-12 * tt3 * 0.20 + -1.020894198687962e+03 * tt4;
    } else {
      h0_RT[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 * 0.50 + -5.684761849244810e-07 * tt1 * OneThird + 1.009704225872734e-10 * tt2 * 0.25 + -6.753354387142974e-15 * tt3 * 0.20 + -9.227966980051905e+02 * tt4;
    };

  };

  template <class Type>
  void getEnthalpiesDerivatives(Type& T, vector<Type>& dh0dT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      dh0dT[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 + 2.015720940000000e-08 * tt2 + -7.376117610000001e-12 * tt3;
    } else {
      dh0dT[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 + -1.795663940000000e-10 * tt2 + 2.002553760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 + 2.300816320000000e-18 * tt2 + -9.277323320000001e-22 * tt3;
    } else {
      dh0dT[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 + -4.735152350000000e-18 * tt2 + 4.981973570000000e-22 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[2] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 + -9.681295090000001e-09 * tt2 + 3.243728370000000e-12 * tt3;
    } else {
      dh0dT[2] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 + 2.094705550000000e-10 * tt2 + -2.167177940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[3] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 + -6.128066240000000e-09 * tt2 + 2.112659710000000e-12 * tt3;
    } else {
      dh0dT[3] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 + -1.001777990000000e-11 * tt2 + 1.228336910000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[4] = 4.125305610000000e+00 + -3.225449390000000e-03 * tt0 + 6.527646910000000e-06 * tt1 + -5.798536430000000e-09 * tt2 + 2.062373790000000e-12 * tt3;
    } else {
      dh0dT[4] = 2.864728860000000e+00 + 1.056504480000000e-03 * tt0 + -2.590827580000000e-07 * tt1 + 3.052186740000000e-11 * tt2 + -1.331958760000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[5] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 + -2.427638940000000e-08 * tt2 + 9.292251240000000e-12 * tt3;
    } else {
      dh0dT[5] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 + 1.142463700000000e-10 * tt2 + -1.079085350000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[6] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 + -2.157708130000000e-08 * tt2 + 8.624543630000000e-12 * tt3;
    } else {
      dh0dT[6] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 + 3.711859860000000e-10 * tt2 + -2.879083050000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[7] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 + -5.487970620000000e-09 * tt2 + 1.771978170000000e-12 * tt3;
    } else {
      dh0dT[7] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 + -9.704198700000000e-11 * tt2 + 1.682009920000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      dh0dT[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };

  };

  template <class Type>
  void getEntropies_R(Type& T, vector<Type>& s0_R) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;
    Type tt4 = 1.0 / T;
    Type tt5 = log(T);

    if(tt0 < 1.0000e+03) {
      s0_R[0] = 2.344331120000000e+00 * tt5 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 * 0.50 + 2.015720940000000e-08 * tt2 * OneThird + -7.376117610000001e-12 * tt3 * 0.25 + 6.830102380000000e-01;
    } else {
      s0_R[0] = 3.337279200000000e+00 * tt5 +  -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 * 0.50 + -1.795663940000000e-10 * tt2 * OneThird + 2.002553760000000e-14 * tt3 * 0.25 + -3.205023310000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[1] = 2.500000000000000e+00 * tt5 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 * 0.50 + 2.300816320000000e-18 * tt2 * OneThird + -9.277323320000001e-22 * tt3 * 0.25 + -4.466828530000000e-01;
    } else {
      s0_R[1] = 2.500000010000000e+00 * tt5 +  -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 * 0.50 + -4.735152350000000e-18 * tt2 * OneThird + 4.981973570000000e-22 * tt3 * 0.25 + -4.466829140000000e-01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[2] = 3.782456360000000e+00 * tt5 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 * 0.50 + -9.681295090000001e-09 * tt2 * OneThird + 3.243728370000000e-12 * tt3 * 0.25 + 3.657675730000000e+00;
    } else {
      s0_R[2] = 3.282537840000000e+00 * tt5 +  1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 * 0.50 + 2.094705550000000e-10 * tt2 * OneThird + -2.167177940000000e-14 * tt3 * 0.25 + 5.453231290000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[3] = 3.168267100000000e+00 * tt5 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 * 0.50 + -6.128066240000000e-09 * tt2 * OneThird + 2.112659710000000e-12 * tt3 * 0.25 + 2.051933460000000e+00;
    } else {
      s0_R[3] = 2.569420780000000e+00 * tt5 +  -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 * 0.50 + -1.001777990000000e-11 * tt2 * OneThird + 1.228336910000000e-15 * tt3 * 0.25 + 4.784338640000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[4] = 3.992015430000000e+00 * tt5 + -2.401317520000000e-03 * tt0 + 4.617938410000000e-06 * tt1 * 0.50 + -3.881133330000000e-09 * tt2 * OneThird + 1.364114700000000e-12 * tt3 * 0.25 + -1.039254580000000e-01;
    } else {
      s0_R[4] = 3.092887670000000e+00 * tt5 +  5.484297160000000e-04 * tt0 + 1.265052280000000e-07 * tt1 * 0.50 + -8.794615559999999e-11 * tt2 * OneThird + 1.174123760000000e-14 * tt3 * 0.25 + 4.476696100000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[5] = 4.301798010000000e+00 * tt5 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 * 0.50 + -2.427638940000000e-08 * tt2 * OneThird + 9.292251240000000e-12 * tt3 * 0.25 + 3.716662450000000e+00;
    } else {
      s0_R[5] = 4.017210900000000e+00 * tt5 +  2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 * 0.50 + 1.142463700000000e-10 * tt2 * OneThird + -1.079085350000000e-14 * tt3 * 0.25 + 3.785102150000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[6] = 4.276112690000000e+00 * tt5 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 * 0.50 + -2.157708130000000e-08 * tt2 * OneThird + 8.624543630000000e-12 * tt3 * 0.25 + 3.435050740000000e+00;
    } else {
      s0_R[6] = 4.165002850000000e+00 * tt5 +  4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 * 0.50 + 3.711859860000000e-10 * tt2 * OneThird + -2.879083050000000e-14 * tt3 * 0.25 + 2.916156620000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[7] = 4.198640560000000e+00 * tt5 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 * 0.50 + -5.487970620000000e-09 * tt2 * OneThird + 1.771978170000000e-12 * tt3 * 0.25 + -8.490322080000000e-01;
    } else {
      s0_R[7] = 3.033992490000000e+00 * tt5 +  2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 * 0.50 + -9.704198700000000e-11 * tt2 * OneThird + 1.682009920000000e-14 * tt3 * 0.25 + 4.966770100000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[8] = 3.298616281368291e+00 * tt5 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 * 0.50 + 5.642920880408571e-09 * tt2 * OneThird + -2.445407041148433e-12 * tt3 * 0.25 + 3.950623591964732e+00;
    } else {
      s0_R[8] = 2.926639911210682e+00 * tt5 +  1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 * 0.50 + 1.009704225872734e-10 * tt2 * OneThird + -6.753354387142974e-15 * tt3 * 0.25 + 5.980528055036107e+00;
    };

  };

  template <class Type>
  void getGibbsFunctions_RT(double& p, Type& T,vector<Type>& g0_RT) {

    vector<Type> h0_RT(kk, 0.0);
    vector<Type> s0_R(kk,  0.0);

    getEnthalpies_RT(T, h0_RT);
    getEntropies_R(T, s0_R);
    for(int k = 0; k < kk; ++k) { g0_RT[k] = h0_RT[k] - s0_R[k]; }

  };

  template <class Type>
  void getEqConstants(double& p, Type& T,vector<Type>& keqs) {

    double       p0 = OneAtm;
    Type         RT = GasConstant * T;
    Type         C0 = p0 / RT;
    vector<Type> g0_RT(kk, 0.0);

    getGibbsFunctions_RT(p, T, g0_RT);
    for(int k = 0; k < kk; ++k) { g0_RT[k] = exp(g0_RT[k]); }

    keqs[0] = (g0_RT[3] * g0_RT[4]);
    keqs[1] = (g0_RT[1] * g0_RT[4]);
    keqs[2] = (g0_RT[1] * g0_RT[7]);
    keqs[3] = (g0_RT[4] * g0_RT[4]);
    keqs[4] = (C0 * g0_RT[0]);
    keqs[5] = (C0 * g0_RT[7]);
    keqs[6] = (C0 * g0_RT[2]);
    keqs[7] = (C0 * g0_RT[4]);
    keqs[8] = (C0 * g0_RT[5]);
    keqs[9] = (C0 * g0_RT[5]);
    keqs[10] = (g0_RT[4] * g0_RT[4]);
    keqs[11] = (g0_RT[0] * g0_RT[2]);
    keqs[12] = (g0_RT[7] * g0_RT[3]);
    keqs[13] = (g0_RT[2] * g0_RT[4]);
    keqs[14] = (g0_RT[7] * g0_RT[2]);
    keqs[15] = (g0_RT[7] * g0_RT[2]);
    keqs[16] = (C0 * g0_RT[6]);
    keqs[17] = (g0_RT[6] * g0_RT[2]);
    keqs[18] = (g0_RT[6] * g0_RT[2]);
    keqs[19] = (g0_RT[0] * g0_RT[5]);
    keqs[20] = (g0_RT[7] * g0_RT[4]);
    keqs[21] = (g0_RT[7] * g0_RT[5]);
    keqs[22] = (g0_RT[7] * g0_RT[5]);
    keqs[23] = (g0_RT[5] * g0_RT[4]);

    keqs[0] /= (g0_RT[1] * g0_RT[2]);
    keqs[1] /= (g0_RT[0] * g0_RT[3]);
    keqs[2] /= (g0_RT[0] * g0_RT[4]);
    keqs[3] /= (g0_RT[7] * g0_RT[3]);
    keqs[4] /= (g0_RT[1] * g0_RT[1]);
    keqs[5] /= (g0_RT[1] * g0_RT[4]);
    keqs[6] /= (g0_RT[3] * g0_RT[3]);
    keqs[7] /= (g0_RT[1] * g0_RT[3]);
    keqs[8] /= (g0_RT[3] * g0_RT[4]);
    keqs[9] /= (g0_RT[1] * g0_RT[2]);
    keqs[10] /= (g0_RT[1] * g0_RT[5]);
    keqs[11] /= (g0_RT[1] * g0_RT[5]);
    keqs[12] /= (g0_RT[1] * g0_RT[5]);
    keqs[13] /= (g0_RT[5] * g0_RT[3]);
    keqs[14] /= (g0_RT[5] * g0_RT[4]);
    keqs[15] /= (g0_RT[5] * g0_RT[4]);
    keqs[16] /= (g0_RT[4] * g0_RT[4]);
    keqs[17] /= (g0_RT[5] * g0_RT[5]);
    keqs[18] /= (g0_RT[5] * g0_RT[5]);
    keqs[19] /= (g0_RT[1] * g0_RT[6]);
    keqs[20] /= (g0_RT[1] * g0_RT[6]);
    keqs[21] /= (g0_RT[6] * g0_RT[4]);
    keqs[22] /= (g0_RT[6] * g0_RT[4]);
    keqs[23] /= (g0_RT[6] * g0_RT[3]);

  };

  template <class Type>
  void getTemperature(double& h, double& Told, vector<Type>& Y, Type& T) {

    double         tol   = 1.0e-08;
    int            niter = 500;
    Type           RT;
    Type           To;
    Type           Tp;
    Type           dT = 1.0;
    Type           FY = 0.0;
    Type           JY = 0.0;
    vector<Type>   hi(kk,    0.0);
    vector<Type>   dhidT(kk, 0.0);
    vector<double> mw = getMolecularWeights();

    To = Told;
    Tp = Told;

    for(int i = 0; i < niter; ++i) { 

      RT = GasConstant * To;
      mech::getEnthalpies_RT(To, hi);
      mech::getEnthalpiesDerivatives(To, dhidT);
      for(int k = 0; k < kk; ++k) { hi[k]    = RT * hi[k] / mw[k]; }
      for(int k = 0; k < kk; ++k) { dhidT[k] = GasConstant * dhidT[k] / mw[k]; }
      for(int k = 0; k < kk; ++k) { FY -= hi[k] * Y[k]; }
      for(int k = 0; k < kk; ++k) { JY -= dhidT[k] * Y[k]; }
      FY = h + FY;
      JY = 1.0 / JY;
      dT = -FY * JY;
      Tp = To + dT;
      To = Tp;

      if( (fabs(dT) < tol)) {
	T = Tp;
	return;
      }

      FY = 0.0;
      JY = 0.0;

    }

    T = Tp;

  };

  /***********************************************************************************/
  /* kinetics                                                                        */
  /***********************************************************************************/
  template <class Type>
  void getFalloffRates(Type& T, vector<Type>& C, vector<Type>& kfwd) {

    int          nfall = 2; 
    Type         tlog  = log(T);
    Type         rt    = 1.0 / T;
    Type         lpr;
    Type         cc;
    Type         nn;
    Type         f1;
    Type         lgf;
    vector<Type> pr(nfall);
    vector<Type> khi(nfall);
    vector<Type> klo(nfall);
    vector<Type> work(nfall);

/*
    khi[0] =  exp(2.226013305655e+01 + 4.4000e-01 * tlog);// kinf(9)
    khi[1] =  exp(2.528239208443e+01 + -2.7000e-01 * tlog);// kinf(16)

    klo[0] =  exp(3.168280606373e+01 + -1.4000e+00 * tlog);// k0(9)
    klo[1] =  exp(4.476434744662e+01 + -3.2000e+00 * tlog);// k0(16)
*/
    khi[0] =  exp( mech::kfwd_coeffs[24].p1 + mech::kfwd_coeffs[24].p2 * tlog + mech::kfwd_coeffs[24].p3 * rt);// kinf for 'H + O2 (+ M) <=> HO2 (+ M)'
    khi[1] =  exp( mech::kfwd_coeffs[25].p1 + mech::kfwd_coeffs[25].p2 * tlog + mech::kfwd_coeffs[25].p3 * rt);// kinf for 'OH + OH (+ M) <=> H2O2 (+ M)'

    klo[0] =  exp( mech::kfwd_coeffs[9].p1 + mech::kfwd_coeffs[9].p2 * tlog + mech::kfwd_coeffs[9].p3 * rt);// k0 for 'H + O2 (+ M) <=> HO2 (+ M)'
    klo[1] =  exp( mech::kfwd_coeffs[16].p1 + mech::kfwd_coeffs[16].p2 * tlog + mech::kfwd_coeffs[16].p3 * rt);// k0 for 'OH + OH (+ M) <=> H2O2 (+ M)'

    pr[0] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.6e+01 * C[7] + 1.0e+00 * C[8]) * klo[0] / khi[0];
    pr[1] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.0e+00 * C[7] + 1.0e+00 * C[8]) * klo[1] / khi[1];

    work[0] = 0.50;
    work[1] = 0.43;
    
    for(int r = 0; r < nfall; ++r) {
      lpr     =  log10(pr[r]);
      cc      = -0.40 - 0.67 * log10(work[r]);
      nn      =  0.75 - 1.27 * log10(work[r]);
      f1      =  (lpr + cc)/(nn - 0.14 * (lpr + cc));
      work[r] =  log10(work[r])/(1 + f1 * f1);
      work[r] =  pow(10.0, work[r]);
      work[r] =  (pr[r] * work[r])/(1 + pr[r]);
    };

    kfwd[9] = khi[0] * work[0];
    kfwd[16] = khi[1] * work[1];

  };


  template <class Type>
  void updateRateConstants(double& p, Type& T, vector<Type>& C,
			   vector<Type>& kfwd, vector<Type>& krev) {

    Type         tlog = log(T);
    Type         rt   = 1.0 / T;
    vector<Type> keqs(ii, 0.0);

    getEqConstants(p, T, keqs);

/*    for (int i=0; i < 24; ++i)
      kfwd[i] = exp( mech::kfwd_coeffs[i].p1 + mech::kfwd_coeffs[i].p2 * tlog + mech::kfwd_coeffs[i].p3 * rt);
*/
    for (int i=0; i < 24; ++i){
    	if(i!=9 && i!=16){ kfwd[i] = exp( mech::kfwd_coeffs[i].p1 + mech::kfwd_coeffs[i].p2 * tlog + mech::kfwd_coeffs[i].p3 * rt);}
	}
/*	for (int i=0; i < 24; ++i) {
		std::cout << "i=" << i << "; " << mech::kfwd_coeffs[i].p1 << "; " << mech::kfwd_coeffs[i].p2 << "; " << mech::kfwd_coeffs[i].p3 << endl;
	}
*/
    /*
      Not used anymore
    kfwd[0] =  exp( + 4.69322933e+00  + 3.805162e+00 * tlog  + 8.327832e+02 * rt);
    kfwd[1] =  exp( + 1.06360705e+01  + 2.131220e+00 * tlog  + 3.023402e+03 * rt);
    kfwd[2] =  exp( + 8.64871289e+00  + 1.018815e+00 * tlog  + 2.590729e+03 * rt);
    kfwd[3] =  exp( + 2.49071942e+01 -3.611042e+00 * tlog -1.540226e+02 * rt);
    kfwd[4] =  exp( + 3.45285820e+00 -1.358729e+00 * tlog  + 5.374278e+01 * rt);
    kfwd[5] =  exp( + 7.24959983e+00  + 2.931154e+00 * tlog -3.649860e+01 * rt);
    kfwd[6] =  exp( + 2.09052088e+01  + 3.963303e+00 * tlog  + 3.091555e+03 * rt);
    kfwd[7] =  exp( + 1.03731787e+01 -1.960075e+00 * tlog  + 7.479079e+02 * rt);
    kfwd[8] =  exp( + 5.90693271e+00 -1.299181e+00 * tlog  + 1.447039e+03 * rt);
    kfwd[10] =  exp( + 2.72170357e+01  + 2.783273e+00 * tlog  + 1.475626e+03 * rt);
    kfwd[11] =  exp( + 1.56058894e+01  + 3.052324e+00 * tlog  + 8.108341e+02 * rt);
    kfwd[12] =  exp( + 6.65441345e+00 -2.389807e+00 * tlog  + 2.790545e+03 * rt);
    kfwd[13] =  exp( + 2.32052605e+01  + 1.258821e-01 * tlog  + 1.659343e+03 * rt);
    kfwd[14] =  exp( + 1.36142808e+01 -2.874468e+00 * tlog  + 3.410635e+03 * rt);
    kfwd[15] =  exp( + 2.55373169e+01 -3.026698e+00 * tlog  + 1.818921e+03 * rt);
    kfwd[17] =  exp( + 4.53929795e+00 -2.655181e+00 * tlog -4.541554e+02 * rt);
    kfwd[18] =  exp( + 1.20825636e+00 -1.026079e-01 * tlog  + 1.317588e+03 * rt);
    kfwd[19] =  exp( + 2.45372403e+01 -1.711131e-01 * tlog  + 1.639791e+03 * rt);
    kfwd[20] =  exp( + 2.14838429e+01  + 1.571397e+00 * tlog  + 2.688730e+03 * rt);
    kfwd[21] =  exp( + 1.74577328e+01 -2.979575e+00 * tlog  + 4.576750e+02 * rt);
    kfwd[22] =  exp( + 1.05185231e+01  + 1.828285e+00 * tlog  + 1.330540e+03 * rt);
    kfwd[23] =  exp( + 1.76863239e+01  + 3.859582e+00 * tlog  + 2.356892e+03 * rt);
    */

    getFalloffRates(T, C, kfwd);

    kfwd[4] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[4];
    kfwd[5] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[5];
    kfwd[6] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[6];
    kfwd[7] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[7];
    kfwd[8] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[8];

    for(int i = 0; i < ii; ++i) { krev[i] = kfwd[i] * keqs[i]; }

  }

  template <class Type>
  void updateSpeciesSourceTerm(double& p, double& h, double& Told,
			       vector<Type>& Y,vector<Type>& FY) {


    Type           rho;
    Type           W;
    Type           T;
    vector<Type>   C(kk,    0.0);
    vector<Type>   kfwd(ii, 0.0);
    vector<Type>   krev(ii, 0.0);
    vector<Type>   Rfwd(ii, 0.0);
    vector<Type>   Rrev(ii, 0.0);
    vector<Type>   Rnet(ii, 0.0);
    vector<Type>   omega(kk,0.0);
    vector<double> mw = getMolecularWeights();

    mech::getTemperature(h, Told, Y, T);

    W   = 0.0;
    for(int k = 0; k < kk; ++k) { W += Y[k] / mw[k]; }
    W   = 1.0/W;
    rho = (p * W)/(GasConstant * T);
    for(int k = 0; k < kk; ++k) { C[k] = rho * Y[k] / mw[k]; }
    mech::updateRateConstants(p, T, C, kfwd, krev);

    Rfwd[0] = kfwd[0] * C[1] * C[2];
    Rfwd[1] = kfwd[1] * C[0] * C[3];
    Rfwd[2] = kfwd[2] * C[0] * C[4];
    Rfwd[3] = kfwd[3] * C[7] * C[3];
    Rfwd[4] = kfwd[4] * C[1] * C[1];
    Rfwd[5] = kfwd[5] * C[1] * C[4];
    Rfwd[6] = kfwd[6] * C[3] * C[3];
    Rfwd[7] = kfwd[7] * C[1] * C[3];
    Rfwd[8] = kfwd[8] * C[3] * C[4];
    Rfwd[9] = kfwd[9] * C[1] * C[2];
    Rfwd[10] = kfwd[10] * C[1] * C[5];
    Rfwd[11] = kfwd[11] * C[1] * C[5];
    Rfwd[12] = kfwd[12] * C[1] * C[5];
    Rfwd[13] = kfwd[13] * C[5] * C[3];
    Rfwd[14] = kfwd[14] * C[5] * C[4];
    Rfwd[15] = kfwd[15] * C[5] * C[4];
    Rfwd[16] = kfwd[16] * C[4] * C[4];
    Rfwd[17] = kfwd[17] * C[5] * C[5];
    Rfwd[18] = kfwd[18] * C[5] * C[5];
    Rfwd[19] = kfwd[19] * C[1] * C[6];
    Rfwd[20] = kfwd[20] * C[1] * C[6];
    Rfwd[21] = kfwd[21] * C[6] * C[4];
    Rfwd[22] = kfwd[22] * C[6] * C[4];
    Rfwd[23] = kfwd[23] * C[6] * C[3];

    Rrev[0] = krev[0] * C[3] * C[4];
    Rrev[1] = krev[1] * C[1] * C[4];
    Rrev[2] = krev[2] * C[1] * C[7];
    Rrev[3] = krev[3] * C[4] * C[4];
    Rrev[4] = krev[4] * C[0];
    Rrev[5] = krev[5] * C[7];
    Rrev[6] = krev[6] * C[2];
    Rrev[7] = krev[7] * C[4];
    Rrev[8] = krev[8] * C[5];
    Rrev[9] = krev[9] * C[5];
    Rrev[10] = krev[10] * C[4] * C[4];
    Rrev[11] = krev[11] * C[0] * C[2];
    Rrev[12] = krev[12] * C[7] * C[3];
    Rrev[13] = krev[13] * C[2] * C[4];
    Rrev[14] = krev[14] * C[7] * C[2];
    Rrev[15] = krev[15] * C[7] * C[2];
    Rrev[16] = krev[16] * C[6];
    Rrev[17] = krev[17] * C[6] * C[2];
    Rrev[18] = krev[18] * C[6] * C[2];
    Rrev[19] = krev[19] * C[0] * C[5];
    Rrev[20] = krev[20] * C[7] * C[4];
    Rrev[21] = krev[21] * C[7] * C[5];
    Rrev[22] = krev[22] * C[7] * C[5];
    Rrev[23] = krev[23] * C[5] * C[4];

    for(int i = 0; i < ii; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }

    omega[0] =  + Rnet[4] + Rnet[11] + Rnet[19] - Rnet[1] - Rnet[2] ;
    omega[1] =  + Rnet[1] + Rnet[2] - Rnet[0] - Rnet[4] - Rnet[4] - Rnet[5] - Rnet[7]
      - Rnet[9] - Rnet[10] - Rnet[11] - Rnet[12] - Rnet[19] - Rnet[20] ;
    omega[2] =  + Rnet[6] + Rnet[11] + Rnet[13] + Rnet[14] + Rnet[15] + Rnet[17] + Rnet[18]
      - Rnet[0] - Rnet[9] ;
    omega[3] =  + Rnet[0] + Rnet[12] - Rnet[1] - Rnet[3] - Rnet[6] - Rnet[6] - Rnet[7]
      - Rnet[8] - Rnet[13] - Rnet[23] ;
    omega[4] =  + Rnet[7] + Rnet[0] + Rnet[1] + Rnet[3] + Rnet[3] + Rnet[10] + Rnet[10]
      + Rnet[13] + Rnet[20] + Rnet[23] - Rnet[2] - Rnet[5] - Rnet[8] - Rnet[14]
      - Rnet[15] - Rnet[16] - Rnet[16] - Rnet[21] - Rnet[22] ;
    omega[5] =  + Rnet[8] + Rnet[9] + Rnet[19] + Rnet[21] + Rnet[22] + Rnet[23] - Rnet[10]
      - Rnet[11] - Rnet[12] - Rnet[13] - Rnet[14] - Rnet[15] - Rnet[17] - Rnet[17]
      - Rnet[18] - Rnet[18] ;
    omega[6] =  + Rnet[16] + Rnet[17] + Rnet[18] - Rnet[19] - Rnet[20] - Rnet[21] - Rnet[22]
      - Rnet[23] ;
    omega[7] =  + Rnet[5] + Rnet[2] + Rnet[12] + Rnet[14] + Rnet[15] + Rnet[20] + Rnet[21]
      + Rnet[22] - Rnet[3] ;

    for(int k = 0; k < kk; ++k) { FY[k] = mw[k] * omega[k] / rho; }

  }

  template <class Type>
  vector<Type> transpose(vector<Type>& mat)
  {

    int nn = mat.size();
    int n  = sqrt(nn);
    
    for (int i = 0; i < n; i++)
      {
	for (int j = i+1; j < n; j++)
	  {
	    Type temp  = mat[i*n+j];
	    mat[i*n+j] = mat[j*n+i];
	    mat[j*n+i] = temp;
	  }
      }
    
    return(mat);

  }
  
} // namespace mech

void nspecies_(int& nspec) {

  nspec = mech::kk;
  
}

void nreactions_(int& nreac) {

  nreac = mech::ii;
  
}

void nelements_(int& nelem) {

  nelem = mech::mm;
  
}

void getmolecularweights_(double* mw) {

  vector<double> mwts = mech::getMolecularWeights();
  for(int k = 0; k < mech::kk; ++k) { mw[k] = mwts[k]; }
  
}

void getmassfractions_(double* x, double* y) {

  double         mmw = 0.0;
  vector<double> mw  = mech::getMolecularWeights();

  for(int k = 0; k < mech::kk; k++) {
    double xk = max(x[k], 0.0);
    mmw  += mw[k] * xk;
  }
  
  mmw = 1.0/mmw;
  for(int k = 0; k < mech::kk; ++k) { y[k] = max(mw[k] * x[k] * mmw, 0.0); }
  
}

void getenthalpymass_(double* T, double* y, double* h) {

  double         RT = mech::GasConstant * (*T);
  vector<double> mw = mech::getMolecularWeights();
  vector<double> hi(mech::kk, 0.0);

  mech::getEnthalpies_RT(*T, hi);
  *h = 0.0;
  for(int k = 0; k < mech::kk; ++k) {
    hi[k]  = RT * hi[k];
    *h    += hi[k] * y[k] / mw[k];
  }
  
}

void getdensity_(double* T, double* p, double* y, double* rho) {

  double         RT  = mech::GasConstant * (*T);
  double         mmw = 0.0;
  vector<double> mw  = mech::getMolecularWeights();
  for(int k = 0; k < mech::kk; ++k) {
    mmw += y[k]/mw[k];
  }
  mmw = 1.0/mmw;

  *rho = (*p) * mmw / RT;
  
}

void gettemperature_(double* h, double* y, double* Told, double* Tnew) {

  double         T;
  vector<double> Y(mech::kk, 0.0);

  for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
  mech::getTemperature(*h, *Told, Y, T);
  *Tnew = T;
  
}

void gettemperaturesource_(double* p, double* h, double* T, double* y, double* fT) {

  double         FT = 0.0;
  double         RT = mech::GasConstant * (*T);
  double         cp;
  double         rho;
  double         mmw = 0.0;
  vector<double> Y(mech::kk);
  vector<double> FY(mech::kk);
  vector<double> hi(mech::kk);
  vector<double> cpi(mech::kk);
  vector<double> mw = mech::getMolecularWeights();

  for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
  for(int k = 0; k < mech::kk; ++k) { mmw += y[k]/mw[k]; }
  mmw = 1.0/mmw;
  rho = (*p) * mmw / RT;

  mech::getSpecificHeats_R(*T, cpi);
  mech::getEnthalpies_RT((*T), hi);
  mech::updateSpeciesSourceTerm(*p, *h, *T, Y, FY);

  cp = 0.0;
  for(int k = 0; k < mech::kk; ++k) {
    hi[k] = RT * hi[k] / mw[k];
    cp   += mech::GasConstant * cpi[k] * Y[k] / mw[k];
  }

  for(int k = 0; k < mech::kk; ++k) { FT -= mw[k] * hi[k] * FY[k] / (cp * rho); }
  
  *fT = FT;
  
}

void getsource_(double* p, double* h, double* Told, double* y, double* fy) {

  /* declarations */
  cout.precision(4);
  cout.setf(ios::scientific);

  vector<double> Y(mech::kk,  0.0);
  vector<double> FY(mech::kk, 0.0);

  for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
  mech::updateSpeciesSourceTerm(*p, *h, *Told, Y, FY);
  copy(FY.begin(), FY.end(), fy);
}

void getjacobian_(double* p, double* h, double* Told, double* y, double* dfy) {

  /* declarations */
  cout.precision(4);
  cout.setf(ios::scientific);

  int                  kk = mech::kk;
  vector<double>       Y(kk, 0.0);
  vector<double>       FY(kk,0.0);
  vector<double>       DFY(kk*kk,0.0);
  vector< AD<double> > a_Y(kk, 0.0);
  vector< AD<double> > a_FY(kk,0.0);

  /* AD stuff */
  for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
  copy(Y.begin(), Y.end(), a_Y.begin());

  /* AD Stuff */
  CppAD::Independent(a_Y);
  mech::updateSpeciesSourceTerm(*p, *h, *Told, a_Y, a_FY);
  CppAD::ADFun<double> tapef(a_Y, a_FY);
  //FY   = tapef.Forward(0, Y);
  
  DFY  = tapef.Jacobian(Y);
  DFY  = mech::transpose(DFY);

  copy(DFY.begin(), DFY.end(), dfy);
   
}

void readreactionrates_(char* filename) {
  ifstream f_in;
  
  mech::kfwd_coeffs.clear();
  correctname(filename); 
  f_in.open(filename, ios::in);
  if(! f_in.good()){std::cout << "the name of the parameters file is invalid" << endl;}
  for(int i=0; i < 26; ++i) {
    mech::Kfentry entry;
    f_in >> entry.p1 >> entry.p2 >> entry.p3;
	  std::cout << "log A: " << entry.p1 <<" n: " << entry.p2 <<" E/r: " << entry.p3 << endl;
    mech::kfwd_coeffs.push_back(entry);
  }
  
  f_in.close();
}
void correctname(char* fname)
{
  char* i = fname;
  char* j = fname;
  while(*j != 0)
  {
    *i = *j++;
    if(*i != ' ')
      i++;
  }
  *i = 0;
}
