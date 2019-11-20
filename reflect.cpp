#include <cmath>
#include <iostream>
#include <stdio.h>
#include <random>
#include <string>
#define SQR(x)  ((x)*(x))  // square of a number

const double degree = M_PI/180.0;

/* TUNABLE PARAMETERS */
const double height = 38; // Height of the observer [km]
const double sineps_min = sin(-40*degree); // Minimum value of sin(epsilon), where epsilon is the observed elevation of the event
const double sineps_max = 0; // Maximum value of sin(epsilon)
constexpr int N_sineps = 41; // Number of bins in sin(epsilon)
constexpr int N_psi_f = 21; // Number of bins in psi_f, the final polarization angle
constexpr int N_P_f = 11; // Number of bins in P_f, the final degree of polarization
const unsigned long int N_pts_MC = 1e7; // Amount of polarization directions to randomly generate
std::uniform_real_distribution<double> dist_psi_i(-M_PI/2., M_PI/2.); // Distribution of incident polarization angles

class Stokes{
public:
  double st_v[3]; // First three components of a normalized Stokes vector describing a *linearly* polarized wave

  Stokes(double P, double psi){
    /*
     * Constructor of a linearly polarized Stokes vector. Arguments:
     *  P   --- degree of polarization
     *  psi --- polarization angle with respect to the horizontal (IN RADIANS)
     */
    setPars(P, psi);
  }

  void setPars(double P, double psi){
    /*
     * Changes the parameters of a linearly polarized Stokes vector. Arguments:
     *  P   --- degree of polarization
     *  psi --- polarization angle with respect to the horizontal (IN RADIANS)
     */
    st_v[0] = 1;
    st_v[1] = P*cos(2*psi);
    st_v[2] = P*sin(2*psi);
  }

  double getReflectedPsi(double n1, double n2, double theta_i){
    /*Returns the reflected polarisation angle. Parameters:
        n1      -- refraction index of the incident medium
        n2      -- refraction index of the reflecting medium
        theta_i -- initial incidence angle **respect to the normal**(IN RADIANS)
    */
    double theta_t = asin(n1/n2 * sin(theta_i));

    // Fresnel coefficients: see eqs.(8-22a)[r_H] and (8-24a)[r_V] in Goldstein
    double r_H = (n1*cos(theta_i) - n2*cos(theta_t))/(n1*cos(theta_i) + n2*cos(theta_t));
    double r_V = (n2*cos(theta_i) - n1*cos(theta_t))/(n2*cos(theta_i) + n1*cos(theta_t));

    // Mueller matrix: see eq.(8.34) in Goldstein
    double M[3][3] = { {(SQR(r_H)+SQR(r_V))/2., (SQR(r_H)-SQR(r_V))/2., 0.},
		       {(SQR(r_H)-SQR(r_V))/2., (SQR(r_H)+SQR(r_V))/2., 0.},
		       {0, 0, r_H*r_V} };
    double ref[3];
    for(int i=0; i<3; ++i){
      ref[i] = 0;
      for(int j=0; j<3; ++j)
	ref[i] += M[i][j]*st_v[j];
    }

    return atan2(ref[2], ref[1])/2;
  }

  double getReflectedP(double n1, double n2, double theta_i){
    /*Returns the reflected polarisation angle. Parameters:
        n1      -- refraction index of the incident medium
        n2      -- refraction index of the reflecting medium
        theta_i -- initial incidence angle **respect to the normal**(IN RADIANS)
    */
    double theta_t = asin(n1/n2 * sin(theta_i));

    // Fresnel coefficients: see eqs.(8-22a)[r_H] and (8-24a)[r_V] in Goldstein
    double r_H = (n1*cos(theta_i) - n2*cos(theta_t))/(n1*cos(theta_i) + n2*cos(theta_t));
    double r_V = (n2*cos(theta_i) - n1*cos(theta_t))/(n2*cos(theta_i) + n1*cos(theta_t));

    // Mueller matrix: see eq.(8.34) in Goldstein
    double M[3][3] = { {(SQR(r_H)+SQR(r_V))/2., (SQR(r_H)-SQR(r_V))/2., 0.},
		       {(SQR(r_H)-SQR(r_V))/2., (SQR(r_H)+SQR(r_V))/2., 0.},
		       {0, 0, r_H*r_V} };
    double ref[3];
    for(int i=0; i<3; ++i){
      ref[i] = 0;
      for(int j=0; j<3; ++j)
	ref[i] += M[i][j]*st_v[j];
    }

    return sqrt(SQR(ref[2]) + SQR(ref[1]))/ref[0];
  }  
};

int main(int argc, char* argv[]){
  if(argc!=2){
    std::cerr<<"Input as an argument the incident degree of polarization"<<std::endl;
    return -1;
  }

  const char* fname = ("Output_" + std::string(argv[1]) + ".dat").c_str();
  FILE * ofile = fopen(fname, "w");
  std::mt19937_64 generator;

  Stokes st(1, 0);
  
  double P = atof(argv[1]);

  unsigned long int N_evts[N_sineps][N_psi_f][N_P_f];
  /* N_evts[i][j][k] is the expected number of events for
   *     sin(epsilon) = sineps_min + (sineps_max - sineps_min) * i/(N_sineps-1)
   *     psi_f = -pi/2 + pi * j/(N_psi_f-1)
   *     P_f = k/(N_P_f-1)
   *
   * epsilon is the observed elevation angle
   * psi_f   ''  ''    ''    polarization angle
   * P_f     ''  ''    ''    degree of polarization
   */
  for(int i=0; i<N_sineps; ++i)
    for(int j=0; j<N_psi_f; ++j)
      for(int k=0; k<N_P_f; ++k)
	N_evts[i][j][k] = 0.;
  
  unsigned long int ctr = 1; // Conunter to compute progress and average

  for(int i_sineps = 0; i_sineps < N_sineps; ++i_sineps){
    double sineps = sineps_min + i_sineps*(sineps_max-sineps_min)/(N_sineps-1.0);
    double epsilon = asin(sineps);
    /* Correct for the curvature of the Earth: Eq.(47) in (1811.00900) 
       The idea is that the events come isotropically in sin(epsilon),
       but the reflection angle theta != epsilon
     */

    /* Check whether we are above or below the horizon */
    bool isReflected;
    double theta = 0.;
    if( (epsilon > 0) || ((height + 6356.75)/6356.75 * cos(epsilon) > 1) ) //We are above the horizon
      isReflected = false;
    else{
      isReflected = true;
      theta = acos((height + 6356.75)/6356.75 * cos(epsilon));
    }

    for(unsigned long int i=0; i<N_pts_MC; ++i){
      /* Generate psi_i randomly */
      double psi_i = dist_psi_i(generator);

      double psi_f, P_f;
      if(!isReflected){
	psi_f = psi_i;
	P_f = P;
      } else{
	/* Reflect psi_i */
	st.setPars(P, psi_i);
	psi_f = st.getReflectedPsi(1, 1.35, 90*degree-theta);
	P_f = st.getReflectedP(1, 1.35, 90*degree-theta);
      }

      /* Compute the indices for N_evts */
      int i_Nevts = round((sin(epsilon) - sineps_min)/ (sineps_max - sineps_min) * (N_sineps-1));
      int j_Nevts = round((psi_f + M_PI/2)/M_PI * (N_psi_f-1));
      int k_Nevts = round(P_f * (N_P_f-1));

      ++N_evts[i_Nevts][j_Nevts][k_Nevts];

      ++ctr;
    }
    std::cout<<"Progress: "<<(1.0*ctr)/N_pts_MC<<"/"<<N_sineps<<std::endl;
  }

  /* Impose periodicity in psi_f
   * The reason for this fix is that ~50% of vertically polarized events will be
   * assigned psi_f = -pi/2 and ~50% of them will be assigned psi_f = pi/2
   * But we know this variable is periodic, so we have to sum the events in both bins
  */
  for(int i=0; i<N_sineps; ++i)
    for(int k=0; k<N_P_f; ++k){
      N_evts[i][0][k] = N_evts[i][0][k] + N_evts[i][N_psi_f-1][k];
      N_evts[i][N_psi_f-1][k] = N_evts[i][0][k];
    }

  fprintf(ofile, "# sin(eps)  psi_f    P_f   phi(sin(eps), psi_f, P_f)\n");
  for(int i=0; i<N_sineps; ++i)
    for(int j=0; j<N_psi_f; ++j)
      for(int k=0; k<N_P_f; ++k){
	double sinepsilon = sineps_min + (sineps_max - sineps_min) * i/(N_sineps-1.0);
	double psi_f = -M_PI/2 + M_PI * j/(N_psi_f-1.0);
	double P_f = k/(N_P_f-1.0);

	fprintf(ofile, " %+.6f %+.6f %.3f %.3e\n", sinepsilon, psi_f, P_f, N_evts[i][j][k]/(1.0*ctr));
      }
}
