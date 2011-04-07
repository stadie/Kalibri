// pt_resolution file from Mikko

#ifndef __ptresolution_h__
#define __ptresolution_h__


#include "TMath.h"

// Switch MC truth or data resolutions
bool _ismcjer = true;
bool _jerkscale = false;

// Hauke's resolutions
const int _nres = 7;//6;
const int _npar = 4;
/* float _pres[_nres][_npar] = */
/* // From CMSSW Fall10_PtResolution_AK5PF (Feb11-Feb25) */
/*   {{/\*0 0.5*\/ 3.96859, 0.183476, 0, 0.626273},// 0.026724 */
/*    {/\*0.5 1*\/ 3.55226, 0.240261, 0, 0.525708},// 0.0254424 */
/*    {/\*1 1.5*\/ 4.54826, 0.226521, 0, 0.589632},// 0.0335984 */
/*    {/\*1.5 2*\/ 4.62622, 0.236643, 0, 0.487379},// 0.0448822 */
/*    {/\*2 2.5*\/ 2.53324, 0.343062, 0, 0.286621},// 0.0517279 */
/*    {/\*2.5 3*\/-3.33814, 0.733604, 0, 0.0826378},// 0.070123 */
/*    {/\*3 9.9*\/ 2.95397, 0.116194, 0, 0.960861}};// 0.0564672 */
/* float _pres_ak5jpt[_nres][_npar] = */
/* // From CMSSW Fall10_PtResolution_AK5JPT (Feb11-Feb25) */
/*   {{/\*0 0.5*\/ 4.76130, 0.224764, 0, 0.546172},// 0.039651 */
/*    {/\*0.5 1*\/ 4.49601, 0.262162, 0, 0.496672},// 0.0188774 */
/*    {/\*1 1.5*\/ 4.91436, 0.298285, 0, 0.497331},// 0.034662 */
/*    {/\*1.5 2*\/ 5.35503, 0.310884, 0, 0.415469},// 0.0416232 */
/*    {/\*2 2.5*\/ 3.96073, 0.411402, 0, 0.362306},// 0.0502304 */
/*    {/\*2.5 3*\/ 0.718871, 0.47073, 0, 0.213428},// 0.0616422 */
/*    {/\*3 9.9*\/ 2.95397, 0.116194, 0, 0.960861}};// 0.0564672 (from PF) */
/* float _pres_ak5calo[_nres][_npar] = */
/* // From CMSSW Fall10_PtResolution_AK5CALO (Feb11-Feb25) */
/*   {{/\*0 0.5*\/ 5.89372, 0.398184, 0, 0.411555},// 0.0325261 */
/*    {/\*0.5 1*\/ 4.70682, 0.586653, 0, 0.277176},// 0.0512329 */
/*    {/\*1 1.5*\/ 5.06249, 0.575817, 0, 0.317870},// 0.0392882 */
/*    {/\*1.5 2*\/ 3.86101, 0.805593, 0, 0.163246},// 0.0673781 */
/*    {/\*2 2.5*\/ 3.05359, 0.487558, 0, 0.243042},// 0.0453157 */
/*    {/\*2.5 3*\/-1.91901, 0.597099, 0, 0.127347},// 0.0708867 */
/*    {/\*3 9.9*\/ 2.81913, 0.0881245, 0, 1.11217}};// 0.0864411 */

float _pres[_nres][_npar] =
 // From CMSSW Spring10_PtResolution_AK5PF
 {{/*0 0.5*/ -0.349206, 0.297831, 0, 0.471121},
  {/*0.5 1*/ -0.499735, 0.336391, 0, 0.430689},
  {/*1 1.5*/ -0.561649, 0.420293, 0, 0.392398},
  {/*1.5 2*/ -1.12329,  0.657891, 0, 0.139595},
  {/*2 2.5*/  1.04792,  0.466763, 0, 0.193137},
  {/*2.5 3*/  2.56933,  0.305802, 0, 0.398929}, // changed since May5?
  {/*2.5 3*/  2.56933,  0.305802, 0, 0.398929}}; // copy
  // Recover old fit for test purposes
  //{ 1.89978, 0.33427, 0.00000, 0.36547}};// 2.5-3.0 (28.3/9), May 5 evening
float _pres_ak5jpt[_nres][_npar] =
 // From CMSSW Spring10_PtResolution_AK5JPT
 {{/*0 0.5*/ 2.51562, 0.316061, 0, 0.445628},
  {/*0.5 1*/ 2.41271, 0.332129, 0, 0.432430},
  {/*1 1.5*/ 2.42631, 0.371536, 0, 0.439477},
  {/*1.5 2*/ 2.67972, 0.511479, 0, 0.257115},
  {/*2 2.5*/ 3.39005, 0.499497, 0, 0.273727},
  {/*2.5 3*/ 2.28498, 0.390866, 0, 0.297617},
  {/*2.5 3*/ 2.28498, 0.390866, 0, 0.297617}}; // copy
float _pres_ak5calo[_nres][_npar] =
 // From CMSSW Spring10_PtResolution_AK5Calo
 {{/*0 0.5*/ 5.0913, 0.511619, 0, 0.325477},
  {/*0.5 1*/ 4.93688, 0.543349, 0, 0.310638},
  {/*1 1.5*/ 4.95729, 0.574202, 0, 0.32344},
  {/*1.5 2*/ 3.37322, 0.96262, 0, 0.0860736},
  {/*2 2.5*/ 4.06368, 0.351561, 0, 0.356167},
  {/*2.5 3*/ 2.59292, 0.357679, 0, 0.321092},
  {/*2.5 3*/ 2.59292, 0.357679, 0, 0.321092}}; // copy

// Hauke's numbers, Feb 11; confirmed Feb25 (in Fall10_PtResolution_AK5*.txt)
float _pconst_ak5pf[_nres] =   {0.027, 0.025, 0.034, 0.045, 0.052, 0.070,0.056};
float _pconst_ak5calo[_nres] = {0.033, 0.051, 0.039, 0.067, 0.045, 0.071,0.086};
float _pconst_ak5jpt[_nres] =  {0.040, 0.019, 0.035, 0.042, 0.050, 0.062,0.062};
// Replace 1.5<|y|<3. numbers with my own for PF to improve ansatz fit (Feb23)
//float _pconst_ak5pf[_nres] =   {0.034, 0.030, 0.032, 0.029, 0.026, 0.043};
//float _pconst_ak5calo[_nres] = {0.033, 0.051, 0.039, 0.067, 0.045, 0.043};
//float _pconst_ak5jpt[_nres] =  {0.040, 0.019, 0.035, 0.042, 0.026, 0.043};
// Replace 1.5<|y|<3. numbers with my own for PF to improve ansatz fit (Feb24)
//float _pconst_ak5pf[_nres] =   {0.037, 0.032, 0.034, 0.028, 0.027, 0.059};
//float _pconst_ak5calo[_nres] = {0.033, 0.051, 0.039, 0.067, 0.045, 0.059};
//float _pconst_ak5jpt[_nres] =  {0.040, 0.019, 0.035, 0.042, 0.027, 0.059};

// k-scale for PFJets from JME-10-014
float _pscale_ak5pf[_nres] = {1.07, 1.07, 1.10, 1.10, 1.07, 1.18};

double ptresolution(float pt, float eta) {
  
  int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  float N = _pres[ieta][0];
  float S = _pres[ieta][1];
  float C = _pres[ieta][2];
  float m = _pres[ieta][3];
  double c = _pconst_ak5pf[ieta];
  double k = _pscale_ak5pf[ieta];
  if (_ismcjer || !_jerkscale) k = 1.;
  if (_ismcjer ||  _jerkscale) c = 0.;

  double res = k*sqrt(TMath::Sign(N*N,N)/(pt*pt) + S*S*pow(pt,m-1) + C*C + c*c);
  //if (ieta==0 && res>0.12) res = 0.12;

  return res;
}

double ptresolution_jpt(float pt, float eta) {
  
  int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  float N = _pres_ak5jpt[ieta][0];
  float S = _pres_ak5jpt[ieta][1];
  float C = _pres_ak5jpt[ieta][2];
  float m = _pres_ak5jpt[ieta][3];
  double c = _pconst_ak5jpt[ieta];
  double k = 1.;
  if (_ismcjer) k = 1.;
  if (_ismcjer) c = 0.;

  double res = k*sqrt(TMath::Sign(N*N,N)/(pt*pt) + S*S*pow(pt,m-1) + C*C + c*c);

  return res;
}

double ptresolution_calo(float pt, float eta) {
  
  int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  float N = _pres_ak5calo[ieta][0];
  float S = _pres_ak5calo[ieta][1];
  float C = _pres_ak5calo[ieta][2];
  float m = _pres_ak5calo[ieta][3];
  double c = _pconst_ak5calo[ieta];
  double k = 1.;
  if (_ismcjer) k = 1.;
  if (_ismcjer) c = 0.;

  double res = k*sqrt(TMath::Sign(N*N,N)/(pt*pt) + S*S*pow(pt,m-1) + C*C + c*c);

  return res;
}

#endif // __ptresolution_h__
