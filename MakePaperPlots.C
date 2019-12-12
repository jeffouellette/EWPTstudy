#ifndef __MakePaperPlots_C__
#define __MakePaperPlots_C__

const Color_t colors[10] = {kBlack, kRed+1, kAzure-1, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1};

double jb_int (const double x, const double y) {
  return x*x*log(1-exp(-sqrt(x*x+y*y)));
}

double jb (const double x) {
  double sum = 0;
  const double dx = 0.005;
  for (int i = 0; i < 2000; i++) {
    sum += jb_int (dx*i, x)*dx;
  }
  return sum;
  //return -pow(pi,4)/45. + pow(pi*x,2)/12. - pi*pow(x,1.5)/6. - pow(x,4)*(log(pow(x,2))-cb)/32.; // x << 1 approximation
}

double jf_int (const double x, const double y) {
  return x*x*log(1+exp(-sqrt(x*x+y*y)));
}

double jf (const double x) {
  double sum = 0;
  const double dx = 0.005;
  for (int i = 0; i < 2000; i++) {
    sum += jf_int (dx*i, x)*dx;
  }
  return sum;
  //return 7*pow(pi,4)/360. - pow(pi*x,2)/24. - pow(x,4)*(log(pow(x,2))-cf)/32.;
}


// top and Higgs mass can be viewed as tuneable parameters in this toy model to exhibit the phase transition at different T
double mtop = 173;
double mh = 125.2;

const double pi = TMath::Pi ();
const double mw = 80.4; // W boson mass
const double mz = 91.2; // Z boson mass
const double v = 246.2; // Higgs condensate value at T=0, measured via coupling constant in Fermi theory
const double D = (2*pow (mw, 2) + pow (mz, 2) + 2*pow (mtop, 2)) / (8*pow (v, 2));
const double E = (2*pow (mw, 3) + pow (mz, 3)) / (4*pi*pow (v, 3));
const double B = (3/(64*pow (pi, 2) * pow (v, 4))) * (2*pow (mw, 4) + pow (mz, 4) - 4*pow (mtop, 4));
const double T0 = TMath::Sqrt ((pow (mh, 2) - 8*B*pow (v, 2)) / (4*D));
const double lambda0 = pow (mh, 2) / (2*pow (v, 2));

const double cb = 1.5 - 2*TMath::EulerGamma () + 2*TMath::Log (2*pi);
const double cf = 1.5 + 2*TMath::Log(pi) - TMath::EulerGamma ();
const double Cb = cb - 1.5;
const double Cf = cf - 1.5;

double lambda (const double temp) {
  if (temp == 0) {
    cout << "error, invalid temperature!" << endl;
  }
  double r = lambda0;
  r -= (3/(16*pi*pi*v*v*v*v))*(2*mw*mw*mw*mw* TMath::Log (mw*mw/(Cb*temp*temp)));
  r -= (3/(16*pi*pi*v*v*v*v))*(mz*mz*mz*mz* TMath::Log (mz*mz/(Cb*temp*temp)));
  r += (3/(16*pi*pi*v*v*v*v))*(4*mtop*mtop*mtop*mtop* TMath::Log (mtop*mtop/(Cb*temp*temp)));
  return r;
}

double potential (const double phic, const double temp) {
  return D*(temp*temp-T0*T0)*phic*phic - E*temp*phic*phic*phic + 0.25*lambda(temp)*phic*phic*phic*phic;
}

double potential_full (const double phic, const double temp) {
  double sum = 0;
  sum += 6*jb (mw*phic/(v*temp));
  sum += 3*jb (mz*phic/(v*temp));
  sum -= 6*jf (mtop*phic/(v*temp));
  sum *= pow (temp, 4) / (2*pi*pi);
  return sum;
}

double scalar_phic (const double _Tc, const double _T) {
  if (_T < 0 || _T > _Tc)
    return 0;
  const double radicand = 0.25*pow (_Tc,2) - 0.25*pow(_T,2);
  if (radicand < 0)
    return 0;
  return sqrt (radicand);
}



void MakePaperPlots () {

  cout << "Input parameters:" << endl;
  cout << "\tm_W = " << mw << " GeV" << endl;
  cout << "\tm_Z = " << mz << " GeV" << endl;
  cout << "\tm_top = " << mtop << " GeV" << endl;
  cout << "\tm_H = " << mh << " GeV" << endl;
  cout << "\tv = " << v << " GeV" << endl;

  cout << "Extracted model parameters:" << endl;
  cout << "\tD = " << D << endl;
  cout << "\tE = " << E << endl;
  cout << "\tT_0 = " << T0 << endl;
  cout << "\tB = " << B << endl;
  const double TC = sqrt ((lambda0*D*T0*T0)/(lambda0*D - E*E));
  cout << "\tT_C ~= " << TC << " GeV" << endl;
  const double T1 = sqrt((-8*lambda0*D*T0*T0)/(9*E*E - 8*lambda0*D));
  cout << "\tT_1 ~= " << T1 << " GeV" << endl;

  const int npoints = 800;

  const int ntemps = 8;
  double temps[ntemps] = {163.2, 163.3, 163.4, 163.5, 163.6, 163.6, 163.7, 163.8};
  double phic_vals[ntemps][npoints] = {};
  double v_vals[ntemps][npoints] = {};
  for (int t = 0; t < ntemps; t++) {
    for (int i = 0; i < npoints; i++) {
      phic_vals[t][i] = (0.4*i/npoints) * TC;
      v_vals[t][i] = potential (phic_vals[t][i], temps[t]) * pow (0.005068, 3);
      //cout << "phi_c = " << phic_vals[t][i] << ", V(phi_c) = " << v_vals[t][i] << endl;
    }
  }

  TCanvas* c1 = new TCanvas ("c1", "", 1000, 1000);
  c1->SetLeftMargin (0.18);
  c1->SetRightMargin (0.03);
  c1->SetTopMargin (0.03);

  TGraph** g = new TGraph*[ntemps];
  for (int t = 0; t < ntemps; t++) {
    g[t] = new TGraph (npoints, phic_vals[t], v_vals[t]);
    g[t]->SetMarkerColor (colors[t+1]);
    g[t]->SetLineColor (colors[t+1]);
    g[t]->SetMarkerSize (0.8);

    g[t]->GetXaxis ()->SetTitle ("#phi_{c} [GeV]");
    g[t]->GetYaxis ()->SetTitle ("V_{eff} (#phi_{c}, T) [GeV / am^{3}]");

    g[t]->GetYaxis ()->SetTitleOffset (1.2*g[t]->GetYaxis ()->GetTitleOffset ());

    g[t]->GetXaxis ()->SetRangeUser (0, 60);
    g[t]->GetYaxis ()->SetRangeUser (-0.005, 0.005);

    g[t]->Draw (t==0 ? "AP" : "P");
  }

  myText (0.25, 0.90, kBlack, "T [GeV]", 0.04);
  myMarkerTextNoLine (0.285, 0.86, colors[1], kFullCircle, Form ("%g", temps[0]),  2.50, 0.04);
  myMarkerTextNoLine (0.285, 0.82, colors[2], kFullCircle, Form ("%g", temps[1]), 2.50, 0.04);
  myMarkerTextNoLine (0.285, 0.78, colors[3], kFullCircle, Form ("%g", temps[2]), 2.50, 0.04);
  myMarkerTextNoLine (0.285, 0.74, colors[4], kFullCircle, Form ("%g", temps[3]), 2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.86, colors[5], kFullCircle, Form ("%g", temps[4]),  2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.82, colors[6], kFullCircle, Form ("%g", temps[5]), 2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.78, colors[7], kFullCircle, Form ("%g", temps[6]), 2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.74, colors[8], kFullCircle, Form ("%g", temps[7]), 2.50, 0.04);

  //myText (0.25, 0.34, kBlack, "m_{h} = 50 GeV", 0.04);
  //myText (0.25, 0.30, kBlack, "m_{top} = 120 GeV", 0.04);

  c1->SaveAs ("sm_effpot.pdf");


  TCanvas* c2 = new TCanvas ("c2", "", 1000, 1000);
  c2->SetLeftMargin (0.10);
  c2->SetRightMargin (0.03);
  c2->SetTopMargin (0.03);

  double xvals[400] = {};
  double jb_vals[400] = {};
  double jf_vals[400] = {};

  for (int i = 0; i < 400; i++) {
    xvals[i] = 5.*(i+1)/400.;
    jb_vals[i] = jb(xvals[i]);
    jf_vals[i] = jf(xvals[i]);
  }

  TGraph* g_jb = new TGraph (400, xvals, jb_vals);
  TGraph* g_jf = new TGraph (400, xvals, jf_vals);

  g_jb->GetXaxis ()->SetTitle ("x = m(#phi_{c})/T");
  g_jf->GetXaxis ()->SetTitle ("x = m(#phi_{c})/T");

  g_jb->GetXaxis ()->SetRangeUser (0, 5);
  g_jf->GetXaxis ()->SetRangeUser (0, 5);
  g_jb->GetYaxis ()->SetRangeUser (-2.8, 2.5);
  g_jf->GetYaxis ()->SetRangeUser (-2.8, 2.5);

  g_jb->SetMarkerColor (colors[1]);
  g_jb->SetLineColor (colors[1]);
  g_jb->SetMarkerSize (0.8);
  g_jf->SetMarkerColor (colors[2]);
  g_jf->SetLineColor (colors[2]);
  g_jf->SetMarkerSize (0.8);

  g_jb->Draw ("AP");
  g_jf->Draw ("P");

  TLine* l0 = new TLine (0, 0, 5, 0);
  l0->SetLineStyle (2);
  l0->SetLineWidth (1);
  l0->Draw ("same");

  myText (0.40, 0.35, colors[1], "J_{B} (x^{2})", 0.04);
  myText (0.35, 0.85, colors[2], "J_{F} (x^{2})", 0.04);

  c2->SaveAs ("thermalfuncs.pdf");


  TCanvas* c3 = new TCanvas ("c3", "", 1000, 1000);
  c3->SetLeftMargin (0.18);
  c3->SetRightMargin (0.03);
  c3->SetTopMargin (0.03);

  temps[0] = 20;
  temps[1] = 50;
  temps[2] = 100;
  temps[3] = 150;
  temps[4] = 200;
  temps[5] = 250;
  temps[6] = 300;
  temps[7] = 400;
  for (int t = 0; t < ntemps; t++) {
    for (int i = 0; i < npoints; i++) {
      phic_vals[t][i] = (0.4*(i+1)/npoints) * TC;
      v_vals[t][i] = potential_full (phic_vals[t][i], temps[t]) * pow (0.005068, 3);
      //cout << "phi_c = " << phic_vals[t][i] << ", V(phi_c) = " << v_vals[t][i] << endl;
    }
  }

  TGraph** g_full = new TGraph*[ntemps];
  for (int t = 0; t < ntemps; t++) {
    g_full[t] = new TGraph (npoints, phic_vals[t], v_vals[t]);
    g_full[t]->SetMarkerColor (colors[t+1]);
    g_full[t]->SetLineColor (colors[t+1]);
    g_full[t]->SetMarkerSize (0.8);

    g_full[t]->GetXaxis ()->SetTitle ("#phi_{c} [GeV]");
    g_full[t]->GetYaxis ()->SetTitle ("V_{eff} (#phi_{c}, T) [GeV / am^{3}]");

    g_full[t]->GetYaxis ()->SetTitleOffset (1.2*g_full[t]->GetYaxis ()->GetTitleOffset ());

    g_full[t]->GetXaxis ()->SetRangeUser (0, 60);
    g_full[t]->GetYaxis ()->SetRangeUser (-0.005, 0.005);

    g_full[t]->Draw (t==0 ? "AP" : "P");
  }

  myText (0.25, 0.90, kBlack, "T [GeV]", 0.04);
  myMarkerTextNoLine (0.285, 0.86, colors[1], kFullCircle, Form ("%g", temps[0]),  2.50, 0.04);
  myMarkerTextNoLine (0.285, 0.82, colors[2], kFullCircle, Form ("%g", temps[1]), 2.50, 0.04);
  myMarkerTextNoLine (0.285, 0.78, colors[3], kFullCircle, Form ("%g", temps[2]), 2.50, 0.04);
  myMarkerTextNoLine (0.285, 0.74, colors[4], kFullCircle, Form ("%g", temps[3]), 2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.86, colors[5], kFullCircle, Form ("%g", temps[4]),  2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.82, colors[6], kFullCircle, Form ("%g", temps[5]), 2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.78, colors[7], kFullCircle, Form ("%g", temps[6]), 2.50, 0.04);
  myMarkerTextNoLine (0.435, 0.74, colors[8], kFullCircle, Form ("%g", temps[7]), 2.50, 0.04);

  //myText (0.25, 0.34, kBlack, "m_{h} = 50 GeV", 0.04);
  //myText (0.25, 0.30, kBlack, "m_{top} = 120 GeV", 0.04);

  c3->SaveAs ("sm_effpot_full.pdf");


  TCanvas* c4 = new TCanvas ("c4", "", 1000, 1000);
  c4->SetLeftMargin (0.18);
  c4->SetRightMargin (0.03);
  c4->SetTopMargin (0.03);
  const int scalar_npoints = 10000;
  double temp_vals[scalar_npoints] = {};
  double scalar_phi_vals[scalar_npoints] = {};
  const double scalar_mu2 = -10000;
  const double scalar_v = 2.88;
  const double scalar_tc = sqrt(24*fabs(scalar_mu2)/scalar_v);
  cout << "scalar theory Tc = " << scalar_tc << " GeV" << endl;
  for (int i = 0; i < scalar_npoints; i++) {
    temp_vals[i] = (1.8*i/scalar_npoints);
    scalar_phi_vals[i] = scalar_phic (scalar_tc, temp_vals[i]*scalar_tc);
  }
  TGraph* g_scalar = new TGraph (scalar_npoints, temp_vals, scalar_phi_vals);
  g_scalar->SetMarkerSize (0.5);
  g_scalar->SetLineWidth(4);
  g_scalar->GetXaxis ()->SetTitle ("T / T_{c}");
  g_scalar->GetYaxis ()->SetTitle ("#LT#phi_{c}#GT(T) [GeV]");
  g_scalar->Draw ("ALP");

  myText (0.6, 0.8, kBlack, "#mu^{2} = -100^{2} GeV^{2}", 0.04);
  myText (0.6, 0.75, kBlack, "v/4! = 0.13", 0.04);
  c4->SaveAs ("scalarFieldExpectationValue.pdf");

}

#endif
