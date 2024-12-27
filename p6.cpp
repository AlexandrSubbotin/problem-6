#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
TRandom3 rnd;
const double m_Ks = 497.611;
const double m_piplus = 139.570;
const double m_piminus = 139.570;
const double E_0 = 1020.0;
const double p_threshold = 40.0; 
const double detector_radius = 30.0; 
const double detector_length = 50.0;
const double m_e =0.511;
const double m_Kl=498.7;
const double t_k = 0.895; //нс
const double t_pi=26;

double generateCosThetaKs() {
  double t;
  TF1 *f= new TF1("f","pow(sin(x),3)",0,M_PI);
  t=f->GetRandom();
  return t;
}
double generatephi(){
    double x = rnd.Uniform(0,1)*2*M_PI;
    return x;
}

double generatedtheta_pi(){
  TF1 *f= new TF1("f","sin(x)",0,M_PI);
  double t=f->GetRandom();
  return t;
}
double fdist (double E,double t,double m){
    double t_m =rnd.Exp(t);
    double b= (E/m);
    double d= (1- 1/(b*b )) *t_m*3*10;
    return d; 
}

double filter(double theta1,double phi1,double theta2,double phi2,double d, double b){
  double zs= d*cos(theta1);
  double xs = d*sin(theta1)*cos(phi1);
  double ys = d*sin(theta1)*sin(phi1);
  double zpi= b*cos(theta2);
  double xpi = b*sin(theta2)*cos(phi2);
  double ypi= b*sin(theta2)*sin(phi2);
  double x0=xs+xpi;
  double y0=ys+ypi;
  double z0=zs+zpi;
  if((x0<=30 && x0>=(-30)) || (y0<=25 && y0>=(-25)) || (z0<=30 && z0>=(-30))){
    return 0;
  }
  else return 1;
}

void func(){
    double E_k;
    TH1F *hist1 = new TH1F("h1"," угол фи пи мезона в системе Ks",100,0,2*M_PI);
    TH1F *hist2 = new TH1F("h2","угол phi пи мезона в лаб системе ",100,0,2*M_PI);
    TH1F *hist3 = new TH1F("h3"," угол theta Ks мезона",100,0,M_PI);
    TH1F *hist4= new TH1F("h4","угол theta пи мезона в лаб системе ",100,0,M_PI);
    TH1F *hist5 = new TH1F("h5"," угол phi пи мезона в системе Ks",100,0, 2*M_PI);
     TH1F *hist6 = new TH1F("h6"," угол theta пи мезона в системе Ks",100,0, 2*M_PI);
    E_k= (1 / (8*E_0))*(2*m_e * m_e + m_Ks*m_Ks+ 4*E_0*E_0- m_Kl*m_Kl);
    double p= sqrt(E_k*E_k-m_Ks*m_Ks);
    double theta;
    double E_pi= m_Ks/2;
    for(int i=0;i<10000;i++){


      theta=generateCosThetaKs();
      double pz=p*cos(theta);
      double phi=generatephi();
      double px=p*sin(theta)*cos(phi);
      double py=p*sin(theta)*sin(phi);
      TLorentzVector p1(px,py,pz,E_k);
      hist3->Fill(theta);


      TLorentzRotation boost(px/E_k,py/E_k,pz/E_k);


      double theta_pi=generatedtheta_pi();
      double phi_pi=generatephi();
      double E_pi= m_Ks/2;
      double p_pi = sqrt(E_pi*E_pi- m_piplus*m_piplus);
      double px_pi=p_pi*sin(theta_pi)*cos(phi_pi);
      double py_pi= p_pi*sin(theta_pi)*sin(phi_pi);
      double pz_pi=p_pi*cos(theta_pi);
      TLorentzVector p2(px_pi,py_pi,pz_pi,E_pi);
      hist5->Fill(phi_pi);


      TLorentzVector pi_lab= boost * p2;
      double phi_lab=std::atan(pi_lab.Py()/pi_lab.Px());
      double theta_lab=std::acos(pi_lab.Pz()/sqrt(pi_lab.Px()*pi_lab.Px()+ pi_lab.Py()*pi_lab.Py()+pi_lab.Pz()*pi_lab.Pz()));
      
      
          hist1->Fill(phi_pi);
          hist2->Fill(phi_lab);
          hist4->Fill(theta_lab);
          hist6->Fill(theta_pi);
      if( sqrt(pi_lab.Px()*pi_lab.Px()+ pi_lab.Py()*pi_lab.Py()+pi_lab.Pz()*pi_lab.Pz()) >40 ){
        double d= fdist(E_k,t_k,m_Ks);
        double b=fdist(pi_lab.E(),t_pi,m_piplus);
        if(filter(theta,phi,theta_lab,phi_lab,d,b)==1 ){
          //hist1->Fill(phi_pi);
          //hist2->Fill(phi_lab);
          //hist4->Fill(theta_lab);
        //}
      }
    }
  TCanvas *c1=  new TCanvas("c1","",800,600);
  c1->Divide(3,3);
  c1->cd(1);
  hist1->Draw();
  c1->cd(2);
  hist2->Draw();
  c1->cd(3);
  hist3->Draw();
  c1->cd(4);
  hist4->Draw();
  c1->cd(5);
  hist5->Draw();
  c1->cd(6);
  hist6->Draw();
}
int main(){
  std::cout<< sqrt(4);
}

