#ifndef __DQ_DB_H__
#define __DQ_DB_H__


//#include "DQcommon.h"
#include <Riostream.h>
#include <TMath.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TProfile2D.h>

#include "DQcommon.h"

#define uint32_t Int_t
#define uint16_t Short_t

/*
 * Object wrapper to store daq and analysis parameters on ttree user info list
 *
 */
class DAQpar : public TObject {

 private:
  uint32_t fFormatID;      // DAQ File Format Identifier
  uint32_t fFormatVersion; // DAQ File Format Version

  uint32_t fMPDs;     // Number of MPD in DAQ
  uint32_t fDaqMode;  // DAQ Mode
  uint32_t fAPVs;     // Number of total front end cards ** version<5 / derived from MPD **

  uint32_t fSamples;      // Number of samples per event ** version<5 / MPD-APV specific **
  uint32_t fCommonNoise;  // Common Noise Subtraction (on/off) ** version<5 / MPD specific **
  uint32_t fEventBuilder; // Event Builder (on/off) ** version<5 / MOD specific **
  uint32_t fBaseline;     // Common Baseline Compensation (on/off) ** version<5 / MPD specific **
  uint32_t fTrigLatency;      // Trigger Latency ** version<5 / MPD specific **
  uint32_t fCalibLatency;     // Calibration Latency ** version<5 / MPD specific **

  uint32_t fAllocatedChs;  // number of allocated channels for the root arrays (usefull to allocate the same arrays to read the root data)

  // following parameters are specific of the decoding/analysis (set by DQdecode)
  Double_t fnRMS;          // applied noise threshold as number of pedestal RMS 
  TString fGeoMapFileName; // fine name of the geometry file
  TString fPedFileName;    // pedestal file name used to subtract pedestal and estimate noise (optional)
  
  // following parameters are specific of the hit detection and track reconstruction (set by DQdisplay)
  Float_t fTimeWindow;        // trigger - signal correlation time window
  Float_t fChargeAsym;        // x/y charge asymmetry
  Float_t fChi2RenThr;        // threshold on the renormalized chi2 of the sample fit
  Float_t fMaxStripInEvent;   // expected maximum number of strips with signal in a single event 
  Float_t fMaxSignalFraction; // threshold on signal as fraction of maximum signal
  Float_t fNumberSigma;       // threshold on enhanced profile as number of rms
  Float_t fAdcThr;            // threshold on single strip adc value

 public:
  DAQpar() { 
    Reset();
  };

  ~DAQpar() {
  };

  // reset all
  void Reset() {
    fChi2RenThr = 0;
    fMaxStripInEvent = 0;
    fMaxSignalFraction = 0;
    fNumberSigma = 0;
    fAdcThr = 0;
    fChargeAsym = 0;
    fTimeWindow = 0;
    SetTrackParAll();

    fnRMS = -99999.; // do not apply threshold
    fGeoMapFileName="_NotDefined_";
    fPedFileName="_NotDefined_";

    fAllocatedChs = 0;
    fFormatID=0;
    fFormatVersion=0;
    fMPDs=0;
    fAPVs=0;
    fSamples=0;
    fDaqMode=0;
    fCommonNoise=0;
    fEventBuilder=0;
    fBaseline=0;
    fTrigLatency=0;
    fCalibLatency=0;
  };

  // get/set (cannot set 0, use reset instead);
  uint32_t FormatID(uint32_t v=0) {
    if (v!=0) {
      fFormatID = v;
    }
    return fFormatID;
  };
     
  uint32_t FormatVersion(uint32_t v=0) {
    if (v!=0) {
      fFormatVersion=v;
    }
    return fFormatVersion;
  };

  uint32_t MPDs(uint32_t v=0) {
    if (v!=0) {
      fMPDs=v;
    }
    return fMPDs;
  };

  uint32_t APVs(uint32_t v=0) {
    if (v!=0) {
      fAPVs=v;
    }
    return fAPVs;
  };

  uint32_t Samples(uint32_t v=0) {
    if (v!=0) {
      fSamples=v;
    }
    return fSamples;
  };

  uint32_t Mode(uint32_t v=0) {
    if (v!=0) {
      fDaqMode=v;
    }
    return fDaqMode;
  };

  uint32_t CommonNoise(uint32_t v=0) {
    if (v!=0) {
      fCommonNoise=v;
    }
    return fCommonNoise;
  };

  uint32_t EventBuilder(uint32_t v=0) {
    if (v!=0) {
      fEventBuilder=v;
    }
    return fEventBuilder;
  };

  uint32_t Baseline(uint32_t v=0) {
    if (v!=0) {
      fBaseline=v;
    }
    return fBaseline;
  };

  uint32_t TrigLatency(uint32_t v=0) {
    if (v!=0) {
      fTrigLatency=v;
    }
    return fTrigLatency;
  };

  uint32_t CalibLatency(uint32_t v=0) {
    if (v!=0) {
      fCalibLatency=v;
    }
    return fCalibLatency;
  };


  uint32_t AllocatedChs(uint32_t v=0) {
    if (v!=0) {
      fAllocatedChs = v;
    }
    return fAllocatedChs;
  };

  // decoding/analysis parameters

  void SetNRMS(Double_t v) {
    fnRMS = v;
  }

  Double_t GetNRMS() { return fnRMS; };

  void SetMapFileName(TString v) {
    fGeoMapFileName=v;
  };
  
  TString MapFileName() { return fGeoMapFileName; };
  
  void SetPedestalFileName(TString v) {
    fPedFileName=v;
  };
  
  TString PedestalFileName() { return fPedFileName; };
  

  // Tracking parameters

  void SetTrackParAll(Float_t adcthr=-20., Float_t maxstrip=15., Float_t chi2ren=1.5, Float_t maxfraction=0.02, Float_t numsigma=1.5, Float_t chargeasy=0.31, Float_t timewindow = 30.) {
    fAdcThr = adcthr;
    fMaxStripInEvent=maxstrip;
    fChi2RenThr=chi2ren;      
    fMaxSignalFraction=maxfraction;
    fNumberSigma=numsigma;       
    fChargeAsym = chargeasy;
    fTimeWindow = timewindow;
  };

  // get/set
  Float_t AdcThr(Float_t v=99999.) { 
    if (v!=99999.) { fAdcThr=v; }
    return fAdcThr;
  };

  Float_t StripInEvent(Float_t v=-1.) { 
    if (v>=0) fMaxStripInEvent=v;
    return fMaxStripInEvent; 
  };

  Float_t Chi2Ren(Float_t v=-1.) { 
    if (v>=0) fChi2RenThr=v;
    return fChi2RenThr;
  };

  Float_t SignalFraction(Float_t v=-1) { 
    if (v>=0) fMaxSignalFraction=v;
    return fMaxSignalFraction; 
  };

  Float_t NumberSigma(Float_t v=-1) {
    if (v>=0) fNumberSigma=v;
    return fNumberSigma;
  };

  Float_t ChargeAsymmetry(Float_t v=-1.) {
    if (v>=0) fChargeAsym = v;
    return fChargeAsym;
  };

  Float_t TimeWindow(Float_t v=-1.) {
    if (v>=0) fTimeWindow = v;
    return fTimeWindow;
  };

  // overidd virtual
  void Print(Option_t *opt="") const {
    cout << "DAQ parameters:" << endl;
    cout << "| File Format ID / Version= 0x" << hex << fFormatID << " / " << fFormatVersion << endl;
    cout << "| Total MPD / APV in DAQ= " << dec << fMPDs << " / " << fAPVs << endl;
    cout << "| DaqMode= 0x " << hex << fDaqMode << endl;
    if (fFormatVersion<5) {
      cout << "| Number of Samples/Evt= " << dec << fSamples << endl;
      cout << "| Flags: Common Noise / Event Builder / Common Baseline= " << dec
	   << " / " << fCommonNoise << " / " << fEventBuilder << " / " << dec << fBaseline << endl;
      cout << "| Latency: Trigger / Calibration= " << dec << fTrigLatency << " / " << fCalibLatency << endl;
    }
    cout << " Analysis parameters:" << endl
	 << " | Geometry File Name = " << fGeoMapFileName << endl
	 << " | Pedestal File Name= " << fPedFileName << endl; 

    cout << " Track Reconstruction Parameters:" << endl
	 << " | AdcThr     = " << fAdcThr << endl 
	 << " | MaxStrip   = " << fMaxStripInEvent << endl
	 << " | Chi2Ren    = " << fChi2RenThr << endl
	 << " | SignalFrac = " << fMaxSignalFraction << endl
	 << " | NumberSigma= " << fNumberSigma << endl
	 << " | ChargeAsymm= " << fChargeAsym << endl
	 << " | TimeWindow = " << fTimeWindow << endl;
  };

  ClassDef (DAQpar, 1) 

};

class APVdb : public TObject {
 private:
  int fI2C;       // I2C address
  int fAdc;       // FIFO (ADC) channel
  int fNumSample; // number of expected samples
  int fLatency;   // trigger latency [x25 ns]
  int fMode;      // see APV user manual

 public:
  APVdb(int gI2C=-1, int gAdc=-1, int gNumSample=-1, int gLatency=-1, int gMode=-1) {
    fI2C=gI2C;
    fAdc=gAdc;
    fNumSample=gNumSample;
    fLatency=gLatency;
    fMode=gMode;
    //  DQ_DUM("allocated with i2c/adc= %d %d\n",fI2C,fAdc);
  };

  ~APVdb() {
  };

  // get/set
  int I2C(int v=-1) { 
    if (v>=0) { fI2C=v; }
    return fI2C;
  };

  int Adc(int v=-1) { 
    if (v>=0) fAdc=v;
    return fAdc; 
  };

  int NumSample(int v=-1) { 
    if (v>=0) fNumSample=v;
    return fNumSample;
  };

  int Latency(int v=-1) { 
    if (v>=0) fLatency=v;
    return fLatency; 
  };

  int Mode(int v=-1) {
    //    if (v==0 || v==1) fMode=v;
    if (v>=0) fMode=v;
    return fMode;
  };

  void Print(Option_t *opt="") const {
    cout << "+ I2C= " << dec << setw(2) << fI2C 
	 << " ADC= " << setw(2) << fAdc
	 << " Sample= " << setw(2) << fNumSample
	 << " Latency= " << setw(3) << fLatency
	 << " Mode= 0x" << hex << fMode << endl;
  };

  ClassDef (APVdb, 1) 

};

class MPDdb : public TObject {
 private:
  int fSlot; // bus slot index
  int fNumTrigger; // number of external triggers
  int fCalibLatency; // calibration latency [x25 ns]
  int fClockPhase[2]; // clock phases value [x25 ns]
  int fBaseline;      // common baseline compensation on/off
  int fCommonMode;   // common mode noise subtraction on/off
  int fEventBuilder; // fpga event builder on/off
  int fLevelZero; // digital zero threshold 
  int fLevelOne; // digital one threshold
  int fApvCount; // number of managed APV boards 
  int fAdcGain[16]; // ADC gain

  TClonesArray fFEC;

 public:

  MPDdb() {
    ApvCount(0);
  };

  MPDdb(int napv) {
    fFEC.SetClass("APVdb",16); // initialize fFEC
    for (int i=0;i<napv;i++) {
      new(fFEC[i]) APVdb();
    }
    ApvCount(napv);
  };

  ~MPDdb() {
    fFEC.Delete();
  };

  APVdb *Apv(int idx) {
    if ((idx<0) || (idx>=fApvCount)) { return 0; };
    return (APVdb *) fFEC[idx];

  };

  // get/set
  int ApvCount(int v=-1) {
    if (v>=0) {fApvCount=v; }
    return fApvCount;
  };
  
  int Slot(int v=-1) {
    if (v>=0) { fSlot = v; }
    return fSlot;
  };

  int NumTrigger(int v=-1) {
    if (v>=0) { fNumTrigger=v; }
    return fNumTrigger;
  };

  int CalibLatency(int v=-1) {
    if (v>=0) { fCalibLatency=v; }
    return fCalibLatency; 
  };

  int ClockPhase(int idx, int v=-1) {
    int i = idx % 2;
    if (v>=0) { fClockPhase[i]=v; }
    return fClockPhase[i];
  };

  int Baseline(int v=-1) { // can be negative ???
    if (v>=0) { fBaseline=v; }
    return fBaseline;
  };

  int CommonMode(int v=-1) {
    if (v>=0) fCommonMode=v;
    return fCommonMode;
  };

  int EventBuilder(int v=-1) {
    if (v>=0) fEventBuilder=v;
    return fEventBuilder;
  };

  int LevelOne(int v=-1) {
    if (v>=0) fLevelOne=v;
    return fLevelOne;
  };

  int LevelZero(int v=-1) {
    if (v>=0) fLevelZero=v;
    return fLevelZero;
  };
  
  int AdcGain(int idx, int v=-1) {
    int i=idx%16;
    if (v>=0) fAdcGain[i]=v;
    return fAdcGain[i];
  };

  // overidd virtual
  void Print(Option_t *opt) const {
    cout << "MPD Slot= " << dec << fSlot << " Total APVs= " << fApvCount << " Int. Triggers= " << fNumTrigger << endl;
    cout << "| Clock Phase 0/1= " << fClockPhase[0] << " / " << fClockPhase[1] << endl;
    cout << "| Common Baseline= " << fBaseline << " Common Mode= " << fCommonMode << " Event Builder= " << fEventBuilder << endl;
    cout << "| Calibr. Latency= " << fCalibLatency << " Thr Level 0/1= " << fLevelZero << " / " << fLevelOne << endl;
    fFEC.Print();
  };

  ClassDef (MPDdb, 1)
 
};


/***
 * Tracking structures (not used ...!)
 *
 * Point is the basic information (strip position, strip charge, sample index)
 * Cluster is a collection of adjacent points, of a given sample
 * Hit is defined by two clusters in x and y on a given plane
 * 
 */

/**

Event
 | nmodule
 | Module[]
 |  | xc, yc, zc // from mapping - should be in the constant info
 |  | Plane[2]
 |  |  | npoint
 |  |  | Point[]  //TClonesArray
 |  |  |  | x
 |  |  |  | nsample
 |  |  |  | Sample[]   //TClonesArray
 |  |  |     | i // sample index (directly related to the sampling time)
 |  |  |     | adc_val
 |  |  |
 |  |  | ncluster
 |  |  | Cluster[] // must be associated to all samples ??? 
 |  |     | x // position centroid of strips
 |  |     | charge
 |  |     | width
 |  |     | time // sample analysis ??
 |  |     | npoint // or size
 |  |     | pPoint[] // pointer to Points
 |  |
 |  | nhit
 |  | Hit[] // hit as combination of two cluster x and y   
 |     | x,y,z // position in lab frame
 |     | pCluster[2] // pointer to cluster x and y (or simply index to them)
 |
 | ntrack
 | Track[]
    | nhit
    | pHit[] // pointer to hits or index to module/hit
    | vx, vy, vz // vertex
    | px, py, pz // momentum 
    | ...

*/

/*
 * strip sample and its index
 */

class Dsample : public TObject {

 private:
  Int_t   fI; // sample index
  Float_t fV; // sample adc value

 public:
  Dsample(Int_t i=0, Float_t v=0.) { // sample index and adc value
    fV=v;
    fI=i;
  };

  ~Dsample() {
  };

  Int_t   I() { return fI; };
  Float_t V() { return fV; }; // return charge (value of sample)
  
  void SetI(Int_t i) { fI=i; };
  void SetV(Float_t v) { fV=v; };
  void Set(Int_t i, Float_t v) { fI=i; fV=v; };

  ClassDef (Dsample, 1)  

};

/*
 * Point is a strip firing in one or more samples
 */

class Dpoint : public TObject {

 private:
  Float_t fX;    // strip index (position)
  Int_t   fN; // number of samples (probably not needed if use TClonesArray)
  //  TClonesArray fS("Dsample",6); //samples
  TClonesArray fS; //samples

  // derived quantities
  Float_t fCharge; // point charge, as maximum of sample value (or integral ?)

  Float_t fT0; // time from trigger (extracted from fit of samples)
  Float_t fTau0, fTau1; // timing constants of fitting samples
  Float_t fNorm; // normalization of fitting function (point charge ??)
  Float_t fChi2r; // Normalized Chi2 of the fit

 public:
  Dpoint(Float_t x=0.) { fX=x; fN=0; Init(); };

  Dpoint(Float_t x, Int_t s, Float_t v) { // position, sample index and adc value
    fX=x;
    fN=1;
    //    new(fS[0]) Dsample(s,v);
    Dsample *ds = (Dsample *)fS.ConstructedAt(0);
    ds->Set(s,v);
    Init();
  };

  ~Dpoint() {
    //   fS.Delete();
    fS.Clear();
  };

  void Init() {
    fCharge = 0;
    fT0=0;
    fTau0 = 0;
    fTau1 = 0;
    fNorm = 0;
    fChi2r = 0;
  }

  Dsample * S(Int_t i) { // return the given sample pointer
    if (i<fN) {
      return (Dsample *) fS[i];
    } else {
      return NULL;
    }
  };

  Float_t X() { // return position
    return fX;
  };

  Int_t N() { // return number of samples
    return fN;
  };

  Float_t V(Int_t i) { // return adc value of sample i
    Dsample *s=NULL;
    if (i<fN) {
      s = (Dsample *) fS[i];
      return s->V();
    };
    return -1;
  };

  Int_t I(Int_t i) { // sample index (basically the time)
    Dsample *s=NULL;
    if (i<fN) {
      s = (Dsample *) fS[i];
      return s->I();
    };
    return -1;
  };

  void Add(Int_t s, Float_t v) {
    //    new(fS[fN]) Dsample(s,v);
    Dsample *ds = (Dsample *)fS.ConstructedAt(fN);
    ds->Set(s,v);
    fN++;
  };

  Float_t Charge() {
    return fCharge;
  };

  Float_t EvalCharge() {
    Float_t max=-1E9;
    for (Int_t i=0;i<N();i++) {
      Float_t y = V(i);
      max = (y > max) ? y : max; 
    }
    
    fCharge = max;
    return fCharge;
  };

  Float_t EvalFit(Float_t xmin, Float_t xmax, // return chi2 
		  Float_t Tau0, Float_t Tau1,  // time constans, if negative, becomes fix parameters 
		  Float_t Norm, Float_t T0) { // fit the samples shape
    TF1 *pulse_fit = new
      TF1("pulse_fit","[2]*(1.-exp(-(x-[3])/[0]))*exp(-(x-[3])/[1])*(TMath::Sign(0.5,x-[3])+0.5)",
	  xmin,xmax);
    if (Tau0<0) {
      Tau0=-Tau0;
      pulse_fit->FixParameter(0,Tau0);
    } else {
      pulse_fit->SetParameter(0,Tau0);
    }

    if (Tau1<0) {
      Tau1=-Tau1;
      pulse_fit->FixParameter(1,Tau1);
    } else {
      pulse_fit->SetParameter(1,Tau1);
    }

    pulse_fit->SetParameter(2,Norm);
    pulse_fit->SetParameter(3,T0);

    TGraph *gph = new TGraph(N());

    Float_t sum=0;
    for (Int_t i=0;i<N();i++) {
      Float_t x,y;
      x = (Float_t) I(i);
      y = V(i);
      gph->SetPoint(i, x, y);
      sum += y;
    }

    gph->Fit(pulse_fit,"QR");

    fT0 = pulse_fit->GetParameter(3);
    fNorm = pulse_fit->GetParameter(2);
    fTau1 = pulse_fit->GetParameter(1);
    fTau0 = pulse_fit->GetParameter(0);
    
    fChi2r = pulse_fit->GetChisquare()/((Float_t) pulse_fit->GetNDF());

    return fChi2r;

  };

  ClassDef (Dpoint, 1)  
  
};

//

class Dcluster : public TObject {

 private:
  Int_t fN; // number of points
  Int_t * fP; //[fN]
  Float_t fX; // centroid position
  Float_t fC; // charge
  Float_t fR; // width of the cluster
  Float_t fT; // time

  Int_t ipcur; // index of point

 public:

  Dcluster(Int_t n, Float_t x=0., Float_t c=0., Float_t r=0., Float_t t=0.) {
    fN=n;
    fP = new Int_t[fN];
    Set(x,c,r,t);
    ipcur = 0;
  };

  ~Dcluster() {
    delete fP;
  };

  Int_t P(int idx) { // index to point idx, -1 if out of range
    if ((idx<0) || (idx>=fN)) {
      return -1;
    } else {
      return fP[idx];
    }
  };

  Int_t N() { return fN; };   // number of points (size of cluster)
  Float_t X() { return fX; }; // position
  Float_t C() { return fC; }; // charge
  Float_t R() { return fR; }; // width of cluster
  Float_t W() { return fR; }; // as before
  Float_t T() { return fT; }; // time

  void SetP(int idx, int ip) { // set point reference (idx=index in Dcluster, ip=point index in Dplane)
    if ((idx<0) || (idx>=fN)) { return; };
    fP[idx]=ip;
  };
  
  void AddP(int ip) {
    if (ipcur<fN) {
      fP[ipcur]=ip;
      ipcur++;
    }
  };

  void Set(Float_t x, Float_t c, Float_t r, Float_t t=0.) { // position, charge, width and time
    fX=x;
    fC=c;
    fR=r;
    fT=t;
  };

  ClassDef (Dcluster, 1)

};

//

class Dplane : public TObject {

 private:
  Int_t fNP; // number of points
  TClonesArray  fP; // points
  Int_t fNC; // number of clusters
  TClonesArray fC; // clusters

 public:

  Dplane() {
    fNP=0;
    fNC=0;
  };

  ~Dplane() {
    if (fNP>0) { fP.Delete(); }
    if (fNC>0) { fC.Delete(); }
  };
  
  Dpoint * AddPoint(Float_t x) {
    new(fP[fNP]) Dpoint(x);
    fNP++;
    return (Dpoint *) fP[fNP-1];
  };

  Dpoint * AddPoint(Float_t x, Int_t s, Float_t v) {
    new(fP[fNP]) Dpoint(x, s, v);
    fNP++;
    return (Dpoint *) fP[fNP-1];
  };

  Dpoint * GetPoint(Int_t ip) { // return point i
    if (ip < fNP) {
      return (Dpoint *) fP[ip];
    };
    return NULL;
  };

  Float_t GetPointPosition(Int_t ip) {
    if (ip < fNP) {    
      return ((Dpoint *) fP[ip])->X();
    };
    return 0.; // @@@
  };

  Float_t GetPointCharge(Int_t ip, Int_t is) { // return charge of sample is of point ip
    if (ip < fNP) {    
      return ((Dpoint *) fP[ip])->V(is);
    };
    return 0.; // @@@
  };

  void AddSampleToPoint(Int_t i, Int_t s, Float_t v) {
    if (i<fNP) {
      ((Dpoint *) fP[i])->Add(s,v);
    }
  };

  Dcluster * AddCluster(Int_t n, Float_t x=0., Float_t c=0., Float_t r=0., Float_t t=0.) {
    new(fC[fNC]) Dcluster(n, x, c, r, t);
    fNC++;
    return (Dcluster *) fC[fNC-1];
  };

  void AddPoint2Cluster(Int_t ic, Int_t ip) {
    if (ic<fNC) {
      ((Dcluster *) fC[ic])->AddP(ip);
    }
  };

  void SetPoint2Cluster(Int_t ic, Int_t idx, Int_t ip) {
    if (ic<fNC) {
      ((Dcluster *) fC[ic])->SetP(idx, ip);
    }
  };

  void EvalCluster(Int_t ic) { // evaluate parameters X,C and R of cluster - points shall be set
    Float_t sum,sum2,mean;
    Float_t rms;
    Float_t norm;
    norm=0.;
    sum=0.;
    sum2=0.;
    Dcluster *cc = (Dcluster *) fC[ic];

    Int_t np = cc->N();
    for (Int_t i=0; i<np; i++) {
      Float_t v = ((Dpoint *) fP[i])->Charge(); // charge
      Float_t x = ((Dpoint *) fP[i])->X(); // position
      sum += v*x;
      sum2 += v*x*x;
      norm += v;
    }

    mean = 0;
    rms=0;

    if (norm!=0) {
      mean = sum / norm;
      rms = TMath::Sqrt((sum2*norm - sum*sum)/norm/norm);
    };

    cc->Set(mean, norm, rms);

  };

  ClassDef (Dplane, 1)

};

//

class Dhit : public TObject {

 private:
  Dcluster *fCX; // cluster X
  Dcluster *fCY; // cluster Y
  Float_t fX;    // x lab coordinate
  Float_t fY;    // y lab coordinate
  Float_t fZ;    // z lab coordinate
  Float_t fT;    // time

 public:

  Dhit() {
  };

  ~Dhit() {
  }

  void SetCluster(Dcluster *cx, Dcluster *cy) {
    fCX = cx;
    fCY = cy;
  };

  void SetPosition(Float_t x, Float_t y, Float_t z) {
    fX = x;
    fY = y;
    fZ = z;
  };

  void SetTime(Float_t t) {
    fT = t;
  };

  Float_t X() {
    return fX;
  };
  Float_t Y() {
    return fY;
  };
  Float_t Z() {
    return fZ;
  };

  Float_t T() {
    return fT;
  };

  Dcluster * XCluster() {
    return fCX;
  };

  Dcluster * YCluster() {
    return fCY;
  };

  ClassDef (Dhit, 1)

};

//

class Dmodule : public TObject {

 private:
  Float_t fX, fY, fZ;
  Dplane fPlane[2];
  Int_t fNH; // number of hits
  TClonesArray fH; // hits

 public:

  Dmodule(Float_t xc, Float_t yc, Float_t zc) {
    fX = xc;
    fY = yc;
    fZ = zc;
  };

  ~Dmodule() {
  };

 void PickPoints(Int_t i, Float_t x, Int_t is, Float_t v) {
  };

 void SearchClusters() {
  };

  ClassDef (Dmodule, 1)
};

/**
 * Simpler DB for track reconstruction
 *
 */

class eTCluster : public TObject {

 private:
  Int_t fN; // number of points
  Int_t fSmax; // sample with max charge in peak
  Int_t fIC; // index ordering the cluster array in correlated doublet on the x and y planes
  Float_t fX; // centroid position
  Float_t fC; // charge
  Float_t fCmax; // charge max
  Float_t fW; // width of the cluster (standard deviation of largest peak sample)

  Float_t fpE; // peak enhanced value

  Float_t ft0; // time
  Float_t fA; // normalization
  Float_t ftau0; // first time constant
  Float_t ftau1; // second time constant
  Float_t fChi2; // normalized chi2 of the samples fit


 public:

  /*
   * Fill a 1D cluster given:
   *  number of hits composing the cluster
   *  hits pointer
   *  number of samples for each hit
   *  time of each sample
   *  values of the strip axis for each hit
   */

  eTCluster() {
    fpE = -1;
    ft0 = 0.;
    fIC = -1;
  };

  ~eTCluster() { };

  virtual Bool_t IsSortable() const {
    return kTRUE; 
  };

  Int_t Compare(const TObject *obj) const { // relative to the charge
    eTCluster * eClu = (eTCluster *) obj;
    if (fC == eClu->C()) { return 0; } // probably not really necessary
    return ((fC > eClu->C()) ? -1 : 1); // sort in descending order
  };

  void SetMainPar(Int_t nhit, Int_t smax, Float_t xcentro, Float_t charge, Float_t mcharge, Float_t width) {
    fN = nhit;
    fSmax = smax;
    fX = xcentro;
    fC = charge;
    fCmax = mcharge;
    fW = width;
  };

  void SetFitPar( Float_t gt0, Float_t gnorm, Float_t gtau0, Float_t gtau1, Float_t gchi2) {
    ft0   = gt0; 
    fA    = gnorm;
    ftau0 = gtau0; // first time constant
    ftau1 = gtau1; // second time constant
    fChi2 = gchi2; // chi2 of the samples fit
  }

  void SetIC( Int_t i) { fIC = i; };

  Int_t N() { return fN; };   // number of points (physical size of cluster)
  Int_t IC() { return fIC; }; // return the index which correlate to the cluster in the other plane
  Float_t X() { return fX; }; // position
  Float_t C() { return fC; }; // charge
  Float_t W() { return fW; }; // width of cluster
  Float_t t0() { return ft0; }; // time
  Float_t norm() { return fA; }; // normalization (charge ??) 
  Float_t tau0() { return ftau0; };
  Float_t tau1() { return ftau1; };
  Float_t chi2() { return fChi2; };
  Float_t peak() { return fpE; };

  ClassDef (eTCluster, 1)

};

/*
 * pedestals and other single channel electronic calibration values (electronics calibration)
 *
 */

class DQeca : public TObject {

 private:
  Int_t c0;         // first channel - absolute channel index for ped and rms
  Int_t nc;         // number of contiguous channels for ped and rms
  Int_t nh;         // number of header indeces
  Double_t *pc; // pedestal values depending on digital header frame configuration
  Double_t *rc;      // rms of pedestals
  Double_t *gc;      // mask of good channel (1=good, 0=to be masked)
  Double_t def_ped;  // default velues
  Double_t def_rms;
  Double_t def_good;

  // ach = absolute channel index
  // hidx = header index
  void setPed(Int_t ach, Int_t hidx, Double_t ped, Double_t rms) {
    Int_t ich = ach-c0; 
    if ((ich>=0) && (ich<nc)) {
      if ((hidx>=0) && (hidx<nh)) {
	pc[ich*nh+hidx]=ped;
	rc[ich*nh+hidx]=rms;
      }
    }
  };
  
  // ich = absolute channel index
  void setMask(Int_t ach, Double_t good) {
    Int_t ich = ach-c0; 
    if ((ich>=0) && (ich<nc)) {
      gc[ich]=good;
    }
  };

 public:

  DQeca(TString ifile, Double_t dped=0., Double_t drms=0., Double_t dgood=1.) {

    Double_t dd1,dd2;
    Int_t cg0;        // first channel - absolute channel index for mask
    Int_t ngc;        // number of contiguous channels for mask

    // init...
    def_ped = dped;
    def_rms = drms;
    def_good = dgood;
    rc = 0;
    pc = 0;
    gc = 0;

    nc = 0;
    ngc = 0;
    nh = 0;
    c0 = 0;
    cg0 = 0;

    TFile *f = new TFile(ifile.Data(),"READ");
    if (f->IsZombie() == 0) {
      DQ_MSG("Channel Calibration data from file %s\n",ifile.Data()); 
      f->ls();
      // try read standard rms
      TGraph *grms = new TGraph();
      Int_t r_rms = grms->Read("rms");

      // try read standard ped
      TGraph *gped = new TGraph();
      Int_t r_ped = gped->Read("profile");

      // try read mask (good/bad channels)
      TGraph *gmask = new TGraph();
      Int_t r_mask = gmask->Read("goodch");

      // try read header dependent pedestal and rms
      TProfile2D *hped = new TProfile2D();
      Int_t r_hped = hped->Read("headped");
      Int_t n_hpedx = 0;
      Int_t n_hpedy = 0;

      DQ_MSG("Got RMS %d, Ped %d, Mask %d, HPed %d (!=0 means that has been read)\n",r_rms, r_ped, r_mask, r_hped);

      if (r_hped !=0) { // load header-dependent pedestal and rms
	r_ped = 0;  // hped replace "standard" pedestal profile
	r_rms = 0;
	n_hpedx = hped->GetNbinsX();
	n_hpedy = hped->GetNbinsY();
	nc = n_hpedx;
	nh = n_hpedy;
	c0 = (Int_t) hped->GetXaxis()->GetBinCenter(1);
	DQ_MSG("header dependent pedestals for %d channels, %d header indeces, first channel = %d\n",nc, nh, c0);
	pc = new Double_t[nc*nh];
	rc = new Double_t[nc*nh];
	for (Int_t i=0;i<nc;i++) {
	  for (Int_t j=0;j<nh;j++) {
	    Int_t ibin = hped->GetBin(i+1,j+1);
	    setPed(i+c0,j,hped->GetBinContent(ibin), hped->GetBinError(ibin));
	  }
	}
      }
      
      if ((r_rms != 0) && (r_ped != 0) && (rc == 0) && (pc == 0)) { // load "standard" ped/rms only if header-ped non loaded (previous block)
	nh = 1;
	nc = grms->GetN();
	c0 = (Int_t) (grms->GetX())[0];
	DQ_MSG("TGraph ped/rms loaded %d channels, first channel = %d\n",nc,c0); 
	rc = new Double_t[nc];
	pc = new Double_t[nc];
	gc = new Double_t[nc];
	for (Int_t i=0;i<nc;i++) {
	  setPed(i+c0,0,(gped->GetY())[i],(grms->GetY())[i]);
	}
      }

      if (r_mask != 0) { // good-bad channel mask
	ngc = gmask->GetN();
	cg0 = (Int_t) (gmask->GetX())[0];
	if ((c0<cg0) || ((c0+nc)>(cg0+ngc))) { // range of channels in pedestal file must be equal or included in mask file
	  DQ_WRN("Channel range of mask data (%d %d) NOT match range of ped. data (%d %d), will use default\n",cg0,cg0+ngc-1,c0,c0+nc-1);
	} 
	if (gc == 0) { gc = new Double_t[nc]; }
	for (Int_t i=0;i<nc;i++) {
	  if (((c0+i-cg0) >= 0) || ((c0+i-cg0)<ngc)) { 
	    dd1 = (gmask->GetY())[c0+i-cg0];
	  } else {
	    dd1 = dgood;
	  }
	  setMask(i+c0, dd1);
	}
	DQ_MSG("Acquired Channel Mask channel limits: %d %d (effective channels= %d)\n",cg0, cg0+ngc-1, nc);
      }

      delete gmask;
      delete gped;
      delete grms;
      delete hped;
      
      f->Close();
      
    }
    
    delete f;

    DQ_MSG("Electronics calibration data initialized in channel range [%d %d)\n", c0, c0+nc-1);

  };
  
  ~DQeca() {
    delete[] rc;
    delete[] gc;
    delete[] pc;
    nc=0;
  };

  Double_t getPed(Int_t ach, Short_t hidx=0) { // absolute channel and header indeces
    if (pc == 0) { return def_ped; }
    Int_t hi = (hidx & 0x1ff)>>1;
    Int_t ich = ach - c0;
    if ((ich>=0) && (ich<nc)) {
      if ((hi>=0) && (hi<nh)) {
	//    cout << ach << " " << ich << " " << hidx << " " << hi << " " << pc[ich*nh+hi] << endl;
	return pc[ich*nh+hi];
      }
    }
    return def_ped;
  };

  Double_t getRMS(Int_t ach, Short_t hidx=0) { // absolute channel and header indeces
    if (rc == 0) { return def_rms; }
    Int_t hi = (hidx & 0x1ff)>>1;
    Int_t ich = ach - c0;
    if ((ich>=0) && (ich<nc)) {
      if ((hi>=0) && (hi<nh)) {
	return rc[ich*nh+hi];
      }
    }
    return def_rms;
  };

  Double_t getGood(Int_t ach) { // absolute channel index
    if (gc == 0) { return def_good; }
    Int_t ich = ach - c0;
    if ((ich>=0) && (ich<nc)) {
      return gc[ich];
    }
    return def_good;
  };

  bool isGood(Int_t ach) {
    if (getGood(ach) > 0.) {
      return true;
    } 
    return false;
  };

  ClassDef (DQeca, 1)

};


#endif


