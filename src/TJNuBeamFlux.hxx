#ifndef JNuBeamFlux_hh_seen
#define JNuBeamFlux_hh_seen

#include <iostream>

#include "TObjString.h"
#include "TObject.h"

///
/// JNuBeam flux pass-through info. Both GRooTracker and NRooTracker inherit
/// from this so as to avoid duplication of code. Currently up to date with
/// 11a version of pass-through info.
///
namespace ND {

const int kNgmax = 12;

// This is the base class for all RooTracker type objects. For now
// it is  empty and just so that an pass round pointers to
// all derived classes (will be usefull for reweighting applications etc...)
class RooTrackerVtxBase : public TObject {
 public:
  RooTrackerVtxBase();
  virtual ~RooTrackerVtxBase() {}
  ClassDef(ND::RooTrackerVtxBase, 1);
};

class JNuBeamFlux : public RooTrackerVtxBase {
 public:
  JNuBeamFlux();
  virtual ~JNuBeamFlux();

  void Reset(void);
  // Using methods from TObject to remove 'hidden' compiler warnings
  using TObject::Copy;
  void Copy(const JNuBeamFlux* flux);

  //
  // Meta information
  //
  long NuFluxEntry;        ///< Corresponding entry in orig flux file
  TObjString* NuFileName;  ///< Original flux file name neutrino parent info
  int NuParentPdg;          ///< Parent hadron pdg code
  int NuParentDecMode;      ///< Parent hadron decay mode
  double NuParentDecP4[4];  ///< Parent hadron 4-momentum at decay
  double NuParentDecX4[4];  ///< Parent hadron 4-position at decay
  float NuCospibm;  ///< Parent hadron cos(theta) relative to beam at decay
  float NuNorm;     ///< Normalization weight
  double NuParentProP4[4];  ///< Parent hadron 4-momentum at production
  double NuParentProX4[4];  ///< Parent hadron 4-position at production
  float
      NuCospi0bm;  ///< Parent hadron cos(theta) relative to beam at production
  float NuRnu;     ///< Neutrino r position at ND5/6 plane
  float NuXnu[2];  ///< Neutrino (x,y) position at ND5/6 plane
  int NuIdfd;      ///< Detector ID primary particle information
  int NuGipart;      ///< Primary particle ID
  float NuGpos0[3];  ///< Primary particle starting point
  float NuGvec0[3];  ///< Primary particle direction at starting point
  float NuGamom0;    ///< Momentum of the primary particle at the starting point Interation History Information
  int NuNg;                ///< Number of interaction steps
  float NuGp[kNgmax][3];   ///< Momentum of ancestor particle
  float NuGcosbm[kNgmax];  ///< Ancestor particle cos(theta) relative to beam
  float NuGv[kNgmax][3];   ///< Position of ancestor particle
  int NuGpid[kNgmax];      ///< PDG ancestor particle ID
  int NuGmec[kNgmax];  ///< Particle production mechanism of ancestor particle Normalization and Transfer Matrix Information
  float NuEnusk;   ///< Neutrino Energy at SK
  float NuNormsk;  ///< Normalization weight at SK
  float NuAnorm;   ///< ND Acceptance Weight Out-of-target Secondary Interactions
  int NuGmat[kNgmax];       ///< Material in which the particle originates
  float NuGdistc[kNgmax];   ///< Distance traveled through carbon
  float NuGdistal[kNgmax];  ///< Distance traveled through aluminum
  float NuGdistti[kNgmax];  ///< Distance traveled through titanium
  float NuGdistfe[kNgmax];  ///< Distance traveled through iron Beam parameter information
  float NuVersion;   ///< Jnubeam flux version
  int NuTuneid;      ///< Beam tune ID #
  int NuNtrig;       ///< Number of triggers
  int NuPint;        ///< Interaction model ID
  float NuBpos[2];   ///< Beam center position
  float NuBtilt[2];  ///< Beam angle
  float NuBrms[2];   ///< Beam RMS width
  float NuEmit[2];   ///< Beam Emittance
  float NuAlpha[2];  ///< Beam Alpha
  float NuHcur[3];   ///< Horn currents
  int NuRand;        ///< Random seed
                     //    int         NuRseed[2];         // Random seed

  ClassDef(ND::JNuBeamFlux, 1);
};

}  // nd280 namespace
#endif

