#include <algorithm>
#include <boost/program_options.hpp>
#include <format>
#include <fstream>
#include <iostream>

char bool_to_char(bool b) { return b ? 'T' : 'F'; }

struct neutrino_induced {
  int nuExp;
  std::string FileNameFlux;
  bool include2pi{true};
};

struct target {
  int Z;
  int A;
};

struct input {
  int numEnsembles = 4000;
  int numTimeSteps = 300;
  double delta_T = 0.1;
  int num_runs_SameEnergy = 10;
  std::string path_to_input;
  std::string version;
};

struct nl_neweN {
  int T;
  int T_2p2h;
};

struct width_Baryon {
  bool mediumSwitch;
};

struct XsectionRatios_input {
  bool flagInMedium;
  int InMediumMode;
  double alpha;
};

const char *const common_input = R"(
&nl_fluxcuts
      Enu_lower_cut    = 0.0
      Enu_upper_cut    = 50.0  
      energylimit_for_Qsrec = T 
/

&initDensity
      densitySwitch=2
/

&initPauli
      pauliSwitch=2             ! 1=dynamic, 2=analytic
/

&nl_neutrinoxsection
      singlePiModel=2           ! 0 for HNV, 1 for MAID, 2 for Bosted-Christy
      invariantMasscut=5.       ! cut events with invariant mass above

      DIScutW=3.0               ! DIS sets in at this mass (in GeV)
      DIScutwidth=0.1

      invariantMasscut=20.      ! mass cuts on Wrec
      invariantMassCut_BG = 20. ! mass cuts on Wrec for background
/

&spectralFunction
      nuclwidth=0.001           ! if mediumSwitch_coll=.false.: this parameter determines the "bare nucleon width", for numerical purpose
      relativistic=.true.       ! if mediumSwitch_coll=.false.: relativistic or non-rel. spectral function
/

&master_2Body
      doScaleResidue = .true.   ! prevents nucleon knock-out with numbers larger than target mass
                                ! requires DetermineResidue=.true.
/

&residue_Input
      DetermineResidue = .true. ! needed for rescaling NN X-sections
/

&baryonPotential
      EQS_Type=5                ! 1=soft mom-dep, 2=hard mom-dep, 3=soft non-mom-dep, 4=hard non-mom-dep, 5=medium mom-dep
      DeltaPot=1                ! 1=2/3 of nucleon potential, 2=100 MeV*rho/rhoNull
/

&insertion
      minimumEnergy=0.005
/

&neutrinoAnalysis
      outputEvents=.true.       ! output list of events and all outgoing particles in each event to the file FinalEvents.dat

      include_W_dist=.false.     ! printout of W-distributions for outgoing channel
	    DoOutChannels = .false.    ! gives cross sections for various final states
	    ZeroPion_analysis=.false.  ! extra cross sections  for events with 0 pions in the final state
	    detailed_diff_output=.false. ! differential xsec; see namelist detailed_diff for max values and bins of the histograms

! now various exp. cuts on final state nucleons and pions; these may affect the FinalEvents file
! should be turned off for own analysis

      kineticEnergyDetectionThreshold_nucleon=0.
      kineticEnergyDetectionThreshold_chargedpion=0.
      kineticEnergyDetectionThreshold_neutralpion=0.

      applyCuts = 2    ! here one may specify cuts, e.g. if particles are bound
                       ! =2: only unbound particles, =0: all particles

      reconstruct_neutrino_energy=.false. ! calculate neutrino energy reconstruction for various final states (see nl_specificEvent)
      specificEvent_analysis=.false.      ! extra cross sections for specific final states; see namelist nl_specificEvent
                                         ! also provides Q2 and energy reconstruction
/

&detailed_diff
      ekinMax=5.0               ! maximum value of Ekin in the output of distribution of hadrons versus kinetic energy
      dEkin=0.02                ! binning of Ekin in the output versus kinetic energy
      fornucleon=.false.         ! produce output versus kinetic energy, angle, ... for this outgoing hadron
      forpion=.false.            ! to the files diff_....
      foreta=.false.
      forkaon=.false.
      forkaonBar=.false.
      forLambda=.false.
      forSigmaResonance=.false.
/

&pythia
      PARP(91)=0.44  ! PYTHIA Parameter tuned to low energies
/

)";

std::string build_neutrino_induced(const neutrino_induced &ni) {
  return std::format(
      R"(
&neutrino_induced
      process_ID      =  2
      flavor_ID       =  2
      nuXsectionMode  = 16
      nuExp           =  {}
      FileNameFlux    = '{}'

      includeQE       = T
      includeDELTA    = T
      includeRES      = T
      include1pi      = T
      includeDIS      = T
      include2p2hQE   = T
      include2pi      = {}

      printAbsorptionXS = F
/
)",
      ni.nuExp, ni.FileNameFlux, bool_to_char(ni.include2pi));
}

std::string build_target(const target &t) {
  return std::format(
      R"(
&target
      Z={}
      A={}
      densitySwitch_Static=2 
      fermiMomentum_input=0.225
      fermiMotion=.true.
      ReAdjustForConstBinding=.true.
      ConstBinding=-0.008 		
/
)",
      t.Z, t.A);
}

std::string build_input(const input &i) {
  return std::format(
      R"(
&input
      eventtype=5 

      numEnsembles={} 
      numTimeSteps={} 

      delta_T={}  

      fullensemble=.true.  
      localEnsemble=.true. 

      num_runs_SameEnergy={}
      num_Energies=1  

      printParticleVectors=.false.
      freezeRealParticles=.true.
      LRF_equals_CALC_frame=.true. !if .false.: no offshelltransport possible

      path_to_input='{}'
    version = {}
/
)",
      i.numEnsembles, i.numTimeSteps, i.delta_T, i.num_runs_SameEnergy,
      i.path_to_input, i.version);
}

std::string build_nl_neweN(const nl_neweN &n) {
  return std::format(
      R"(
&nl_neweN
      new_eN = T
!     new_eNres = F

      T = {}
      T2p2h = {}
      ME_ODW = 2
/
)",
      n.T, n.T_2p2h);
}

std::string build_width_Baryon(const width_Baryon &w) {
  return std::format(
      R"(
&width_Baryon
      mediumSwitch = {0}
      mediumSwitch_Delta={0}
      mediumSwitch_coll=.false.
/
)",
      bool_to_char(w.mediumSwitch));
}

std::string build_XsectionRatios_input(const XsectionRatios_input &x) {
  return std::format(
      R"(
&XsectionRatios_input
      flagScreen = .false.      ! turns off any in-medium changes of nucleon-nucleon cross sections  
      ScreenMode = 1            ! in-medium screening of NN X-sections, following Li and Machleidt

      flagInMedium = {}         ! turns off any in-medium changes in deciding to keep an event
      InMediumMode = {}          ! BB -> BB* is reduced following Song,Ko, PRC91, 014901 (2015), NN elastic reduced according to Machleidt-Li
      alpha = {}               ! parameter in the Song-Ko exponential suppression of NN->NDelta Xsection
/
)",
      bool_to_char(x.flagInMedium), x.InMediumMode, x.alpha);
}

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  // clang-format off
      desc.add_options()
                  ("path-to-input", po::value<std::string>()->required(), "Path to BUUInput")
                  ("output-file",   po::value<std::string>()->required(), "Output JobCard file")
                  ("experiment",    po::value<std::string>()->required(), "minerva or microboone")
                  ("T",             po::value<int>()->default_value(1),         "T value")
                  ("T2p2h",         po::value<int>()->default_value(0),         "T value for 2p2h process")
                  ("2pibg",         po::value<bool>()->required(),              "Include 2pi BG")
                  ("mediumSwitch",  po::value<bool>()->required(),              "mediumSwitch")
                  ("flagInMedium",  po::value<bool>()->required(),              "NN Cross Section (Machleidt-Li and Song-Ko)")
                  ("InMediumMode",  po::value<int>()->default_value(2),         "Oset width modification for Delta")
                  ("alpha",         po::value<double>()->default_value(1.2),      "alpha for Song-Ko")
                  ("version",         po::value<std::string>()->default_value("2025"),      "GiBUU version");
  // clang-format on
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << e.what() << '\n';
    std::cout << desc << '\n';
    return 1;
  }
  std::string path_to_input = vm["path-to-input"].as<std::string>();
  std::string output_file = vm["output-file"].as<std::string>();
  std::string experiment = vm["experiment"].as<std::string>();
  int T = vm["T"].as<int>();
  int T2p2h = vm["T2p2h"].as<int>();
  bool include2pi = vm["2pibg"].as<bool>();
  bool mediumSwitch = vm["mediumSwitch"].as<bool>();
  bool flagInMedium = vm["flagInMedium"].as<bool>();
  int InMediumMode = vm["InMediumMode"].as<int>();
  double alpha = vm["alpha"].as<double>();

  // upper to lower for experiment
  std::ranges::transform(experiment, experiment.begin(),
                         [](unsigned char c) { return std::tolower(c); });

  std::ofstream out(output_file, std::ios::trunc);
  auto version = vm["version"].as<std::string>();

  if (experiment == "minerva") {
    out << build_neutrino_induced(
               {.nuExp = 25, .FileNameFlux = "", .include2pi = include2pi})
        << build_target({.Z = 6, .A = 12})
        << build_input({.path_to_input = path_to_input, .version = version});
  } else if (experiment == "t2k") {
    out << build_neutrino_induced(
               {.nuExp = 9, .FileNameFlux = "", .include2pi = include2pi})
        << build_target({.Z = 6, .A = 12})
        << build_input({.path_to_input = path_to_input, .version = version});
  } else if (experiment == "microboone") {
    out << build_neutrino_induced(
               {.nuExp = 99,
                .FileNameFlux = BASEPATH "/jobcard/microboone.dat",
                .include2pi = include2pi})
        << build_target({.Z = 18, .A = 40})
        << build_input({.numEnsembles = 3000,
                        .path_to_input = path_to_input,
                        .version = version});
  } else {
    std::cerr << "Unknown experiment: " << experiment << '\n';
    return 1;
  }

  out << common_input << build_width_Baryon({.mediumSwitch = mediumSwitch})
      << build_nl_neweN({.T = T, .T_2p2h = T2p2h})
      << build_XsectionRatios_input({.flagInMedium = flagInMedium,
                                     .InMediumMode = InMediumMode,
                                     .alpha = alpha})
      << '\n';

  return 0;
}