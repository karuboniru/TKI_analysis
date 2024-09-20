#pragma once
#include "tkigeneral.h"

#include <TLorentzVector.h>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

class NeutrinoEvent;

// TLorentzVector get_full_hadron_MINERvAPi0(const TKIEvent &e);
// TLorentzVector get_full_hadron_MINERvA0PI(const TKIEvent &e);

ROOT::RDF::RNode DoTKICut_MINERvA_pi0(ROOT::RDF::RNode df) ;
ROOT::RDF::RNode DoTKICut_MINERvA0PI(ROOT::RDF::RNode df) ;



ROOT::RDF::RNode CommonVariableDefinePI0(ROOT::RDF::RNode df);
ROOT::RDF::RNode CommonVariableDefine0PI(ROOT::RDF::RNode df);
