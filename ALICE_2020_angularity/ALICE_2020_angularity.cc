// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2020_angularity : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2020_angularity);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const ChargedFinalState fs(Cuts::abseta < 0.9);
      //const ALICE::Primary primary(Cuts::abseta < 1.0 && Cuts::abscharge > 0);

      declare(FastJets(fs, FastJets::ANTIKT, 0.2), "JetsAK2");
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");

      // Book histograms
      // specify custom binning
      book(_h["ang_R02_B10"], "ang_R02_B10", 20, 0.0, 1.0);
      book(_h["ang_R02_B15"], "ang_R02_B15", 20, 0.0, 1.0);
      book(_h["ang_R02_B20"], "ang_R02_B20", 20, 0.0, 1.0);
      book(_h["ang_R02_B30"], "ang_R02_B30", 20, 0.0, 1.0);

      book(_h["ang_R04_B10"], "ang_R04_B10", 20, 0.0, 1.0);
      book(_h["ang_R04_B15"], "ang_R04_B15", 20, 0.0, 1.0);
      book(_h["ang_R04_B20"], "ang_R04_B20", 20, 0.0, 1.0);
      book(_h["ang_R04_B30"], "ang_R04_B30", 20, 0.0, 1.0);
        
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      // book(_h["angR02"], 1, 1, 1);
      // book(_h["angR04"], 2, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
      // Note: Event weight does not need to be included anymore.

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets fjR02 = apply<FastJets>(event, "JetsAK2").jetsByPt(Cuts::pT > 5*GeV);
      Jets fjR04 = apply<FastJets>(event, "JetsAK4").jetsByPt(Cuts::pT > 5*GeV);

      //skip the event if there are no useful jets inside
      if (fjR02.empty() && fjR04.empty()) vetoEvent;
        
      fill_jet_histograms(fjR02, "ang_R02_B10", 0.2, 1.0, 5, 40);
      fill_jet_histograms(fjR02, "ang_R02_B15", 0.2, 1.5, 5, 40);
      fill_jet_histograms(fjR02, "ang_R02_B20", 0.2, 2.0, 5, 40);
      fill_jet_histograms(fjR02, "ang_R02_B30", 0.2, 3.0, 5, 40);

      fill_jet_histograms(fjR04, "ang_R04_B10", 0.4, 1.0, 5, 40);
      fill_jet_histograms(fjR04, "ang_R04_B15", 0.4, 1.5, 5, 40);
      fill_jet_histograms(fjR04, "ang_R04_B20", 0.4, 2.0, 5, 40);
      fill_jet_histograms(fjR04, "ang_R04_B30", 0.4, 3.0, 5, 40);
        
    }
      
    ///Select good jets and fill histograms
    void fill_jet_histograms(const Jets& jets, std::string hname, double R, double beta, double min_pt, double max_pt) {

      for(const Jet& jet : jets)
      {
        if ((jet.abseta()<_jetetamax-R) && (jet.perp()>min_pt) && (jet.perp()<max_pt))
        {
          _h[hname]->fill(angularity(jet, R, beta));
        }
      }
    }
        
    ///Compute angularity
    double angularity(const Jet& jet, double R, double beta) {
      double lambda = 0;
      for (auto p : jet.constituents()) {
        double lambda_i = p.perp() * pow(deltaR(p, jet)/R, beta);
        lambda += lambda_i/jet.perp();
      }
      return lambda;
    }
    
    double deltaR(PseudoJet j1, PseudoJet j2) {
      double deta = j1.eta() - j2.eta();
      double dphi = j1.delta_phi_to(j2);
      return sqrt(deta*deta + dphi*dphi);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize to unity
      normalize(_h["ang_R02_B10"]);
      normalize(_h["ang_R02_B15"]);
      normalize(_h["ang_R02_B20"]);
      normalize(_h["ang_R02_B30"]);

      normalize(_h["ang_R04_B10"]);
      normalize(_h["ang_R04_B15"]);
      normalize(_h["ang_R04_B20"]);
      normalize(_h["ang_R04_B30"]);
        
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}

  private:
    float _jetetamax = 0.9;
  };

  DECLARE_RIVET_PLUGIN(ALICE_2020_angularity);

}
