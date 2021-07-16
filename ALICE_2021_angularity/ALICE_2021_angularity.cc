// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/contrib/SoftDrop.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class ALICE_2021_angularity : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2021_angularity);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
        
      // Measurement parameters to loop over
      alpha_list = {1.0, 1.5, 2.0, 3.0};
      R_list = {0.2, 0.4};
      pt_bins = {20, 40, 60, 80, 100};

      // Initialise and register projections
      // The basic final-state projection: all final-state particles within the given eta acceptance
      const ChargedFinalState fs(Cuts::abseta < 0.9);
      //const ALICE::Primary primary(Cuts::abseta < 1.0 && Cuts::abscharge > 0);

      declare(FastJets(fs, FastJets::ANTIKT, 0.2), "JetsAK2");
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");

      // Book histograms
      for (int pt_bin=0; pt_bin<n_pt_bins; pt_bin++) {
        for (auto R : R_list) {
          for (auto alpha : alpha_list) {
            string hname = hname_string(pt_bin, R, alpha);
            string hname_SD = hname + "_SD";
            book(_h[hname], hname, 20, 0.0, 1.0);
            book(_h[hname_SD], hname_SD, 20, 0.0, 1.0);
              
            // TODO: take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
            // book(_h["angR02"], 1, 1, 1);
            // book(_h["angR04"], 2, 1, 1);
              
            // book(_h_Table1, 1,1,1);
            // book(_h_Table2, 2,1,1);
            // book(_h_Table3, 3,1,1);
          }
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
      // Note: Event weight does not need to be included anymore.

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets fjR02 = apply<FastJets>(event, "JetsAK2").jetsByPt(Cuts::pT > 20*GeV);
      Jets fjR04 = apply<FastJets>(event, "JetsAK4").jetsByPt(Cuts::pT > 20*GeV);

      //skip the event if there are no useful jets inside
      if (fjR02.empty() && fjR04.empty()) vetoEvent;
        
      for (int pt_bin=0; pt_bin<n_pt_bins; pt_bin++) {
        for (auto alpha : alpha_list) {
          fill_jet_histograms(fjR02, 0.2, alpha, pt_bin);
          fill_jet_histograms(fjR04, 0.4, alpha, pt_bin);
        }
      }
        
    }
      
    ///Select good jets and fill histograms
    void fill_jet_histograms(const Jets& jets, float R, float alpha, int pt_bin) {
        
      string hname = hname_string(pt_bin, R, alpha);
      string hname_SD = hname + "_SD";
        
      float min_pt = pt_bins.at(pt_bin);
      float max_pt = pt_bins.at(pt_bin+1);

      for(const Jet& jet : jets)
      {
        if ((jet.abseta()<_jetetamax-R) && (jet.perp()>min_pt) && (jet.perp()<max_pt))
        {
          // Ungroomed angularity
          _h[hname]->fill(angularity(jet, R, alpha));
            
          // Groomed angularity
          float zcut = 0.2;
          float beta = 0.;
          ClusterSequence cs_CA(jet.constituents(), JetDefinition(fastjet::cambridge_algorithm, R));
          PseudoJets jets_CA = sorted_by_pt(cs_CA.inclusive_jets());
          if (jets_CA.size() > 0) {
            fastjet::contrib::SoftDrop sd(beta, zcut);
            PseudoJet jet_SD = sd(jets_CA[0]);
            _h[hname_SD]->fill(angularity(jet_SD, R, alpha));
          }
        }
      }
    }
        
    ///Compute angularity
    float angularity(const Jet& jet, float R, float alpha) {
      float lambda = 0;
      for (auto p : jet.constituents()) {
        float lambda_i = p.perp() * pow(deltaR(p, jet)/R, alpha);
        lambda += lambda_i/jet.perp();
      }
      return lambda;
    }
    
    float deltaR(PseudoJet j1, PseudoJet j2) {
      float deta = j1.eta() - j2.eta();
      float dphi = j1.delta_phi_to(j2);
      return sqrt(deta*deta + dphi*dphi);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalize ungroomed angularities to unity
      // TODO: Normalize groomed angularities to SD tagged fraction
      for (int pt_bin=0; pt_bin<n_pt_bins; pt_bin++) {
        for (auto R : R_list) {
          for (auto alpha : alpha_list) {
            string hname = hname_string(pt_bin, R, alpha);
            string hname_SD = hname + "_SD";
            normalize(_h[hname]);
            normalize(_h[hname_SD]);
          }
        }
      }

    }
      
    /// Get string label for each unique jet setting
    string hname_string(int pt_bin, float R, float alpha) {
        return "ang_R" + to_string(R) + "_alpha" + to_string(alpha) + "_pt" + to_string(pt_bin);
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}

  private:
    float _jetetamax = 0.9;
    int n_pt_bins = 4;
    vector<float> pt_bins;
    vector<float> R_list;
    vector<float> alpha_list;
  };

  DECLARE_RIVET_PLUGIN(ALICE_2021_angularity);

}
