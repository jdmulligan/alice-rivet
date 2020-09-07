// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2020_I1755387 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2020_I1755387);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs(Cuts::abseta < 0.9);
      
      declare(FastJets(fs, FastJets::ANTIKT, 0.1), "JetsAK1");
      declare(FastJets(fs, FastJets::ANTIKT, 0.2), "JetsAK2");
      declare(FastJets(fs, FastJets::ANTIKT, 0.3), "JetsAK3");
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, FastJets::ANTIKT, 0.5), "JetsAK5");
      declare(FastJets(fs, FastJets::ANTIKT, 0.6), "JetsAK6");

      // Book histograms
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
      // Inclusive single jet spectra
      book(_h["jetsR01"], 1, 1, 1);
      book(_h["jetsR02"], 2, 1, 1);
      book(_h["jetsR03"], 3, 1, 1);
      book(_h["jetsR04"], 4, 1, 1);
      book(_h["jetsR05"], 5, 1, 1);
      book(_h["jetsR06"], 6, 1, 1);

      // Jet cross section ratios
      book(_h_ratio["jetsRatioR02R03"], 18, 1, 1);
      book(_h_ratio["jetsRatioR02R04"], 19, 1, 1);
      book(_h_ratio["jetsRatioR02R05"], 20, 1, 1);
      book(_h_ratio["jetsRatioR02R06"], 21, 1, 1);
      
      book(_h_ratio["jetsRatioR01R02"], 13, 1, 1);
      book(_h_ratio["jetsRatioR01R03"], 14, 1, 1);
      book(_h_ratio["jetsRatioR01R04"], 15, 1, 1);
      book(_h_ratio["jetsRatioR01R05"], 16, 1, 1);
      book(_h_ratio["jetsRatioR01R06"], 17, 1, 1);
      
      // Temporary histograms for ratio
      book(_h["jetsR01_ratio_binning"], "jetsR01_ratio_binning", {20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0});
      book(_h["jetsR02_ratio_binning"], "jetsR02_ratio_binning", {20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      // Note: Event weight does not need to be included anymore.
         
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets fjR01 = apply<FastJets>(event, "JetsAK1").jetsByPt(Cuts::pT > 20*GeV);
      Jets fjR02 = apply<FastJets>(event, "JetsAK2").jetsByPt(Cuts::pT > 20*GeV);
      Jets fjR03 = apply<FastJets>(event, "JetsAK3").jetsByPt(Cuts::pT > 20*GeV);
      Jets fjR04 = apply<FastJets>(event, "JetsAK4").jetsByPt(Cuts::pT > 20*GeV);
      Jets fjR05 = apply<FastJets>(event, "JetsAK5").jetsByPt(Cuts::pT > 20*GeV);
      Jets fjR06 = apply<FastJets>(event, "JetsAK6").jetsByPt(Cuts::pT > 20*GeV);
      
      //skip the event if there are no useful jets inside
      if (fjR01.empty() && fjR02.empty() && fjR03.empty() && fjR04.empty() && fjR05.empty() && fjR06.empty()) vetoEvent;
      
      fill_jet_histograms(fjR01, "jetsR01", _jet_r01, 20, 140);
      fill_jet_histograms(fjR02, "jetsR02", _jet_r02, 20, 140);
      fill_jet_histograms(fjR03, "jetsR03", _jet_r03, 20, 140);
      fill_jet_histograms(fjR04, "jetsR04", _jet_r04, 20, 140);
      fill_jet_histograms(fjR05, "jetsR05", _jet_r05, 20, 140);
      fill_jet_histograms(fjR06, "jetsR06", _jet_r06, 20, 100);
      
      // Fill also temporary histograms for ratio
      fill_jet_histograms(fjR01, "jetsR01_ratio_binning", _jet_r01, 20, 100);
      fill_jet_histograms(fjR02, "jetsR02_ratio_binning", _jet_r02, 20, 100);

      
    }
      
    ///Select good jets and fill histograms
    void fill_jet_histograms(const Jets& jets, std::string hname, double R, double min_pt, double max_pt) {

      for(const Jet& jet : jets)
      {
        if ((jet.abseta()<_jetetamax-R) && (jet.perp()>min_pt) && (jet.perp()<max_pt))
        {
          _h[hname]->fill(jet.perp()/GeV);
        }
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      cout << "CrossSect = " << crossSection() << " NbOfEvents = " << sumOfWeights() << endl;;
  
      double normFactor = crossSection() / (millibarn * sumW() * 2 * pi);
      
      //to account for the fact that we reconstruct jets in eta=2*(0.7-R)
      double normFactor01 = normFactor/(2 * (_jetetamax-_jet_r01));
      double normFactor02 = normFactor/(2 * (_jetetamax-_jet_r02));
      double normFactor03 = normFactor/(2 * (_jetetamax-_jet_r03));
      double normFactor04 = normFactor/(2 * (_jetetamax-_jet_r04));
      double normFactor05 = normFactor/(2 * (_jetetamax-_jet_r05));
      double normFactor06 = normFactor/(2 * (_jetetamax-_jet_r06));
      
      scale(_h["jetsR01"], normFactor01 );
      scale(_h["jetsR02"], normFactor02 );
      scale(_h["jetsR03"], normFactor03 );
      scale(_h["jetsR04"], normFactor04 );
      scale(_h["jetsR05"], normFactor05 );
      scale(_h["jetsR06"], normFactor06 );
      
      scale(_h["jetsR01_ratio_binning"], normFactor01 );
      scale(_h["jetsR02_ratio_binning"], normFactor02 );

      //cross section ratios
      divide(_h["jetsR02"], _h["jetsR03"], _h_ratio["jetsRatioR02R03"]);
      divide(_h["jetsR02"], _h["jetsR04"], _h_ratio["jetsRatioR02R04"]);
      divide(_h["jetsR02"], _h["jetsR05"], _h_ratio["jetsRatioR02R05"]);
      divide(_h["jetsR02_ratio_binning"], _h["jetsR06"], _h_ratio["jetsRatioR02R06"]);

      divide(_h["jetsR01"], _h["jetsR02"], _h_ratio["jetsRatioR01R02"]);
      divide(_h["jetsR01"], _h["jetsR03"], _h_ratio["jetsRatioR01R03"]);
      divide(_h["jetsR01"], _h["jetsR04"], _h_ratio["jetsRatioR01R04"]);
      divide(_h["jetsR01"], _h["jetsR05"], _h_ratio["jetsRatioR01R05"]);
      divide(_h["jetsR01_ratio_binning"], _h["jetsR06"], _h_ratio["jetsRatioR01R06"]);
      
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _h_ratio;
    ///@}

    
  private:
    float _jetetamax = 0.7; //EMCal acceptance
    float _jet_r01 = 0.1;
    float _jet_r02 = 0.2;
    float _jet_r03 = 0.3;
    float _jet_r04 = 0.4;
    float _jet_r05 = 0.5;
    float _jet_r06 = 0.6;
  };

  DECLARE_RIVET_PLUGIN(ALICE_2020_I1755387);

}
