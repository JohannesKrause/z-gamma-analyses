// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {

   /*Cuts:  
    * 2 myons, opposite charge, pT> 10 GeV, |eta| < 2.5 
    * hardest myon: pT> 20 GeV 
    * photonen: pT>5GeV , |eta| < 1.4
    *      -> highest pT-photon selected that does not come from a decay 
    *         and has a dR in (0.05, 3) to the next myon
    * myon pair: invariant mass in (30, 87)
    * 
    *  */
   

  /// @brief Add a short analysis description here
  class MC_Z_FSR : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Z_FSR);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      Cut c_photons = Cuts::pT >= 5.0*GeV && (Cuts::etaIn(-1.4, 1.4));
      IdentifiedFinalState photons(c_photons);
      photons.acceptId(PID::PHOTON);
      declare(photons, "PHOTFS");

      Cut c_muons   = Cuts::pT > 10*GeV && Cuts::abseta < 2.4;
      IdentifiedFinalState muons(c_muons);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "MUFS");


      _hist_pho_et           = bookHisto1D("photon_et_all", 95, 5, 250);  // photon transverse energy
      _hist_pho_et_wide      = bookHisto1D("photon_et_wide", 95, 5, 250);  // photon transverse energy (0.5 < dr < 3.0)
      _hist_pho_et_close     = bookHisto1D("photon_et_close", 95, 5, 250);  // photon transverse energy (0.05 < dr < 0.5)
      _hist_pho_et_lqt       = bookHisto1D("photon_et_lqt", 95, 5, 250);  // photon transverse energy (q_T < 10)
      _hist_pho_et_hqt       = bookHisto1D("photon_et_hqt", 95, 5, 250);  // photon transverse energy (q_T > 50)
      _hist_pho_dr           = bookHisto1D("dr_all",60,0,5);  // delta_R
      _hist_pho_dr_lqt       = bookHisto1D("dr_lqt",60,0,5);  // delta_R (q_T < 10)
      _hist_pho_dr_hqt       = bookHisto1D("dr_hqt",60,0,5);  // delta_R  (q_T > 50)
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    
      const Particles muons = apply<IdentifiedFinalState>(event, "MUFS").particlesByPt();

      if (muons.size() < 2) vetoEvent;
      if (muons[0].pT()/GeV < 20) vetoEvent;
      if (muons[0].charge()*muons[1].charge() > 0) vetoEvent;
      const double mZ = (muons[0].momentum() + muons[1].momentum()).mass();
      if (!inRange(mZ, 30*GeV, 87*GeV)) vetoEvent;

      const Particles photons = apply<IdentifiedFinalState>(event, "PHOTFS").particlesByPt();
      // We want the photon with the highest pT that does not come from a decay
      foreach(const Particle& p, photons) {
        if (p.fromDecay() || !p.isStable()) continue;

        const double dR = std::min(deltaR(p, muons[0]), deltaR(p, muons[1]) );
        if (!inRange(dR, 0.05, 3.0)) continue;

        // Calculate the three-body (mu,mu,gamma) transverse momentum
        const double qT = (muons[0].mom() + muons[1].mom() + p.mom()).pT();

        // Fill the analysis histograms
        _hist_pho_et->fill(p.pT()/GeV, event.weight());
        _hist_pho_dr->fill(dR, event.weight());

        (dR <= 0.5 ? _hist_pho_et_close : _hist_pho_et_wide)->fill(p.pT()/GeV, event.weight());

        if (qT / GeV < 10.) {
          _hist_pho_et_lqt->fill(p.pT()/GeV, event.weight());
          _hist_pho_dr_lqt->fill(dR, event.weight());
        }

        if (qT / GeV > 50.) {
          _hist_pho_et_hqt->fill(p.pT()/GeV, event.weight());
          _hist_pho_dr_hqt->fill(dR, event.weight());
        }

        break; // Exit the loop since we found the highest pT lepton already
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_hist_pho_et,       crossSection() / sumOfWeights());
      scale(_hist_pho_et_wide,  crossSection() / sumOfWeights());
      scale(_hist_pho_et_close, crossSection() / sumOfWeights());
      scale(_hist_pho_et_lqt,   crossSection() / sumOfWeights());
      scale(_hist_pho_et_hqt,   crossSection() / sumOfWeights());
      scale(_hist_pho_dr,       crossSection() / sumOfWeights());
      scale(_hist_pho_dr_lqt,   crossSection() / sumOfWeights());
      scale(_hist_pho_dr_hqt,   crossSection() / sumOfWeights());

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _hist_pho_et;
    Histo1DPtr _hist_pho_et_wide, _hist_pho_et_close;
    Histo1DPtr _hist_pho_et_lqt,  _hist_pho_et_hqt;
    Histo1DPtr _hist_pho_dr;
    Histo1DPtr _hist_pho_dr_lqt, _hist_pho_dr_hqt;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Z_FSR);


}
