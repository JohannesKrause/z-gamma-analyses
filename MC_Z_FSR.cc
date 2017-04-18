// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {

   /*Cuts:  
    * 2 leptons, opposite charge, pT> 15 GeV, |eta| < 2.5 
    * photonen: pT>5GeV , |eta| < 2.5
    *      -> highest pT-photon selected that does not come from a decay 
    *         and has a dR in (0.05, 5) to the next myon
    * lepton pair: invariant mass in (30, 87)
    * 
    *  */
   

  /// @brief Add a short analysis description here
  class MC_Z_FSR : public Analysis {
  public:

    /// Constructor
    MC_Z_FSR(string name="MC_Z_FSR")
      : Analysis(name)
    {
      // mode 1: electron, 2: muon
      _mode = 2;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      Cut c_photons = Cuts::pT >= 5.0*GeV && Cuts::abseta < 2.5;
      IdentifiedFinalState photons(c_photons);
      photons.acceptId(PID::PHOTON);
      declare(photons, "PHOTFS");

      Cut c_leptons   = Cuts::pT > 15*GeV && Cuts::abseta < 2.5;
      IdentifiedFinalState leptonfs(c_leptons);
      if (_mode==1)  {
		   leptonfs.acceptIdPair(PID::ELECTRON);
	   }  
      else leptonfs.acceptIdPair(PID::MUON);
      declare(leptonfs, "LEPFS");


      _hist_pho_et           = bookHisto1D("photon_et_all", 100, 5, 250);  // photon transverse energy
      _hist_pho_et_wide      = bookHisto1D("photon_et_wide", 100, 5, 250);  // photon transverse energy (0.5 < dr < 3.0)
      _hist_pho_et_close     = bookHisto1D("photon_et_close", 100, 5, 250);  // photon transverse energy (0.05 < dr < 0.5)
      _hist_pho_dr           = bookHisto1D("dr_all",80,0.05,5);  // delta_R to the closest lepton
      _hist_pho_dphi         = bookHisto1D("dphi_all", 80, 0.05, 2); //delta_Phi to the closest lepton
      _hist_pho_deta         = bookHisto1D("deta_all", 80, 0.05, 3); //delta_Phi to the closest lepton
      _hist_mZgamma          = bookHisto1D("mZgamma", 100,30, 200 ); //dSigma / dmZgamma
      _hist_m_gl_next        = bookHisto1D("m_gamma_lepton_next", 80, 10, 250 ); //dSigma / dm(gamma,  closest lepton)
      _hist_qt               = bookHisto1D("qT", 80, 1, 200 ); //dSigma / dm(gamma,  closest lepton)


    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    
      const Particles leptonen = apply<IdentifiedFinalState>(event, "LEPFS").particlesByPt();

      if (leptonen.size() < 2) vetoEvent;
      if (leptonen[0].charge()*leptonen[1].charge() > 0) vetoEvent;
      const double mZ = ( leptonen[0].momentum() + leptonen[1].momentum()).mass();
      if (!inRange(mZ, 30*GeV, 87*GeV)) vetoEvent;

      const Particles photons = apply<IdentifiedFinalState>(event, "PHOTFS").particlesByPt();
      // We want the photon with the highest pT that does not come from a decay
      foreach(const Particle& p, photons) {
        if (p.fromDecay() || !p.isStable()) continue;
        const double dR = std::min(deltaR(p,leptonen[0]), deltaR(p,leptonen[1]) );
        if (!inRange(dR, 0.05, 5.0)) continue;

        const double dphi = std::min(deltaPhi(p,leptonen[0]), deltaPhi(p,leptonen[1]) );
        const double deta = std::min(deltaEta(p,leptonen[0]), deltaEta(p,leptonen[1]) );
        const double qT = (leptonen[0].mom() + leptonen[1].mom() + p.mom()).pT();
        const double mZgamma = (leptonen[0].mom() + leptonen[1].mom() + p.mom()).mass() * GeV; 

				const Particle next_lepton = (deltaR(leptonen[0], p) > deltaR(leptonen[1], p) ? leptonen[0] : leptonen[1]);

   			 // Fill the analysis histograms
				_hist_pho_et->fill(p.pT()/GeV, event.weight());
				_hist_pho_dr->fill(dR, event.weight());
				_hist_pho_dphi->fill(dphi, event.weight());
				_hist_pho_deta->fill(deta, event.weight());
				_hist_mZgamma->fill(mZgamma, event.weight());
				_hist_m_gl_next->fill((next_lepton.mom() + p.mom()).mass()*GeV, event.weight());
				_hist_qt->fill(qT*GeV, event.weight());

        (dR <= 0.5 ? _hist_pho_et_close : _hist_pho_et_wide)->fill(p.pT()/GeV, event.weight());

        break; // Exit the loop since we found the highest pT lepton already
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_hist_pho_et,       crossSection() / sumOfWeights());
      scale(_hist_pho_et_wide,  crossSection() / sumOfWeights());
      scale(_hist_pho_et_close, crossSection() / sumOfWeights());
      scale(_hist_pho_dr,       crossSection() / sumOfWeights());
      scale(_hist_m_gl_next, crossSection() / sumOfWeights());
      scale(_hist_mZgamma, crossSection() / sumOfWeights());
      scale(_hist_pho_deta, crossSection() / sumOfWeights());
      scale(_hist_pho_dphi, crossSection() / sumOfWeights());
      scale(_hist_qt, crossSection() / sumOfWeights());      

    }

    //@}

  protected:
      size_t _mode;


  private:
   

    /// @name Histograms
    //@{
    Histo1DPtr _hist_pho_et;
    Histo1DPtr _hist_pho_et_wide, _hist_pho_et_close;
    Histo1DPtr _hist_pho_dr, _hist_pho_dphi, _hist_pho_deta;
    Histo1DPtr _hist_mZgamma , _hist_m_gl_next;
    Histo1DPtr _hist_qt;
    //@}


  };


  class MC_Z_FSR_E : public MC_Z_FSR {
  public:
    MC_Z_FSR_E()
      :MC_Z_FSR("MC_Z_FSR_E")
    {
      _mode = 1;
    }
  };
  
  
    class MC_Z_FSR_MU : public MC_Z_FSR {
  public:
    MC_Z_FSR_MU()
      :MC_Z_FSR("MC_Z_FSR_MU")
    {
      _mode = 2;
    }
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Z_FSR);
  DECLARE_RIVET_PLUGIN(MC_Z_FSR_E);
  DECLARE_RIVET_PLUGIN(MC_Z_FSR_MU);





}
