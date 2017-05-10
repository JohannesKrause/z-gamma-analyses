// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {

   /*Cuts:  
    * 2 leptons, opposite charge, pT> 15 GeV, |eta| < 2.5 
    * photonen: pT>5GeV , |eta| < 2.5
    *      -> leading photon selected:
    *  (removed)  not from a decay 
    *           dR in (0.05, 5) to the next lepton
    *           isolated from everything else
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
      //all photons
      Cut c_photons = Cuts::pT >= 5.0*GeV && Cuts::abseta < 2.5;
      LeadingParticlesFinalState photons(c_photons);
      photons.addParticleId(PID::PHOTON);
      declare(photons, "PHOTFS");
      
      //all leptons
      Cut c_leptons   = Cuts::pT > 15*GeV && Cuts::abseta < 2.5;
      IdentifiedFinalState leptonfs(c_leptons);
      if (_mode==1)  {
		   leptonfs.acceptIdPair(PID::ELECTRON);
	   }  
      else leptonfs.acceptIdPair(PID::MUON);
      declare(leptonfs, "LEPFS");

      //all particles
      FinalState fs;
      declare(fs, "FS");  

      // book histograms      
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

      // lepton selection      
      const Particles leptonen = apply<IdentifiedFinalState>(event, "LEPFS").particlesByPt();
      if (leptonen.size() < 2) vetoEvent;
      if (leptonen[0].charge()*leptonen[1].charge() > 0) vetoEvent;
      const double mZ = ( leptonen[0].momentum() + leptonen[1].momentum()).mass();
      if (!inRange(mZ, 30*GeV, 87*GeV)) vetoEvent;

      const Particles photons = apply<LeadingParticlesFinalState>(event, "PHOTFS").particles();
      if (photons.empty()) vetoEvent;
      const Particle & photon = photons[0];
      Particles fs = apply<FinalState>(event, "FS").particles();
      // check if the leading photon passes all cuts
     // if (photon.fromDecay() ) vetoEvent;
      const double dR = std::min(deltaR(photon,leptonen[0]), deltaR(photon,leptonen[1]) );
      if (!inRange(dR, 0.05, 5.0)) vetoEvent;
      // check photon isolation
      double coneEnergy(0.0);
      foreach(const Particle& p, fs) {
            if ( deltaR(photon, p) < 0.4 )  coneEnergy += p.E();
      }
      if ( (coneEnergy - photon.E() - leptonen[0].E() - leptonen[1].E() ) >= 0.5*photon.E() )  vetoEvent;

      const double dphi = std::min(deltaPhi(photon,leptonen[0]), deltaPhi(photon,leptonen[1]) );
      const double deta = std::min(deltaEta(photon,leptonen[0]), deltaEta(photon,leptonen[1]) );
      const double qT = (leptonen[0].mom() + leptonen[1].mom() + photon.mom()).pT();
      const double mZgamma = (leptonen[0].mom() + leptonen[1].mom() + photon.mom()).mass() * GeV; 
			const Particle next_lepton = (deltaR(leptonen[0], photon) > deltaR(leptonen[1], photon) ? leptonen[0] : leptonen[1]);
   	  // Fill the analysis histograms
			_hist_pho_et->fill(photon.pT()/GeV, event.weight());
			_hist_pho_dr->fill(dR, event.weight());
			_hist_pho_dphi->fill(dphi, event.weight());
			_hist_pho_deta->fill(deta, event.weight());
			_hist_mZgamma->fill(mZgamma, event.weight());
			_hist_m_gl_next->fill((next_lepton.mom() + photon.mom()).mass()*GeV, event.weight());
			_hist_qt->fill(qT*GeV, event.weight());

      (dR <= 0.5 ? _hist_pho_et_close : _hist_pho_et_wide)->fill(photon.pT()/GeV, event.weight());


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
