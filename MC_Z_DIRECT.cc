// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {
   /* Cuts:
    * 	1 bare Z-Boson with mass >40GeV
    *   2 leptons: |eta| < 2.47 && pT> 25GeV
    *   photons: |eta| < 2.37 && pT> _ptcut
    *   photon isolation: coneEnergy(dr=0.4) / leadingPhoton.E() >= 0.5
    *   deltaR(leadingPhoton, lepton) < 0.4 
    * 
    * */

  /// @brief Add a short analysis description here
  class MC_Z_DIRECT : public Analysis {
  public:

    /// Constructor
    MC_Z_DIRECT(string name="MC_Z_DIRECT")
      : Analysis(name)
    {
      // the myon mode is used by default
      _mode = 3;
      _ptcut = 15;
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      declare(fs, "FS");

      Cut cuts = Cuts::abseta < 2.47 && Cuts::pT > 25*GeV;

      // Z finder
      ZFinder zf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 40.0*GeV, 1000.0*GeV, 0.0);
      declare(zf, "ZF");

      // leading photon
      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 2.37 && Cuts::pT > _ptcut*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZF"));
      jet_fs.addVetoOnThisFinalState(getProjection<LeadingParticlesFinalState>("LeadingPhoton"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "Jets");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "isolatedFS");


      // Book histograms
      _hist_EgammaT_incl   = bookHisto1D(_name("et_incl",_ptcut, _mode), logspace(80,_ptcut, 1500) ); // dSigma / dE^gamma_T for Njet >= 0
      _hist_EgammaT_excl   = bookHisto1D(_name("et_excl",_ptcut, _mode), logspace(80,_ptcut, 1500) ); // dSigma / dE^gamma_T for Njet = 0
      _hist_Njet_incl       = bookHisto1D(_name("njet_incl",_ptcut, _mode), 5, 0,5 ); //dSigma / dNJET with NJET>= Njet
      _hist_Njet_excl       = bookHisto1D(_name("njet_excl",_ptcut, _mode), 5, 0,5 ); //dSigma / dNJET with NJET= Njet
      _hist_mZgamma        = bookHisto1D(_name("mZgamma",_ptcut, _mode), logspace(80,_ptcut, 1500) ); //dSigma / dmZgamma
      _hist_dR_gl_min         = bookHisto1D(_name("dR_gamma_lepton_min",_ptcut, _mode), 80, 0.4, 5 ); //dSigma / dR(gamma, lepton1)
      _hist_dR_gl_max         = bookHisto1D(_name("dR_gamma_lepton_max",_ptcut, _mode), 80, 0.4, 5 ); //dSigma / dR(gamma, lepton2)
      _hist_m_gl_min          = bookHisto1D(_name("m_gamma_lepton_min",_ptcut, _mode), logspace(80, 40, 1500) ); //dSigma / dm(gamma, lepton1)
      _hist_m_gl_max          = bookHisto1D(_name("m_gamma_lepton_max",_ptcut, _mode), logspace(80, 40, 1500) ); //dSigma / dm(gamma, lepton2)

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // retrieve leading photon
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1)  vetoEvent;
      const Particle& leadingPhoton = photons[0];

      // check photon isolation
      double coneEnergy(0.0);
      Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
      foreach(const Particle& p, fs) {
        if ( deltaR(leadingPhoton, p) < 0.4 )  coneEnergy += p.E();
      }
      if ( coneEnergy / leadingPhoton.E() >= 0.5 )  vetoEvent;
      
      
      // retrieve Z boson candidate
      const ZFinder& zf = apply<ZFinder>(event, "ZF");
      if ( zf.bosons().size() != 1 )  vetoEvent; // only one Z boson candidate
      const Particle& Zboson  = zf.boson();
      if ( !(Zboson.mass() > 40.0*GeV) )  vetoEvent;
      
      // check charge of constituent leptons
      const ParticleVector& leptons = zf.constituents();
      if (leptons.size() != 2 || leptons[0].charge() * leptons[1].charge() > 0.)  vetoEvent;

      // check photon-lepton overlap
      foreach(const Particle& p, leptons) {
        if ( deltaR(leadingPhoton, p) < 0.4)  vetoEvent;
      }

      // count jets
      const FastJets& jetfs = apply<FastJets>(event, "Jets");
      Jets jets = jetfs.jets(cmpMomByEt);
      int goodJets = 0;
      foreach (const Jet& j, jets) {
        if ( !(j.Et() > 30.0*GeV) )  break;
        if ( (j.abseta() < 4.4) && \
             (deltaR(leadingPhoton, j) > 0.3) &&    \
             (deltaR(leptons[0],    j) > 0.3) &&            \
             (deltaR(leptons[1],    j) > 0.3) )  ++goodJets;
      }

      double Njets = double(goodJets) + 0.5;
      double photonEt = leadingPhoton.Et()*GeV;
      double mZgamma = (Zboson.momentum() + leadingPhoton.momentum()).mass() * GeV; 
      double dr_g_lep1= deltaR(leptons[0].momentum(), leadingPhoton.momentum() );
      double dr_g_lep2= deltaR(leptons[1].momentum(), leadingPhoton.momentum() );
      double m_g_lep1 = (leptons[0].momentum() + leadingPhoton.momentum()).mass(); 
      double m_g_lep2 = (leptons[1].momentum() + leadingPhoton.momentum()).mass(); 

      //ET histogramms
      _hist_EgammaT_incl->fill(photonEt, weight);
      if (!goodJets) _hist_EgammaT_excl->fill(photonEt, weight);
      
      // Multiplicities
      _hist_Njet_excl->fill(Njets, weight);
      for (size_t i = 0; i < 5; ++i) {
        if (Njets > i) {
          _hist_Njet_incl->fill(i+0.5, weight);
        }
      }
      
      // masses
      _hist_mZgamma->fill(mZgamma, weight);
      _hist_m_gl_min->fill(min(m_g_lep1, m_g_lep2), weight);
      _hist_m_gl_max->fill(max(m_g_lep1, m_g_lep2), weight);
      
      // deltaR
      _hist_dR_gl_min->fill(min(dr_g_lep1, dr_g_lep2), weight);
      _hist_dR_gl_max->fill(max(dr_g_lep1, dr_g_lep2), weight);


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_EgammaT_incl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_EgammaT_excl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_Njet_incl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_Njet_excl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_mZgamma, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_dR_gl_min, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_dR_gl_max, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_m_gl_min, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_m_gl_max, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }

    //@}

  protected:
     size_t _mode;
     double _ptcut;


  private:
    string _name(const string &basename, const double &pt, const int &mode){
       return basename + "_" +  to_string(int(pt)) + "_" + to_string(mode);
    } 

    /// @name Histograms
    //@{
    Histo1DPtr _hist_EgammaT_incl;
	Histo1DPtr _hist_EgammaT_excl;
	Histo1DPtr _hist_Njet_incl;
	Histo1DPtr _hist_Njet_excl;
	Histo1DPtr _hist_mZgamma;
	Histo1DPtr _hist_dR_gl_min;
	Histo1DPtr _hist_dR_gl_max;
	Histo1DPtr _hist_m_gl_min;
	Histo1DPtr _hist_m_gl_max;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Z_DIRECT);


}
