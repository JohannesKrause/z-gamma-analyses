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
    *   2 leptons: |eta| < 2.5 && pT> 25GeV
    *   photons: |eta| < 2.5 && pT> _ptcut
    *   photon isolation: coneEnergy(dr=0.4) / leadingPhoton.E() >= 0.5
    *   deltaR(leadingPhoton, lepton) < 0.4 
    * 
    *   jets:  |eta| < 4.4
    *           deltaR (jet, {gamma, lepton1, lepton2} ) < 0.3
    * 
    *  lepton1: high pT, lepton2: low pT
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

      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;

      // Z finder
      ZFinder zf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 40.0*GeV, 1000.0*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zf, "ZF");


      // leading photon
      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 2.5 && Cuts::pT > _ptcut*GeV));
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
      _hist_EgammaT_incl      = bookHisto1D(_name("et_incl",_ptcut, _mode), logspace(80,_ptcut, 1500) ); // dSigma / dE^gamma_T for Njet >= 0
      _hist_EgammaT_excl      = bookHisto1D(_name("et_excl",_ptcut, _mode), logspace(80,_ptcut, 1500) ); // dSigma / dE^gamma_T for Njet = 0
      _hist_eta_gamma_incl    = bookHisto1D(_name("eta_incl",_ptcut, _mode), 80, -2.5, 2.5 ); // dSigma / dE^gamma_T for Njet >= 0
      _hist_eta_gamma_excl    = bookHisto1D(_name("eta_excl",_ptcut, _mode), 80, -2.5, 2.5 ); // dSigma / dE^gamma_T for Njet = 0
      _hist_Njet_incl         = bookHisto1D(_name("njet_incl",_ptcut, _mode), 5, -0.5, 4.5 ); //dSigma / dNJET with NJET>= Njet
      _hist_Njet_excl         = bookHisto1D(_name("njet_excl",_ptcut, _mode), 5, -0.5, 4.5 ); //dSigma / dNJET with NJET= Njet
      _hist_mZgamma           = bookHisto1D(_name("mZgamma",_ptcut, _mode), logspace(80, 40, 5000) ); //dSigma / dmZgamma
      _hist_dR_gl_min         = bookHisto1D(_name("dR_gamma_lepton_min",_ptcut, _mode), 80, 0.4, 5 ); //dSigma / dR(gamma, lepton1)
      _hist_dR_gl_max         = bookHisto1D(_name("dR_gamma_lepton_max",_ptcut, _mode), 80, 0.4, 5 ); //dSigma / dR(gamma, lepton2)
      _hist_m_gl_min          = bookHisto1D(_name("m_gamma_lepton_min",_ptcut, _mode), logspace(80, 20, 1500) ); //dSigma / dm(gamma, lepton1)
      _hist_m_gl_max          = bookHisto1D(_name("m_gamma_lepton_max",_ptcut, _mode), logspace(80, 20, 1500) ); //dSigma / dm(gamma, lepton2)
      _hist_pt_jet1           = bookHisto1D(_name("pt_leadingjet",_ptcut, _mode), logspace(80, 30, 1500) ); 
      _hist_pt_jet2           = bookHisto1D(_name("pt_subleadingjet",_ptcut, _mode), logspace(80, 30, 1500) ); 
      _hist_pt_lep1           = bookHisto1D(_name("pt_leadinglepton",_ptcut, _mode), logspace(80, 25, 1500) ); 
      _hist_pt_lep2           = bookHisto1D(_name("pt_subleadinglepton",_ptcut, _mode), logspace(80, 25, 1500) ); 
      _hist_ht	      		  = bookHisto1D(_name("HT",_ptcut, _mode), logspace(80,30, 5000) ); // HT
      
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
      //leptons are orderd by pt, starting with lowest pt
      const ParticleVector& leptons = zf.constituents(cmpMomByAscPt);
      if (leptons.size() != 2 || leptons[0].charge() * leptons[1].charge() > 0.)  vetoEvent;

      // check photon-lepton overlap
      foreach(const Particle& p, leptons) {
        if ( deltaR(leadingPhoton, p) < 0.4)  vetoEvent;
      }

      // count jets
      const FastJets& jetfs = apply<FastJets>(event, "Jets");
      Jets jets = jetfs.jets(cmpMomByEt);
      Jets goodJets;
      int num_goodJets = 0;
      foreach (const Jet& j, jets) {
        if ( !(j.Et() > 30.0*GeV) )  break;
        if ( (j.abseta() < 4.4) && \
             (deltaR(leadingPhoton, j) > 0.3) &&    \
             (deltaR(leptons[0],    j) > 0.3) &&            \
             (deltaR(leptons[1],    j) > 0.3) ) {
					++num_goodJets;
					goodJets.push_back(j);
		}
      }

      double Njets = double(num_goodJets);
      double photonEt = leadingPhoton.Et()*GeV;
      double photonEta = leadingPhoton.eta();
      double mZgamma = (Zboson.momentum() + leadingPhoton.momentum()).mass() * GeV; 
      double dr_g_lep2= deltaR(leptons[0].momentum(), leadingPhoton.momentum() );
      double dr_g_lep1= deltaR(leptons[1].momentum(), leadingPhoton.momentum() );
      double m_g_lep1 = (leptons[1].momentum() + leadingPhoton.momentum()).mass(); 
      double m_g_lep2 = (leptons[0].momentum() + leadingPhoton.momentum()).mass(); 

      //gamma histogramms
      _hist_EgammaT_incl->fill(photonEt, weight);
      _hist_eta_gamma_incl->fill(photonEta, weight);
      if (!num_goodJets) {
		  _hist_EgammaT_excl->fill(photonEt, weight);
		  _hist_eta_gamma_excl->fill(photonEta, weight);
	  }	
      // Multiplicities
      _hist_Njet_excl->fill(Njets, weight);
      for (size_t i = 0; i < 5; ++i) {
        if (Njets >= i) {
          _hist_Njet_incl->fill(i, weight);
        }
      }
      // HT
      double HT=0;
      foreach (const Jet& gj, goodJets)  HT+=gj.pt();
	  _hist_ht->fill(HT, weight);
            
      // masses
      _hist_mZgamma->fill(mZgamma, weight);
      _hist_m_gl_min->fill(min(m_g_lep1, m_g_lep2), weight);
      _hist_m_gl_max->fill(max(m_g_lep1, m_g_lep2), weight);
      
      // deltaR
      _hist_dR_gl_min->fill(min(dr_g_lep1, dr_g_lep2), weight);
      _hist_dR_gl_max->fill(max(dr_g_lep1, dr_g_lep2), weight);
      
      //jets
      if(goodJets.size()>0) _hist_pt_jet1->fill(goodJets[0].pT(), weight);
      if(goodJets.size()>1) _hist_pt_jet2->fill(goodJets[1].pT(), weight);

	  //leptons
	  _hist_pt_lep1->fill(leptons[1].pT(), weight);
	  _hist_pt_lep2->fill(leptons[0].pT(), weight);
	  
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_EgammaT_incl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_EgammaT_excl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_eta_gamma_incl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_eta_gamma_excl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_Njet_incl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_Njet_excl, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_mZgamma, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_dR_gl_min, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_dR_gl_max, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_m_gl_min, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_m_gl_max, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_pt_jet1, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_pt_jet2, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_pt_lep1, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_pt_lep2, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_hist_ht, crossSection()/picobarn/sumOfWeights()); // norm to cross section
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
    Histo1DPtr _hist_EgammaT_incl, _hist_EgammaT_excl;
	Histo1DPtr _hist_eta_gamma_incl, _hist_eta_gamma_excl;
	Histo1DPtr _hist_Njet_incl;
	Histo1DPtr _hist_Njet_excl;
	Histo1DPtr _hist_mZgamma;
	Histo1DPtr _hist_dR_gl_min;
	Histo1DPtr _hist_dR_gl_max;
	Histo1DPtr _hist_m_gl_min;
	Histo1DPtr _hist_m_gl_max;
	Histo1DPtr _hist_pt_jet1, _hist_pt_lep1;
	Histo1DPtr _hist_pt_jet2, _hist_pt_lep2;
	Histo1DPtr _hist_ht;

    //@}


  };

  class MC_Z_DIRECT_60 : public MC_Z_DIRECT{
  public:
		MC_Z_DIRECT_60():
		    MC_Z_DIRECT("MC_Z_DIRECT_60")
		{ _mode = 3;
		  _ptcut = 60; 
		}
  };

  class MC_Z_DIRECT_E : public MC_Z_DIRECT{
  public:
		MC_Z_DIRECT_E():
		     MC_Z_DIRECT("MC_Z_DIRECT_E")
		{ _mode = 1;
		  _ptcut = 15; 
		}
  };
  class MC_Z_DIRECT_E_60 : public MC_Z_DIRECT{
  public:
		MC_Z_DIRECT_E_60():
		     MC_Z_DIRECT("MC_Z_DIRECT_E_60")
		{ _mode = 1;
		  _ptcut = 60; 
		}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Z_DIRECT);
  DECLARE_RIVET_PLUGIN(MC_Z_DIRECT_E);
  DECLARE_RIVET_PLUGIN(MC_Z_DIRECT_60);
  DECLARE_RIVET_PLUGIN(MC_Z_DIRECT_E_60);


}
