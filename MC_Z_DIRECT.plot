# BEGIN PLOT /MC_Z_DIRECT/*
LogY=1
RatioPlotYMax=1.5
RatioPlotYMin=0.5
ErrorBars=1
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/et_excl_*
LogX=1
Title= $\pT^\gamma$, exclusive
XLabel=$E^\gamma_\text{T}$ [GeV]
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d}E^\gamma_\text{T}$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/et_incl_*
LogX=1
Title= $\pT ^\gamma$, inclusive
XLabel=$E^\gamma_\text{T}$ [GeV]
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d}E^\gamma_\text{T}$ [fb $\text{GeV}^{-1}$]
# END PLOT


# BEGIN PLOT /MC_Z_DIRECT/njet_excl_*
Title=Exclusive Jet multiplicity
XLabel=Jet multiplicity 
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} N_\text{jet}$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/njet_incl_*
Title=Inclusive Jet multiplicity
XLabel=Jet multiplicity 
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} N_\text{jet}$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/mZgamma_*
LogX=1
Title = Invariant mass of $Z\gamma$
XLabel= invariant mass($Z\gamma$)
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} m^{Z\gamma}$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/dR_gamma_lepton_min_*
Title= (lepton $\gamma$) pair with low $\Delta R$
XLabel=$ Delta R$($\gamma$, lepton)
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} \Delta R$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/dR_gamma_lepton_max_*
Title= (lepton $\gamma$) pair with high $\Delta  R$
XLabel= $Delta R$($\gamma$, lepton)
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} \Delta R$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/m_gamma_lepton_min_*
LogX=1
Title= (lepton $\gamma$ )pair with low mass
XLabel= invariant mass($\gamma$ lepton)
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} m^{\gamma\text{lepton}}$ [fb $\text{GeV}^{-1}$]
# END PLOT

# BEGIN PLOT /MC_Z_DIRECT/m_gamma_lepton_max_*
LogX=1
Title= (lepton $\gamma$) pair with high mass
XLabel= invariant mass($\gamma$ lepton)
YLabel=$\text{d}\sigma_{Z\gamma} / \text{d} m^{\gamma\text{lepton}}$ [fb $\text{GeV}^{-1}$]
# END PLOT
