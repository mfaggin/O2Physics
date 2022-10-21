// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskSigmaC.cxx
/// \brief Task for Σc0,++ → Λc+(→pK-π+) π-,+ analysis
/// \note Σc0,++ candidates built in O2Physics/PWGHF/TableProducer/HFCandidateCreatorSigmaCZeroPlusPlus.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN PADOVA

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_sc;

struct TaskSigmaC{

    /// analysis histograms
    HistogramRegistry registry{
        "registry",
        {/// Σc0
        {"RecoData/hPtSigmaCZero", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hDeltaMassSigmaCZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Σc++
        {"RecoData/hPtSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hDeltaMassSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Σc0,++
        {"RecoData/hPtSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hDeltaMassSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc0
        {"RecoData/hPtLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hDeltaMassLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc++
        {"RecoData/hPtLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hDeltaMassLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc0,++
        {"RecoData/hPtLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}
    };

    /// @brief init function, to define the additional analysis histograms
    /// @param  
    void init(InitContext&) {
        /// TO DO: add histograms for MC, when required
        /// [...]
    }; /// end init

    /// @brief process function to fill the histograms needed in analysis (data)
    /// @param candidatesSigmaC are the reconstructed candidate Σc0,++
    /// @param 
    void process(const aod::HfCandSigmaC& candidatesSigmaC,
    soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate> const&, const soa::Join<aod::Tracks, aod::TracksDCA>&) {

        /// loop over the candidate Σc0,++
        for(auto& candSigmaC : candidatesSigmaC) {
            
            const int chargeSigmaC = candSigmaC.charge();   // either Σc0 or Σc++

            /// get the candidate Λc+ used to build the candidate Σc0,++
            /// and understand which mass hypotheses are possible
            const auto& candLambdaC = candSigmaC.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>>();
            const int isCandLambdaCpKpi = (candLambdaC.isSelLcpKpi() >= 1) && candSigmaC.statusSpreadLcMinvpKpiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
            const int isCandLambdaCpiKp = (candLambdaC.isSelLcpiKp() >= 1) && candSigmaC.statusSpreadLcMinvpiKpFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
            double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.), ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
            if(isCandLambdaCpKpi) {
                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                massSigmaC = InvMassScRecoLcpKpi(candSigmaC);
                massLambdaC = InvMassLcpKpi(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("RecoData/hPtSigmaCZero"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZero"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("RecoData/hPtSigmaCPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → pK-π+ (and charge conjugate)
            if(isCandLambdaCpiKp) {
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                massSigmaC = InvMassScRecoLcpiKp(candSigmaC);
                massLambdaC = InvMassLcpiKp(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("RecoData/hPtSigmaCZero"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZero"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("RecoData/hPtSigmaCPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → π+K-p (and charge conjugate)
        } /// end loop over the candidate Σc0,++
    };  /// end process



    /// @brief process function to fill the histograms needed in analysis (MC)
    /// @param candidatesSigmaC are the reconstructed candidate Σc0,++ with MC info
    /// @param particlesMC are the generated particles with flags wheter they are Σc0,++ or not
    /// @param 
    void processMC(const soa::Join<aod::HfCandSigmaC, aod::HfCandSigmaCMCRec>& candidatesSigmaC,
    soa::Join<aod::McParticles, aod::HfCandSigmaCMCGen> const& particlesMC,
    soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate, aod::HfCandProng3MCRec> const&, const soa::Join<aod::Tracks, aod::TracksDCA>&) {

        /// MC generated particles
        for(auto& particle : particlesMC) {

            /// look for the generated Σc0,++
            /// reject immediately different particles
            bool isSigmaCZeroGen = (std::abs(particle.flagMCMatchGen()) == (1 << aod::hf_cand_sc::DecayType::SigmaCZeroToPKPiPi));
            bool isSigmaCPlusPlusGen = (std::abs(particle.flagMCMatchGen()) == (1 << aod::hf_cand_sc::DecayType::SigmaCPlusPlusToPKPiPi));
            if(!isSigmaCZeroGen && !isSigmaCPlusPlusGen)
                continue;
            
            /// check for generated particles in GenLimAcc, GenAccMother, GenAcc ....
            /// [...]
            

        } /// end loop over generated particles
    };

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskSigmaC>(cfgc, TaskName{"hf-task-sigmac"})};
}