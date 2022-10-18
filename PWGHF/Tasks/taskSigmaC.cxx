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
    HistogramRegistry registry{};

    /// @brief init function, to define the analysis histograms
    /// @param  
    void init(InitContext&) {
        /// Σc0
        registry.add("hDeltaMassSigmaCZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {20, 0., 20.}}});
        /// Σc++
        registry.add("hDeltaMassSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {20, 0., 20.}}});
        /// Σc0,++
        registry.add("hDeltaMassSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {20, 0., 20.}}});
        /// Λc+ ← Σc0
        registry.add("hDeltaMassLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {20, 0., 20.}}});
        /// Λc+ ← Σc++
        registry.add("hDeltaMassLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {20, 0., 20.}}});
        /// Λc+ ← Σc0,++
        registry.add("hDeltaMassLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {20, 0., 20.}}});
    }; /// end init

    /// @brief process function to fill the histograms needed in analysis
    /// @param collision is the reconstruction collision
    /// @param candidatesSigmaC is the candidate SigmaC
    /// @param 
    void process(const aod::Collision& collision, const aod::HfCandScBase& candidatesSigmaC,
    soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate> const&, const soa::Join<aod::Tracks, aod::TracksDCA>&) {

        /// loop over the candidate Σc0,++
        for(auto& candSigmaC : candidatesSigmaC) {
            
            const int chargeSigmaC = candSigmaC.charge();   // either Σc0 or Σc++

            /// get the candidate Λc+ used to build the candidate Σc0,++
            /// and understand which mass hypotheses are possible
            const auto& candLambdaC = candSigmaC.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>>();
            const int isCandLambdaCpKpi = candLambdaC.isSelLcpKpi();
            const int isCandLambdaCpiKp = candLambdaC.isSelLcpiKp();
            double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.), ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
            if(isCandLambdaCpKpi >= 1 && candSigmaC.statusSpreadLcMinvpKpiFromPDG()) {
                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                massSigmaC = InvMassScRecoLcpKpi(candSigmaC);
                massLambdaC = InvMassLcpKpi(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → pK-π+ (and charge conjugate)
            if(isCandLambdaCpiKp >=1 && candSigmaC.statusSpreadLcMinvpiKpFromPDG()) {
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                massSigmaC = InvMassScRecoLcpiKp(candSigmaC);
                massLambdaC = InvMassLcpiKp(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → π+K-p (and charge conjugate)
        } /// end loop over the candidate Σc0,++
    };  /// end process

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskSigmaC>(cfgc, TaskName{"hf-task-sigmac"})};
}