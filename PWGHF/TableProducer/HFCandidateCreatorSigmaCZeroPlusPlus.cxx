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

/// \file HFCandidateCreatorSigmaCZeroPlusPlus.cxx
/// \brief Σc0,++ → Λc+(→pK-π+) π-,+ analysis task
/// \note Λc± candidates selected from the HFLcCandidateSelector.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN PADOVA

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong3;

struct HFCandidateCreatorSigmaCZeroPlusPlus {

    /// Table with Σc0,++ info
    Produces<aod::HfCandScBase> rowScCandidates;

    /// Selection of candidates Λc+
    Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
    Configurable<double> cutYCandLcMax{"cutYCandLcMax", -1., "max. candLc. Lc rapidity"};
    Configurable<double> cutMinvpKpiCandLcMax{"cutMinvpKpiCandLcMax", 0.03, "max. spread (abs. value9) between PDG(Lc) and Minv(pKpi)"};
    Configurable<double> cutMinvpK0sCandLcMax{"cutMinvpK0sCandLcMax", 0.03, "max. spread (abs. value9) between PDG(Lc) and Minv(pK0s)"};

    /// Selections on candidate soft π-,+
    Configurable<float> softPiEta{"softPiEta", 0.9f, "Soft pion max value for pseudorapidity (abs vale)"};
    Configurable<int> softPiITSHitMap{"softPiITSHitMap", 127, "Soft pion ITS hitmap"};
    Configurable<int> softPiITSHitsMin{"softPiITSHitsMin", 1, "Minimum number of ITS layers crossed by the soft pion among those in \"softPiITSHitMap\""};
    Configurable<float> softPidcaXYmax{"softPidcaXYmax", 0.065, "Soft pion max dcaXY (cm)"};
    Configurable<float> softPidcaZmax{"softPidcaZmax", 0.065, "Soft pion max dcaZ (cm)"};

    /// Filter the candidate Λc+ used for the Σc0,++ creation
    Filter filterSelectCandidateLc = (aod::hf_selcandidate_lc::isSelLcpKpi >= selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= selectionFlagLc);

    /// Cut selection object for soft π-,+
    TrackSelection softPiCuts;

    /// @brief init function, to define the soft pion selections and histograms
    /// @param  
    void init(InitContext&) {

        ////////////////////////////////////////
        /// set the selections for soft pion ///
        ////////////////////////////////////////
        softPiCuts.SetEtaRange(-softPiEta, softPiEta);  // eta
        softPiCuts.SetMaxDcaXY(softPidcaXYmax); // dcaXY
        softPiCuts.SetMaxDcaZ(softPidcaZmax);   // dcaZ
        // ITS hitmap
        std::set<uint8_t> set_softPiITSHitMap; // = {};
        for (int ITSlayerId = 0; ITSlayerId < 7; ITSlayerId++) {
            if ((softPiITSHitMap & (1 << ITSlayerId)) > 0) {
              set_softPiITSHitMap.insert(static_cast<uint8_t>(ITSlayerId));
            }
        }
        LOG(info) << "### ITS hitmap for soft pion";
        LOG(info) << "    >>> set_softPiITSHitMap.size(): " << set_softPiITSHitMap.size();
        LOG(info) << "    >>> Custom ITS hitmap checked: ";
        for (std::set<uint8_t>::iterator it = set_softPiITSHitMap.begin(); it != set_softPiITSHitMap.end(); it++) {
          LOG(info) << "        Layer " << (int)(*it) << " ";
        }
        LOG(info) << "############";
        softPiCuts.SetRequireITSRefit();
        softPiCuts.SetRequireHitsInITSLayers(softPiITSHitsMin, set_softPiITSHitMap);
        
    }

    /// @brief process function for Σc0,++ → Λc+(→pK-π+) π- candidate reconstruction
    /// @param collision is a o2::aod::Collision
    /// @param tracks are the tracks (with dcaXY, dcaZ information) in the collision → soft-pion candidate tracks
    /// @param candidates are 3-prong candidates satisfying the analysis selections for Λc+ → pK-π+ (and charge conj.)
    void process(const o2::aod::Collision& collision, const soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>> const& candidates) {

      /// loop over Λc+ → pK-π+ (and charge conj.) candidates
      for(auto& candLc : candidates) {

        /// keep only the candidates flagged as possible Λc+ (and charge conj.) decaying into a charged pion, kaon and proton
        /// if not selected, skip it and go to the next one
        if (!(candLc.hfflag() & 1 << DecayType::LcToPKPi)) {
          continue;
        }
        /// keep only the candidates Λc+ (and charge conj.) within the desired rapidity
        /// if not selected, skip it and go to the next one
        if (cutYCandLcMax >= 0. && std::abs(YLc(candLc)) > cutYCandLcMax) {
          continue;
        }

        /// selection on the Λc+ inv. mass window we want to consider for Σc0,++ candidate creation
        auto statusSpreadMinvpKpiFromPDG = 0;
        auto statusSpreadMinvpiKpFromPDG = 0;
        if(candLc.isSelLcpKpi() >= 1 && std::abs( InvMassLcpKpi(candLc) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus) ) <= statusSpreadMinvpKpiFromPDG) {
          statusSpreadMinvpKpiFromPDG = 1;
        }
        if(candLc.isSelLcpiKp() >= 1 && std::abs( InvMassLcpiKp(candLc) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus) ) <= statusSpreadMinvpiKpFromPDG) {
          statusSpreadMinvpiKpFromPDG = 1;
        }
        if (statusSpreadMinvpKpiFromPDG == 0 && statusSpreadMinvpiKpFromPDG == 0){
          /// none of the two possibilities are satisfied, therefore this candidate Lc can be skipped
          continue;
        }

        /// loop over tracks
        for(auto& trackSoftPi : tracks) {

          /////////////////////////////////////////////////////////////////////////////////
          ///                       Σc0,++ candidate creation                           ///
          ///                                                                           ///
          /// For each candidate Λc, let's loop over all the candidate soft-pion tracks ///
          /////////////////////////////////////////////////////////////////////////////////

          /// keep only soft-pion candidate tracks
          /// if not selected, skip it and go to the next one
          if( !softPiCuts.IsSelected(trackSoftPi) ) {
            continue;
          }

          /// Exclude the current candidate soft pion if it corresponds already to a candidate Lc prong
          int indexProng0 = candLc.index0_as<aod::Tracks>().globalIndex();
          int indexProng1 = candLc.index1_as<aod::Tracks>().globalIndex();
          int indexProng2 = candLc.index2_as<aod::Tracks>().globalIndex();
          int indexSoftPi = trackSoftPi.globalIndex();
          if (indexSoftPi == indexProng0 || indexSoftPi == indexProng1 || indexSoftPi == indexProng2)
            continue;

          /// determine the Σc candidate charge
          int chargeLc = candLc.index0_as<aod::Tracks>().sign() + candLc.index1_as<aod::Tracks>().sign() + candLc.index2_as<aod::Tracks>().sign();
          int chargeSoftPi = trackSoftPi.sign();
          int chargeSc = chargeLc + chargeSoftPi;
          if(chargeSc != 0 && chargeSc != 2) {
            /// this shall never happen
            LOG(fatal) << ">>> SigmaC candidate with charge +1 built, not possible! Charge Lc: " << chargeLc << ", charge soft pion: " << chargeSoftPi;
            continue;
          }

          ///
          /// TODO: fill the table
          ///
          rowScCandidates(/* general columns */
                          candLc.collisionId(),
                          /* 2-prong specific columns */
                          candLc.px(), candLc.py(), candLc.pz(),
                          trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz(),
                          candLc.globalIndex(), trackSoftPi.globalIndex(),
                          candLc.hfflag(),
                          /* Σc0,++ specific columns */
                          chargeSc,
                          statusSpreadMinvpKpiFromPDG, statusSpreadMinvpiKpFromPDG);


          /// THE FOLLOWING THINGS CAN BE MOVED TO THE ANALYSIS TASK!!!
          /// HERE USE THE FUNCTIONS FROM hf_can_sc !!! They are templated and calculated on-the-fly (dynamic columns)

          /// kinematic properties for Σc0,++
          //double energySoftPi = RecoDecay::e(trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz(), RecoDecay::getMassPDG(kPiPlus));
          //double ptLc = candLc.pt();
          //double ptSc = RecoDecay::pt(candLc.px()+trackSoftPi.px(), candLc.py()+trackSoftPi.py());
          //double ySc = RecoDecay::y(std::array{candLc.px()+trackSoftPi.px(), candLc.py()+trackSoftPi.py(), candLc.pz()+trackSoftPi.pz()}, RecoDecay::getMassPDG(pdgSc));
          //if(candLc.isSelLcpKpi() >= selectionFlagLc) { /// first prong: proton; third prong: pion
          //  massSc = RecoDecay::m( array{ array{candLc.px(), candLc.py(), candLc.pz()}, array{trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz()} }, array{(double)InvMassLcpKpi(candLc), RecoDecay::getMassPDG(kPiPlus)});
          //  //
          //  //  TODO: fill the table with ptSc and massSc
          //  //
          //}
          //if(candLc.isSelLcpiKp() >= selectionFlagLc) { /// first prong: pion; third prong: proton
          //  massSc = RecoDecay::m( array{ array{candLc.px(), candLc.py(), candLc.pz()}, array{trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz()} }, array{(double)InvMassLcpiKp(candLc), RecoDecay::getMassPDG(kPiPlus)});
          //  //
          //  //  TODO: fill the table with ptSc and massSc
          //  //
          //}
          

        } /// end loop over tracks

      } /// end loop over candidtes

    } /// end process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFCandidateCreatorSigmaCZeroPlusPlus>(cfgc, TaskName{"hf-sigmac-candidate-creator"})};
}