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

//
// Task to add a table of track parameters propagated to the primary vertex
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"

// The Run 3 AO2D stores the tracks at the point of innermost update. For a track with ITS this is the innermost (or second innermost)
// ITS layer. For a track without ITS, this is the TPC inner wall or for loopers in the TPC even a radius beyond that.
// In order to use the track parameters, the tracks have to be propagated to the collision vertex which is done by this task.
// The task consumes the TracksIU and TracksCovIU tables and produces Tracks and TracksCov to which then the user analysis can subscribe.
//
// This task is not needed for Run 2 converted data.
// There are two versions of the task (see process flags), one producing also the covariance matrix and the other only the tracks table.

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

struct TrackPropagation {
  Produces<aod::StoredTracks> tracksParPropagated;
  Produces<aod::TracksExtension> tracksParExtensionPropagated;

  Produces<aod::StoredTracksCov> tracksParCovPropagated;
  Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;

  Produces<aod::TracksDCA> tracksDCA;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  bool fillTracksDCA = false;
  int runNumber = -1;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  const o2::dataformats::MeanVertexObject* mVtx = nullptr;
  o2::parameters::GRPObject* grpo = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  HistogramRegistry histos;

  void init(o2::framework::InitContext& initContext)
  {
    if (doprocessCovariance == true && doprocessStandard == true) {
      LOGF(fatal, "Cannot enable processStandard and processCovariance at the same time. Please choose one.");
    }

    // Checking if the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "TracksDCA") {
          fillTracksDCA = true;
        }
      }
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    auto hTrkParsBug = histos.add<TH1>("hTrkParsBug", "", kTH1D, {{24, -1.5, 22.5}});
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(-1), "Tracks");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(0), "Bugged tracks");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(1), "Bug x");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(2), "Bug y");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(3), "Bug z");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(4), "Bug snp");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(5), "Bug tgl");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(6), "Bug signed1Pt");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(7), "Bug alpha");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(8), "Bug cYY");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(9), "Bug cZY");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(10), "Bug cZZ");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(11), "Bug cSnpY");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(12), "Bug cSnpZ");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(13), "Bug cSnpSnp");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(14), "Bug cTglY");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(15), "Bug cTglZ");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(16), "Bug cTglSnp");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(17), "Bug cTglTgl");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(18), "Bug c1PtY");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(19), "Bug c1PtZ");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(20), "Bug c1PtSnp");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(21), "Bug c1PtTgl");
    hTrkParsBug->GetXaxis()->SetBinLabel(hTrkParsBug->FindBin(22), "Bug c1Pt21Pt2");
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, bc.timestamp());
    LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object", grpo->getNominalL3Field(), bc.runNumber());
    o2::base::Propagator::initFieldFromGRP(grpo);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  template <typename TTrack, typename TTrackPar>
  void FillTracksPar(TTrack& track, TTrackPar& trackPar)
  {
    tracksParPropagated(track.collisionId(), track.trackType(), trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
    tracksParExtensionPropagated(trackPar.getPt(), trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
  }

  /// Protection against NaN trackIU parameters (no cov. matrix)
  void ProtectTrackPar(o2::track::TrackParF& trackPar)
  {

    histos.fill(HIST("hTrkParsBug"), -1); // all tracks

    bool isTrkBug = false;
    if (isnan(trackPar.getX())) {
      trackPar.setX(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 1);
    }
    if (isnan(trackPar.getY())) {
      trackPar.setY(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 2);
    }
    if (isnan(trackPar.getZ())) {
      trackPar.setZ(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 3);
    }
    if (isnan(trackPar.getSnp())) {
      trackPar.setSnp(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 4);
    }
    if (isnan(trackPar.getTgl())) {
      trackPar.setTgl(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 5);
    }
    if (isnan(trackPar.getQ2Pt())) {
      trackPar.setQ2Pt(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 6);
    }
    if (isnan(trackPar.getAlpha())) {
      trackPar.setAlpha(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 7);
    }

    if (isTrkBug) {
      histos.fill(HIST("hTrkParsBug"), 0); // all bugged tracks
    }
  }

  /// Protection against NaN trackIU parameters (with cov. matrix)
  void ProtectTrackParCov(o2::track::TrackParCovF& trackParCov)
  {

    histos.fill(HIST("hTrkParsBug"), -1); // all tracks

    bool isTrkBug = false;
    if (isnan(trackParCov.getX())) {
      trackParCov.setX(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 1);
    }
    if (isnan(trackParCov.getY())) {
      trackParCov.setY(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 2);
    }
    if (isnan(trackParCov.getZ())) {
      trackParCov.setZ(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 3);
    }
    if (isnan(trackParCov.getSnp())) {
      trackParCov.setSnp(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 4);
    }
    if (isnan(trackParCov.getTgl())) {
      trackParCov.setTgl(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 5);
    }
    if (isnan(trackParCov.getQ2Pt())) {
      trackParCov.setQ2Pt(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 6);
    }
    if (isnan(trackParCov.getAlpha())) {
      trackParCov.setAlpha(999.);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 7);
    }
    if (isnan(trackParCov.getSigmaY2())) {
      trackParCov.setCov(999., o2::track::kSigY2);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 8);
    }
    if (isnan(trackParCov.getSigmaZY())) {
      trackParCov.setCov(999., o2::track::kSigZY);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 9);
    }
    if (isnan(trackParCov.getSigmaZ2())) {
      trackParCov.setCov(999., o2::track::kSigZ2);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 10);
    }
    if (isnan(trackParCov.getSigmaSnpY())) {
      trackParCov.setCov(999., o2::track::kSigSnpY);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 11);
    }
    if (isnan(trackParCov.getSigmaSnpZ())) {
      trackParCov.setCov(999., o2::track::kSigSnpZ);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 12);
    }
    if (isnan(trackParCov.getSigmaSnp2())) {
      trackParCov.setCov(999., o2::track::kSigSnp2);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 13);
    }
    if (isnan(trackParCov.getSigmaTglY())) {
      trackParCov.setCov(999., o2::track::kSigTglY);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 14);
    }
    if (isnan(trackParCov.getSigmaTglZ())) {
      trackParCov.setCov(999., o2::track::kSigTglZ);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 15);
    }
    if (isnan(trackParCov.getSigmaTglSnp())) {
      trackParCov.setCov(999., o2::track::kSigTglSnp);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 16);
    }
    if (isnan(trackParCov.getSigmaTgl2())) {
      trackParCov.setCov(999., o2::track::kSigTgl2);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 17);
    }
    if (isnan(trackParCov.getSigma1PtY())) {
      trackParCov.setCov(999., o2::track::kSigQ2PtY);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 18);
    }
    if (isnan(trackParCov.getSigma1PtZ())) {
      trackParCov.setCov(999., o2::track::kSigQ2PtZ);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 19);
    }
    if (isnan(trackParCov.getSigma1PtSnp())) {
      trackParCov.setCov(999., o2::track::kSigQ2PtSnp);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 20);
    }
    if (isnan(trackParCov.getSigma1PtTgl())) {
      trackParCov.setCov(999., o2::track::kSigQ2PtTgl);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 21);
    }
    if (isnan(trackParCov.getSigma1Pt2())) {
      trackParCov.setCov(999., o2::track::kSigQ2Pt2);
      isTrkBug = true;
      histos.fill(HIST("hTrkParsBug"), 21);
    }

    if (isTrkBug) {
      histos.fill(HIST("hTrkParsBug"), 0); // all bugged tracks
    }
  }

  void processStandard(aod::StoredTracksIU const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    gpu::gpustd::array<float, 2> dcaInfo;

    for (auto& track : tracks) {
      dcaInfo[0] = 999;
      dcaInfo[1] = 999;
      auto trackPar = getTrackPar(track);
      ProtectTrackPar(trackPar);
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.x() < o2::constants::geom::XTPCInnerRef + 0.1) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        } else {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({mVtx->getX(), mVtx->getY(), mVtx->getZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        }
      }
      FillTracksPar(track, trackPar);
      if (fillTracksDCA) {
        tracksDCA(dcaInfo[0], dcaInfo[1]);
      }
    }
  }
  PROCESS_SWITCH(TrackPropagation, processStandard, "Process without covariance", true);

  void processCovariance(soa::Join<aod::StoredTracksIU, aod::TracksCovIU> const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;

    for (auto& track : tracks) {
      dcaInfoCov.set(999, 999, 999, 999, 999);
      auto trackParCov = getTrackParCov(track);
      ProtectTrackParCov(trackParCov);
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.x() < o2::constants::geom::XTPCInnerRef + 0.1) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          vtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          vtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, 2.f, matCorr, &dcaInfoCov);
        } else {
          vtx.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
          vtx.setCov(mVtx->getSigmaX() * mVtx->getSigmaX(), 0.0f, mVtx->getSigmaY() * mVtx->getSigmaY(), 0.0f, 0.0f, mVtx->getSigmaZ() * mVtx->getSigmaZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, 2.f, matCorr, &dcaInfoCov);
        }
      }
      FillTracksPar(track, trackParCov);
      if (fillTracksDCA) {
        tracksDCA(dcaInfoCov.getY(), dcaInfoCov.getZ());
      }
      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCovPropagated(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                             std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtensionPropagated(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                                      trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                                      trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                                      trackParCov.getSigma1Pt2());
    }
  }
  PROCESS_SWITCH(TrackPropagation, processCovariance, "Process with covariance", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackPropagation>(cfgc)};
  return workflow;
}
