//
// Created by mikhail on 6/16/20.
//

#include "analysis_task.h"

namespace AnalysisTree {
void AnalysisTask::Init(std::map<std::string, void *> &branch_map) {
  // linking pointers with branch fields
  event_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
  mdc_vtx_tracks_ = static_cast<Particles *>(branch_map.at("mdc_vtx_tracks"));
  meta_hits_ = static_cast<HitDetector *>(branch_map.at("meta_hits"));
  mdc_meta_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2meta_hits"));

  // getting branch configurations, which store information about fields in branches
  auto event_header_config = config_->GetBranchConfig("event_header");
  auto mdc_vtx_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");
  auto meta_hits_config = config_->GetBranchConfig("meta_hits");

  // linking necessary for analysis fields with enumerator for fast access to them
  fields_id_.insert(std::make_pair(FIELDS::HITS_TOF, event_header_config.GetFieldId("selected_tof_hits")));

  // initializing histograms
  tof_multiplicity_distribution_ = new TH1F( "tof_multiplicity", ";TOF hits;counts", 100, 0, 100 );
  pT_distribution_ = new TH1F( "pT_distribution", ";p_{T}, [GeV/c];counts", 250, 0, 2.5 );
}

void AnalysisTask::Exec() {
  auto hits_tof = event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_TOF)); // getting multiplicity from event header
  tof_multiplicity_distribution_->Fill(hits_tof); // filling histogram

  int n_tracks = mdc_vtx_tracks_->GetNumberOfChannels(); // number of tracks in current event
  for (int i = 0; i < n_tracks; ++i) { // loop over all tracks if current event
	auto track = mdc_vtx_tracks_->GetChannel(i); // getting track from track detector
	int match_meta_hit = mdc_meta_matching_->GetMatchDirect(i); // getting index of matched with track TOF-system hit
	auto hit = meta_hits_->GetChannel(i); // getting matched with track hit in TOF-system
	auto pT = track.GetPt(); // getting transverse momentum
	auto eta = track.GetEta(); // getting pseudorapidity
	auto p = track.GetP(); // getting absolute value of momentum
	// filling distributions
	pT_distribution_->Fill(pT);
  }
}

void AnalysisTask::Finish() {
  // Writing histograms to file
  tof_multiplicity_distribution_->Write();
  pT_distribution_->Write();
}
} // namespace AnalysisTree