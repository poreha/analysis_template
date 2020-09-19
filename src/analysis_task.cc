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
  fields_id_.insert(std::make_pair(FIELDS::CHARGE, mdc_vtx_tracks_config.GetFieldId("charge")));
  fields_id_.insert(std::make_pair(FIELDS::TOF_DE_DX, meta_hits_config.GetFieldId("dEdx")));
  fields_id_.insert(std::make_pair(FIELDS::MASS_2, meta_hits_config.GetFieldId("mass2")));
  fields_id_.insert(std::make_pair(FIELDS::BETA, meta_hits_config.GetFieldId("beta")));

  // initializing histograms
  tof_multiplicity_distribution_ = new TH1F( "tof_multiplicity", ";TOF hits;counts", 100, 0, 100 );
  pT_distribution_ = new TH1F( "pT_distribution", ";p_{T}, [GeV/c];counts", 250, 0, 2.5 );
  pT_vs_pseudorapidity_distribution_ = new TH2F( "pT_vs_pseudorapidity", ";#eta;p_{T}, [GeV/c]", 250, 0, 2.5, 250, 0.0, 2.5 );
  pT_vs_mass2_distribution_ = new TH2F( "pT_vs_mass2", ";p_{T} / z, [GeV/c];m^{2}, [GeV^{2}/c^{4}]", 200, -2.0, 2.0, 200, 0.0, 10.0 );
  p_vs_beta_distribution_ = new TH2F( "p_vs_beta", ";p / z, [GeV/c];#beta", 250, -1.5, 3.5, 240, 0.0, 1.2 );
  p_vs_dEdx_tof_distribution_ = new TH2F( "p_vs_dEdx_tof", ";p / z, [GeV/c];dE/dx", 250, -1.5, 3.5, 200, 0.0, 15.0 );
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
	auto charge = track.GetField<int>( fields_id_.at(FIELDS::CHARGE) ); // getting charge
	auto mass2 = hit.GetField<float>( fields_id_.at(FIELDS::MASS_2) ); // getting squared mass
	auto beta = hit.GetField<float>( fields_id_.at(FIELDS::BETA) ); // getting v/c
	auto dEdx = hit.GetField<float>( fields_id_.at(FIELDS::TOF_DE_DX) ); // getting energy loss in TOF-system
	// filling distributions
	pT_distribution_->Fill(pT);
	pT_vs_pseudorapidity_distribution_->Fill(eta, pT);
	pT_vs_mass2_distribution_->Fill( pT/charge, mass2 );
	p_vs_beta_distribution_->Fill( p/charge, beta );
	p_vs_dEdx_tof_distribution_->Fill( p/charge, dEdx );
  }
}

void AnalysisTask::Finish() {
  // Writing histograms to file
  tof_multiplicity_distribution_->Write();
  pT_distribution_->Write();
  pT_vs_pseudorapidity_distribution_->Write();
  pT_vs_mass2_distribution_->Write();
  p_vs_beta_distribution_->Write();
  p_vs_dEdx_tof_distribution_->Write();
}
} // namespace AnalysisTree