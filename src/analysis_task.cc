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
  fields_id_.insert(std::make_pair(FIELDS::BETA, meta_hits_config.GetFieldId("beta"))); //new edit looking for beta
  fields_id_.insert(std::make_pair(FIELDS::M2, meta_hits_config.GetFieldId("mass2")));
  // initializing histograms
  tof_multiplicity_distribution_ = new TH1F( "tof_multiplicity", ";TOF hits;counts", 100, 0, 100 );
  pT_distribution_ = new TH1F( "pT_distribution", ";p_{T}, [GeV/c];counts", 256, 0, 2.5 );
  phi_distribution_ = new TH1F("phi_distribution", ";phi; grad", 256, -4, 4);
  beta_distribution_ = new TH1F("beta", "beta; counts", 256, -1, 10);
  mass2_distribution_ = new TH1F("mass_squarred", "mass^2, [GeV]", 1024, -1, 13);
  p_distribution_ = new TH1F("p_distribution", "momentum, [GeV/c]", 1024, -50, 300);

  PTvRAPIDITY_ = new TH2F("Rapidity_vs_pT", ";Rapidity;p_{T}, GeV/c", 250, 0, 2, 250, 0, 1.8);
  PTvPSEUDORAPIDITY_ = new TH2F("Pseudo-Rapidity_vs_pT", ";Psudo-Rapidity;p_{T}, GeV/c", 250, 0, 2, 250, 0, 1.8);
  PHIvPT_ = new TH2F("Phi_vs_pT", ";Phi;p_{T}, GeV", 256, -3.5, 3.5, 250, 0, 5);
  MOMENTUMvBETA_ = new TH2F("Momentum_vs_beta", ";p, GeV/c; beta", 256, -3, 6, 256, 0, 1);
  M2vMOMENTUM_ = new TH2F("Momentum_vs_Mass2", ";p, GeV/c;mass^2, GeV", 256, 0, 6, 512, -10, 60); // check it
  RAPIDITYvPHI_ = new TH2F("Rapidity_vs_Phi", ";Phi ;Rapidity", 256, -3.5, 3.5, 256, 0, 2);
}

void AnalysisTask::Exec() {
  auto hits_tof = event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_TOF)); // getting multiplicity from event header
  tof_multiplicity_distribution_->Fill(hits_tof); // filling histogram
  
  int n_tracks = mdc_vtx_tracks_->GetNumberOfChannels(); // number of tracks in current event
  for (int i = 0; i < n_tracks; ++i) { // loop over all tracks if current event
	  auto track = mdc_vtx_tracks_->GetChannel(i); // getting track from track detector
	  
    int match_meta_hit = mdc_meta_matching_->GetMatchDirect(i); // getting index of matched with track TOF-system hit
    auto hit = meta_hits_->GetChannel(i); // getting matched with track hit in TOF-system
	  auto charge = hit.GetField<int>(0);
    auto pT = track.GetPt(); // getting transverse momentum
	  auto eta = track.GetEta(); // getting pseudorapidity
    auto rapidity = track.GetRapidity(); // getting rapidity
	  auto p = track.GetP(); // getting absolute value of momentum
    auto phi = track.GetPhi(); // getting phi

    auto beta = hit.GetField<float>(fields_id_.at(FIELDS::BETA)); // getting beta from meta_hits
    auto m2 = hit.GetField<float>(fields_id_.at(FIELDS::M2)); // getting mass squarred from meta_hits
	  m2 *= charge;
    // filling distributions
    if (pid == 211)//pie-mesons
    {
    PHIvPT_->Fill(phi, pT);
    PTvRAPIDITY_->Fill(rapidity, pT);
    PTvPSEUDORAPIDITY_->Fill(eta, pT);
    MOMENTUMvBETA_->Fill(p*charge, beta); // *charge to see negatively charged
    M2vMOMENTUM_->Fill(p, m2);
    RAPIDITYvPHI_->Fill(phi, rapidity);
    }
  }
}

void AnalysisTask::Finish() {
  // Writing histograms to file
  
  PTvRAPIDITY_->Write();
  PTvPSEUDORAPIDITY_->Write();
  PHIvPT_->Write();
  MOMENTUMvBETA_->Write();
  M2vMOMENTUM_->Write();
  RAPIDITYvPHI_->Write();
  //tof_multiplicity_distribution_->Write();
  //pT_distribution_->Write();
  //beta_distribution_->Write();
  //mass2_distribution_->Write();
  //p_distribution_->Write();
}
} // namespace AnalysisTree