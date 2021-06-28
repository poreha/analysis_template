//
// Created by mikhail on 6/16/20.
// Edited by alexey on 2/25/21
//

#include "analysis_task.h"
#include <TCanvas.h>
#include <math.h>

namespace AnalysisTree
{
  void AnalysisTask::Init(std::map<std::string, void *> &branch_map)
  {
    // linking pointers with branch fields
    event_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
    mdc_vtx_tracks_ = static_cast<Particles *>(branch_map.at("mdc_vtx_tracks"));
    meta_hits_ = static_cast<HitDetector *>(branch_map.at("meta_hits"));
    mdc_meta_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2meta_hits"));

    // for simulated
    sim_header_ = static_cast<EventHeader *>(branch_map.at("sim_header"));
    sim_tracks_ = static_cast<Particles *>(branch_map.at("sim_tracks"));
    track_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2sim_tracks"));

    // getting branch configurations, which store information about fields in branches
    auto event_header_config = config_->GetBranchConfig("event_header");
    auto mdc_vtx_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");
    auto meta_hits_config = config_->GetBranchConfig("meta_hits");

    auto sim_header_config = config_->GetBranchConfig("sim_header");
    auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");

    // linking necessary for analysis fields with enumerator for fast access to them
    fields_id_.insert(std::make_pair(FIELDS::HITS_TOF, event_header_config.GetFieldId("selected_tof_hits")));
    fields_id_.insert(std::make_pair(FIELDS::BETA, meta_hits_config.GetFieldId("beta")));
    fields_id_.insert(std::make_pair(FIELDS::M2, meta_hits_config.GetFieldId("mass2")));
    fields_id_.insert(std::make_pair(FIELDS::CHARGE, meta_hits_config.GetFieldId("charge")));
    fields_id_.insert(std::make_pair(FIELDS::PT2, event_header_config.GetFieldId("physical_trigger_2")));
    fields_id_.insert(std::make_pair(FIELDS::TOF, meta_hits_config.GetFieldId("time_of_flight")));
    fields_id_.insert(std::make_pair(FIELDS::PATH_LEN, meta_hits_config.GetFieldId("path_length")));

    // initializing histograms
    tof_multiplicity_distribution_ = new TH1F("tof_multiplicity", ";TOF hits;counts", 100, 0, 100);

    //pT_distribution_ = new TH1F( "pT_distribution", ";p_{T}, [GeV/c];counts", 256, 0, 2.5 );
    //phi_distribution_ = new TH1F("phi_distribution", ";phi; grad", 256, -4, 4);
    beta_distribution_ = new TH1F("beta", "beta; counts", 256, -1, 10);
    mass2_distribution_ = new TH1F("mass_squarred", "mass^2, [GeV]", 1024, -6, 6);
    mass2_distribution_branch_ = new TH1F("mass_squarred_branch", "mass^2, [GeV]", 1024, -6, 6);
    mass_comparrison_ = new TH1F("mass_comparrison", "delta_mass^2, [what]", 1024, -6, 6);
    //p_distribution_ = new TH1F("p_distribution", "momentum, [GeV/c]", 1024, -50, 300);
    impact_parameter_ = new TH1F("impact parameter", "?", 1024, 0, 9);
    selected_tof_rpc_centrality_ = new TH1F("centrality_tofrpc", "?", 1024, 0, 10);

    IMPvPT_ = new TH2F("Impact_parameter_vs_Pt","IMPvPt;impact parameter;p_{T}",512,0,9,250,0,2.);

    PTvRAPIDITY_ = new TH2F("Rapidity_vs_pT", ";Rapidity;p_{T}, GeV/c", 250, 0, 2.5, 250, 0, 2.);
    PTvPSEUDORAPIDITY_ = new TH2F("Pseudo-Rapidity_vs_pT", ";Pseudo-Rapidity;p_{T}, GeV/c", 250, 0, 2.5, 250, 0, 2.);

    // use legend
    PID_RECO_ = new TH2F("PID_reco", "All reconstructed protons;Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);
    PID_PRIM_ = new TH2F("PID_prim", "Reconstructed primary protons matched with MC-primary particles;Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);
    PID_SEC_ = new TH2F("PID_sec", "Reconstructed secondary protons matched with MC-secondary particles;Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);

    PID_PRIM_PDG_ = new TH2F("PID_prim_pdg", "Reconstructed primary protons matched with MC-primary particles (NO CUTS);Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);
    PID_SEC_PDG_ = new TH2F("PID_sec_pdg", "Reconstructed secondary protons matched with MC-secondary particles (NO CUTS);Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);
    PID_PRIM_CUTS_ = new TH2F("PID_prim_cuts", "Reconstructed primary protons matched with MC-primary particles (ONLY CUTS);Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);
    PID_SEC_CUTS_ = new TH2F("PID_sec_cuts", "Reconstructed secondary protons matched with MC-secondary particles (ONLY CUTS);Y;p_{T}, [GeV/c]", 120, 0., 2.5, 120, 0., 2.5);

    PDG_PRIM_ = new TH2F("PDG_prim_reco", "Particles matched with MC-primary proton;Y;p_{T}, [GeV/c]", 120, 0., 2.5, 250, 0., 2.5);
    PDG_SEC_ = new TH2F("PDG_sec_reco", "Particles matched with MC-secondary proton;Y;p_{T}, [GeV/c]", 120, 0., 2.5, 250, 0., 2.5);

    MISMATCH_ = new TH2F("Mismatch", "Protons matched with all MC-particles except protons;Y;p_{T}, [GeV/c]", 120, 0, 2.5, 120, 0, 2.5);
    LOSS_ = new TH2F("Loss", "All MC-particles, except protons, matched with MC-protons;Y;p_{T}, [GeV/c]", 120, 0, 2.5, 120, 0, 2.5);

    GEN_ = new TH2F("PDG", "All MC-protons;Y;p_{T}, [GeV/c]", 120, -1.5, 2.5, 120, 0, 2.5);
    GEN_PRIM_ = new TH2F("PDG_prim", "All primary MC-protons;Y;p_{T}, [GeV/c]", 120, -1.5, 2.5, 120, 0., 2.5);
    GEN_SEC_ = new TH2F("PDG_sec", "All secondary MC-protons;Y;p_{T}, [GeV/c]", 120, -1.5, 2.5, 120, 0, 2.5);

    PHIvPT_ = new TH2F("Phi_vs_pT", ";Phi;p_{T}, GeV", 256, -3.5, 3.5, 250, 0, 5);
    MOMENTUMvBETA_ = new TH2F("Momentum_vs_beta", ";p/q, GeV/c; beta", 256, -3, 6.2, 256, 0, 1.3);
    M2vMOMENTUM_ = new TH2F("Momentum_vs_Mass2", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60); // check it
    RAPIDITYvPHI_ = new TH2F("Rapidity_vs_Phi", ";Phi ;Rapidity", 256, -3.5, 3.5, 256, 0, 2);

    M2vP_PID_PRIM_ = new TH2F("pid_prim", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_PID_SEC_ = new TH2F("pid_sec", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_PDG_PRIM_ = new TH2F("pdg_prim", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_PDG_SEC_ = new TH2F("pdg_sec", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_MISMATCH_ = new TH2F("mismatch", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_MISMATCH_SIM_ = new TH2F("mismatch_sim", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_LOSS_ = new TH2F("loss", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_LOSS_SIM_ = new TH2F("loss_sim", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);

    M2vP_GEN_ = new TH2F("gen", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_GEN_PRIM_ = new TH2F("gen_prim", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);
    M2vP_GEN_SEC_ = new TH2F("gen_sec", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60);

    PvTOF_ = new TH2F("Momentum_vs_TimeOfFlight", ";p/q, GeV/c; t", 256, -6, 6, 512, 0, 100);
    PvPathLen_ = new TH2F("Momentum_vs_PathLength", ";p/q, GeV/c; path", 256, -6, 6, 512, 1000, 3000);
  }

  void AnalysisTask::Exec()
  {
    int proton = 2212; //int particle ID
    //auto hits_tof = event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_TOF)); // getting multiplicity from event header
    //tof_multiplicity_distribution_->Fill(hits_tof); // filling histogram

//auto im_param = sim_header_->GetField<float>(1); // impact parameter // посмотреть по событиям, а не по трекам

    int n_tracks = mdc_vtx_tracks_->GetNumberOfChannels(); // number of tracks in current event
    for (int i = 0; i < n_tracks; ++i)
    {                                              // loop over all tracks if current event
      auto track = mdc_vtx_tracks_->GetChannel(i); // getting track from track detector
      int pid = track.GetPid();

      int matched_track = track_matching_->GetMatchDirect(i); // matching simulated track's ID to current reconstructed track

      //  для no_match(id = -999) ввести обработку
      if (matched_track == -999)
        continue;

      auto sim_track = sim_tracks_->GetChannel(matched_track);    // getting matched simulated track // changed i->matched_track
      int sim_pid = sim_track.GetPid();                           // getting PID of matched simulated track
      int match_meta_hit = mdc_meta_matching_->GetMatchDirect(i); // getting index of matched with track TOF-system hit

      auto hit = meta_hits_->GetChannel(match_meta_hit); // getting matched with track hit in TOF-system

      auto pT = track.GetPt();         // getting transverse momentum
      auto sim_pT = sim_track.GetPt(); // getting transverse momentum of matched simulated track

      auto eta = track.GetEta();                   // getting pseudorapidity
      auto rapidity = track.GetRapidity();         // getting rapidity
      auto sim_rapidity = sim_track.GetRapidity(); // getting rapidity of matched simulated track

      auto sim_p = sim_track.GetP();
      auto sim_m2 = sim_track.GetMass();

      auto p = track.GetP(); // getting absolute value of momentum

      auto phi = track.GetPhi(); // getting phi

      auto charge = hit.GetField<int>(fields_id_.at(FIELDS::CHARGE)); // getting charge from meta_hits
      auto beta = hit.GetField<float>(fields_id_.at(FIELDS::BETA));   // getting beta from meta_hits
      auto m2 = hit.GetField<float>(fields_id_.at(FIELDS::M2));       // getting mass squarred from meta_hits

      auto tof = hit.GetField<float>(fields_id_.at(FIELDS::TOF));           // getting time of flight from meta_hits
      auto path_len = hit.GetField<float>(fields_id_.at(FIELDS::PATH_LEN)); // check if zero

      p = p / charge;

      // play w filters and watch how phi-distr reacts

      bool pt2 = event_header_->GetField<bool>(5);             // physical_trigger_2 (id=5)
      bool kGoodVertexCand = event_header_->GetField<bool>(1); // good_vertex_candidate (id=1)
      bool kGoodTrigger = event_header_->GetField<bool>(13);   // good_trigger (id=13)
      bool kNoVeto = event_header_->GetField<bool>(15);        // no_veto (id=15)
      bool kNoPileUpMeta = event_header_->GetField<bool>(7);

      bool kGoodStart = event_header_->GetField<bool>(2);
      bool kGoodStartMeta = event_header_->GetField<bool>(17);

      float dca_xy = track.GetField<float>(3);
      float dca_z = track.GetField<float>(4);
      float chi2 = track.GetField<float>(0);
      // filling distributions
      /*
      filter: pie+ 211
              pie0 111
              proton 2212
              K0 311
              K+ 321
              e- 11
    */

      bool is_prim = sim_track.GetField<bool>(0); // sim_particle is a primary

      int particle = proton; // proton
      auto selected_centrality = event_header_->GetField<float>(0);
      selected_tof_rpc_centrality_->Fill(selected_centrality);
      // Filling histograms
      //mass2_distribution_branch_->Fill(m2);
      if (kNoVeto && kGoodVertexCand && kGoodStartMeta){
      if (pid == particle)
        M2vMOMENTUM_->Fill(p, m2);

      if (pid == particle)
      {
        PID_RECO_->Fill(rapidity, pT); // all reconstructed certain particles
        // PID
        if (dca_xy < 15. && dca_xy > -15. && dca_z < 15. && chi2 < 100.)
        {
          PID_PRIM_CUTS_->Fill(rapidity, pT);
          // if passes the cut as primary
          if (is_prim)
            {                                
              PID_PRIM_->Fill(rapidity, pT); // any matched certain particle is primary
              M2vP_PID_PRIM_->Fill(p, m2);
            }
        }
        else
        {
          PID_SEC_CUTS_->Fill(rapidity, pT);
          // if didn't pass the cut
          if (!is_prim)
          {
            PID_SEC_->Fill(rapidity, pT); // any matched certain particle is secondary
            M2vP_PID_SEC_->Fill(p, m2);
          }
        }

        if (is_prim) PID_PRIM_PDG_->Fill(rapidity, pT);
        else PID_SEC_PDG_->Fill(rapidity, pT);
      }
      // MISMATCH if pid == particle , sim_pid != particle
      else
      {
        MISMATCH_->Fill(rapidity, pT);
        M2vP_MISMATCH_->Fill(p, m2);
        M2vP_MISMATCH_SIM_->Fill(sim_p, sim_m2);
      }
      //PTvRAPIDITY_->Fill(rapidity, pT); // PTvRAPIDITY == PID_RECO
      PTvPSEUDORAPIDITY_->Fill(eta, pT);

      // LOSS
      if (pid != particle && sim_pid == particle)
      {
        LOSS_->Fill(rapidity, pT);
        M2vP_LOSS_->Fill(p, m2);
        M2vP_LOSS_SIM_->Fill(sim_p, sim_m2);
      }
      }
      // PDG
      if (sim_pid == particle)
      {
        if (is_prim)
        {
          PDG_PRIM_->Fill(sim_rapidity, sim_pT); // matched particle is primary // changed to sim_
          M2vP_PDG_PRIM_->Fill(p, m2);
        }
        else
        {
          PDG_SEC_->Fill(sim_rapidity, sim_pT); // matched particle is secondary
          M2vP_PDG_SEC_->Fill(p, m2);
        }
      }

      //PHIvPT_->Fill(phi, pT);
      //MOMENTUMvBETA_->Fill(p, beta); // charge to see negatively charged
      //M2vMOMENTUM_->Fill(p, m2);
      //RAPIDITYvPHI_->Fill(phi, rapidity);

      beta_distribution_->Fill(beta);
    }

    // GEN
    int sim_ntracks = sim_tracks_->GetNumberOfChannels();
    for (int sim_i = 0; sim_i < sim_ntracks; sim_i++)
    {
      auto sim_track = sim_tracks_->GetChannel(sim_i); // getting matched simulated track // changed i->matched_track
      int sim_pid = sim_track.GetPid();                // getting PID of matched simulated track

      auto sim_pT = sim_track.GetPt();             // getting transverse momentum of matched simulated track
      auto sim_rapidity = sim_track.GetRapidity(); // getting rapidity of matched simulated track

      auto sim_p = sim_track.GetP();
      auto sim_m2 = sim_track.GetMass();

      bool is_prim = sim_track.GetField<bool>(0); // sim_particle is a primary

      if (sim_pid == 2212)
      {

        GEN_->Fill(sim_rapidity, sim_pT);
        M2vP_GEN_->Fill(sim_p, sim_m2);
        if (is_prim)
        {
          GEN_PRIM_->Fill(sim_rapidity, sim_pT);
          M2vP_GEN_PRIM_->Fill(sim_p, sim_m2);
        }
        else
        {
          GEN_SEC_->Fill(sim_rapidity, sim_pT);
          M2vP_GEN_SEC_->Fill(sim_p, sim_m2);
        }
      }
    }
    
  }

  void AnalysisTask::Finish()
  {
    GEN_PRIM_->Sumw2();
    PDG_PRIM_->Sumw2();
    // Writing histograms to file
    selected_tof_rpc_centrality_->Write();
    
    PID_PRIM_PDG_->Write();
    PID_SEC_PDG_->Write();
    PID_PRIM_CUTS_->Write();
    PID_SEC_CUTS_->Write();

    PID_RECO_->Write();
    PID_PRIM_->Write();
    PID_SEC_->Write();
    PDG_PRIM_->Write();
    PDG_SEC_->Write();
    MISMATCH_->Write();
    LOSS_->Write();
    GEN_->Write();
    GEN_PRIM_->Write();
    GEN_SEC_->Write();

    //PTvRAPIDITY_->Write();
    PTvPSEUDORAPIDITY_->Write();
    M2vP_PID_PRIM_->Write();
    M2vP_PID_SEC_->Write();
    M2vP_PDG_PRIM_->Write();
    M2vP_PDG_SEC_->Write();
    M2vP_MISMATCH_->Write();
    M2vP_MISMATCH_SIM_->Write();
    M2vP_LOSS_->Write();
    M2vP_LOSS_SIM_->Write();

    M2vP_GEN_->Write();
    M2vP_GEN_PRIM_->Write();
    M2vP_GEN_SEC_->Write();

    //mass2_distribution_->Write();
    //mass2_distribution_branch_->Write();

    //PvTOF_->Write();
    //PvPathLen_->Write();

    //PHIvPT_->Write();
    //MOMENTUMvBETA_->Write();
    M2vMOMENTUM_->Write();
    //RAPIDITYvPHI_->Write();
    //tof_multiplicity_distribution_->Write();
    //pT_distribution_->Write();
    //beta_distribution_->Write();
    //mass2_distribution_->Write();
    //p_distribution_->Write();
  }
} // namespace AnalysisTree

/*
    TO DO
    make a Factory to create multiple histos
*/