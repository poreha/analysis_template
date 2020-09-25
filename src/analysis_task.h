//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <TChain.h>
#include <TH3F.h>
#include <TProfile2D.h>

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Matching.hpp>

namespace AnalysisTree {
class AnalysisTask : public FillTask{
public:
 AnalysisTask() = default;
  ~AnalysisTask() override = default;
  void Init( std::map<std::string, void*>& branch_map ) override;
  void Exec() override;
  void Finish() override;
private:
 enum class FIELDS{ // enumerator to fast access to detectors' fields
   HITS_TOF, 	// Hits in TOF-system
  };
  std::map<FIELDS, int> fields_id_; // map to match detectors' fields with enumerator

  /* pointers to link tree's branches with */
  EventHeader* event_header_{nullptr}; 		// event info
  Particles* mdc_vtx_tracks_{nullptr}; 		// tracks
  HitDetector* meta_hits_{nullptr}; 		// TOF-system
  Matching* mdc_meta_matching_{nullptr}; 	// matching between tracking system and TOF-system

  TH1F* tof_multiplicity_distribution_; 	// 1D histogram of sum multiplicity in TOF system
  TH1F* pT_distribution_; 					// 1D histogram of track transverse momentum distribution
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
