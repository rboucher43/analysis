#include <fun4all/Fun4AllServer.h>

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
//#include <globalvertex/GlobalVertex.h>
//#include <globalvertex/GlobalVertexMap.h>
#pragma GCC diagnostic pop

R__LOAD_LIBRARY(libkfparticle_sphenix.so)

namespace HeavyFlavorReco
{
   // https://wiki.bnl.gov/sPHENIX/index.php/KFParticle
   string decayDescriptor = "[Upsilon -> e^- e^+]cc";  //See twiki on how to set this 
   string reconstructionName = "myTestReco";         //Used for naming output folder, file and node
   string outputRecoFile;
   string outputHFEffFile;
   bool runTruthTrigger = false;  //Decay Finder
   bool runTrackEff = false;  //HF track efficiency
   bool getTruthInfo = false;      //Add truth matching to output file
   bool getCaloInfo = false;
   bool runTracking = false;  //Run tracking on DSTs
   bool buildTruthTable = false;
   bool runQA = false;        //Run QA, needs set up
   int VERBOSITY = INT_MAX;
  //TH1D* vtx = new TH1D*("h", "Zvtx", "100", 0, 100);
};  // namespace HeavyFlavorReco

using namespace HeavyFlavorReco;

void myHeavyFlavorReco()
{
  Fun4AllServer *se = Fun4AllServer::instance();

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(reconstructionName);
  kfparticle->Verbosity(INT_MAX);

  kfparticle->setDecayDescriptor(decayDescriptor);

  if (runTrackEff) kfparticle->setTrackMapNodeName("HFSelected_SvtxTrackMap");

  kfparticle->doTruthMatching(getTruthInfo);
  kfparticle->getDetectorInfo(false);
  kfparticle->getCaloInfo(getCaloInfo);
  kfparticle->getAllPVInfo(false);
  kfparticle->allowZeroMassTracks(true);
  kfparticle->saveDST(false);
  kfparticle->saveParticleContainer(false);

  bool fixToPV = true;
  bool useFakePV = false;

  if (useFakePV)
  {
    fixToPV = false;  //Constraining to a fake PV results in some gibberish variables
    kfparticle->useFakePrimaryVertex(true);
  }

  if (fixToPV)
  {
    kfparticle->constrainToPrimaryVertex(true);
    kfparticle->setMotherIPchi2(FLT_MAX);
    kfparticle->setFlightDistancechi2(-1.0);
    kfparticle->setMinDIRA(-1.1);
  }

  //Track parameters
  kfparticle->setMinimumTrackPT(0);
  //kfparticle->setMinimumTrackIPchi2(0.0);
  kfparticle->setMinimumTrackIP(-0.01);
  kfparticle->setMaximumTrackchi2nDOF(10.0);

  //Vertex parameters
  kfparticle->setMaximumVertexchi2nDOF(10.0);
  kfparticle->setMaximumDaughterDCA(10.0);

  //Parent parameters
  kfparticle->setMotherPT(0);
  kfparticle->setMinimumMass(6);
  kfparticle->setMaximumMass(10);
  kfparticle->setDecayLengthRange(-1*FLT_MAX, FLT_MAX);
  kfparticle->setMaximumMotherVertexVolume(10.0);
  //Intermediate parameters
  /*std::vector<std::pair<float, float>> intermediate_mass_range;
  intermediate_mass_range.push_back(make_pair(0.4, 0.6));
  kfparticle->setIntermediateMassRange(intermediate_mass_range);
  std::vector<float> intermediate_min_pt = {0.5};
  kfparticle->setIntermediateMinPT(intermediate_min_pt);
  std::vector<std::pair<float, float>> intermediate_IP_range;
  intermediate_IP_range.push_back(make_pair(0., 1.));
  kfparticle->setIntermediateIPRange(intermediate_IP_range);
  std::vector<std::pair<float, float>> intermediate_IPchi2_range;
  intermediate_IPchi2_range.push_back(make_pair(0., 100.));
  kfparticle->setIntermediateIPchi2Range(intermediate_IPchi2_range);
  std::vector<float> intermediate_min_dira = {0.5};
  kfparticle->setIntermediateMinDIRA(intermediate_min_dira);
  std::vector<float> intermediate_min_FDchi2 = {0.};
  kfparticle->setIntermediateMinFDchi2(intermediate_min_FDchi2);*/

  kfparticle->setOutputName(outputRecoFile);
  se->registerSubsystem(kfparticle);

}
