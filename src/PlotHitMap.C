#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>

// Variables to hold branch data
vector<unsigned int>   *m_apvID = nullptr;
vector<unsigned int>   *m_apvCH = nullptr;
vector<unsigned int>   *m_mmStrip = nullptr;
vector<short>   *m_apvQMax = nullptr;
vector<short>   *m_apvTBQmax = nullptr;

TTree *tree_raw = NULL;
TTree *tree_data = NULL;

void OpenFile(std::string fname) {
  // Open the ROOT file
  TFile *file = TFile::Open(fname.c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error opening file!" << std::endl;
      return;
  }

  // Get the TTree from the file
  tree_raw = (TTree*)file->Get("raw");
  if (!tree_raw) {
      std::cerr << "Error: TTree not found!" << std::endl;
      return;
  }
  // Get the TTree from the file
  tree_data = (TTree*)file->Get("data");
  if (!tree_data) {
      std::cerr << "Error: TTree not found!" << std::endl;
      return;
  }

  // Set branch addresses
  tree_raw->SetBranchAddress("apv_id", &m_apvID);
  tree_raw->SetBranchAddress("apv_ch", &m_apvCH);
  tree_raw->SetBranchAddress("mm_strip", &m_mmStrip);
  tree_data->SetBranchAddress("apv_qmax", &m_apvQMax);
  tree_data->SetBranchAddress("apv_tbqmax", &m_apvTBQmax);

}

void PlotHitMap()
{

  if (tree_raw == NULL || tree_data == NULL) {return;}
  // Make some Histograms
  TH1F* hStripsX = new TH1F("hStripsX", "hStripsX",256, 0, 255);
  TH1F* hStripsY = new TH1F("hStripsY", "hStripsY",256, 0, 255);
  TH1F* hChX = new TH1F("hChX", "hChX",256, 0, 255);
  TH1F* hChY = new TH1F("hChY", "hChY",256, 0, 255);
  TH1F* hHitX = new TH1F("hHitX", "hHitX",256, 0, 255);
  TH1F* hHitY = new TH1F("hHitY", "hHitY",256, 0, 255);
  TH2F* hHitMap = new TH2F("hHitMap", "hHitMap",256, 0, 255, 256, 0, 255);
  
  // Loop over entries
  Long64_t nentries = tree_raw->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree_raw->GetEntry(i);
    tree_data->GetEntry(i);
    
    float chargeCenterX = 0;
    int cntX = 0;
    float chargeCenterY = 0;
    int cntY = 0;
    for (unsigned int j = 0; j < m_apvID->size(); j++){
      if (m_apvID->at(j) == 0 || m_apvID->at(j) == 1){
        hStripsX->Fill(m_mmStrip->at(j));
        if (m_apvID->at(j) == 0){
          hChX->Fill(m_apvCH->at(j));
          chargeCenterX += m_apvCH->at(j);
          cntX ++;
        } else {
          hChX->Fill(m_apvCH->at(j)+128);
          chargeCenterX += m_apvCH->at(j)+128;
          cntX ++;
        }
      } else {
        hStripsY->Fill(m_mmStrip->at(j));
        if (m_apvID->at(j) == 3){
          hChY->Fill(m_apvCH->at(j));
          chargeCenterY += m_apvCH->at(j);
          cntY ++;
        } else {
          hChY->Fill(m_apvCH->at(j)+128);
          chargeCenterY += m_apvCH->at(j)+128;
          cntY ++;
        }
      }
    }
    chargeCenterX /= cntX;
    chargeCenterY /= cntY;
    hHitX->Fill(chargeCenterX);
    hHitY->Fill(chargeCenterY);
    hHitMap->Fill(chargeCenterX, chargeCenterY);
    
  }
  TCanvas *c = new TCanvas();
  hStripsX->Draw();
  c = new TCanvas();
  hStripsY->Draw();
  c = new TCanvas();
  hChX->Draw();
  c = new TCanvas();
  hChY->Draw();
  c = new TCanvas();
  hHitX->Draw();
  c = new TCanvas();
  hHitY->Draw();
  c = new TCanvas();
  hHitMap->Draw();
  
}

void PlotHitMap(std::string fName)
{
  OpenFile(fName);
  PlotHitMap();
}
