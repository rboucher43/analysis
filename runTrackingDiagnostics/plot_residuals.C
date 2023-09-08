#include <cmath>
#include <TFile.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TCut.h>
#include <TStyle.h>
#include <Eigen/Dense>
#include <map>
#include <TLatex.h>
#include <TCanvas.h>
#include <iostream>


int getSector(float phi, float z)
{
  int angleRange[13] = {180,150,120,90,60,30,0,-30,-60,-90,-120,-150,-180};
  float phioffset    = 15;
  int sector         = 0;

  for(int i=0; i<12; ++i)
    {
      if(phi<(angleRange[i]-phioffset) && phi>(angleRange[i+1]-phioffset))
	{
	  //std::cout << "Warning: only see this once continuously" << std::endl;
	  if(z>0)
	    {
	      sector = i;
	      if(abs(phi) > 165){sector = 11;}
	    }
	  if(z<0)
	    {
	      sector = i+12;
	      if(abs(phi) > 165){sector=23;}
	    }
	}
    }
      
  return sector;

}

void plot_residuals(const std::string inputFile = "aligned_99.root", const std::string outputFile = "residuals_plots.root")
{
  std::cout << "Initializing " << std::endl;

  bool plotBool       = true;  // Determines whether plots are displayed or not
  bool plotSector     = true;
  bool savePlot       = false;  // Determines whether histograms are saved to output file or not
  bool sectorVariance = false;
  float ymax          = 0.1;  // y limits for mean and variance sector plots

  gStyle->SetOptStat(1);
  gStyle->SetTitleSize(0.1,"t");

  TFile fin(inputFile.c_str());

  TNtuple *ntuple;
  fin.GetObject("ntp_residuals",ntuple);

  float seed_id;
  float layer;
  float dphi;
  float dz;
  float x;
  float y;
  float z;
  float pt;
  float px;
  float py;
  float pz;
  float crossing;
  float isSilicon;
  float isTpc;

  ntuple->SetBranchAddress("seed_id",&seed_id);
  ntuple->SetBranchAddress("layer",&layer);
  ntuple->SetBranchAddress("dphi",&dphi);
  ntuple->SetBranchAddress("dz",&dz);
  ntuple->SetBranchAddress("x",&x);
  ntuple->SetBranchAddress("y",&y);
  ntuple->SetBranchAddress("z",&z);
  ntuple->SetBranchAddress("pt",&pt);
  ntuple->SetBranchAddress("px",&px);
  ntuple->SetBranchAddress("py",&py);
  ntuple->SetBranchAddress("pz",&pz);
  ntuple->SetBranchAddress("crossing",&crossing);
  ntuple->SetBranchAddress("isSilicon",&isSilicon);
  ntuple->SetBranchAddress("isTpc",&isTpc);

  int entries = ntuple->GetEntries();
 
  // initializing histogram for non-sector plots
  TH2D *xy            = new TH2D("xy","cluster x vs. y",5000,-80,80,5000,-80,80);
  TH2D *tpc_xy        = new TH2D("tpc_xy","tpc cluster x vs. y",5000,-80,80,5000,-80,80);
  TH2D *si_xy         = new TH2D("si_xy","si cluster x vs. y",5000,-10,10,5000,-10,10);
  TH2D *dphiLayer     = new TH2D("dphiLayer","dphi vs. layer",5000,-0.02,0.02,56,0,56);
  TH2D *dzLayer       = new TH2D("dzLayer","dz vs. layer",5000,-1,1,5000,0,56);
  TH2D *resLayer      = new TH2D("resLayer","residual vs. layer",5000,0,1,5000,0,56);
  TH2D *dphiPt        = new TH2D("dphiPt","dphi vs. pt",5000,-0.02,0.02,5000,0,7);
  TH2D *dphiPz        = new TH2D("dphiPz","dphi vs. pz",5000,-0.02,0.02,5000,0,5);
  TH2D *rmsdphiLayer  = new TH2D("rmsdphiLayer","rms dphi vs. layer",5000,-0.02,0.02,5000,0,56);
  TH2D *rmsdphiNhits  = new TH2D("rmsdphiNhits","rms dphi vs. nHits",5000,0.0,0.02,5000,0,350);
  TH2D *meanPhi       = new TH2D("meanPhi","Average dphi vs. sector",5000,-10.0,10.0,5000,-10,10);
  TH2D *rdphiLayer    = new TH2D("rdphiLayer","Rdphi vs. Layer",5000,-1.0,1.0,5000,0,56);


  std::multimap<int,float> nhitsMap;
  std::multimap<int,float> layerMap;
  std::multimap<int,float> sector0Map;
  std::multimap<int,float> sector1Map;
  std::multimap<int,float> sector2Map;
  std::multimap<int,float> sector3Map;
  std::multimap<int,float> sector4Map;
  std::multimap<int,float> sector5Map;
  std::multimap<int,float> sector6Map;
  std::multimap<int,float> sector7Map;
  std::multimap<int,float> sector8Map;
  std::multimap<int,float> sector9Map;
  std::multimap<int,float> sector10Map;
  std::multimap<int,float> sector11Map;
  std::multimap<int,float> sector12Map;
  std::multimap<int,float> sector13Map;
  std::multimap<int,float> sector14Map;
  std::multimap<int,float> sector15Map;
  std::multimap<int,float> sector16Map;
  std::multimap<int,float> sector17Map;
  std::multimap<int,float> sector18Map;
  std::multimap<int,float> sector19Map;
  std::multimap<int,float> sector20Map;
  std::multimap<int,float> sector21Map;
  std::multimap<int,float> sector22Map;
  std::multimap<int,float> sector23Map;


  std::multimap<int,float> sectormaps[24] = {sector0Map,sector1Map,sector2Map,sector3Map,sector4Map,sector5Map,sector6Map,sector7Map,sector8Map,sector9Map,sector10Map,sector11Map,sector12Map,sector13Map,sector14Map,sector15Map,sector16Map,sector17Map,sector18Map,sector19Map,sector20Map,sector21Map,sector22Map,sector23Map};
  
  TH1D *meanhists[24];
  TH1D *varhists[24];
  TLatex *lt[24];
  
  for(int i=0; i<24; ++i) //initialize histograms for mean and variance plots
    {
      char hname[100];
      sprintf(hname,"mean%i",i);
      char htitle[100];
      if(i<12){sprintf(htitle,"Mean North Sector %i",i);}
      else{sprintf(htitle,"Mean South Sector %i",i);}
      meanhists[i] = new TH1D(hname, htitle, 48,7,54);
      meanhists[i]->SetMaximum(ymax);
      meanhists[i]->SetMinimum(-ymax);
      meanhists[i]->GetYaxis()->SetLabelSize(0.06);
      meanhists[i]->GetXaxis()->SetLabelSize(0.06);

      char varname[100];
      sprintf(varname,"var%i",i);
      char vtitle[100];
      if(i<12){sprintf(vtitle,"Variance North Sector %i",i);}
      else{sprintf(vtitle,"Variance South Sector %i",i);}
      varhists[i] = new TH1D(varname, vtitle, 48,7,54);
      varhists[i]->SetMaximum(1.0);
      varhists[i]->SetMinimum(0);
      varhists[i]->GetYaxis()->SetLabelSize(0.06);
      varhists[i]->GetXaxis()->SetLabelSize(0.06);
    }

  std::cout << "Filling Maps and Histograms" << std::endl;
  for(int i=0; i<entries; ++i)
    {
      ntuple->GetEntry(i);
     
      float sectorPhi = atan2(y,x)*180/M_PI;
      if(isnan(sectorPhi))
	{
	  //std::cout << "WARNING: instance of sector phi angle as NaN in plot_residuals.C line 174 "<<std::endl;
	  continue;
	}
      //if(sectorPhi<0){std::cout << "negative phi"<< sectorPhi<<std::endl;}
      if(isnan(dphi))
	{
	  std::cout << "WARNING: instance of residual phi as NaN plot_residuals.C line 179"<< std::endl;
	  std::cout << "sectorphi " <<sectorPhi<<std::endl;
	  continue;
	}

      float r = sqrt(x*x+y*y);
      //std::cout << " radius "<< r<<std::endl;
      int sector = getSector(sectorPhi,z); //find sector of entry
      sectormaps[sector].insert({layer,r*dphi}); 
      //std::cout<<" dphi "<< dphi<<std::endl;
      //std::cout<<" r*dphi "<< r*dphi<<std::endl;
      //if(!isSilicon)
      //	{
      rdphiLayer->Fill(r*dphi,layer);
      //}
      if(plotBool)
	{
	  nhitsMap.insert({seed_id,dphi*dphi});     
	  layerMap.insert({layer,dphi*dphi});
	  float r        = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	  float residual = sqrt(pow(r*dphi,2) + pow(dz,2));

	  resLayer->Fill(residual,layer);
	  dphiLayer->Fill(dphi,layer);
	  dzLayer->Fill(dz,layer);
	  xy->Fill(x,y);
	  dphiPt->Fill(dphi,pt);
	  dphiPz->Fill(dphi,pz);

	  if(isSilicon == 0){ tpc_xy->Fill(x,y); }
	  if(isTpc == 0){ si_xy->Fill(x,y); }
	}
    } 

  if(plotBool)
    {      
      //plot nhits v rmsdphi  
      for(multimap<int,float>::iterator it = nhitsMap.begin(), end = nhitsMap.end(); it != end; it = nhitsMap.upper_bound(it->first))
	{
	  int nhits             = nhitsMap.count(it->first);
	  auto second_it        = nhitsMap.equal_range(it->first);
	  double dphisquaredsum = 0.0;
	  for (auto itr = second_it.first; itr != second_it.second; ++itr) // loops over all entries with the same key it->first 
	    {
	      dphisquaredsum += itr->second;
	    }
	  float rmsdphi = sqrt(dphisquaredsum/nhitsMap.count(it->first));
	  rmsdphiNhits->Fill(rmsdphi,it->first);
      
	}
      //sector mapplots rmsdphi v layer
      std::cout<<"intializing layermap"<<std::endl;
      for(multimap<int,float>::iterator it = layerMap.begin(), end = layerMap.end(); it != end; it = layerMap.upper_bound(it->first))
	{      
	  auto second_it        = layerMap.equal_range(it->first);
	  double dphisquaredsum = 0.0;
	  for (auto itr = second_it.first; itr != second_it.second; ++itr) // loops over all entries with the same key it->first 
	    {
	      dphisquaredsum += itr->second;
	    }
	  float rmsdphi = sqrt(dphisquaredsum/layerMap.count(it->first));
	  rmsdphiLayer->Fill(rmsdphi,it->first);
	}
    }


  //make loop over layers instead of multimap
  std::cout << "Creating sector plots..."<< std::endl;
  for(int i = 0; i < 24; ++i) // loop over all 24 sectors and create plots
    {
      std::cout<<"sector "<< i << std::endl;

      float meanSum = 0;
      float n       = 0;

      for(int layer=7; layer<54; ++layer)
	{
	  if(layer< 7 or layer>54){continue;}

	  auto second_it = sectormaps[i].equal_range(layer);

	  double dphiSum = 0.0;
	  for (auto itr = second_it.first; itr != second_it.second; ++itr) // loops over all entries with the same key(layer) it->first 
	    {
	      if(isnan(itr->second)){std::cout <<"itr-second NaN "<< " itr->first "<<itr->first<<" sector "<<i<<std::endl;}
	      if(isnan(itr->second)){continue;} //undertsand why some residuals from itr->second are NaN 

	      dphiSum += itr->second; 
	    }
	    
	  if(sectormaps[i].count(layer) != 0)
	    {

	      double meanphi = dphiSum/sectormaps[i].count(layer);

	      // case of layerphi with NaN
	      // if(layer==53 && sectormaps[i].count(layer)==2102)
	      // 	{
	      // 	  std::cout << "SPECIAL CASE mean: " << meanphi << " dphisum "<<dphiSum<< " layer " <<layer<<std::endl; 
	      // 	}
	      meanhists[i]->SetBinContent(layer-6,meanphi);
	      meanSum += meanphi;
	      n++;
	      
	      //Calculate variance
	      double varianceSum = 0;
	      for (auto itr = second_it.first; itr != second_it.second; ++itr) // loops over all entries with the same key(layer) it->first 
		{
		  //varianceSum += pow((itr->second - meanphi),2);
		  varianceSum += pow((itr->second),2);
		  //if(itr->second>0.5){std::cout<<"Large Value "<< " dphi " << itr->second << " layer "<< layer <<" sector "<<i<<std::endl;}
		}
	      float variancephi = sqrt(varianceSum/(sectormaps[i].count(layer)));
	      if(sectormaps[i].count(layer)==1){variancephi = 0;}
	      if(!isnan(variancephi)){varhists[i]->SetBinContent(layer-6,variancephi);}
	      if(isnan(variancephi)){std::cout<<"Warning: Instance of variancephi as NaN"<<std::endl<<"sector "<<i<<" layer "<<layer <<" variance sum "<<varianceSum<<" layer_count "<<sectormaps[i].count(layer)<<std::endl;} 
	    }
	}

      float averageMean = meanSum/n; // add this to meanhists[i]
      //std::cout << "average mean: " << averageMean << " sector " << i  << std::endl;
      char meanchar[100];
      sprintf(meanchar,"mean: %.2e",averageMean);
      lt[i] = new TLatex(0.25,0.8,meanchar);
      lt[i]->SetNDC(1);
      lt[i]->SetTextSize(0.065);
    }

  std::cout<<"plotting"<<std::endl;
  if(plotBool)
    {
      TCanvas *c1 = new TCanvas("c1","",10,10,800,800);
      xy->DrawCopy();
      TCanvas *c2 = new TCanvas("c2","",10,10,800,800);
      tpc_xy->DrawCopy();
      TCanvas *c3 = new TCanvas("c3","",10,10,600,600); 
      si_xy->DrawCopy();
      TCanvas *c4 = new TCanvas("c4","",10,10,800,800);
      dphiLayer->DrawCopy();
      TCanvas *c5 = new TCanvas("c5","",10,10,800,800); 
      dzLayer->DrawCopy();
      TCanvas *c6 = new TCanvas("c6","",10,10,800,800); 
      resLayer->DrawCopy();
      TCanvas *c7 = new TCanvas("c7","",10,10,800,800); 
      dphiPt->DrawCopy();
      TCanvas *c8 = new TCanvas("c8","",10,10,800,800); 
      dphiPz->DrawCopy();
      TCanvas *c9 = new TCanvas("c9","",10,10,800,800); 
      rmsdphiLayer->DrawCopy();
      TCanvas *c10 = new TCanvas("c10","",10,10,800,800);
      rmsdphiNhits->DrawCopy();
    }

  TCanvas *rp = new TCanvas("rp","",10,10,800,800);
  rdphiLayer->DrawCopy();

  TCanvas *project = new TCanvas("project","",10,10,800,800);
  int binlow  = rdphiLayer->GetYaxis()->FindBin(2.5);
  int binhigh = rdphiLayer->GetYaxis()->FindBin(6.5);

  rdphiLayer->ProjectionX("px",binlow,binhigh)->DrawCopy();


  if(plotSector)
    {
      TCanvas *canvasArr[24];

      int count = 0;
      for(int i=0; i<6; ++i)
	{
	  std::cout << "saving plots" <<std::endl;
  
	  char canvasname[100];
	  sprintf(canvasname,"cArr%i",i);
	  canvasArr[i] = new TCanvas(canvasname,canvasname,10,10,1600,400); //last two arguments are pixels in x,y try 1600 400 for rectangle
	  canvasArr[i]->Divide(4,2);
	  int cdCount = 0;
	  for(int j=0; j<4; j++)
	    {
	      cdCount++;
	      canvasArr[i]->cd(cdCount);
	      gPad->SetLeftMargin(0.15);
	      char sectorChar[100];
	      sprintf(sectorChar,"mean%count",count);
	      meanhists[count]->SetMarkerStyle(20);
	      meanhists[count]->SetMarkerSize(0.5);
	      meanhists[count]->DrawCopy("P");
	      lt[count]->Draw();
	      canvasArr[i]->cd(cdCount+4);
	      gPad->SetLeftMargin(0.15);

	      char varChar[100];
	      sprintf(varChar,"var%count",count);
	      varhists[count]->SetMarkerStyle(20);
	      varhists[count]->SetMarkerSize(0.5);
	      varhists[count]->DrawCopy("P");
	      count++;
	    }
	}
      
    }  
  


  if(savePlot && plotBool)
    {
      std::cout << "saving plots" <<std::endl;
      std::unique_ptr<TFile> residualsFile(TFile::Open(outputFile.c_str(), "recreate"));
      residualsFile->WriteObject(xy,"cluster_xvy");
      residualsFile->WriteObject(tpc_xy,"tpc_cluster_xvy ");
      residualsFile->WriteObject(si_xy,"si_cluster_xvy");
      residualsFile->WriteObject(dphiLayer,"dphi_v_layer");
      residualsFile->WriteObject(dzLayer,"dz_v_layer");
      residualsFile->WriteObject(resLayer,"residual_v_layer");
      residualsFile->WriteObject(dphiPt,"dphi_v_pt");
      residualsFile->WriteObject(dphiPz,"dphi_v_pz");
      residualsFile->WriteObject(rmsdphiNhits,"rmsdphi_v_Nhits");
      //residualsFile->WriteObject(mean11,"mean1dphi_v_layer");
    }

 if(savePlot && plotSector)
    {
      std::cout << "saving plots" <<std::endl;

    }



}



