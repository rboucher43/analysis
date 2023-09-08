
void plot_D0()
{

  char filename[200];

  int file_number = 500;

  float D0_mass;
  float D0_pT;
  float D0_phi;
  float D0_theta;
  float D0_rapidity;
  float DIRA;
  float D0_decayLength;
  int event_number;

  TH1D *file_hist = new TH1D("file_hist","",100000,0,100000);


  TFile *D0_reco = new TFile("D0_reco.root","recreate");

  TNtuple *ntp_D0_reco = new TNtuple("ntp_D0_reco", "D0_pairs", "file_number:event_number:D01_mass:D01_pT:D01_phi:D01_theta:D01_rapidity:D01_DIRA:D01_decayLength:D02_mass:D02_pT:D02_phi:D02_theta:D02_rapidity:D02_DIRA:D02_decayLength");




  for(int i=0; i<file_number; ++i)
    {
      sprintf(filename,"charm_output_%i.root",i);
      TChain *tree = new TChain("DecayTree"); 
      tree->Add(filename);
      tree->SetBranchAddress("D0_mass",&D0_mass);
      tree->SetBranchAddress("D0_pT",&D0_pT);
      tree->SetBranchAddress("D0_phi",&D0_phi);
      tree->SetBranchAddress("D0_theta",&D0_theta);
      tree->SetBranchAddress("D0_rapidity",&D0_rapidity);
      tree->SetBranchAddress("D0_DIRA",&DIRA);
      tree->SetBranchAddress("D0_decayLength",&D0_decayLength);
      tree->SetBranchAddress("eventNumber",&event_number);
      
      int entries = tree->GetEntries();
      //std::cout << "entries: "<< entries<<std::endl;
      //std::cout <<"entering loop" << std::endl;

      std::multimap<int, int> bin_map;
      std::set<int> entry_set; 

      for(int j=1; j<entries; ++j)
	{
	  tree->GetEntry(j);
	  int index = i*200 + event_number;
	  
	  if(DIRA<0.95){continue;}
	  if(D0_decayLength<0.02){continue;}
	  if(D0_mass<1.83){continue;}
	  if(D0_mass>1.91){continue;}
	  //std::cout << "D0_mass " << D0_mass << std::endl;
	  bin_map.insert(std::make_pair(event_number, j));
	  file_hist->Fill(index);
	  //std::cout << "event number "<< event_number<<std::endl;
	  //std::cout << "index        "<< index<<std::endl;
	}
      
      for(int event=0; event<200;event++)
	{
	  if(bin_map.count(event)>1) //check 161 later in the last file
	    {

	      float data[16]; 
	      auto pair_map = bin_map.equal_range(event);
	      int count =0;
	      for(auto itr = pair_map.first; itr != pair_map.second; ++itr)
		{
		  count++;
		  std::cout << "Key: " << itr->first << ", Value: " << itr->second << std::endl;
	  	
		  tree->GetEntry(itr->second);
		  std::cout << "event number: "<<event_number << std::endl;
		  std::cout << "D0_mass  " <<D0_mass << std::endl;
		  std::cout << "D0_theta " <<D0_theta << std::endl;
		  std::cout << "D0_phi   " <<D0_phi << std::endl;

		  //TNtuple ntp_D0_reco = new TNtuple("ntp_D0_reco", "D0_pairs", "file_number:event_number:D01_mass:D01_pT:D01_phi:D01_theta:D01_rapidity:D01_DIRA:D01_decayLength:D02_mass:D02_pT:D02_phi:D02_theta:D02_rapidity:D02_DIRA:D02_decayLength");

		  if(count==1)
		    {
		      data[0] = i;
		      data[1] = event_number;
		      data[2] = D0_mass;
		      data[3] = D0_pT;
		      data[4] = D0_phi;
		      data[5] = D0_theta;
		      data[6] = D0_rapidity;
		      data[7] = DIRA;
		      data[8] = D0_decayLength;
		    }
		  if(count==2)
		    {
		      data[9]  = D0_mass;
		      data[10] = D0_pT;
		      data[11] = D0_phi;
		      data[12] = D0_theta;
		      data[13] = D0_rapidity;
		      data[14] = DIRA;
		      data[15] = D0_decayLength;
		    }
		  
		 
		  //fill ntuple here wiht event numebr and all stuff for d01 an d02 
		  
		}
	      ntp_D0_reco->Fill(data);

	    }
	}
    
       
    }

  D0_reco->cd();
  ntp_D0_reco->Write();
  D0_reco->Close();
  file_hist->DrawCopy();

}
