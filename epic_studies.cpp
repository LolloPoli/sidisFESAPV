#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <set>
#include <vector>
#include <TMath.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TLorentzVector.h"
#include <cmath>
#include <random>
#include <fstream>
//#include <TROOT.h>

// lancia con: root 'esempio1.cpp("lista.txt")' -l -b -q
// lista.txt dovrà avere nella prima riga il nome del file di output (con relativo percorso se necessario) e poi i file in input
// per non avere le canvas a schermo in ambiente locale
//gROOT->SetBatch(kTRUE);

std::vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    std::vector<double> bin_edges(nbins + 1);
    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = std::pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}


void epic_studies(const char* fileList){
    std::ifstream inputFile(fileList);
    if (!inputFile.is_open()) {
        std::cerr << "Errore: impossibile aprire il file di lista " << fileList << std::endl;
        return;
    }
    std::string outputFile;
    std::vector<std::string> fileNames;
    std::string line;
    bool firstLine = true;
    while (std::getline(inputFile, line)) {
        if (firstLine) {
            outputFile = line;  // Prima riga = file di output
            firstLine = false;
        } else {
            fileNames.push_back(line);  // Altre righe = file di input
        }
    }
    inputFile.close();
    // Creiamo la TChain e aggiungiamo i file
    TChain *mychain = new TChain("events");
    for (const auto& file : fileNames) {
        mychain->Add(file.c_str());
    }
    // Controllo se la catena ha eventi
    if (mychain->GetEntries() == 0) {
        std::cerr << "Errore: Nessun evento trovato nella catena!" << std::endl;
        return;
    }
    // Apri il file di output
    TFile *ofile = TFile::Open(outputFile.c_str(), "RECREATE");

    TTreeReader tree_reader(mychain);
    if (mychain->GetEntries() == 0) {
        std::cerr << "Errore: Nessun evento nella chain!" << std::endl;
        return;
    }
    
    // Get Particle Information
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<int> parentsIndex(tree_reader, "_MCParticles_parents.index");
    TTreeReaderArray<int> daughterIndex(tree_reader, "_MCParticles_daughters.index");
    TTreeReaderArray<unsigned int> par(tree_reader, "MCParticles.parents_end");
    //
    TTreeReaderArray<float> dRICHx(tree_reader, "_DRICHAerogelTracks_points.position.x");
    TTreeReaderArray<float> dRICHy(tree_reader, "_DRICHAerogelTracks_points.position.y");
    TTreeReaderArray<float> dRICHz(tree_reader, "_DRICHAerogelTracks_points.position.z");
    TTreeReaderArray<float> dRICH_momX(tree_reader, "_DRICHAerogelTracks_points.momentum.x");
    TTreeReaderArray<float> dRICH_momY(tree_reader, "_DRICHAerogelTracks_points.momentum.y");
    TTreeReaderArray<float> dRICH_momZ(tree_reader, "_DRICHAerogelTracks_points.momentum.z");
    TTreeReaderArray<float> dRICH_theta(tree_reader, "_DRICHAerogelTracks_points.theta");
    TTreeReaderArray<float> dRICH_phi(tree_reader, "_DRICHAerogelTracks_points.phi");


    // Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedRealPIDParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedRealPIDParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedRealPIDParticles.momentum.z");
    TTreeReaderArray<int> recPdg(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> goodnessOfPID(tree_reader, "ReconstructedChargedRealPIDParticles.goodnessOfPID");
    TTreeReaderArray<float> chi2(tree_reader, "CentralCKFTracks.chi2");


    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> scatElAssoc(tree_reader, "MCScatteredElectronAssociations_objIdx.collectionID"); // trovato questo ma non so come funzioni ancora
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

    // GRAFICI 
    const int nbins = 120;
    const double xmin_xB = 1e-4, xmax_xB = 1;
    const double xmin_Q2 = 1, xmax_Q2 = 100.;
    const double xmin_a = 1, xmax_a = 10000.;
    auto make_bins = [](int bins, double min, double max) {
        return CreateLogBinning(bins, min, max);
      };
    const auto log_bins_Q2 = make_bins(nbins, xmin_Q2, xmax_Q2);
    const auto log_bins_xB = make_bins(nbins, xmin_xB, xmax_xB);
    const auto log_bins_pr = make_bins(nbins, xmin_a, xmax_a);
    // DIstributions
    
    // _______________________________________________________________________________________________________________________________________________



    double hadron_mom_mc, hadron_Q2_mc, hadron_xB_mc, hadron_xF_mc, hadron_z_mc, hadron_PhT_mc;
    double hadron_Phi_h_mc, hadron_Phi_s_mc, hadron_Phi_lab_mc, hadron_Theta_mc, hadron_eta_mc, hadron_y_mc, hadron_W_mc, hadron_Mx_mc;
    double hel_mc, eps_mc, hadron_px_mc, hadron_py_mc, hadron_pz_mc;
    double el_px_mc, el_py_mc, el_pz_mc, el_theta_mc, el_phi_mc, el_eta_mc, el_mom_mc;
    double pr_mom_mc, pr_px_mc, pr_py_mc, pr_pz_mc, pr_phi_mc, pr_theta_mc, pr_eta_mc;
    double hadron_mom, hadron_Q2, hadron_xB, hadron_xF, hadron_z, hadron_PhT;
    double hadron_Phi_h, hadron_Phi_s, hadron_Phi_lab, hadron_Theta, hadron_eta, hadron_y, hadron_W, hadron_Mx;
    double hel, eps, hadron_px, hadron_py, hadron_pz, hadron_pdg_mc;
    double el_px, el_py, el_pz, el_theta, el_phi, el_eta, el_mom;
    double pr_mom, pr_px, pr_py, pr_pz, pr_phi, pr_theta, pr_eta;
    double hadron_goodPID, hadron_pdg, hadron_index;
    double el_y_mc;
    double el_pdg;
    int hadron_count;

    TTree ElectronTreeMC("ElectronTree_MC", "");
    ElectronTreeMC.Branch("el_px_mc", &el_px_mc, "el_px_mc/D");
    ElectronTreeMC.Branch("el_py_mc", &el_py_mc, "el_py_mc/D");
    ElectronTreeMC.Branch("el_pz_mc", &el_pz_mc, "el_pz_mc/D");
    ElectronTreeMC.Branch("el_mom_mc", &el_mom_mc, "el_mom_mc/D");
    ElectronTreeMC.Branch("el_y_mc", &el_y_mc, "el_y_mc/D");

    TTree HadronTreeMC("HadronTree_MC", "");
    HadronTreeMC.Branch("hadron_px_mc", &hadron_px_mc, "hadron_px_mc/D");
    HadronTreeMC.Branch("hadron_py_mc", &hadron_py_mc, "hadron_py_mc/D");
    HadronTreeMC.Branch("hadron_pz_mc", &hadron_pz_mc, "hadron_pz_mc/D");

    TTree ElectronTreeRECO("ElectronTree_RECO", "");
    ElectronTreeRECO.Branch("el_px", &el_px, "el_px/D");
    ElectronTreeRECO.Branch("el_py", &el_py, "el_py/D");
    ElectronTreeRECO.Branch("el_pz", &el_pz, "el_pz/D");
    ElectronTreeRECO.Branch("el_mom", &el_mom, "el_mom/D");
    ElectronTreeRECO.Branch("el_pdg", &el_pdg, "el_pdg/D");

    TTree HadronTreeRECO("HadronTree_RECO", "");
    HadronTreeRECO.Branch("hadron_px", &hadron_px, "hadron_px/D");
    HadronTreeRECO.Branch("hadron_py", &hadron_py, "hadron_py/D");
    HadronTreeRECO.Branch("hadron_pz", &hadron_pz, "hadron_pz/D");
    HadronTreeRECO.Branch("hadron_pdg", &hadron_pdg, "hadron_pdg/D");
    HadronTreeRECO.Branch("hadron_y", &hadron_y, "hadron_y/D");

    
    double count3 = 0;
    double event = -1;

    std::ofstream csv("output.csv");
    csv << "event , status, parents, reco_pdg, true_pdg, px_reco, py_reco, pz_reco, mom_reco, px_true, py_true, pz_true, mom_true\n";  

    // __________________________________________________________________________________________________________________________

    while(tree_reader.Next()) { // Loop over events
        // ALCUNI VETTORI UTILI
        std::vector<TVector3> recScatElectron;
        TVector3 ElBeam(0.,0.,-18.); 
        const double m_p = 0.9382720813; // GeV
        const double m_e = 0.00051099895; // GeV
        double pz_e = -18.0;
        double p_e = std::abs(pz_e);
        double E_e = sqrt(p_e*p_e + m_e*m_e);
        double pz_p = 275.0;
        double p_p = std::abs(pz_p);
        double E_p = sqrt(p_p*p_p + m_p*m_p);
        TLorentzVector ElectronScattered;
        TLorentzVector ElectronScattered_mc;
        TLorentzVector ProtonScattered;
        TLorentzVector ProtonScattered_mc;
        //TLorentzVector ElectronBeam(0, 0, -18, 18);
        //TLorentzVector ProtonBeam(0, 0, 275,275);
        TLorentzVector ElectronBeam(0., 0., pz_e, E_e);
        TLorentzVector ProtonBeam(0., 0., pz_p, E_p);
        TLorentzVector Lab = ElectronBeam + ProtonBeam;
        TLorentzVector MC_ProtonBeam;
        TLorentzVector MC_ElectronBeam;
        double currentPhi = 0;
        double currentMom = 0;
        // per la costruzione del pione
        TVector3 ScatElectron_mc;
        double count2 = 0;
        count3++;
        double count = 0;
        //if(event >= 19) continue;
        //std::cout << event << std::endl;
        for(unsigned int i=0; i<partGenStat.GetSize(); i++){ // Loop over thrown particles
            int pdg = (partPdg[i]);
            int pdg2 = (std::abs(partPdg[i]));
            //int reco_pdg = recPdg[i];
            // status = 4 is the beam (ref 1767) in HepMC
            // status = 1 stable particle in the final state
            if(partGenStat[i] != -1){    // stable
                if(parentsIndex[i] > -1){
                    if (partGenStat[i] == 4 && pdg == 2212){
                        TVector3 Mom(partMomX[i],partMomY[i],partMomZ[i]);
                        double E = sqrt(0.9382720813*0.9382720813 + partMomX[i]*partMomX[i] + partMomY[i]*partMomY[i] + partMomZ[i]*partMomZ[i]);
                        MC_ProtonBeam.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], E);
                        event++;
                        //csv << "----------------------------- event " << event << " ----------------------------" << endl;
                    }
                    if(partGenStat[i] == 4 && pdg == 11){
                        TVector3 Mom(partMomX[i],partMomY[i],partMomZ[i]);
                        double E = sqrt(0.00051099895*0.00051099895 + partMomX[i]*partMomX[i] + partMomY[i]*partMomY[i] + partMomZ[i]*partMomZ[i]);
                        MC_ElectronBeam.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], E);
                    }
                    if(pdg2 == 11 || pdg2 == 211 || pdg2 == 321 || pdg2 == 2212){
                        for(unsigned int j=0; j<recPdg.GetSize(); j++){
                            if(simuAssoc[j] == i){
                                TVector3 ElMom(partMomX[i],partMomY[i],partMomZ[i]);
                                TVector3 MomReco(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                                csv << event << ",  " << partGenStat[i] << ",  "  << parentsIndex[i] << ",  " << recPdg[j] << ",  " << pdg << ",  " << MomReco.X() << ",  " << MomReco.Y() << ",  " << MomReco.Z()
                                << ",  " << MomReco.Mag() << ",  " << ElMom.X() << ",  " << ElMom.Y()  << ",  " << ElMom.Z() << ",  " << ElMom.Mag()  << "\n";
                                count++;
                            }
                        }
                    }
                }
                if(pdg == 11 || pdg2 == 211 || pdg2 == 321 || pdg2 == 2212){
                    count2++;
                    if(parentsIndex[i] > -1){
                        for(unsigned int j=0; j<recPdg.GetSize(); j++){
                            if(simuAssoc[j] == i){ // Find association index matching the index of the thrown particle we are looking at
                                if (pdg == 11 ){
                                    TVector3 ElMom(partMomX[i],partMomY[i],partMomZ[i]);
                                    el_eta_mc = ElMom.PseudoRapidity();
                                    double angleR = ElMom.Theta();
                                    double angle = angleR * (180.0 / TMath::Pi());
                                    el_px_mc = ElMom.X();
                                    el_py_mc = ElMom.Y();
                                    el_pz_mc = ElMom.Z();
                                    el_theta_mc = ElMom.Theta();
                                    el_phi_mc = ElMom.Phi();
                                    el_mom_mc = ElMom.Mag();
                                    //if (el_mom_mc < 1.5) continue; // to skip the background electron
                                    currentPhi = ElMom.Theta();
                                    currentMom = ElMom.Mag();
                                    el_pdg = recPdg[j];
                                    ScatElectron_mc.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                                    ElectronScattered_mc.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], el_mom_mc);
                                    TLorentzVector el_q_mc = MC_ElectronBeam - ElectronScattered_mc;
                                    el_y_mc = (MC_ProtonBeam.Dot(el_q_mc))/(MC_ProtonBeam.Dot(MC_ElectronBeam));
                                    if(el_eta_mc >= -3.5 && el_eta_mc <= 3.5 && el_y_mc >= 0.01 && el_y_mc <= 0.99){
                                        int recpdg = (recPdg[j]);
                                        TVector3 recElmom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                                        double eta = recElmom.PseudoRapidity();
                                        el_px = recElmom.X(), el_py = recElmom.Y(), el_pz = recElmom.Z();
                                        double momE = recElmom.Mag();
                                        el_mom = momE;
                                        double theta = recElmom.Theta();
                                        double phi = recElmom.Phi();
                                        //if (el_y_mc < 0.01 || el_y_mc > 0.99) std::cout << eta_mc << std::endl;
                                        //
                                        ElectronTreeMC.Fill();
                                        ElectronTreeRECO.Fill();
                                    }
                                } else if (partPdg[i] == recPdg[j]){
                                    TVector3 particle(partMomX[i],partMomY[i],partMomZ[i]);
                                    TVector3 momentum_vect(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                                    double eta_particle = particle.PseudoRapidity();
                                    double eta = momentum_vect.PseudoRapidity();
                                    if(eta >= -3.5 && eta <= 3.5 && el_y_mc >= 0.01 && el_y_mc <= 0.99){
                                        hadron_px_mc = particle.X(), hadron_py_mc = particle.Y(), hadron_pz_mc = particle.Z(); 
                                        hadron_pdg = pdg; 
                                        int recpdg = (recPdg[j]);
                                        double goodPID_cont = goodnessOfPID[j];
                                        double mom = particle.Mag();
                                        TLorentzVector photon_hadron_noBoost = ElectronBeam - ElectronScattered_mc;
                                        double E_hadron_mc = sqrt(mom*mom + 0.13957*0.13957);
                                        TLorentzVector hadron_noBoost(particle.X(), particle.Y(), particle.Z(), E_hadron_mc);
                                        hadron_y = el_y_mc;
                                        hadron_px = particle.X(); hadron_py = particle.Y(); hadron_pz = particle.Z();
                                        hadron_index++;
                                        hadron_pdg = pdg;
                                        hadron_mom = particle.Mag();
                                        hadron_Q2 = -photon_hadron_noBoost.M2();  
                                        hadron_xB = hadron_Q2 / (2 * ProtonBeam.Dot(photon_hadron_noBoost));
                                        hadron_eta = eta;
                                        hadron_Theta = particle.Theta(); hadron_Phi_lab = particle.Phi();
                                        hadron_z = (ProtonBeam * hadron_noBoost) / (ProtonBeam * photon_hadron_noBoost);
                                        hadron_goodPID = goodPID_cont;
                                        //
                                        TLorentzVector hadron_mc = hadron_noBoost;
                                        TLorentzVector photon_hadron_mc = photon_hadron_noBoost;
                                        // boost gamma*N - necessario per il calcolo di P_hT e gli angolo Phi_h e Phi_s direi
                                        TLorentzVector gammaN_mc = photon_hadron_noBoost + ProtonBeam;
                                        TVector3 boost_gammaN_mc = -gammaN_mc.BoostVector();
                                        hadron_mc.Boost(boost_gammaN_mc);
                                        photon_hadron_mc.Boost(boost_gammaN_mc);
                                        TVector3 zAxis_mc = photon_hadron_mc.Vect().Unit();
                                        hadron_PhT = hadron_mc.Perp(zAxis_mc);
                                        // 
                                        HadronTreeMC.Fill();
                                        HadronTreeRECO.Fill();
                                    } 
                                }
                            }
                        }
                    }    
                }
            } 
        }
    } 


    //ElectronTreeMC.Write();
    //HadronTreeMC.Write();
    //ElectronTreeRECO.Write();
    //HadronTreeRECO.Write();
      
    std::cout << " " << std::endl;
    std::cout << "ROOT output file: " << outputFile << endl;
    std::cout << "______________________________________________________________________________________" << std::endl;
    std::cout << " " << std::endl;

    
    ofile->Write(); // Write histograms 
    ofile->Close(); 

    mychain->Delete();
    ofile->Delete();
  }
