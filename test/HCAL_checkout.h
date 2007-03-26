#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"
#include "TCut.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TF1.h"
#include "TStyle.h"
#include "Rtypes.h"
#include "TText.h"
#include "TLine.h"

double et2e(int eta);

void initStyle(TStyle *sty);

void SetupTowerDisplay(TH2F *hist);

int check_out(char* input_file, char* output_file);
