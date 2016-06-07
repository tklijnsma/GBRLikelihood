#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"

#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooLinearVar.h"

#include "../interface/RooHybridBDTAutoPdf.h"
#include "../interface/HybridGBRForest.h"
#include "../interface/RooDoubleCBFast.h"

#include "../interface/HybridGBRForestFlex.h"

//#include "../interface/HybridGBRForestD.h"

using namespace RooFit;
 
// For debugging purposes
#include <typeinfo>
#include <cstring>


//#######################################
// effsigma function from Chris
//#######################################

Double_t effSigma(TH1 * hist){

    TAxis *xaxis = hist->GetXaxis();
    Int_t nb = xaxis->GetNbins();
    if(nb < 10) {
        cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
        return 0.;
        }

    Double_t bwid = xaxis->GetBinWidth(1);
    if(bwid == 0) {
        cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
        return 0.;
        }

    Double_t xmax = xaxis->GetXmax();
    Double_t xmin = xaxis->GetXmin();
    Double_t ave = hist->GetMean();
    Double_t rms = hist->GetRMS();

    Double_t total=0.;
    for(Int_t i=0; i<nb+2; i++) {
        total+=hist->GetBinContent(i);
        }
    //   if(total < 100.) {
    //     cout << "effsigma: Too few entries " << total << endl;
    //     return 0.;
    //   }
    Int_t ierr=0;
    Int_t ismin=999;

    Double_t rlim=0.683*total;
    Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
    if(nrms > nb/10) nrms=nb/10; // Could be tuned...

    Double_t widmin=9999999.;
    
    for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
        
        Int_t ibm=(ave-xmin)/bwid+1+iscan;
        Double_t x=(ibm-0.5)*bwid+xmin;
        Double_t xj=x;
        Double_t xk=x;
        Int_t jbm=ibm;
        Int_t kbm=ibm;
        Double_t bin=hist->GetBinContent(ibm);
        total=bin;
        
        for(Int_t j=1;j<nb;j++){
            
            if(jbm < nb) {
                jbm++;
                xj+=bwid;
                bin=hist->GetBinContent(jbm);
                total+=bin;
                if(total > rlim) break;
                }
            else ierr=1;
            
            if(kbm > 0) {
                kbm--;
                xk-=bwid;
                bin=hist->GetBinContent(kbm);
                total+=bin;
                if(total > rlim) break;
                }
            else ierr=1;

            }

        Double_t dxf=(total-rlim)*bwid/bin;
        Double_t wid=(xj-xk+bwid-dxf)*0.5;
        
        if(wid < widmin) {
            widmin=wid;
            ismin=iscan;
            }   

        }
    
    if(ismin == nrms || ismin == -nrms) ierr=3;
    if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;

    return widmin;
    }


//#######################################
// Main
//#######################################

void eregtest(
    // T: Disable input options
    //bool dobarrel, bool doele
    ) {

    // There should be a switch between barrel and endcap
    // bool dobarrel = true;
    // bool dobarrel = true;

    // TCut selcut;
    // if (dobarrel)
    //     selcut = "scIsEB";
    // else
    //     // selcut = "gen_pt>25. && !isEB";
    //     selcut = "!scIsEB";


    bool dobarrel;
    TCut selcut;

    if ( getenv("REGION")==NULL ) {
        cout << "Unclear whether endcap or barrel should be run" << endl;
        return;
        }
    else if ( (string)getenv("REGION")=="EB" ){
        selcut = "scIsEB";
        dobarrel = true;
        }
    else if ( (string)getenv("REGION")=="EE" ) {
        selcut = "!scIsEB";
        dobarrel = false;
        }

    // TString dirname = "/afs/cern.ch/work/t/tklijnsm/public/CMSSW_8_0_4/src/RegressionTraining/" ;
    TString dirname = getenv("PLOTDIR_FULLPATH") ;
    cout << dirname << endl;
    gSystem->mkdir(dirname,true);
    gSystem->cd(dirname);


    //=======================================
    // Specify what to run over

    // TString fname = "default_config_results.root";
    // // TString fname = "FullTraining_ComparableBranches_30May.root";
    // // TString fname = "Lastyears_Regression.root";
    // // TString fname = "FullTraining_30May_a1is3_forEB.root";

    // TString Ntup_path  = "/afs/cern.ch/work/t/tklijnsm/public/CMSSW_8_0_4/src/NTuples/";
    // TString Ntup_fname = "Ntup_01June_DoublePhoton.root";
    // TString Ntup       = Ntup_path + Ntup_fname;

    // TString tree_name  = "PhotonTree";

    TString Ntup      = getenv("FLATNTUPLE");
    TString tree_name = getenv("NTUPLETREE");
    TString fname     = getenv("TRAININGOUTPUT");

    cout << "In macro" << endl;
    cout << getenv("FLATNTUPLE") << endl;
    cout << getenv("NTUPLETREE") << endl;
    cout << getenv("TRAININGOUTPUT") << endl;


    //=======================================
    // Regression output - opening the workspace

    TString infile = TString::Format( fname.Data() );

    // Open the RooWorkspace
    TFile *fws = TFile::Open(infile); 

    // if (dobarrel)     
    //     RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg_eb");
    // else
    //     RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg_ee");

    RooWorkspace *ws;
    if (dobarrel)     
        ws = (RooWorkspace*)fws->Get("wereg_eb");  
    else
        ws = (RooWorkspace*)fws->Get("wereg_ee");  
    ws->Print();

    // Get Target var
    RooRealVar *tgtvar = ws->var("targetvar");
    // cout << tgtvar << endl;
    // cout << typeid( tgtvar ).name() << endl;
    // tgtvar->Print();

    // Get GBRFunction
    RooGBRFunctionFlex *func;
    if (dobarrel) 
        func = static_cast<RooGBRFunctionFlex*>( ws->arg("sigmeantfuncEB") );
    else
        func = static_cast<RooGBRFunctionFlex*>( ws->arg("sigmeantfuncEE") );

    // Some test statements
    cout << func << endl;
    cout << typeid( func ).name() << endl;
    func->Print();

    // Strange type only needed to get a list of variables from
    RooGBRTargetFlex *tgtFlex;
    if (dobarrel) 
        tgtFlex = static_cast<RooGBRTargetFlex*>( ws->arg("sigmeantEB") );
    else
        tgtFlex = static_cast<RooGBRTargetFlex*>( ws->arg("sigmeantEE") );

    // Set it in the RooArgList
    RooArgList vars;
    vars.add( tgtFlex->FuncVars() );
    vars.add( *tgtvar );


    //=======================================
    // Open the NTuple

    // Not sure what this does
    RooRealVar weightvar("weightvar","",1.);

    // TFile *fdin = TFile::Open("/afs/cern.ch/work/t/tklijnsm/public/CMSSW_8_0_4/src/NTuples/Ntup_01June_DoubleElectron.root");
    // TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("een_analyzer");
    // TTree *dtree = (TTree*)ddir->Get("ElectronTree");

    TFile *fdin = TFile::Open( Ntup );
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("een_analyzer");
    TTree *dtree = (TTree*)ddir->Get(tree_name);


    //#######################################
    // Start with cuts, count trees
    //#######################################

    // Cut on event number
    // TCut eventcut = "eventNumber%2==1 && eventNumber>400 && eventNumber<=1000";
    TCut eventcut = "eventNumber%2==1";

    if ( (string)getenv("TESTRUN")=="Y" ){
        // Use only part of sample
        TCut NtupIDcut = "NtupID>400&&NtupID<800";
        eventcut *= NtupIDcut;
        }

    weightvar.SetTitle( eventcut * selcut );

        // Other cuts possible (original)

        // //  TCut selcut = "ph.pt>25. && ph.isbarrel && ph.ispromptgen && abs(ph.sceta)<1.0"; 
        // //TCut selcut = "ph.pt>25. && ph.isbarrel && (ph.scrawe/ph.gene)>0. && (ph.scrawe/ph.gene)<2. && ph.ispromptgen";
        // //TCut selcut = "ph.pt>25. && ph.isbarrel && (ph.gene/ph.scrawe)>0. && (ph.gene/ph.scrawe)<2.";

        // TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
        // TCut prescale10 = "(evt%10==0)";
        // TCut prescale10alt = "(evt%10==1)";
        // TCut prescale25 = "(evt%25==0)";
        // TCut prescale100 = "(evt%100==0)";  
        // TCut prescale1000 = "(evt%1000==0)";  
        // TCut evenevents = "(evt%2==0)";
        // TCut oddevents = "(evt%2==1)";
        // TCut prescale100alt = "(evt%100==1)";
        // TCut prescale1000alt = "(evt%1000==1)";
        // TCut prescale50alt = "(evt%50==1)";
        // //TCut oddevents = prescale100;

        // if (doele) 
        //     weightvar.SetTitle(prescale100alt*selcut);
        // else
        //     weightvar.SetTitle(selcut);

        // RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   

        // if (doele) 
        // weightvar.SetTitle(prescale1000alt*selcut);
        // else
        // weightvar.SetTitle(prescale10alt*selcut);


    cout << "Creating dataset" << endl;
    //RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet( "hdatasmall", dtree, vars, weightvar );
    RooDataSet *hdata = RooTreeConvert::CreateDataSet( "hdata", dtree, vars, weightvar );

    // const HybridGBRForest *forest = func->Forest();
    const HybridGBRForestFlex *forest = func->Forest();

    // const auto trees = forest->Trees();
    // cout << typeid( trees[1] ).name() << endl;
    // cout << &trees << endl;
    // cout << endl;

    // HybridGBRForest.h:
    //     std::vector<std::vector<HybridGBRTree> > &Trees() { return fTrees; }
    //     const std::vector<std::vector<HybridGBRTree> > &Trees() const { return fTrees; }
    //     So vector of vector of trees

    // HybridGBRForestFlex.h:
    //     std::vector<HybridGBRTreeD> &Trees() { return fTrees; }
    //     const std::vector<HybridGBRTreeD> &Trees() const { return fTrees; }
    //     So just a vector of treeDs

    // --> For now assume n_tgts = 1 then, and that forest->Trees returns a vector of trees


    int ntrees = 0;
    for ( unsigned int itree=0; itree < forest->Trees().size(); ++itree) {

        if ( forest->Trees()[itree].Responses().size() > 1 )
            ntrees++;

        }

    cout << "ntrees = " << ntrees << endl;

    // Old code:
        // for ( unsigned int itgt=0; itgt < forest->Trees().size(); ++itgt) {
            
        //     cout << "target " << itgt << endl;

        //     int ntrees = 0;

        //     for (unsigned int itree = 0; itree < forest->Trees().at(itgt).size(); ++itree) {
        //         if ( forest->Trees()[itgt][itree].Responses().size() > 1 )
        //             ++ntrees;
        //         }

        //     printf("itgt = %i, ntrees = %i\n", int(itgt),ntrees);

        //     }


    //#######################################
    // Prepare plotting variables
    //#######################################

    RooRealVar *scetavar = ws->var("var_2");
    
    RooAbsPdf *sigpdf;
    RooAbsReal *sigmeanlim;
    RooAbsReal *sigwidthlim;
    RooAbsReal *signlim;
    RooAbsReal *sign2lim;

    if (dobarrel){
        sigpdf = ws->pdf("sigpdfEB");
        sigmeanlim = ws->function("sigmeanlimEB");
        sigwidthlim = ws->function("sigwidthlimEB");
        signlim = ws->function("signlimEB");
        sign2lim = ws->function("sign2limEB");
        }
    else {
        sigpdf = ws->pdf("sigpdfEE");
        sigmeanlim = ws->function("sigmeanlimEE");
        sigwidthlim = ws->function("sigwidthlimEE");
        signlim = ws->function("signlimEE");
        sign2lim = ws->function("sign2limEE");
        }

    // No clear alternative for these 2 in my output root file
    // RooAbsReal *alphalim = ws->function("alpha1");
    // RooAbsReal *alpha2lim = ws->function("alpha2");  

    RooFormulaVar ecor("ecor","","1./(@0*@1)",RooArgList(*tgtvar,*sigmeanlim));
    //RooFormulaVar ecor("ecor","","@1/@0",RooArgList(*tgtvar,*sigmeanlim));
    //RooFormulaVar ecor("ecor","","@0/@1",RooArgList(*tgtvar,*sigmeanlim));
    //RooFormulaVar ecor("ecor","","@0",RooArgList(*tgtvar));
    RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
    ecorvar->setRange(0.,2.);
    ecorvar->setBins(800);

    RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
    RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
    rawvar->setRange(0.,2.);
    rawvar->setBins(800);

    /*  RooFormulaVar eraw("eraw","","@0",RooArgList(*tgtvar));
    RooRealVar *erawvar = (RooRealVar*)hdatasig->addColumn(eraw);
    erawvar->setRange(0.,2.);
    erawvar->setBins(400); */ 

    RooDataSet *hdataclone = new RooDataSet( *hdata, "hdataclone" );
    RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
    RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
    RooRealVar *nvar = (RooRealVar*)hdataclone->addColumn(*signlim);
    RooRealVar *n2var = (RooRealVar*)hdataclone->addColumn(*sign2lim);

    // No clear alternative for these 2 in my output root file
    // RooRealVar *alphavar = (RooRealVar*)hdataclone->addColumn(*alphalim);
    // RooRealVar *alpha2var = (RooRealVar*)hdataclone->addColumn(*alpha2lim);

    // // hdataclone = (RooDataSet*)hdataclone->reduce("sigwidthlim>0.017");


    //#######################################
    // Plotting
    //#######################################

    TCanvas *craw = new TCanvas;
    //RooPlot *plot = tgtvar->frame(0.6,1.2,100);
    RooPlot *plot = tgtvar->frame(0.6,2.0,100);
    hdata->plotOn(plot);
    sigpdf->plotOn(plot,ProjWData(*hdata));
    plot->Draw();
    craw->SaveAs("RawE.eps");
    craw->SaveAs("RawE.png");

    craw->SetLogy();
    plot->SetMinimum(0.1);
    craw->SaveAs("RawElog.eps");
    craw->SaveAs("RawElog.png");

    /*  new TCanvas;
    RooPlot *plotsig = tgtvar->frame(0.6,1.2,100);
    hdatasig->plotOn(plotsig);
    sigpdf.plotOn(plotsig,ProjWData(*hdatasig));
    plotsig->Draw(); */ 

    TCanvas *cmean = new TCanvas;
    RooPlot *plotmean = meanvar->frame(0.8,2.0,100);
    hdataclone->plotOn(plotmean);
    plotmean->Draw();
    cmean->SaveAs("mean.eps");
    cmean->SaveAs("mean.png");

    TCanvas *cwidth = new TCanvas;
    RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
    hdataclone->plotOn(plotwidth);
    plotwidth->Draw();
    cwidth->SaveAs("width.eps");
    cwidth->SaveAs("width.png");

    TCanvas *cn = new TCanvas;
    RooPlot *plotn = nvar->frame(0.,111.,200);
    hdataclone->plotOn(plotn);
    plotn->Draw();
    cn->SaveAs("n.eps");
    cn->SaveAs("n.png");

    TCanvas *cn2 = new TCanvas;
    RooPlot *plotn2 = n2var->frame(0.,111.,100);
    hdataclone->plotOn(plotn2);
    plotn2->Draw();
    cn2->SaveAs("n2.eps");
    cn2->SaveAs("n2.png");

    // TCanvas *calpha = new TCanvas;
    // RooPlot *plotalpha = alphavar->frame(0.,6.,200);
    // hdataclone->plotOn(plotalpha);
    // plotalpha->Draw();    
    // calpha->SaveAs("alpha.eps");
    // calpha->SetLogy();
    // plotalpha->SetMinimum(0.1);

    // TCanvas *calpha2 = new TCanvas;
    // RooPlot *plotalpha2 = alpha2var->frame(0.,6.,200);
    // hdataclone->plotOn(plotalpha2);
    // plotalpha2->Draw();      
    // calpha2->SaveAs("alpha2.eps");

    TCanvas *ceta = new TCanvas;
    RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
    hdataclone->plotOn(ploteta);
    ploteta->Draw();      
    ceta->SaveAs("eta.eps");  
    ceta->SaveAs("eta.png");  

    //TH1 *heold = hdatasigtest->createHistogram("heold",testvar);
    //TH1 *heraw = hdata->createHistogram("heraw",*tgtvar,Binning(800,0.,2.));
    TH1 *heraw = hdata->createHistogram("hraw",*rawvar,Binning(800,0.,2.));
    TH1 *hecor = hdata->createHistogram("hecor",*ecorvar);

    //heold->SetLineColor(kRed);
    hecor->SetLineColor(kBlue);
    heraw->SetLineColor(kMagenta);

    hecor->GetXaxis()->SetRangeUser(0.6,1.2);
    //heold->GetXaxis()->SetRangeUser(0.6,1.2);

    TCanvas *cresponse = new TCanvas;

    hecor->Draw("HIST");
    //heold->Draw("HISTSAME");
    heraw->Draw("HISTSAME");
    cresponse->SaveAs("response.eps");
    cresponse->SaveAs("response.png");
    cresponse->SetLogy();
    cresponse->SaveAs("responselog.eps");
    cresponse->SaveAs("responselog.png");


    printf("make fine histogram\n");
    TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(20e3,0.,2.));

    printf("calc effsigma\n");

    double effsigma = effSigma(hecorfine);

    printf("effsigma = %5f\n",effsigma);

    /*  new TCanvas;
    RooPlot *ploteold = testvar.frame(0.6,1.2,100);
    hdatasigtest->plotOn(ploteold);
    ploteold->Draw();    

    new TCanvas;
    RooPlot *plotecor = ecorvar->frame(0.6,1.2,100);
    hdatasig->plotOn(plotecor);
    plotecor->Draw(); */   
  
    }


void eregtestThomas(){
    
    printf( "\nCalling the eregtest() function\n" );

    eregtest();

    printf( "\nProgram finished.\n" );    

    }
