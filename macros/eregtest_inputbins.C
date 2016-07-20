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

#include "TColor.h"

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




void SingleVarPlot( RooDataSet *hdata, RooRealVar *somevar, Double_t x_min, Double_t x_max, Int_t n_bins, TString title ){

    TCanvas *c1 = new TCanvas;
    RooPlot *plot = somevar->frame( x_min, x_max, n_bins );
    hdata->plotOn( plot );
    plot->Draw();
    // c1->SaveAs( title + ".eps" );
    c1->SaveAs( "singlevar_" + title + ".png" );

    }

void PlotSlice_2x2(
        RooRealVar* yvar,
        TString ytitle,
        Double_t ymin,
        Double_t ymax,
        std::vector<RooRealVar*> slice_vars,
        std::vector<TString>     slice_vartitles,
        RooDataSet *hdata,
        TString filename_prefix
        ){
    
    //Define Canvas properties
    Double_t canvas_width(1300);
    Double_t canvas_height(800);
    
    TCanvas *canvas = new TCanvas( ytitle, ytitle, canvas_width, canvas_height );
    canvas->Divide(2,2);

    int n_pads = 4;

    for( int i_pad = 1; i_pad <= n_pads; i_pad++ ){

        RooRealVar *slice_var = slice_vars[i_pad-1];
        TString slice_vartitle = slice_vartitles[i_pad-1];

        canvas->cd(i_pad);
        gPad->SetGridx();
        gPad->SetGridy();
        TH2* histo2 = hdata->createHistogram( *slice_var, *yvar );
        // TProfile *hprof = histo2->ProfileX();
        TProfile *hprof = histo2->ProfileX( "ProfX" + ytitle + slice_vartitle , 1, -1, "X" );
        
        hprof->SetErrorOption("s");
        hprof->GetYaxis()->SetRangeUser( ymin, ymax );
        hprof->GetXaxis()->SetTitle( slice_vartitles[i_pad-1] );
        hprof->GetYaxis()->SetTitle( ytitle );
        hprof->Draw();

        }

    canvas->SaveAs( filename_prefix + ytitle + ".png" );

    }


void PlotSlice_2x2_2d(
        RooRealVar* yvar,
        TString ytitle,
        Double_t ymin,
        Double_t ymax,
        std::vector<RooRealVar*> slice_vars,
        std::vector<TString>     slice_vartitles,
        RooDataSet *hdata,
        TString filename_prefix
        ){
    
    //Define Canvas properties
    Double_t canvas_width(1300);
    Double_t canvas_height(800);
    
    TCanvas *canvas = new TCanvas( ytitle, ytitle, canvas_width, canvas_height );
    canvas->Divide(2,2);

    int n_pads = 4;

    for( int i_pad = 1; i_pad <= n_pads; i_pad++ ){

        RooRealVar *slice_var = slice_vars[i_pad-1];
        TString slice_vartitle = slice_vartitles[i_pad-1];

        canvas->cd(i_pad);
        gPad->SetGridx();
        gPad->SetGridy();
        TH2F* histo2 = hdata->createHistogram( *slice_var, *yvar );
        // TProfile *hprof = histo2->ProfileX();
        // TProfile *hprof = histo2->ProfileX( "someprofile", 1, -1, "X" );
        
        // hprof->SetErrorOption("s");
        // hprof->GetYaxis()->SetRangeUser( ymin, ymax );
        // hprof->GetXaxis()->SetTitle( slice_vartitles[i_pad-1] );
        // hprof->GetYaxis()->SetTitle( ytitle );
        // hprof->Draw();

        histo2->SetContour( 103 );
        histo2->GetZaxis()->SetRangeUser(0.0, 2.0);

        histo2->Draw("COLZ");

        }

    canvas->SaveAs( filename_prefix + ytitle + ".png" );

    }

void PlotSlice_2x3(
        RooRealVar* yvar,
        TString ytitle,
        Double_t ymin,
        Double_t ymax,
        std::vector<RooRealVar*> slice_vars,
        std::vector<TString>     slice_vartitles,
        RooDataSet *hdata,
        TString filename_prefix
        ){
    
    //Define Canvas properties
    Double_t canvas_width(800);
    Double_t canvas_height(1300);
    
    TCanvas *canvas = new TCanvas( ytitle, ytitle, canvas_width, canvas_height );
    canvas->Divide(2,3);

    int n_pads = 6;

    for( int i_pad = 1; i_pad <= n_pads; i_pad++ ){

        RooRealVar *slice_var = slice_vars[i_pad-1];
        TString slice_vartitle = slice_vartitles[i_pad-1];

        canvas->cd(i_pad);
        gPad->SetGridx();
        gPad->SetGridy();
        TH2* histo2 = hdata->createHistogram( *slice_var, *yvar );
        // TProfile *hprof = histo2->ProfileX();
        TProfile *hprof = histo2->ProfileX( "someprofile", 1, -1, "X" );

        hprof->SetErrorOption("s");
        hprof->GetYaxis()->SetRangeUser( ymin, ymax );
        hprof->GetXaxis()->SetTitle( slice_vartitles[i_pad-1] );
        hprof->GetYaxis()->SetTitle( ytitle );
        hprof->Draw();

        }

    canvas->SaveAs( filename_prefix + ytitle + ".png" );

    }


void PlotSlice_2x2_overlap(
        RooRealVar* yvar1,
        RooRealVar* yvar2,
        TString y1title,
        TString y2title,
        Double_t ymin,
        Double_t ymax,
        std::vector<RooRealVar*> slice_vars,
        std::vector<TString>     slice_vartitles,
        RooDataSet *hdata,
        TString filename_prefix
        ){
    
    //Define Canvas properties
    Double_t canvas_width(1300);
    Double_t canvas_height(800);
    
    TCanvas *canvas = new TCanvas( y1title, y1title, canvas_width, canvas_height );
    canvas->Divide(2,2);

    int n_pads = 4;

    for( int i_pad = 1; i_pad <= n_pads; i_pad++ ){

        RooRealVar *slice_var = slice_vars[i_pad-1];
        TString slice_vartitle = slice_vartitles[i_pad-1];

        canvas->cd(i_pad);
        gPad->SetGridx();
        gPad->SetGridy();
        TH2* histo2 = hdata->createHistogram( *slice_var, *yvar1 );
        // TProfile *hprof = histo2->ProfileX();
        TProfile *hprof = histo2->ProfileX( "ProfX_" + y1title + slice_vartitle , 1, -1, "X" );

        hprof->SetErrorOption("s");
        hprof->GetYaxis()->SetRangeUser( ymin, ymax );
        hprof->GetXaxis()->SetTitle( slice_vartitle );
        hprof->GetYaxis()->SetTitle( y1title + " " + y2title );
        hprof->Draw();

        TH2* histo22 = hdata->createHistogram( *slice_var, *yvar2 );
        // TProfile *hprofs = histo2->ProfileX();
        TProfile *hprof2 = histo22->ProfileX( "ProfX_" + y2title + slice_vartitle , 1, -1, "X" );

        hprof2->SetErrorOption("s");
        hprof2->SetLineColor(2);
        hprof2->SetMarkerColor(2);
        // hprof2->GetYaxis()->SetRangeUser( ymin, ymax );
        // hprof2->GetXaxis()->SetTitle( slice_vartitles[i_pad-1] );
        // hprof2->GetYaxis()->SetTitle( ytitle );
        hprof2->Draw("SAME");

        }

    canvas->SaveAs( filename_prefix + y1title + ".png" );

    }


//#######################################
// Main
//#######################################

void eregtest(){

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

    RooWorkspace *ws;
    if (dobarrel)     
        ws = (RooWorkspace*)fws->Get("wereg_eb");  
    else
        ws = (RooWorkspace*)fws->Get("wereg_ee");  
    ws->Print();

    // Get Target var
    RooRealVar *tgtvar = ws->var("targetvar");

    // Get GBRFunction
    RooGBRFunctionFlex *func;
    if (dobarrel) 
        func = static_cast<RooGBRFunctionFlex*>( ws->arg("sigmeantfuncEB") );
    else
        func = static_cast<RooGBRFunctionFlex*>( ws->arg("sigmeantfuncEE") );

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

    // Add more vars
    RooRealVar pt_p( "pt", "pt", 1., 0.0, 300.0 );
    RooRealVar *ptvar = &pt_p;
    vars.add( *ptvar );


    //=======================================
    // Open the NTuple

    // This probably means all events are weighted equally now
    RooRealVar weightvar("weightvar","",1.);

    TFile *fdin = TFile::Open( Ntup );
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("een_analyzer");
    TTree *dtree = (TTree*)ddir->Get(tree_name);


    //#######################################
    // Start with cuts, count trees
    //#######################################

    // Cut on event number, should be orthogonal to training cuts
    TCut eventcut = "eventNumber%2==1";

    if ( (string)getenv("TESTRUN")=="Y" ){
        // Use only part of sample
        // TCut NtupIDcut = "(NtupID>400&&NtupID<800) || (NtupID>12000&&NtupID<12400)";
        TCut NtupIDcut = "eventNumber%20==1";
        eventcut *= NtupIDcut;
        }

    weightvar.SetTitle( eventcut * selcut );

    cout << "Creating dataset" << endl;
    RooDataSet *hdata = RooTreeConvert::CreateDataSet( "hdata", dtree, vars, weightvar );

    // const HybridGBRForest *forest = func->Forest();
    const HybridGBRForestFlex *forest = func->Forest();

    int ntrees = 0;
    for ( unsigned int itree=0; itree < forest->Trees().size(); ++itree) {
        if ( forest->Trees()[itree].Responses().size() > 1 ) ntrees++;
        }
    cout << "ntrees = " << ntrees << endl;


    /*******************************************************************************************
    ********************************************************************************************
       FROM HERE ON PLOTTING
    ********************************************************************************************
    ********************************************************************************************/


    //#######################################
    // Prepare plotting variables
    //#######################################

    // Simple variables

    // Access training variables from workspace
    RooRealVar *nVtx                               = ws->var("var_0");
    RooRealVar *scRawEnergy                        = ws->var("var_1");
    RooRealVar *scEta                              = ws->var("var_2");
    RooRealVar *scPhi                              = ws->var("var_3");
    RooRealVar *scEtaWidth                         = ws->var("var_4");
    RooRealVar *scPhiWidth                         = ws->var("var_5");
    RooRealVar *scSeedR9                           = ws->var("var_6");
    RooRealVar *scSeedRawEnergy_scRawEnergy        = ws->var("var_7");
    RooRealVar *scSeedEmax                         = ws->var("var_8");
    RooRealVar *scSeedE2nd                         = ws->var("var_9");
    RooRealVar *scSeedLeftRightAsym                = ws->var("var_10");
    RooRealVar *scSeedTopBottomAsym                = ws->var("var_11");
    RooRealVar *scSeedSigmaIetaIeta                = ws->var("var_12");
    RooRealVar *scSeedSigmaIetaIphi                = ws->var("var_13");
    RooRealVar *scSeedSigmaIphiIphi                = ws->var("var_14");
    RooRealVar *N_ECALClusters                     = ws->var("var_15");
    RooRealVar *clusterMaxDR                       = ws->var("var_16");
    RooRealVar *clusterMaxDRDPhi                   = ws->var("var_17");
    RooRealVar *clusterMaxDRDEta                   = ws->var("var_18");
    RooRealVar *clusterMaxDRRawEnergy_scRawEnergy  = ws->var("var_19");

    scEta           ->setRange( -2.5 ,2.5  );
    scPhi           ->setRange( -4.  ,4.  );
    scSeedR9        ->setRange( 0.65 ,1.  );
    scEtaWidth      ->setRange( 0.   ,0.025  );
    scPhiWidth      ->setRange( 0.   ,0.2  );
    N_ECALClusters  ->setRange( 0.   ,20.  );
    nVtx            ->setRange( 0.   ,50.  );
    // hoveretower  ->setRange( 0.   ,0.2  );
    // rho          ->setRange( 0.   ,40.  );

   
    RooAbsPdf  *sigpdf;
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


    // RooFormulaVar ecor("ecor","","1./(@0*@1)",RooArgList(*tgtvar,*sigmeanlim));
    RooFormulaVar ecor("ecor","","(@1/@0)",RooArgList(*tgtvar,*sigmeanlim));
    RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
    ecorvar->setRange(0.,2.);
    ecorvar->setBins(800);

    RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
    RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
    rawvar->setRange(0.,2.);
    rawvar->setBins(800);

    // RooDataSet *hdataclone = new RooDataSet( *hdata, "hdataclone" );
    RooRealVar *meanvar  = (RooRealVar*)hdata->addColumn(*sigmeanlim);
    RooRealVar *widthvar = (RooRealVar*)hdata->addColumn(*sigwidthlim);
    RooRealVar *nvar     = (RooRealVar*)hdata->addColumn(*signlim);
    RooRealVar *n2var    = (RooRealVar*)hdata->addColumn(*sign2lim);

    meanvar->setRange(0.,2.);
    meanvar->setBins(800);

    // What does this do, and should it be necessary at some point?
    // hdataclone = (RooDataSet*)hdataclone->reduce("sigwidthlim>0.017");

    // Really not sure if this is necessary
    // RooFormulaVar pt_formula("pt","","@0",RooArgList(*pt));
    // RooRealVar *ptvar = (RooRealVar*)hdata->addColumn(pt_formula);
    // ptvar->setRange( 0., 300. );
    // ptvar->setBins(200);


    //#######################################
    // Plotting standard plots
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

    // Simple single variable plots
    SingleVarPlot( hdata, meanvar,  0.8, 2.0, 100, "mean" );
    SingleVarPlot( hdata, widthvar, 0.0, 0.05, 100, "width" );
    SingleVarPlot( hdata, nvar,     0.0, 111.0, 200 , "n" );
    SingleVarPlot( hdata, n2var,    0.0, 111.0, 100 , "n2" );
    SingleVarPlot( hdata, scEta,    -2.6, 2.6, 200 , "scEta" );
    SingleVarPlot( hdata, ptvar,       0.0, 350.0, 200 , "pt" );


    //#######################################
    // Response curve plot
    //#######################################

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

    // Does it make sense to do this for the whole distribution?
    // printf("make fine histogram\n");
    TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(20e3,0.,2.));
    // printf("calc effsigma\n");
    double effsigma = effSigma(hecorfine);
    printf("effsigma = %5f\n",effsigma);


    //#######################################
    // eta-phi plot
    //#######################################

    TCanvas *canvas_phi_eta = new TCanvas("phi_eta_Histo","Histograms of phi and eta variables", 800, 600 );
    TPad *padTop2 = new TPad("padTop2","padTop",0.05,0.5,0.95,0.95,256);
    TPad *padBottom2 = new TPad("padBottom2","padBottom",0.05,0.05,0.95,0.45,256);
    padTop2->Divide(3,1);
    padBottom2->Divide(3,1);

    padTop2->Draw();
    padBottom2->Draw();

    padTop2->cd(1);
    dtree->Draw("scEta", selcut );
    padTop2->cd(2);
    dtree->Draw("scEtaWidth", selcut );  
    padTop2->cd(3);
    dtree->Draw("scSeedCryIetaV2", selcut );

    padBottom2->cd(1);
    dtree->Draw("scPhi", selcut );
    padBottom2->cd(2);
    dtree->Draw("scPhiWidth", selcut );
    padBottom2->cd(3);
    dtree->Draw("scSeedCryIphiV2", selcut );
    
    canvas_phi_eta->SaveAs( "phi_eta_KH.png" );


    //#######################################
    // Plotting ecor and eraw binned in variables
    //#######################################

    // ======================================
    // Prepare vectors of variables for slices

    // Length 4 vector for 2x2 plots
    std::vector<RooRealVar*> slice4_vars;
    std::vector<TString>     slice4_vartitles;
    
    std::vector<RooRealVar*> slice4_vars2;
    std::vector<TString>     slice4_vartitles2;

    slice4_vars.push_back( scEta       );
    slice4_vars.push_back( scPhi       );
    slice4_vars.push_back( scEtaWidth  );
    slice4_vars.push_back( scPhiWidth  );
    
    slice4_vars2.push_back( ptvar       );
    slice4_vars2.push_back( nVtx        );
    slice4_vars2.push_back( scSeedR9    );
    slice4_vars2.push_back( N_ECALClusters );

    slice4_vartitles.push_back( "scEta"       );
    slice4_vartitles.push_back( "scPhi"       );
    slice4_vartitles.push_back( "scEtaWidth"  );
    slice4_vartitles.push_back( "scPhiWidth"  );

    slice4_vartitles2.push_back( "ptvar"       );
    slice4_vartitles2.push_back( "nVtx"        );
    slice4_vartitles2.push_back( "r9"          );
    slice4_vartitles2.push_back( "N_ECALClusters" );


    // Length 6 vector for 2x3 plots
    // std::vector<RooRealVar*> slice6_vars;
    // std::vector<TString>     slice6_vartitles;

    // slice6_vars.clear();
    // slice6_vartitles.clear();

    // slice6_vars.push_back( nVtx        );
    // slice6_vars.push_back( scSeedR9          );
    // slice6_vars.push_back( scEtaWidth  );
    // slice6_vars.push_back( scPhiWidth  );
    // slice6_vars.push_back( N_ECALClusters );
    // slice6_vars.push_back( ptvar       );

    // slice6_vartitles.push_back( "nVtx"        );
    // slice6_vartitles.push_back( "r9"          );
    // slice6_vartitles.push_back( "scEtaWidth"  );
    // slice6_vartitles.push_back( "scPhiWidth"  );
    // slice6_vartitles.push_back( "N_ECALClusters" );
    // slice6_vartitles.push_back( "ptvar"       );


    // Actual plotting

    // Format:
    // void PlotSlice(
    //         RooRealVar* yvar,
    //         Double_t ymin,
    //         Double_t ymax,
    //         TString ytitle,
    //         std::vector<RooRealVar*> slice_vars,
    //         std::vector<TString>     slice_vartitles,
    //         RooDataSet *hdata
    //         ){

    // PlotSlice_2x2(
    //     ecorvar, "ecor", 0.0, 1.1,
    //     slice4_vars, slice4_vartitles, hdata,
    //     "slice_etaphi_"
    //     );
    // // Same thing, now for raw    
    // PlotSlice_2x2(
    //     rawvar, "rawvar", 0.0, 1.1,
    //     slice4_vars, slice4_vartitles, hdata,
    //     "slice_etaphi_"
    //     );
    // PlotSlice_2x2(
    //     meanvar, "meanvar", 0.0, 2.0,
    //     slice4_vars, slice4_vartitles, hdata,
    //     "slice_etaphi_"
    //     );

    // PlotSlice_2x2(
    //     ecorvar, "ecor", 0.0, 1.1,
    //     slice4_vars2, slice4_vartitles2, hdata,
    //     "slice_other_"
    //     );
    // // Same thing, now for raw    
    // PlotSlice_2x2(
    //     rawvar, "rawvar", 0.0, 1.1, 
    //     slice4_vars2, slice4_vartitles2, hdata,
    //     "slice_other_"
    //     );
    // PlotSlice_2x2(
    //     meanvar, "meanvar", 0.0, 2.0, 
    //     slice4_vars2, slice4_vartitles2, hdata,
    //     "slice_other_"
    //     );


    // PlotSlice_2x2_overlap(
    //     rawvar, ecorvar, "rawvar", "ecor",
    //     0.0, 2.0,
    //     slice4_vars, slice4_vartitles, hdata,
    //     "slice_overlap_etaphi_"
    //     );
    // PlotSlice_2x2_overlap(
    //     rawvar, ecorvar, "rawvar", "ecor",
    //     0.0, 2.0,
    //     slice4_vars2, slice4_vartitles2, hdata,
    //     "slice_overlap_other_"
    //     );


    // PlotSlice_2x2_2d(
    //     meanvar, "meanvar", 0.0, 2.0,
    //     slice4_vars, slice4_vartitles, hdata,
    //     "slice_2d_etaphi_"
    //     );

    // PlotSlice_2x2_2d(
    //     meanvar, "meanvar", 0.0, 2.0,
    //     slice4_vars2, slice4_vartitles2, hdata,
    //     "slice_2d_other_"
    //     );


    // Later!
    // std::vector<Double_t> gen_pt_bins = { 0., 100., 300., 1000., 3000. };
    // int n_pt_bins = sizeof( gen_pt_bins );

    // Double_t pt_left, pt_right;
    // RooDataSet hdataclone;

    // for ( int ibin = 0; ibin < n_pt_bins-1; ibin++ ){
    //     pt_left = gen_pt_bins[ibin]; pt_right = gen_pt_bins[ibin+1];

    //     //*hdataclone = new RooDataSet( *hdata, "hdataclone" );
    //     }


    // PlotSlice_2x3(
    //     ecorvar, 0.0, 1.1, "ecor",
    //     slice6_vars, slice6_vartitles, hdata
    //     );
    // // Same thing, now for raw    
    // PlotSlice_2x3(
    //     rawvar, 0.0, 1.1, "rawvar",
    //     slice6_vars, slice6_vartitles, hdata
    //     );
    // PlotSlice_2x3(
    //     meanvar, 0.0, 2.0, "meanvar",
    //     slice6_vars, slice6_vartitles, hdata
    //     );


    }


void eregtest_inputbins(){

    // # Color map
    // config['NRGBs'] = 7
    // config['NCont'] = 103
    // config['stops'] = [ 0.00, 0.05, 0.25, 0.36, 0.58, 0.88, 1.00 ]
    // config['red']   = [ 1.00, 1.00, 0.00, 0.40, 0.87, 1.00, 0.85 ]
    // config['green'] = [ 1.00, 1.00, 0.71, 0.81, 1.00, 0.20, 0.00 ]
    // config['blue']  = [ 1.00, 1.00, 1.00, 1.00, 0.18, 0.00, 0.00 ]

    // ROOT.TColor.CreateGradientColorTable( self.config['NRGBs'],
    //     array('d', self.config['stops']),
    //     array('d', self.config['red']) ,
    //     array('d', self.config['green']) ,
    //     array('d', self.config['blue']) ,
    //     self.config['NCont'] )

    const int NRGBs             = 7;
    Double_t Red[NRGBs]    = { 0.00, 0.05, 0.25, 0.36, 0.58, 0.88, 1.00 };
    Double_t Green[NRGBs]  = { 1.00, 1.00, 0.00, 0.40, 0.87, 1.00, 0.85 };
    Double_t Blue[NRGBs]   = { 1.00, 1.00, 0.71, 0.81, 1.00, 0.20, 0.00 };
    Double_t Length[NRGBs] = { 1.00, 1.00, 1.00, 1.00, 0.18, 0.00, 0.00 };

    TColor::CreateGradientColorTable(
        7,
        Red,
        Green,
        Blue,
        Length,
        103
        );
    
    printf( "\nCalling the eregtest() function\n" );

    eregtest();

    printf( "\nProgram finished.\n" );    

    }
