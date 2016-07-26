#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
// #include "RooHybridBDTAutoPdf.h"
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
#include "TLatex.h"
// #include "HybridGBRForest.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
// #include "RooDoubleCBFast.h"
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
#include "TH2.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TGraphErrors.h"


#include "../interface/RooHybridBDTAutoPdf.h"
#include "../interface/HybridGBRForest.h"
#include "../interface/RooDoubleCBFast.h"
#include "../interface/HybridGBRForestFlex.h"


using namespace RooFit;


#include "BinPlot.h"


//#######################################
// MAIN
//#######################################

void DrawBinPlots(
    // bool dobarrel=true, bool weighted=false
    ) {
    
    gStyle->SetLineScalePS(1);

    // bool SILENCE = false;
    bool SILENCE = true;

    if (SILENCE){
        gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
        RooMsgService::instance().setSilentMode(true);
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Eval );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Generation );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Minimization );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Plotting );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Fitting );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Integration );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::LinkStateMgmt );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Caching );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Optimization );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::ObjectHandling );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::InputArguments );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Tracing );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Contents );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::DataHandling );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::NumIntegration );
        RooMsgService::instance().getStream(1).removeTopic( RooFit::Eval );
        }

    //#######################################
    // Determine particle and region
    //#######################################

    bool dobarrel;
    bool weighted=false;
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

    bool isElectron, isPhoton ;
    if ( getenv("PARTICLE")==NULL ) {
        cout << "Unclear whether particle is photon or electron" << endl;
        return;
        }
    else if ( (string)getenv("PARTICLE") == "electron" ){
        isElectron = true;
        isPhoton  = false;
        }
    else if ( (string)getenv("PARTICLE") == "photon" ){
        isElectron = false;
        isPhoton  = true;
        }

    bool highpt = false;
    if ( getenv("HIGHPT")==NULL )                highpt = false;
    else if ( (string)getenv("HIGHPT") == "N" )  highpt = false;
    else if ( (string)getenv("HIGHPT") == "Y" )  highpt = true ;

    bool useglobalptbins = false;
    if ( getenv("USEGLOBALPTBINS")==NULL )                useglobalptbins = false;
    else if ( (string)getenv("USEGLOBALPTBINS") == "N" )  useglobalptbins = false;
    else if ( (string)getenv("USEGLOBALPTBINS") == "Y" )  useglobalptbins = true ;
    else                                                  useglobalptbins = false;

    bool testcombination = true;


    //#######################################
    // Event selection
    //#######################################

    // Cut on event number, should be orthogonal to training cuts
    // TCut eventcut = "eventNumber%2==1";
    TCut eventcut = "";
    
    TCut NtupIDcut;
    if ( (string)getenv("TESTRUN")=="Y" ){
        // Use only part of sample
        // TCut NtupIDcut = "(NtupID>400&&NtupID<800) || (NtupID>12000&&NtupID<12400)";
        // TCut NtupIDcut = "NtupID%20==1||NtupID%21==1";
        TCut NtupIDcut = "eventNumber%20==1||eventNumber%20==0";;
        // if(highpt) NtupIDcut = "eventNumber%400==1" ;
        // else       NtupIDcut = "eventNumber%20==1" ;
        eventcut *= NtupIDcut;
        }

    RooRealVar weightvar("weightvar","",1.);
    TCut selweight;
    if(weighted)
        selweight= "(weight)";
    else
        selweight= "(1.)";
    
    weightvar.SetTitle( eventcut * selcut * selweight );


    //#######################################
    // Read IO stuff from environment
    //#######################################

    cout << "Getting Ntuple" << endl;

    //output dir
    TString dirname = getenv("PLOTDIR_FULLPATH") ;
    gSystem->mkdir(dirname,true);
    gSystem->cd(dirname);    

    // Read the Ntuple
    TString Ntup      = getenv("FLATNTUPLE");
    TString tree_name = getenv("NTUPLETREE");
    TString fname     = getenv("TRAININGOUTPUT");
    
    TFile *fdin = TFile::Open( Ntup );
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("een_analyzer");
    TTree *dtree = (TTree*)ddir->Get(tree_name);

    cout << "Finished getting Ntuple" << endl;


    //#######################################
    // Read variables
    //#######################################

    // ======================================
    // Workspace 

    TString infile = TString::Format("./%s",fname.Data());
    TFile *fws = TFile::Open(infile); 

    RooWorkspace *ws;
    if (dobarrel)     
        ws = (RooWorkspace*)fws->Get("wereg_eb");  
    else
        ws = (RooWorkspace*)fws->Get("wereg_ee");  
    // ws->Print();
  
    RooWorkspace *ws_comb = (RooWorkspace*)fws->Get("wereg_comb");


    // Read variables from workspace
    RooGBRTargetFlex *meantgt;
    if (dobarrel) 
        meantgt = static_cast<RooGBRTargetFlex*>( ws->arg("sigmeantEB") );
    else
        meantgt = static_cast<RooGBRTargetFlex*>( ws->arg("sigmeantEE") );

    RooRealVar *tgtvar = ws->var("targetvar");


    // ======================================
    // Define all the input variables

    cout << "Defining input variables" << endl;

    RooRealVar *scRawEnergy = new RooRealVar( "scRawEnergy", "scRawEnergy", 0.);
    RooRealVar *r9          = new RooRealVar( "r9", "r9", 0.);
    RooRealVar *nVtx        = new RooRealVar( "nVtx", "nVtx", 0.);
    RooRealVar *pt          = new RooRealVar( "pt", "pt", 0.);
    RooRealVar *genEta      = new RooRealVar( "genEta", "genEta", 0.);
    RooRealVar *genE        = new RooRealVar( "genEnergy", "genEnergy", 0.);
    RooRealVar *genPt       = new RooRealVar( "genPt", "genPt", 0.);

    // Set variables ranges -- Should be big enough!!!

    r9->setRange( 0., 1.2);

    if (highpt){
        pt     ->setRange( 0., 10000.);
        genE   ->setRange( 0., 10000. );
        genPt  ->setRange( 0., 10000. );
        }
    else {
        pt     ->setRange( 0., 300.);
        genE   ->setRange( 0., 1000. );
        genPt  ->setRange( 0., 300. );
        }

    if(dobarrel) 
        genEta->setRange( -1.5,1.5 );
    else
        genEta->setRange( -3,3 );

    RooArgList vars;
    vars.add(meantgt->FuncVars());
    vars.add(*tgtvar);

    vars.add(*genE);
    vars.add(*genPt);
    vars.add(*genEta);
    vars.add(*scRawEnergy);
    vars.add(*r9);
    vars.add(*pt);
    vars.add(*nVtx);

    // Add the old regression results if particle is electron
    RooRealVar *cor74E      = new RooRealVar( "corrEnergy74X",      "corrEnergy74X", 0. );
    RooRealVar *cor74Eerror = new RooRealVar( "corrEnergy74XError", "corrEnergy74XError", 0. );
    vars.add( *cor74E );
    vars.add( *cor74Eerror );


    // ======================================
    // Create the testing dataset

    RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);
    cout << "In total " << hdata->numEntries() << " entries in Ntuple" << endl;


    // ======================================
    // Define the variables for the regression output
    
    RooAbsPdf  *sigpdf;
    RooAbsReal *sigmeanlim;
    RooAbsReal *sigwidthlim;
    RooAbsReal *signlim;
    RooAbsReal *sign2lim;
    if (dobarrel){
        sigpdf      = ws->pdf("sigpdfEB");
        sigmeanlim  = ws->function("sigmeanlimEB");
        sigwidthlim = ws->function("sigwidthlimEB");
        signlim     = ws->function("signlimEB");
        sign2lim    = ws->function("sign2limEB");
        }
    else {
        sigpdf      = ws->pdf("sigpdfEE");
        sigmeanlim  = ws->function("sigmeanlimEE");
        sigwidthlim = ws->function("sigwidthlimEE");
        signlim     = ws->function("signlimEE");
        sign2lim    = ws->function("sign2limEE");
        }

    //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
    RooFormulaVar ecor("ecor","","(@1/@0)",RooArgList(*tgtvar,*sigmeanlim));
    RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
    ecorvar->setRange(0.,2);
    ecorvar->setBins(800);
    
    //formula for raw energy/true energy
    RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
    RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
    rawvar->setRange(0.,2.);
    rawvar->setBins(800);

    RooFormulaVar ecor74( "ecor74", "@0/@1", RooArgList( *cor74E, *genE ) );
    RooRealVar *ecor74var = (RooRealVar*)hdata->addColumn(ecor74);
    ecor74var->setRange(0.,2.);
    ecor74var->setBins(800);


    // if (testcombination):

    //     //formula for target of the EP combination
    //     RooFormulaVar EPcomb("EPcomb", "", "(@1/@0)", RooArgList(*tgtvar,*sigmeanlim) );
    //     RooRealVar *EPcombvar = (RooRealVar*)hdata->addColumn(EPcomb);
    //     EPcombvar->setRange(0.,2);
    //     EPcombvar->setBins(800);



    cout << "Finished reading variables" << endl;


    //#######################################
    // Load into classes and plot
    //#######################################

    Double_t global_ymin = 0.95;
    Double_t global_ymax = 1.03;



    //#######################################
    // Low pt plots, no further binning in pt
    //#######################################
    
    if (!useglobalptbins) {
    if (!highpt) {

        // ======================================
        // genpt slices

        BinPlot L_genpt_plot = BinPlot();

        const Int_t L_genpt_nbounds = 21;
        Double_t L_genpt_xbins[L_genpt_nbounds] = {
            0.,   10.,  20.,  30.,  40.,
            50.,  60.,  70.,  80.,  90.,
            100., 120., 140., 160., 180.,
            200., 220., 240., 260., 280.,
            300.
            };
        L_genpt_plot.Set_xbins( L_genpt_nbounds , L_genpt_xbins );

        // Load variables into class
        L_genpt_plot.hdata_          = hdata ;
        L_genpt_plot.rawvar_         = rawvar ;
        L_genpt_plot.ecorvar_        = ecorvar ;
        L_genpt_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            L_genpt_plot.draw_old_regression_ = true ;
            L_genpt_plot.ecor74var_ = ecor74var ;
            // }

        L_genpt_plot.slicevar_       = genPt ;
        L_genpt_plot.slicevartitle_  = "genPt"  ;
        L_genpt_plot.slicevarname_   = "genPt"  ;

        // These two number multiplied should be > nbins_
        L_genpt_plot.n_columns_      = 5  ;
        L_genpt_plot.n_rows_         = 5  ;
        
        L_genpt_plot.ymin_           = global_ymin ;
        L_genpt_plot.ymax_           = global_ymax ;

        L_genpt_plot.MakeSlicePlot();


        // ======================================
        // genEnergy slices

        BinPlot L_genE_plot = BinPlot();

        const Int_t L_genE_nbounds = 18;
        Double_t L_genE_xbins[L_genE_nbounds] = {
            0,    20,   40,   60,   80,
            100,  120,  140,  160,  180,
            200,  220,  240,  260,  280,
            300,  500,  800,
            };
        L_genE_plot.Set_xbins( L_genE_nbounds , L_genE_xbins );

        // Load variables into class
        L_genE_plot.hdata_          = hdata ;
        L_genE_plot.rawvar_         = rawvar ;
        L_genE_plot.ecorvar_        = ecorvar ;
        L_genE_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            L_genE_plot.draw_old_regression_ = true ;
            L_genE_plot.ecor74var_ = ecor74var ;
            // }

        L_genE_plot.slicevar_       =  genE ;
        L_genE_plot.slicevartitle_  = "genEnergy"  ;
        L_genE_plot.slicevarname_   = "genEnergy"  ;

        // These two number multiplied should be > nbins_
        L_genE_plot.n_columns_      = 5  ;
        L_genE_plot.n_rows_         = 4  ;
        
        L_genE_plot.ymin_           = global_ymin ;
        L_genE_plot.ymax_           = global_ymax ;

        L_genE_plot.MakeSlicePlot();


        // ======================================
        // genEta slices

        BinPlot L_genEta_plot = BinPlot();

        const Int_t L_genEta_nbounds_EB = 16;
        Double_t L_genEta_xbins_EB[L_genEta_nbounds_EB] = {
            0.0,  0.1,  0.2,  0.3,  0.4,
            0.5,  0.6,  0.7,  0.8,  0.9,
            1.0,  1.1,  1.2,  1.3,  1.4,
            1.5
            };

        const Int_t L_genEta_nbounds_EE = 11;
        Double_t L_genEta_xbins_EE[L_genEta_nbounds_EE] = {
            1.5, 1.6, 1.7, 1.8, 1.9,
            2.0, 2.2, 2.4, 2.6, 2.8,
            3.0
            };

        if (dobarrel) L_genEta_plot.Set_xbins( L_genEta_nbounds_EB , L_genEta_xbins_EB );
        else          L_genEta_plot.Set_xbins( L_genEta_nbounds_EE , L_genEta_xbins_EE );

        // Load variables into class
        L_genEta_plot.hdata_          = hdata ;
        L_genEta_plot.rawvar_         = rawvar ;
        L_genEta_plot.ecorvar_        = ecorvar ;
        L_genEta_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            L_genEta_plot.draw_old_regression_ = true ;
            L_genEta_plot.ecor74var_ = ecor74var ;
            // }

        L_genEta_plot.slicevar_       =  genEta ;
        L_genEta_plot.slicevartitle_  = "genEta"  ;
        L_genEta_plot.slicevarname_   = "genEta"  ;

        // These two number multiplied should be > nbins_
        L_genEta_plot.n_columns_      = 5  ;
        L_genEta_plot.n_rows_         = 4  ;
        
        L_genEta_plot.ymin_           = global_ymin ;
        L_genEta_plot.ymax_           = global_ymax ;

        L_genEta_plot.MakeSlicePlot();


        // ======================================
        // r9 slices

        BinPlot L_r9_plot = BinPlot();

        const Int_t L_r9_nbounds = 10;
        Double_t L_r9_xbins[L_r9_nbounds] = {
            0.80,  0.85,  0.90,  0.92,  0.94,
            0.95,  0.96,  0.97,  0.98,  1.02,
            };
        L_r9_plot.Set_xbins( L_r9_nbounds , L_r9_xbins );

        // Load variables into class
        L_r9_plot.hdata_          = hdata ;
        L_r9_plot.rawvar_         = rawvar ;
        L_r9_plot.ecorvar_        = ecorvar ;
        L_r9_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            L_r9_plot.draw_old_regression_ = true ;
            L_r9_plot.ecor74var_ = ecor74var ;
            // }

        L_r9_plot.slicevar_       =  r9 ;
        // L_r9_plot.slicevartitle_  = "scSeedR9"  ;
        // L_r9_plot.slicevarname_   = "scSeedR9"  ;
        L_r9_plot.slicevartitle_  = "r9"  ;
        L_r9_plot.slicevarname_   = "r9"  ;

        // These two number multiplied should be > nbins_
        L_r9_plot.n_columns_      = 4  ;
        L_r9_plot.n_rows_         = 3  ;
        
        L_r9_plot.ymin_           = global_ymin ;
        L_r9_plot.ymax_           = global_ymax ;

        L_r9_plot.MakeSlicePlot();


        }


    //#######################################
    // High pt plots, no further binning in pt
    //#######################################
    
    else if (highpt) {

        // ======================================
        // genpt slices

        BinPlot H_genpt_plot = BinPlot();

        const Int_t H_genpt_nbounds = 18;
        Double_t H_genpt_xbins[H_genpt_nbounds] = {
            0.,    20.,   40.,   60.,   80.,
            100.,  150.,  200.,  300.,  500.,
            700.,  1000., 1500., 2000., 3000.,
            4000., 5000., 6500.
            };
        
        H_genpt_plot.Set_xbins( H_genpt_nbounds , H_genpt_xbins );

        // Load variables into class
        H_genpt_plot.hdata_          = hdata ;
        H_genpt_plot.rawvar_         = rawvar ;
        H_genpt_plot.ecorvar_        = ecorvar ;
        H_genpt_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            H_genpt_plot.draw_old_regression_ = true ;
            H_genpt_plot.ecor74var_ = ecor74var ;
            // }

        H_genpt_plot.slicevar_       = genPt ;
        H_genpt_plot.slicevartitle_  = "genPt"  ;
        H_genpt_plot.slicevarname_   = "genPt"  ;

        // These two number multiplied should be > nbins_
        H_genpt_plot.n_columns_      = 5  ;
        H_genpt_plot.n_rows_         = 5  ;
        
        H_genpt_plot.ymin_           = global_ymin ;
        H_genpt_plot.ymax_           = global_ymax ;

        H_genpt_plot.MakeSlicePlot();


        // ======================================
        // genEta slices

        BinPlot H_genEta_plot = BinPlot();

        const Int_t H_genEta_nbounds_EB = 16;
        Double_t H_genEta_xbins_EB[H_genEta_nbounds_EB] = {
            0.0,  0.1,  0.2,  0.3,  0.4,
            0.5,  0.6,  0.7,  0.8,  0.9,
            1.0,  1.1,  1.2,  1.3,  1.4,
            1.5
            };

        const Int_t H_genEta_nbounds_EE = 11;
        Double_t H_genEta_xbins_EE[H_genEta_nbounds_EE] = {
            1.5, 1.6, 1.7, 1.8, 1.9,
            2.0, 2.2, 2.4, 2.6, 2.8,
            3.0
            };

        if (dobarrel) H_genEta_plot.Set_xbins( H_genEta_nbounds_EB , H_genEta_xbins_EB );
        else          H_genEta_plot.Set_xbins( H_genEta_nbounds_EE , H_genEta_xbins_EE );

        // Load variables into class
        H_genEta_plot.hdata_          = hdata ;
        H_genEta_plot.rawvar_         = rawvar ;
        H_genEta_plot.ecorvar_        = ecorvar ;
        H_genEta_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            H_genEta_plot.draw_old_regression_ = true ;
            H_genEta_plot.ecor74var_ = ecor74var ;
            // }

        H_genEta_plot.slicevar_       =  genEta ;
        H_genEta_plot.slicevartitle_  = "genEta"  ;
        H_genEta_plot.slicevarname_   = "genEta"  ;

        // These two number multiplied should be > nbins_
        H_genEta_plot.n_columns_      = 5  ;
        H_genEta_plot.n_rows_         = 4  ;
        
        H_genEta_plot.ymin_           = global_ymin ;
        H_genEta_plot.ymax_           = global_ymax ;

        H_genEta_plot.MakeSlicePlot();


        // ======================================
        // r9 slices

        BinPlot H_r9_plot = BinPlot();

        const Int_t H_r9_nbounds = 10;
        Double_t H_r9_xbins[H_r9_nbounds] = {
            0.80,  0.85,  0.90,  0.92,  0.94,
            0.95,  0.96,  0.97,  0.98,  1.02,
            };
        H_r9_plot.Set_xbins( H_r9_nbounds , H_r9_xbins );

        // Load variables into class
        H_r9_plot.hdata_          = hdata ;
        H_r9_plot.rawvar_         = rawvar ;
        H_r9_plot.ecorvar_        = ecorvar ;
        H_r9_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
            H_r9_plot.draw_old_regression_ = true ;
            H_r9_plot.ecor74var_ = ecor74var ;
            // }

        H_r9_plot.slicevar_       =  r9 ;
        H_r9_plot.slicevartitle_  = "r9"  ;
        H_r9_plot.slicevarname_   = "r9"  ;

        // These two number multiplied should be > nbins_
        H_r9_plot.n_columns_      = 4  ;
        H_r9_plot.n_rows_         = 3  ;
        
        H_r9_plot.ymin_           = global_ymin ;
        H_r9_plot.ymax_           = global_ymax ;

        H_r9_plot.MakeSlicePlot();


        }
        }//End of if(!useglobalptbins)
    

    //#######################################
    // High pt plots, plots redone for separate pt ranges
    //#######################################

    else if (useglobalptbins) {

        const Int_t n_globalptbins = 4;
        Double_t globalptbins[n_globalptbins+1] = {
            0.,    100.,  500.,  2500., 6500.
            };

        // These are the boundaries used for the genpt bin plots
        const Int_t n_globalptbounds = 25;
        Double_t globalptbounds[n_globalptbounds] = {
            0.,    20.,   30.,   40.,   50.,
            60.,   80.,   100.,
            150.,  200.,  250.,  300.,  400.,
            500.,  
            750.,  1000., 1250., 1500., 2000.,
            2500.,
            3000., 3500., 4000., 5000., 6500.
            };

        TString min_pt_str;
        TString max_pt_str;
        TString id_str;

        // Loop over the set global pt bins; make plots for every bin
        for( Int_t i_ptbin=0; i_ptbin < n_globalptbins; i_ptbin++ ) {

            Double_t min_pt = globalptbins[i_ptbin];
            Double_t max_pt = globalptbins[i_ptbin+1];

            min_pt_str.Form( "%.2f", min_pt );
            max_pt_str.Form( "%.2f", max_pt );
            TString sel_str = "(genPt>" + min_pt_str + "&&genPt<" + max_pt_str + ")" ;
            id_str.Form( "GENPT%.0f-%.0f_", min_pt, max_pt );

            // Create a new dataset specific for this pt bin
            RooDataSet *hdata_ptbin = (RooDataSet*)hdata->reduce(sel_str);

            cout << endl << "Processing genPt from " << min_pt_str << " to " << max_pt_str
                << " (" << id_str << ")" << endl;
            cout << "  Number of entries in this genPt selection = " << hdata_ptbin->numEntries() << endl;

            // ======================================
            // genpt slices

            // Bins are a bit tricky; extract the right bins from globalptbounds
            Int_t index_minpt_bound, index_maxpt_bound;
            Double_t current_minpt_bound=0., current_maxpt_bound=0.;
            for( Int_t i_ptbound=1; i_ptbound < n_globalptbounds; i_ptbound++ ) {
                
                if( globalptbounds[i_ptbound] > min_pt && globalptbounds[i_ptbound-1] <= min_pt )
                    index_minpt_bound = i_ptbound-1;
                
                if( globalptbounds[i_ptbound] >= max_pt && globalptbounds[i_ptbound-1] < max_pt )
                    index_maxpt_bound = i_ptbound;
                }

            // cout << "Using the following bins in the genPt plot:" << endl;
            // for( Int_t i_ptbound=index_minpt_bound; i_ptbound <= index_maxpt_bound; i_ptbound++ ) {
            //     cout << globalptbounds[i_ptbound] << ", " ;
            //     }
            // cout << endl << endl;

            // Load the found bins in a properly sized Double_t array
            Int_t H_genpt_nbounds = index_maxpt_bound - index_minpt_bound + 1;
            Double_t *H_genpt_xbins;
            H_genpt_xbins = (Double_t*)malloc( sizeof(Double_t)*(H_genpt_nbounds) );
            for ( Int_t ibin = 0; ibin < H_genpt_nbounds; ibin++ ){
                H_genpt_xbins[ibin] = globalptbounds[ index_minpt_bound + ibin ];
                // cout << "Loaded value " << globalptbounds[ index_minpt_bound + ibin ] << " at index " << ibin << endl;
                }

            BinPlot H_genpt_plot = BinPlot();
            
            H_genpt_plot.Append_filename( id_str );
            H_genpt_plot.Set_xbins( H_genpt_nbounds , H_genpt_xbins );

            // Load variables into class
            H_genpt_plot.hdata_          = hdata_ptbin ;
            H_genpt_plot.rawvar_         = rawvar ;
            H_genpt_plot.ecorvar_        = ecorvar ;
            H_genpt_plot.tgtvar_         = tgtvar ;

            H_genpt_plot.draw_old_regression_ = true ;
            H_genpt_plot.ecor74var_ = ecor74var ;

            H_genpt_plot.slicevar_       = genPt ;
            H_genpt_plot.slicevartitle_  = "genPt"  ;
            H_genpt_plot.slicevarname_   = "genPt"  ;

            // These two number multiplied should be > nbins_
            H_genpt_plot.n_columns_      = 3  ;
            H_genpt_plot.n_rows_         = 3  ;
            
            H_genpt_plot.ymin_           = global_ymin ;
            H_genpt_plot.ymax_           = global_ymax ;

            H_genpt_plot.MakeSlicePlot();


            // ======================================
            // genEta slices

            BinPlot H_genEta_plot = BinPlot();

            const Int_t H_genEta_nbounds_EB = 16;
            Double_t H_genEta_xbins_EB[H_genEta_nbounds_EB] = {
                0.0,  0.1,  0.2,  0.3,  0.4,
                0.5,  0.6,  0.7,  0.8,  0.9,
                1.0,  1.1,  1.2,  1.3,  1.4,
                1.5
                };

            const Int_t H_genEta_nbounds_EE = 11;
            Double_t H_genEta_xbins_EE[H_genEta_nbounds_EE] = {
                1.5, 1.6, 1.7, 1.8, 1.9,
                2.0, 2.2, 2.4, 2.6, 2.8,
                3.0
                };

            if (dobarrel) H_genEta_plot.Set_xbins( H_genEta_nbounds_EB , H_genEta_xbins_EB );
            else          H_genEta_plot.Set_xbins( H_genEta_nbounds_EE , H_genEta_xbins_EE );

            H_genEta_plot.Append_filename( id_str );

            // Load variables into class
            H_genEta_plot.hdata_          = hdata_ptbin ;
            H_genEta_plot.rawvar_         = rawvar ;
            H_genEta_plot.ecorvar_        = ecorvar ;
            H_genEta_plot.tgtvar_         = tgtvar ;

            H_genEta_plot.draw_old_regression_ = true ;
            H_genEta_plot.ecor74var_ = ecor74var ;

            H_genEta_plot.slicevar_       =  genEta ;
            H_genEta_plot.slicevartitle_  = "genEta"  ;
            H_genEta_plot.slicevarname_   = "genEta"  ;

            // These two number multiplied should be > nbins_
            H_genEta_plot.n_columns_      = 5  ;
            H_genEta_plot.n_rows_         = 4  ;
            
            H_genEta_plot.ymin_           = global_ymin ;
            H_genEta_plot.ymax_           = global_ymax ;

            H_genEta_plot.MakeSlicePlot();


            // ======================================
            // r9 slices

            BinPlot H_r9_plot = BinPlot();

            const Int_t H_r9_nbounds = 10;
            Double_t H_r9_xbins[H_r9_nbounds] = {
                0.80,  0.85,  0.90,  0.92,  0.94,
                0.95,  0.96,  0.97,  0.98,  1.02,
                };
            H_r9_plot.Set_xbins( H_r9_nbounds , H_r9_xbins );
            H_r9_plot.Append_filename( id_str );

            // Load variables into class
            H_r9_plot.hdata_          = hdata_ptbin ;
            H_r9_plot.rawvar_         = rawvar ;
            H_r9_plot.ecorvar_        = ecorvar ;
            H_r9_plot.tgtvar_         = tgtvar ;

            H_r9_plot.draw_old_regression_ = true ;
            H_r9_plot.ecor74var_ = ecor74var ;

            H_r9_plot.slicevar_       =  r9 ;
            H_r9_plot.slicevartitle_  = "r9"  ;
            H_r9_plot.slicevarname_   = "r9"  ;

            // These two number multiplied should be > nbins_
            H_r9_plot.n_columns_      = 4  ;
            H_r9_plot.n_rows_         = 3  ;
            
            H_r9_plot.ymin_           = global_ymin ;
            H_r9_plot.ymax_           = global_ymax ;

            H_r9_plot.MakeSlicePlot();


            }
        }







    }