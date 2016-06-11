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


//#######################################
// Calculates effective sigma
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
// Container class that contains the plotting script
//#######################################

class BinPlot {

    public:

    RooDataSet* hdata_;

    RooRealVar* rawvar_;
    RooRealVar* ecorvar_;
    RooRealVar* tgtvar_;

    // If old regression result is available it can be drawn as well
    bool draw_old_regression_ = false ;
    RooRealVar* ecor74var_;

    RooRealVar* slicevar_;
    TString slicevarname_;
    TString slicevartitle_;

    Int_t nbins_;
    Double_t binwidth_;
    Double_t binoffset_;
    Double_t* xbins_ = 0;

    Double_t fitxmin_ = 0.8, fitxmax_ = 1.1;

    // For plotting purposes
    Double_t ymin_ = 0.9, ymax_ = 1.1 ;

    Double_t ymin_sigma_ = 0.0, ymax_sigma_ = 0.1 ;

    // Dimensions of main slice canvas
    Int_t c_width_ = 2000, c_height_ = 1000;

    // Dimensions (in pixels) of 1 subplot of the perbin plots
    Int_t perbin_width_ = 500, perbin_height_ = 400;
    Int_t n_rows_, n_columns_;

    void FitOneSlice(
        RooDataSet *hdata_reduced,
        TString sel_str,
        RooRealVar *var,
        TString varname,
        TCanvas *c,
        Int_t ibin,

        Double_t *p_mean,
        Double_t *p_error_on_mean,
        Double_t *p_sigma,
        Double_t *p_error_on_sigma,
        Double_t *p_effsigma
        );

    void MakeSlicePlot();

    };


void BinPlot::FitOneSlice(
        // Input
        RooDataSet *hdata_reduced,
        TString sel_str,
        RooRealVar *var,
        TString varname,
        TCanvas *c,
        Int_t ibin,
        // Output , probably pointers!
        Double_t *p_mean,
        Double_t *p_error_on_mean,
        Double_t *p_sigma,
        Double_t *p_error_on_sigma,
        Double_t *p_effsigma
        ) {

    RooDataSet *hdata_var = (RooDataSet*)hdata_reduced->reduce(*var);

    RooRealVar mean( "mean_" + varname, "mean_" + varname, 1.,0.9,1.1);
    RooRealVar sig(  "sig_" + varname,  "sig_" + varname,  0.01,0.0002,0.8);
    RooRealVar a1(   "a1_" + varname,   "a1_" + varname,   3,0.05,10);
    RooRealVar a2(   "a2_" + varname,   "a2_" + varname,   3,0.05,10);
    RooRealVar n1(   "n1_" + varname,   "n1_" + varname,   3,1.01,500);
    RooRealVar n2(   "n2_" + varname,   "n2_" + varname,   3,1.01,500);

    RooDoubleCBFast pdfCB(
        "pdfCB_" + varname, "pdfCB_" + varname,
        *var, 
        mean, sig, a1, n1, a2, n2
        );

    pdfCB.fitTo( *hdata_var, Range( fitxmin_, fitxmax_ ) );

    // h_raw->AddBinContent( ibin+1, mean.getVal() );
    // h_raw->SetBinError(   ibin+1, mean_raw.getError() );

    // Translate pdfCB_raw to a TH1
    // cout << "Turning CB_raw into a histogram (var=" + slicevarname_ + ")" << endl;
    TH1* h_CB = pdfCB.createHistogram( "h_CB", *var, Binning(1000) );
    // cout << "Getting effective sigma for the histogram" << endl;
    effsigma = effSigma( h_CB );

    *p_effsigma         = effsigma;
    *p_mean             = mean.getVal();
    *p_error_on_mean    = mean.getError();
    *p_sigma            = sig.getVal();
    *p_error_on_sigma   = sig.getError();


    // ======================================
    // Drawing

    c->cd(ibin+1);

    RooPlot *var_datapoints = var->frame(0.,2, 250);
    
    var_datapoints->SetTitle( "Distribution of raw/true " + sel_str );
    var_datapoints->GetXaxis()->SetTitle("raw/true");
    gPad->SetGridx();
    gPad->SetGridy();
    
    hdata_reduced->plotOn( var_datapoints, Name( "datapoints" + varname ), MarkerSize(0.02));
    pdfCB.plotOn( var_datapoints, Name( "fit" + varname ), LineColor(kRed));
    var_datapoints->Draw();
    
    TLegend *legend = new TLegend( 0.6, 0.75, 0.9, 0.9 );
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry( "datapoints" + varname, "datapoints " + varname, "lep");
    legend->AddEntry( "fit" + varname, "fit " + varname, "l");
    legend->Draw();

    }


void BinPlot::MakeSlicePlot(){


    //#######################################
    // Prepare objects
    //#######################################

    // Make xbins axis
    xbins_ = (Double_t*)malloc( sizeof(Double_t)*(nbins_+1) );
    for( int ibin = 0; ibin <= nbins_; ibin++ ) {
        xbins_[ibin] = (Double_t)( binoffset_ + ibin * binwidth_ );
        }

    // Main histograms, binned in xbins
    TH1F* h_raw   = new TH1F( "h_raw", "", nbins_, xbins_ );
    TH1F* h_cor   = new TH1F( "h_cor", "", nbins_, xbins_ );
    TH1F* h_cor74 = new TH1F( "h_cor74", "", nbins_, xbins_ );

    TH1F* h_raw_sigma   = new TH1F( "h_raw_sigma", "", nbins_, xbins_ );
    TH1F* h_cor_sigma   = new TH1F( "h_cor_sigma", "", nbins_, xbins_ );
    TH1F* h_cor74_sigma = new TH1F( "h_cor74_sigma", "", nbins_, xbins_ );

    // Canvasses for datapoints + fit plots
    TCanvas *craw = new TCanvas("craw","craw", perbin_width_ * n_columns_, perbin_height_ * n_rows_ );
    craw->Divide( n_columns_, n_rows_ );
    TCanvas *ccor = new TCanvas("ccor","ccor", perbin_width_ * n_columns_, perbin_height_ * n_rows_ );
    ccor->Divide( n_columns_, n_rows_ );
    TCanvas *ccor74 = new TCanvas("ccor74","ccor74", perbin_width_ * n_columns_, perbin_height_ * n_rows_ );
    ccor74->Divide( n_columns_, n_rows_ );


    //#######################################
    // Fit DSCB in slices; also makes plots per bin
    //#######################################

    // Loop variables
    Double_t x_left, x_right;
    TString x_left_str, x_right_str;
    TString sel_str;

    RooDataSet *hdata_reduced;

    std::vector<Double_t> errors_on_mean_raw;
    std::vector<Double_t> errors_on_mean_cor;
    std::vector<Double_t> errors_on_mean_cor74;

    // Double_t effsigma;


    // RooRealVar &a_mean ;
    // RooRealVar &a_sigma ;
    // Double_t   &a_sigma_eff ;
    
    Double_t mean, error_on_mean, sigma, error_on_sigma, effsigma;

    for(Int_t ibin=0; ibin<nbins_; ibin++) {
        
        x_left = xbins_[ibin]; x_right = xbins_[ibin+1];
        x_left_str.Form(  "%.2f", x_left );
        x_right_str.Form( "%.2f", x_right );

        sel_str = slicevarname_ + ">" + x_left_str + "&&" + slicevarname_ + "<" + x_right_str;

        // Get the dataset for this bin
        hdata_reduced = (RooDataSet*)hdata_->reduce(sel_str);

        FitOneSlice(
            // Input
            hdata_reduced,
            sel_str,
            rawvar_,
            "raw",
            craw,
            ibin,
            // Output
            &mean,
            &error_on_mean,
            &sigma,
            &error_on_sigma,
            &effsigma
            );

        h_raw->AddBinContent( ibin+1, mean );
        h_raw->SetBinError(   ibin+1, effsigma );
        h_raw_sigma->AddBinContent( ibin+1, effsigma        );
        h_raw_sigma->SetBinError(   ibin+1, error_on_sigma  );
        errors_on_mean_raw.push_back( error_on_mean );


        FitOneSlice(
            // Input
            hdata_reduced,
            sel_str,
            ecorvar_,
            "cor",
            ccor,
            ibin,
            // Output
            &mean,
            &error_on_mean,
            &sigma,
            &error_on_sigma,
            &effsigma
            );

        h_cor->AddBinContent( ibin+1, mean );
        h_cor->SetBinError(   ibin+1, effsigma );
        h_cor_sigma->AddBinContent( ibin+1, effsigma        );
        h_cor_sigma->SetBinError(   ibin+1, error_on_sigma  );
        errors_on_mean_cor.push_back( error_on_mean );


        if (draw_old_regression_){

            FitOneSlice(
                // Input
                hdata_reduced,
                sel_str,
                ecor74var_,
                "cor74",
                ccor74,
                ibin,
                // Output
                &mean,
                &error_on_mean,
                &sigma,
                &error_on_sigma,
                &effsigma
                );

            h_cor74->AddBinContent( ibin+1, mean );
            h_cor74->SetBinError(   ibin+1, effsigma );
            h_cor74_sigma->AddBinContent( ibin+1, effsigma        );
            h_cor74_sigma->SetBinError(   ibin+1, error_on_sigma  );
            errors_on_mean_cor74.push_back( error_on_mean );

            }

        }

        // h_raw->AddBinContent( ibin+1, mean_raw.getVal() );
        // h_raw->SetBinError( ibin+1, effsigma );

        // h_raw_sigma->AddBinContent( ibin+1, effsigma           );
        // h_raw_sigma->SetBinError(   ibin+1, sig_raw.getError() );

        // errors_on_mean_raw.push_back( mean_raw.getError() );


        // // ======================================
        // // cor histogram

        // // Reduce further to the range of ecorvar    
        // RooDataSet *hdata_cor = (RooDataSet*)hdata_reduced->reduce(*ecorvar_);

        // // Set DSCB variables with some initial ranges and values
        // RooRealVar mean_cor( "mean_cor","mean_cor",1,0.9,1.1);
        // RooRealVar sig_cor(  "sig_cor","sig_cor",0.01,0.0002,0.8);
        // RooRealVar a1_cor(   "a1_cor","a1_cor",3,0.05,10);
        // RooRealVar a2_cor(   "a2_cor","a2_cor",3,0.05,10);
        // RooRealVar n1_cor(   "n1_cor","n1_cor",3,1.01,500);
        // RooRealVar n2_cor(   "n2_cor","n2_cor",3,1.01,500);

        // // Construct function
        // RooDoubleCBFast pdfCB_cor(
        //     "pdfCB_cor", "pdfCB_cor", *ecorvar_,
        //     mean_cor, sig_cor, a1_cor, n1_cor, a2_cor, n2_cor);

        // pdfCB_cor.fitTo( *hdata_cor, Range( fitxmin_, fitxmax_ ) );
        
        // h_cor->AddBinContent( ibin+1, mean_cor.getVal()   );

        // // Translate pdfCB_raw to a TH1
        // cout << "Turning CB_cor into a histogram (var=" + slicevarname_ + ")" << endl;
        // TH1* h_CB_cor = pdfCB_cor.createHistogram( "h_CB_cor", *ecorvar_, Binning(1000) );

        // cout << "Getting effective sigma for the histogram" << endl;
        // effsigma = effSigma( h_CB_cor );

        // h_cor->SetBinError( ibin+1, effsigma );

        // // Save errors for later use
        // errors_on_mean_cor.push_back( mean_cor.getError() );

        // // h_cor_sigma->AddBinContent( ibin+1, sig_cor.getVal()   );
        // h_cor_sigma->AddBinContent( ibin+1, effsigma           );
        // h_cor_sigma->SetBinError(   ibin+1, sig_cor.getError() );



        // // ======================================
        // // raw histogram

        // RooDataSet *hdata_raw = (RooDataSet*)hdata_reduced->reduce(*rawvar_);

        // RooRealVar mean_raw("mean_raw","mean_raw",1.,0.9,1.1);
        // RooRealVar sig_raw("sig_raw","sig_raw",0.01,0.0002,0.8);
        // RooRealVar a1_raw("a1_raw","a1_raw",3,0.05,10);
        // RooRealVar a2_raw("a2_raw","a2_raw",3,0.05,10);
        // RooRealVar n1_raw("n1_raw","n1_raw",3,1.01,500);
        // RooRealVar n2_raw("n2_raw","n2_raw",3,1.01,500);

        // RooDoubleCBFast pdfCB_raw(
        //     "pdfCB_raw", "pdfCB_raw", *rawvar_, 
        //     mean_raw, sig_raw, a1_raw, n1_raw, a2_raw, n2_raw );

        // pdfCB_raw.fitTo( *hdata_raw, Range( fitxmin_, fitxmax_ ) );

        // h_raw->AddBinContent( ibin+1, mean_raw.getVal() );
        // // h_raw->SetBinError(   ibin+1, mean_raw.getError() );

        // // Translate pdfCB_raw to a TH1
        // cout << "Turning CB_raw into a histogram (var=" + slicevarname_ + ")" << endl;
        // TH1* h_CB_raw = pdfCB_raw.createHistogram( "h_CB_raw", *rawvar_, Binning(1000) );

        // cout << "Getting effective sigma for the histogram" << endl;
        // effsigma = effSigma( h_CB_raw );

        // h_raw->SetBinError( ibin+1, effsigma );

        // // Save errors for later use
        // errors_on_mean_raw.push_back( mean_raw.getError() );

        // // h_raw_sigma->AddBinContent( ibin+1, sig_raw.getVal()   );
        // h_raw_sigma->AddBinContent( ibin+1, effsigma           );
        // h_raw_sigma->SetBinError(   ibin+1, sig_raw.getError() );


        // ======================================
        // Draw slicevar datapoints and fit

        // ccor->cd(ibin+1);
        
        // RooPlot *cor_datapoints = ecorvar_->frame( 0., 2., 250 );
        
        // cor_datapoints->SetTitle("Distribution of cor/true, " + sel_str );
        // cor_datapoints->GetXaxis()->SetTitle("cor/true");
        // gPad->SetGridx();
        // gPad->SetGridy();
        
        // hdata_reduced->plotOn( cor_datapoints, Name("cor_datapoints"), MarkerSize(0.02));
        // pdfCB_cor.plotOn(      cor_datapoints, Name("cor_fit"), LineColor(kRed));
        // cor_datapoints->Draw();

        // TLegend *cor_legend = new TLegend( 0.6, 0.75, 0.9, 0.9 );
        // cor_legend->SetFillStyle(0);
        // cor_legend->SetBorderSize(0);
        // cor_legend->AddEntry( "cor_datapoints", "cor_datapoints", "lep" );
        // cor_legend->AddEntry( "cor_fit", "cor_fit", "l" );
        // cor_legend->Draw();
        

        // craw->cd(ibin+1);

        // RooPlot *raw_datapoints = rawvar_->frame(0.,2, 250);
        
        // raw_datapoints->SetTitle( "Distribution of raw/true " + sel_str );
        // raw_datapoints->GetXaxis()->SetTitle("raw/true");
        // gPad->SetGridx();
        // gPad->SetGridy();
        
        // hdata_reduced->plotOn( raw_datapoints, Name("raw_datapoints"), MarkerSize(0.02));
        // pdfCB_raw.plotOn(      raw_datapoints, Name("raw_fit"), LineColor(kRed));
        // raw_datapoints->Draw();
        
        // TLegend *raw_legend = new TLegend( 0.6, 0.75, 0.9, 0.9 );
        // raw_legend->SetFillStyle(0);
        // raw_legend->SetBorderSize(0);
        // raw_legend->AddEntry( "raw_datapoints", "raw_datapoints", "lep");
        // raw_legend->AddEntry( "raw_fit", "raw_fit", "l");
        // raw_legend->Draw();

        // }

    craw->SaveAs( slicevarname_ + "_PerBinFits_raw.png" );
    ccor->SaveAs( slicevarname_ + "_PerBinFits_cor.png" );

    if (draw_old_regression_) ccor74->SaveAs( slicevarname_ + "_PerBinFits_cor74.png" );

    }

    /*
    //#######################################
    // Make main slice plot
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* c = new TCanvas( "c_" + slicevarname_, "c_" + slicevarname_, c_width_, c_height_ );
    c->cd();

    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);


    // ======================================
    // Draw both profiles ('averages')


    cout << "Getting Profiles for " + slicevarname_ << endl;

    TH2 *h2D_cor, *h2D_raw;
    h2D_cor = hdata_->createHistogram( *slicevar_, *ecorvar_ ); 
    h2D_raw = hdata_->createHistogram( *slicevar_, *rawvar_ );

    // Get the profile, rebinned
    TProfile *ProfX_cor_wrongbins = h2D_cor->ProfileX( "ProfX_raw_" + slicevarname_ + "_wrongbins" , 1, -1, "S" );
    ProfX_cor_wrongbins->Rebin( nbins_, "ProfX_cor_" + slicevarname_ , xbins_ );
    TProfile *ProfX_cor = (TProfile*)gDirectory->Get( "ProfX_cor_" + slicevarname_ );

    TProfile *ProfX_raw_wrongbins = h2D_raw->ProfileX( "ProfX_raw_" + slicevarname_ + "_wrongbins" , 1, -1, "S" );
    ProfX_raw_wrongbins->Rebin( nbins_, "ProfX_raw_" + slicevarname_ , xbins_ );
    TProfile *ProfX_raw = (TProfile*)gDirectory->Get( "ProfX_raw_" + slicevarname_ );


    ProfX_raw->SetMarkerStyle(26);
    ProfX_raw->SetMarkerColor(kRed);
    ProfX_raw->SetLineColor(kRed);
    ProfX_raw->SetName("ProfX_raw");

    ProfX_raw->GetYaxis()->SetRangeUser( ymin_, ymax_ );
    ProfX_raw->GetYaxis()->SetTitle("mean");
    ProfX_raw->GetYaxis()->SetTitleOffset(1.4);

    ProfX_raw->GetXaxis()->SetRangeUser( xbins_[0], xbins_[nbins_] );
    ProfX_raw->GetXaxis()->SetTitle( slicevartitle_ );

    ProfX_raw->SetFillColor(kMagenta);
    ProfX_raw->SetFillStyle(3003);
    
    // ProfX_raw->Draw("HISTSAMEPE3");
    ProfX_raw->Draw("HISTSAMEP");


    ProfX_cor->SetMarkerStyle(26);
    ProfX_cor->SetMarkerColor( kBlue );
    ProfX_cor->SetLineColor(   kBlue );
    ProfX_cor->SetName("ProfX_cor");

    ProfX_cor->SetFillColor(kBlue);
    ProfX_cor->SetFillStyle(3005);
    
    // ProfX_cor->Draw("HISTSAMEPE3");
    ProfX_cor->Draw("HISTSAMEP");


    // ======================================
    // Draw the means of DSCB

    h_raw->SetMarkerStyle(22);
    h_raw->SetMarkerColor(kRed);
    h_raw->SetLineColor(kRed);
    h_raw->SetName("h_raw");
    h_raw->Draw("e p same");

    h_cor->SetMarkerStyle(8);
    h_cor->SetMarkerColor(kBlue);
    h_cor->SetLineColor(kBlue);
    h_cor->SetName("h_cor");
    h_cor->Draw("e p same");    


    // ======================================
    // Legend
    
    TLegend *sliceplot_legend = new TLegend(0.70,0.7,0.87,0.9);
    sliceplot_legend->SetFillStyle(0);
    sliceplot_legend->SetBorderSize(0);
    sliceplot_legend->AddEntry( "h_raw",      "E_{raw}/E_{true}:  #mu_{CB} #pm #sigma_{eff}", "pe" );
    sliceplot_legend->AddEntry( "h_cor",      "E_{corr}/E_{true}:  #mu_{CB} #pm #sigma_{eff}", "pe" );
    // sliceplot_legend->AddEntry( "ProfX_raw",  "E_{raw}/E_{true}:  mean #pm RMS", "pf" );
    // sliceplot_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean #pm RMS", "pf" );
    sliceplot_legend->AddEntry(  "ProfX_raw",  "E_{raw}/E_{true}:  mean", "p" );
    sliceplot_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean", "p" );
    
    sliceplot_legend->Draw("same");

    // Output to file
    c->SaveAs( slicevarname_ + "_Sliced.png" );


    //#######################################
    // Same type of plot, but now with errors on the mean
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* c2 = new TCanvas( "c2_" + slicevarname_, "c2_" + slicevarname_, c_width_, c_height_ );
    c2->cd();

    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);

    // ProfX_raw->Draw("HISTSAMEPE3");
    // ProfX_cor->Draw("HISTSAMEPE3");
    ProfX_raw->Draw("HISTSAMEP");
    ProfX_cor->Draw("HISTSAMEP");

    // Overwrite bin errors with error on mean
    for(Int_t ibin=0; ibin<nbins_; ibin++) {
        h_raw->SetBinError( ibin+1, errors_on_mean_raw[ibin] );
        h_cor->SetBinError( ibin+1, errors_on_mean_cor[ibin] );
        }

    h_raw->Draw("e p same");
    h_cor->Draw("e p same");    

    // ======================================
    // Legend
    
    TLegend *sliceplot_errorsonmean_legend = new TLegend(0.70,0.7,0.87,0.9);
    sliceplot_errorsonmean_legend->SetFillStyle(0);
    sliceplot_errorsonmean_legend->SetBorderSize(0);
    sliceplot_errorsonmean_legend->AddEntry( "h_raw",      "E_{raw}/E_{true}:  #mu_{CB} #pm #sigma_{mean}", "pe" );
    sliceplot_errorsonmean_legend->AddEntry( "h_cor",      "E_{corr}/E_{true}:  #mu_{CB} #pm #sigma_{mean}", "pe" );
    // sliceplot_errorsonmean_legend->AddEntry( "ProfX_raw",  "E_{raw}/E_{true}:  mean #pm RMS", "pf" );
    // sliceplot_errorsonmean_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean #pm RMS", "pf" );
    sliceplot_errorsonmean_legend->AddEntry(  "ProfX_raw",  "E_{raw}/E_{true}:  mean", "p" );
    sliceplot_errorsonmean_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean", "p" );
    
    sliceplot_errorsonmean_legend->Draw("same");

    // Output to file
    c2->SaveAs( slicevarname_ + "_Sliced_errorsonmean.png" );


    //#######################################
    // Same type of plot, now plotting sigma
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* csigma = new TCanvas( "csigma_" + slicevarname_, "csigma_" + slicevarname_, c_width_, c_height_ );
    csigma->cd();
    
    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);

    h_raw_sigma->SetMarkerStyle(22);
    h_raw_sigma->SetMarkerColor(kRed);
    h_raw_sigma->SetLineColor(kRed);
    h_raw_sigma->SetName("h_raw_sigma");
    h_raw_sigma->Draw("e p same");
    
    h_raw_sigma->GetYaxis()->SetRangeUser( ymin_sigma_, ymax_sigma_ );
    h_raw_sigma->GetYaxis()->SetTitle("#sigma_{eff}");
    h_raw_sigma->GetYaxis()->SetTitleOffset(1.4);

    h_raw_sigma->GetXaxis()->SetRangeUser( xbins_[0], xbins_[nbins_] );
    h_raw_sigma->GetXaxis()->SetTitle( slicevartitle_ );

    h_cor_sigma->SetMarkerStyle(8);
    h_cor_sigma->SetMarkerColor(kBlue);
    h_cor_sigma->SetLineColor(kBlue);
    h_cor_sigma->SetName("h_cor_sigma");
    h_cor_sigma->Draw("e p same");    

    // Plot RMSs
    TH1F* h_raw_RMS = new TH1F( "h_raw_RMS", "", nbins_, xbins_ );
    TH1F* h_cor_RMS = new TH1F( "h_cor_RMS", "", nbins_, xbins_ );

    // Set the RMS values in the histograms
    for(Int_t ibin=0; ibin<nbins_; ibin++) {
        h_raw_RMS->AddBinContent( ibin+1, ProfX_raw->GetBinError( ibin+1 )  );
        h_raw_RMS->SetBinError(   ibin+1, 0.00001 );
        h_cor_RMS->AddBinContent( ibin+1, ProfX_cor->GetBinError( ibin+1 )  );
        h_cor_RMS->SetBinError(   ibin+1, 0.00001 );
        }

    h_raw_RMS->SetMarkerStyle(26);
    h_cor_RMS->SetMarkerStyle(26);
    h_raw_RMS->SetMarkerColor(kRed );
    h_raw_RMS->SetLineColor(  kRed );
    h_cor_RMS->SetMarkerColor(kBlue);
    h_cor_RMS->SetLineColor(  kBlue);

    h_raw_RMS->Draw( "SAMEPE" );
    h_cor_RMS->Draw( "SAMEPE" );


    // ======================================
    // Legend
    
    TLegend *sliceplot_sigma_legend = new TLegend(0.70,0.7,0.87,0.9);
    sliceplot_sigma_legend->SetFillStyle(0);
    sliceplot_sigma_legend->SetBorderSize(0);
    sliceplot_sigma_legend->AddEntry( "h_raw_sigma", "E_{raw}/E_{true}:  #sigma_{eff}", "pe" );
    sliceplot_sigma_legend->AddEntry( "h_cor_sigma", "E_{corr}/E_{true}:  #sigma_{eff}", "pe" );
    sliceplot_sigma_legend->AddEntry( "h_raw_RMS",  "E_{raw}/E_{true}:  RMS", "pe" );
    sliceplot_sigma_legend->AddEntry( "h_cor_RMS", "E_{corr}/E_{true}:  RMS", "pe" );
    sliceplot_sigma_legend->Draw("same");

    // Output to file
    csigma->SaveAs( slicevarname_ + "_Sliced_sigma.png" );

    }
    */


//#######################################
// MAIN
//#######################################

void scale_fitMean_RawCor(
    // bool dobarrel=true, bool weighted=false
    ) {

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


    //#######################################
    // Event selection
    //#######################################

    // Cut on event number, should be orthogonal to training cuts
    TCut eventcut = "eventNumber%2==1";
    if ( (string)getenv("TESTRUN")=="Y" ){
        // Use only part of sample
        TCut NtupIDcut = "(NtupID>400&&NtupID<800) || (NtupID>12000&&NtupID<12400)";
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

    //output dir
    //TString dirname = "Barrel_Log_sig3_alpha2-3_evts15/";
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
    ws->Print();
  
    // Read variables from workspace
    RooGBRTargetFlex *meantgt;
    if (dobarrel) 
        meantgt = static_cast<RooGBRTargetFlex*>( ws->arg("sigmeantEB") );
    else
        meantgt = static_cast<RooGBRTargetFlex*>( ws->arg("sigmeantEE") );

    RooRealVar *tgtvar = ws->var("targetvar");


    // ======================================
    // Read simple variables

    //RooRealVar *etrue = new RooRealVar("etrue","etrue",0.);
    // RooRealVar *etrue       = new RooRealVar("genEnergy","genEnergy",0.);
    RooRealVar *scRawEnergy = new RooRealVar( "scRawEnergy", "scRawEnergy", 0.);
    RooRealVar *r9          = new RooRealVar( "scSeedR9", "scSeedR9", 0.);
    RooRealVar *nVtx        = new RooRealVar( "nVtx", "nVtx", 0.);
    RooRealVar *pt          = new RooRealVar( "pt", "pt", 0.);
    RooRealVar *genE        = new RooRealVar( "genEnergy", "genEnergy", 0.);
    RooRealVar *genPt       = new RooRealVar( "genPt", "genPt", 0.);
    //RooRealVar *maxEnergyXtal = new RooRealVar("maxEnergyXtal","maxEnergyXtal",0.);

    RooArgList vars;
    vars.add(meantgt->FuncVars());
    vars.add(*tgtvar);

    vars.add(*genE);
    vars.add(*genPt);
    vars.add(*scRawEnergy);
    vars.add(*r9);
    vars.add(*pt);
    vars.add(*nVtx);
    //vars.add(*maxEnergyXtal);

    // Add the old regression results if particle is electron
    RooRealVar *CorrEcalE  = new RooRealVar( "CorrectedEcalEnergy", "CorrectedEcalEnergy", 0. );
    RooRealVar *CorrEcalEerror = new RooRealVar( "CorrectedEcalEnergyError", "CorrectedEcalEnergyError", 0. );
    if ( isElectron ){
        vars.add(*CorrEcalE);
        vars.add(*CorrEcalEerror);
        }
    

    //create the testing dataset
    RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);
    
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


    r9->setRange(0.,1.2);
    RooRealVar *scetavar = ws->var("var_2");
    //r9->setRange(0.75,1.);
    if(dobarrel)
        scetavar->setRange(-1.5,1.5);
    else
        scetavar->setRange(-3,3);
    
    pt->setRange(0.,300.);
    
    genE->setRange( 0., 1500. );
    genPt->setRange( 0., 300. );


    //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
    //RooFormulaVar ecor("ecor","","1./exp(@0)*exp(@1)",RooArgList(*tgtvar,*sigmeanlim));
    // RooFormulaVar ecor("ecor","","1./exp(@0)*exp(@1)",RooArgList(*tgtvar,*sigmeanlim)); // <--
    RooFormulaVar ecor("ecor","","(@1/@0)",RooArgList(*tgtvar,*sigmeanlim));
    RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
    ecorvar->setRange(0.,2);
    ecorvar->setBins(800);
    
    //formula for raw energy/true energy
    //RooFormulaVar raw("raw","","1./exp(@0)",RooArgList(*tgtvar));
    // RooFormulaVar raw("raw","","1./exp(@0)",RooArgList(*tgtvar));
    RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
    RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
    rawvar->setRange(0.,2.);
    rawvar->setBins(800);

    RooRealVar *ecor74var;
    if ( isElectron ){
        RooFormulaVar ecor74( "ecor74", "@0/@1", RooArgList( *CorrEcalE, *genE ) );
        ecor74var = (RooRealVar*)hdata->addColumn(ecor74);
        ecor74var->setRange(0.,2.);
        ecor74var->setBins(800);
        }


    //#######################################
    // Load into classes and plot
    //#######################################

    // ======================================
    // pt slices

    pt_plot = BinPlot();

    // Load variables into class
    pt_plot.hdata_          = hdata ;
    pt_plot.rawvar_         = rawvar ;
    pt_plot.ecorvar_        = ecorvar ;
    pt_plot.tgtvar_         = tgtvar ;

    if (isElectron){
        pt_plot.draw_old_regression_ = true ;
        pt_plot.ecor74var_ = ecor74var ;
        }

    pt_plot.slicevar_       = pt ;
    pt_plot.slicevartitle_  = "pt"  ;
    pt_plot.slicevarname_  = "pt"  ;

    pt_plot.nbins_          = 15  ;
    pt_plot.binwidth_       = 20.  ;
    pt_plot.binoffset_      = 0.   ;

    // These two number multiplied should be > nbins_
    pt_plot.n_columns_      = 5  ;
    pt_plot.n_rows_         = 3  ;
    
    pt_plot.ymin_           = 0.9 ;
    pt_plot.ymax_           = 1.1 ;

    pt_plot.MakeSlicePlot();


    // ======================================
    // genpt slices

    genpt_plot = BinPlot();

    // Load variables into class
    genpt_plot.hdata_          = hdata ;
    genpt_plot.rawvar_         = rawvar ;
    genpt_plot.ecorvar_        = ecorvar ;
    genpt_plot.tgtvar_         = tgtvar ;

    genpt_plot.slicevar_       = genPt ;
    genpt_plot.slicevartitle_  = "genPt"  ;
    genpt_plot.slicevarname_   = "genPt"  ;

    genpt_plot.nbins_          = 15  ;
    genpt_plot.binwidth_       = 20.  ;
    genpt_plot.binoffset_      = 0.   ;

    // These two number multiplied should be > nbins_
    genpt_plot.n_columns_      = 5  ;
    genpt_plot.n_rows_         = 3  ;
    
    genpt_plot.ymin_           = 0.9 ;
    genpt_plot.ymax_           = 1.1 ;

    genpt_plot.MakeSlicePlot();


    // ======================================
    // r9 slices

    r9_plot = BinPlot();

    // Load variables into class
    r9_plot.hdata_          = hdata ;
    r9_plot.rawvar_         = rawvar ;
    r9_plot.ecorvar_        = ecorvar ;
    r9_plot.tgtvar_         = tgtvar ;

    r9_plot.slicevar_       = r9 ;
    r9_plot.slicevartitle_  = "r9"  ;
    r9_plot.slicevarname_   = "scSeedR9"  ;

    r9_plot.nbins_          = 12    ;
    r9_plot.binwidth_       = 0.02  ;
    r9_plot.binoffset_      = 0.76   ;

    // These two number multiplied should be > nbins_
    r9_plot.n_columns_      = 4  ;
    r9_plot.n_rows_         = 3  ;
    
    r9_plot.ymin_           = 0.9 ;
    r9_plot.ymax_           = 1.1 ;

    r9_plot.MakeSlicePlot();


    // ======================================
    // genE slices

    genE_plot = BinPlot();

    // Load variables into class
    genE_plot.hdata_          = hdata ;
    genE_plot.rawvar_         = rawvar ;
    genE_plot.ecorvar_        = ecorvar ;
    genE_plot.tgtvar_         = tgtvar ;

    genE_plot.slicevar_       = genE ;
    genE_plot.slicevartitle_  = "genEnergy"  ;
    genE_plot.slicevarname_   = "genEnergy"  ;

    genE_plot.nbins_          = 6    ;
    genE_plot.binwidth_       = 150.  ;
    genE_plot.binoffset_      = 0.0   ;

    // These two number multiplied should be > nbins_
    genE_plot.n_columns_      = 3  ;
    genE_plot.n_rows_         = 2  ;
    
    genE_plot.ymin_           = 0.9 ;
    genE_plot.ymax_           = 1.1 ;

    genE_plot.MakeSlicePlot();








    }