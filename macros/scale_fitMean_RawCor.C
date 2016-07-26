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

    // Draws RMS in sigma plots and means in mean plots if set to true
    bool draw_meanRMS_ = false;
    bool draw_pdf_ = true;

    RooRealVar* slicevar_;
    TString slicevarname_;
    TString slicevartitle_;

    bool Build_xbins_ = true;
    Int_t nbins_;
    Double_t binwidth_;
    Double_t binoffset_;
    Double_t* xbins_ = 0;

    Double_t fitxmin_ = 0.8, fitxmax_ = 1.1;

    // For plotting purposes
    Double_t ymin_ = 0.9, ymax_ = 1.1 ;
    Double_t ymin_sigma_ = 0.0, ymax_sigma_ = 0.1 ;
    Double_t sliceplot_hmargin_ = 0.07;

    // Dimensions of main slice canvas
    Int_t c_width_ = 2000, c_height_ = 800;

    // Dimensions (in pixels) of 1 subplot of the perbin plots
    Int_t perbin_width_ = 1000, perbin_height_ = 800;
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

    void Set_xbins( Int_t nbounds, Double_t xbins[] );

    };


void BinPlot::Set_xbins( Int_t nbounds, Double_t xbins[] ){

    // Instruct the class that xbins_ is filled by the user
    Build_xbins_ = false;

    // Fill in nbins_
    nbins_ = nbounds - 1;
    
    // Copy in the list
    xbins_ = (Double_t*)malloc( sizeof(Double_t)*(nbounds) );
    for ( Int_t ibin = 0; ibin < nbounds; ibin++ ){
        xbins_[ibin] = (Double_t)xbins[ibin];
        cout << xbins_[ibin] << endl;
        }
    }


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
    
    // var_datapoints->SetTitle( "Distribution of " + varname + "/true " + sel_str );
    var_datapoints->SetTitle( varname + "/true " + sel_str );
    var_datapoints->GetXaxis()->SetTitle( varname + "/true" );
    // gPad->SetGridx();
    // gPad->SetGridy();
    
    hdata_reduced->plotOn( var_datapoints, Name( "datapoints" + varname ), MarkerSize(0.02));
    pdfCB.plotOn( var_datapoints, Name( "fit" + varname ), LineColor(kRed));
    var_datapoints->Draw();


    // --------------------------------
    // Labels

    TLatex *l = new TLatex();
    l->SetNDC();
    l->SetTextSize(0.05);

    Double_t text_x     = 0.72;
    Double_t textheight = 0.65;
    Double_t nextline   = 0.06;

    TString TString_effsigma        ;
    TString TString_mean            ;
    TString TString_error_on_mean   ;
    TString TString_sigma           ;
    TString TString_error_on_sigma  ;

    TString_effsigma        .Form( "%.4f", effsigma        );   
    TString_mean            .Form( "%.4f", mean.getVal()   );   
    TString_error_on_mean   .Form( "%.4f", mean.getError() );   
    TString_sigma           .Form( "%.4f", sig.getVal()    );   
    TString_error_on_sigma  .Form( "%.4f", sig.getError()  );   

    l->DrawLatex( text_x - 0.1, textheight, "#sigma_{eff}"         );
    l->DrawLatex( text_x      , textheight, TString_effsigma       ); textheight -= nextline;
    l->DrawLatex( text_x - 0.1, textheight, "#mu_{CB}"             );
    l->DrawLatex( text_x      , textheight, TString_mean           ); textheight -= nextline;
    l->DrawLatex( text_x - 0.1, textheight, "#Delta#mu_{CB}"       );
    l->DrawLatex( text_x      , textheight, TString_error_on_mean  ); textheight -= nextline;
    l->DrawLatex( text_x - 0.1, textheight, "#sigma_{CB}"          );
    l->DrawLatex( text_x      , textheight, TString_sigma          ); textheight -= nextline;
    l->DrawLatex( text_x - 0.1, textheight, "#Delta#sigma_{CB}"    );
    l->DrawLatex( text_x      , textheight, TString_error_on_sigma ); textheight -= nextline;

    
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

    // Make xbins axis if not already built
    if ( Build_xbins_ ){
        xbins_ = (Double_t*)malloc( sizeof(Double_t)*(nbins_+1) );
        for( int ibin = 0; ibin <= nbins_; ibin++ ) {
            xbins_[ibin] = (Double_t)( binoffset_ + ibin * binwidth_ );
            }
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
    
    Double_t mean, error_on_mean, sigma, error_on_sigma, effsigma;

    for(Int_t ibin=0; ibin<nbins_; ibin++) {
        
        x_left = xbins_[ibin]; x_right = xbins_[ibin+1];
        x_left_str.Form(  "%.2f", x_left );
        x_right_str.Form( "%.2f", x_right );

        sel_str = slicevarname_ + ">" + x_left_str + "&&" + slicevarname_ + "<" + x_right_str;

        // Get the dataset for this bin
        hdata_reduced = (RooDataSet*)hdata_->reduce(sel_str);


        // Fitting bins for the raw energy
        FitOneSlice(
            // Input
            hdata_reduced, sel_str, rawvar_, "raw", craw, ibin,
            // Output
            &mean, &error_on_mean, &sigma, &error_on_sigma, &effsigma
            );

        h_raw->AddBinContent( ibin+1, mean );
        h_raw->SetBinError(   ibin+1, effsigma );
        h_raw_sigma->AddBinContent( ibin+1, effsigma        );
        h_raw_sigma->SetBinError(   ibin+1, error_on_sigma  );
        errors_on_mean_raw.push_back( error_on_mean );


        // Fitting bins for the cor energy
        FitOneSlice(
            // Input
            hdata_reduced, sel_str, ecorvar_, "cor", ccor, ibin,
            // Output
            &mean, &error_on_mean, &sigma, &error_on_sigma, &effsigma
            );

        h_cor->AddBinContent( ibin+1, mean );
        h_cor->SetBinError(   ibin+1, effsigma );
        h_cor_sigma->AddBinContent( ibin+1, effsigma        );
        h_cor_sigma->SetBinError(   ibin+1, error_on_sigma  );
        errors_on_mean_cor.push_back( error_on_mean );


        // Fitting bins for the old cor energy
        if (draw_old_regression_){
            FitOneSlice(
                // Input
                hdata_reduced, sel_str, ecor74var_, "cor74", ccor74, ibin,
                // Output
                &mean, &error_on_mean, &sigma, &error_on_sigma, &effsigma
                );

            h_cor74->AddBinContent( ibin+1, mean );
            h_cor74->SetBinError(   ibin+1, effsigma );
            h_cor74_sigma->AddBinContent( ibin+1, effsigma        );
            h_cor74_sigma->SetBinError(   ibin+1, error_on_sigma  );
            errors_on_mean_cor74.push_back( error_on_mean );
            }

        }

    craw->SaveAs( slicevarname_ + "_PerBinFits_raw.png" );
    ccor->SaveAs( slicevarname_ + "_PerBinFits_cor.png" );
    if (draw_old_regression_) ccor74->SaveAs( slicevarname_ + "_PerBinFits_cor74.png" );

    if (draw_pdf_){
        craw->SaveAs( slicevarname_ + "_PerBinFits_raw.pdf" );
        ccor->SaveAs( slicevarname_ + "_PerBinFits_cor.pdf" );
        if (draw_old_regression_) ccor74->SaveAs( slicevarname_ + "_PerBinFits_cor74.pdf" );
        }


    //#######################################
    // Make main slice plot
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* c = new TCanvas( "c_" + slicevarname_, "c_" + slicevarname_, c_width_, c_height_ );
    c->cd();

    c->SetLeftMargin( sliceplot_hmargin_);
    c->SetRightMargin(sliceplot_hmargin_);

    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);


    // ======================================
    // Draw the means of DSCB

    h_raw->SetMarkerStyle(22);
    h_raw->SetMarkerColor(kRed);
    h_raw->SetLineColor(kRed);
    h_raw->SetName("h_raw");
    h_raw->Draw("e p same");

    // Set axis properties here
    h_raw->GetYaxis()->SetRangeUser( ymin_, ymax_ );
    h_raw->GetYaxis()->SetTitle("#mu_{CB}");
    h_raw->GetYaxis()->SetTitleSize(0.043);
    h_raw->GetYaxis()->SetTitleOffset(0.62);

    h_raw->GetXaxis()->SetRangeUser( xbins_[0], xbins_[nbins_] );
    h_raw->GetXaxis()->SetTitle( slicevartitle_ );
    h_raw->GetXaxis()->SetTitleSize(0.043);


    h_cor->SetMarkerStyle(8);
    h_cor->SetMarkerColor(kBlue);
    h_cor->SetLineColor(kBlue);
    h_cor->SetName("h_cor");
    h_cor->Draw("e p same");

    if (draw_old_regression_){
        h_cor74->SetMarkerStyle(8);
        h_cor74->SetMarkerColor(kGreen);
        h_cor74->SetLineColor(kGreen);
        h_cor74->SetName("h_cor74");
        h_cor74->Draw("e p same");
        }


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

    if (draw_meanRMS_){

        ProfX_raw->SetMarkerStyle(26);
        ProfX_raw->SetMarkerColor(kRed);
        ProfX_raw->SetLineColor(kRed);
        ProfX_raw->SetName("ProfX_raw");

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

        }


    // for(Int_t ibin=0; ibin<nbins_; ibin++) {
    //     cout << h_raw->GetBinContent( ibin+1 ) << endl;
    //     cout << h_cor->GetBinContent( ibin+1 ) << endl;
    //     if (draw_old_regression_) cout << h_cor74->GetBinContent( ibin+1 ) << endl;
    //     cout << endl;
    //     }


    // ======================================
    // Legend
    
    TLegend *sliceplot_legend = new TLegend( 0.15, 0.81,  0.85, 0.89 );
    sliceplot_legend->SetNColumns(3);
    sliceplot_legend->SetFillStyle(0);
    sliceplot_legend->SetBorderSize(0);
    
    sliceplot_legend->AddEntry( "h_raw",      "E_{raw}/E_{true}:  #mu_{CB} #pm #sigma_{eff}", "pe" );
    sliceplot_legend->AddEntry( "h_cor",      "E_{corr}/E_{true}:  #mu_{CB} #pm #sigma_{eff}", "pe" );
    
    if (draw_old_regression_) sliceplot_legend->AddEntry( "h_cor74", "E_{corr}/E_{true}:  #mu_{CB} #pm #sigma_{eff} (74X)", "pe" ) ;

    // sliceplot_legend->AddEntry( "ProfX_raw",  "E_{raw}/E_{true}:  mean #pm RMS", "pf" );
    // sliceplot_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean #pm RMS", "pf" );
    if (draw_meanRMS_) sliceplot_legend->AddEntry(  "ProfX_raw",  "E_{raw}/E_{true}:  mean", "p" );
    if (draw_meanRMS_) sliceplot_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean", "p" );
    
    sliceplot_legend->Draw("same");

    // Output to file

    c->SaveAs( slicevarname_ + "_Sliced.png" );
    if (draw_pdf_) c->SaveAs( slicevarname_ + "_Sliced.pdf" );


    //#######################################
    // Same type of plot, but now with errors on the mean
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* c2 = new TCanvas( "c2_" + slicevarname_, "c2_" + slicevarname_, c_width_, c_height_ );
    c2->cd();

    c2->SetLeftMargin( sliceplot_hmargin_);
    c2->SetRightMargin(sliceplot_hmargin_);


    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);

    // ProfX_raw->Draw("HISTSAMEPE3");
    // ProfX_cor->Draw("HISTSAMEPE3");
    if (draw_meanRMS_) ProfX_raw->Draw("HISTSAMEP");
    if (draw_meanRMS_) ProfX_cor->Draw("HISTSAMEP");

    // Overwrite bin errors with error on mean
    for(Int_t ibin=0; ibin<nbins_; ibin++) {
        h_raw->SetBinError( ibin+1, errors_on_mean_raw[ibin] );
        h_cor->SetBinError( ibin+1, errors_on_mean_cor[ibin] );
        if (draw_old_regression_) h_cor74->SetBinError( ibin+1, errors_on_mean_cor74[ibin] );
        }

    h_raw->Draw("e p same");
    h_cor->Draw("e p same");
    if (draw_old_regression_) h_cor74->Draw("e p same");

    // ======================================
    // Legend
    
    TLegend *sliceplot_errorsonmean_legend = new TLegend( 0.15, 0.81,  0.85, 0.89 );
    sliceplot_errorsonmean_legend->SetNColumns(3);
    sliceplot_errorsonmean_legend->SetFillStyle(0);
    sliceplot_errorsonmean_legend->SetBorderSize(0);
    sliceplot_errorsonmean_legend->AddEntry( "h_raw",      "E_{raw}/E_{true}:  #mu_{CB} #pm #Delta#mu_{CB}", "pe" );
    sliceplot_errorsonmean_legend->AddEntry( "h_cor",      "E_{corr}/E_{true}:  #mu_{CB} #pm #Delta#mu_{CB}", "pe" );
    if (draw_old_regression_) sliceplot_errorsonmean_legend->AddEntry( "h_cor74", "E_{corr}/E_{true}:  #mu_{CB} #pm #Delta#mu_{CB} (74X)", "pe" ) ;

    // sliceplot_errorsonmean_legend->AddEntry( "ProfX_raw",  "E_{raw}/E_{true}:  mean #pm RMS", "pf" );
    // sliceplot_errorsonmean_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean #pm RMS", "pf" );
    if (draw_meanRMS_) sliceplot_errorsonmean_legend->AddEntry(  "ProfX_raw",  "E_{raw}/E_{true}:  mean", "p" );
    if (draw_meanRMS_) sliceplot_errorsonmean_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean", "p" );
    
    sliceplot_errorsonmean_legend->Draw("same");

    // Output to file
    c2->SaveAs( slicevarname_ + "_Sliced_errorsonmean.png" );
    if (draw_pdf_) c2->SaveAs( slicevarname_ + "_Sliced_errorsonmean.pdf" );


    //#######################################
    // Same type of plot, now plotting sigma
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* csigma = new TCanvas( "csigma_" + slicevarname_, "csigma_" + slicevarname_, c_width_, c_height_ );
    csigma->cd();

    csigma->SetLeftMargin( sliceplot_hmargin_);
    csigma->SetRightMargin(sliceplot_hmargin_);

    
    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);

    h_raw_sigma->SetMarkerStyle(22);
    h_raw_sigma->SetMarkerColor(kRed);
    h_raw_sigma->SetLineColor(kRed);
    h_raw_sigma->SetName("h_raw_sigma");
    h_raw_sigma->Draw("e p same");
    
    h_raw_sigma->GetYaxis()->SetRangeUser( ymin_sigma_, ymax_sigma_ );
    h_raw_sigma->GetYaxis()->SetTitle( "#sigma_{eff}" );
    h_raw_sigma->GetYaxis()->SetTitleSize(0.043);
    h_raw_sigma->GetYaxis()->SetTitleOffset(0.62);

    h_raw_sigma->GetXaxis()->SetRangeUser( xbins_[0], xbins_[nbins_] );
    h_raw_sigma->GetXaxis()->SetTitle( slicevartitle_ );
    h_raw_sigma->GetXaxis()->SetTitleSize(0.043);

    h_cor_sigma->SetMarkerStyle(8);
    h_cor_sigma->SetMarkerColor(kBlue);
    h_cor_sigma->SetLineColor(kBlue);
    h_cor_sigma->SetName("h_cor_sigma");
    h_cor_sigma->Draw("e p same");    

    if (draw_old_regression_){
        h_cor74_sigma->SetMarkerStyle(8);
        h_cor74_sigma->SetMarkerColor(kGreen);
        h_cor74_sigma->SetLineColor(kGreen);
        h_cor74_sigma->SetName("h_cor74_sigma");
        h_cor74_sigma->Draw("e p same");
        }

    // Plot RMSs
    if (draw_meanRMS_){
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

        }


    // ======================================
    // Legend
    
    TLegend *sliceplot_sigma_legend = new TLegend( 0.15, 0.81,  0.85, 0.89 );
    sliceplot_sigma_legend->SetNColumns(3);
    sliceplot_sigma_legend->SetFillStyle(0);
    sliceplot_sigma_legend->SetBorderSize(0);
    sliceplot_sigma_legend->AddEntry( "h_raw_sigma", "E_{raw}/E_{true}:  #sigma_{eff} #pm #Delta#sigma_{eff}", "pe" );
    sliceplot_sigma_legend->AddEntry( "h_cor_sigma", "E_{corr}/E_{true}:  #sigma_{eff} #pm #Delta#sigma_{eff}", "pe" );
    if (draw_old_regression_) sliceplot_sigma_legend->AddEntry( "h_cor74_sigma", "E_{corr}/E_{true}:  #sigma_{eff} #pm #Delta#sigma_{eff} (74X)", "pe" );
    if (draw_meanRMS_)        sliceplot_sigma_legend->AddEntry( "h_raw_RMS",  "E_{raw}/E_{true}:  RMS", "pe" );
    if (draw_meanRMS_)        sliceplot_sigma_legend->AddEntry( "h_cor_RMS", "E_{corr}/E_{true}:  RMS", "pe" );
    sliceplot_sigma_legend->Draw("same");

    // Output to file
    csigma->SaveAs( slicevarname_ + "_Sliced_sigma.png" );
    if (draw_pdf_) csigma->SaveAs( slicevarname_ + "_Sliced_sigma.pdf" );

    }


//#######################################
// MAIN
//#######################################

void scale_fitMean_RawCor(
    // bool dobarrel=true, bool weighted=false
    ) {

    
    gStyle->SetLineScalePS(1);

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
        // TCut NtupIDcut = "NtupID%200==1";
        // TCut NtupIDcut = "eventNumber%20==1";;
        // if(highpt) NtupIDcut = "eventNumber%10==1" ;
        // else       NtupIDcut = "eventNumber%10==1" ;
        NtupIDcut = "eventNumber%20==0||eventNumber%20==1" ;
        eventcut *= NtupIDcut;
        }

    // Apply additional cut on pt if the up-to-2TeV training is used
    if ( !(getenv("USE2TEVCUT")==NULL) ){
        if ( (string)getenv("USE2TEVCUT")=="Y" ){
            TCut pt2TeV = "genPt<2000" ;
            eventcut *= pt2TeV;
            }
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

    cout << "Reading variables" << endl;

    //RooRealVar *etrue = new RooRealVar("etrue","etrue",0.);
    // RooRealVar *etrue       = new RooRealVar("genEnergy","genEnergy",0.);
    RooRealVar *scRawEnergy = new RooRealVar( "scRawEnergy", "scRawEnergy", 0.);
    // RooRealVar *r9          = new RooRealVar( "scSeedR9", "scSeedR9", 0.);
    RooRealVar *r9          = new RooRealVar( "r9", "r9", 0.);
    RooRealVar *nVtx        = new RooRealVar( "nVtx", "nVtx", 0.);
    RooRealVar *pt          = new RooRealVar( "pt", "pt", 0.);
    RooRealVar *genEta      = new RooRealVar( "genEta", "genEta", 0.);
    RooRealVar *genE        = new RooRealVar( "genEnergy", "genEnergy", 0.);
    RooRealVar *genPt       = new RooRealVar( "genPt", "genPt", 0.);
    //RooRealVar *maxEnergyXtal = new RooRealVar("maxEnergyXtal","maxEnergyXtal",0.);

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
    //vars.add(*maxEnergyXtal);

    // // Add the old regression results if particle is electron
    // RooRealVar *CorrEcalE  = new RooRealVar( "CorrectedEcalEnergy", "CorrectedEcalEnergy", 0. );
    // RooRealVar *CorrEcalEerror = new RooRealVar( "CorrectedEcalEnergyError", "CorrectedEcalEnergyError", 0. );
    // RooRealVar *regression1E =  new RooRealVar( "regression1Energy", "regression1Energy", 0. );
    // RooRealVar *regression1Eerror =  new RooRealVar( "regression1EnergyError", "regression1EnergyError", 0. );

    // if ( isElectron ){
    //     vars.add(*CorrEcalE);
    //     vars.add(*CorrEcalEerror);
    //     }
    // else if ( isPhoton ){
    //     vars.add( *regression1E );
    //     vars.add( *regression1Eerror );
    //     }

    // Add the old regression results if particle is electron
    RooRealVar *cor74E  = new RooRealVar(     "corrEnergy74X", "corrEnergy74X", 0. );
    RooRealVar *cor74Eerror = new RooRealVar( "corrEnergy74XError", "corrEnergy74XError", 0. );
    vars.add( *cor74E );
    vars.add( *cor74Eerror );


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


    // ======================================
    // Set variables ranges -- Should be big enough!!!
    
    // RooRealVar *scetavar = ws->var("var_2");
    //r9->setRange(0.75,1.);
    // if(dobarrel)
    //     scetavar->setRange(-1.5,1.5);
    // else
    //     scetavar->setRange(-3,3);
    // pt->setRange(0.,300.);
    

    r9     ->setRange( 0., 1.2);

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
    // if ( isElectron ){
    //     RooFormulaVar ecor74( "ecor74", "@0/@1", RooArgList( *CorrEcalE, *genE ) );
    //     ecor74var = (RooRealVar*)hdata->addColumn(ecor74);
    //     ecor74var->setRange(0.,2.);
    //     ecor74var->setBins(800);
    //     }
    // else if ( isPhoton ){
    //     RooFormulaVar ecor74( "ecor74", "@0/@1", RooArgList( *regression1E, *genE ) );
    //     ecor74var = (RooRealVar*)hdata->addColumn(ecor74);
    //     ecor74var->setRange(0.,2.);
    //     ecor74var->setBins(800);
    //     }
    RooFormulaVar ecor74( "ecor74", "@0/@1", RooArgList( *cor74E, *genE ) );
    ecor74var = (RooRealVar*)hdata->addColumn(ecor74);
    ecor74var->setRange(0.,2.);
    ecor74var->setBins(800);



    cout << "Finished reading variables" << endl;


    //#######################################
    // Load into classes and plot
    //#######################################

    Double_t global_ymin = 0.95;
    Double_t global_ymax = 1.03;


    if (!highpt) {

        // // ======================================
        // // pt slices

        // BinPlot L_pt_plot = BinPlot();

        // const Int_t L_pt_nbounds = 21;
        // Double_t L_pt_xbins[L_pt_nbounds] = {
        //     0.,   10.,  20.,  30.,  40.,
        //     50.,  60.,  70.,  80.,  90.,
        //     100., 120., 140., 160., 180.,
        //     200., 220., 240., 260., 280.,
        //     300.
        //     };
        // L_pt_plot.Set_xbins( L_pt_nbounds , L_pt_xbins );

        // // Load variables into class
        // L_pt_plot.hdata_          = hdata ;
        // L_pt_plot.rawvar_         = rawvar ;
        // L_pt_plot.ecorvar_        = ecorvar ;
        // L_pt_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
        //     L_pt_plot.draw_old_regression_ = true ;
        //     L_pt_plot.ecor74var_ = ecor74var ;
        //     }

        // L_pt_plot.slicevar_       = pt ;
        // L_pt_plot.slicevartitle_  = "pt"  ;
        // L_pt_plot.slicevarname_   = "pt"  ;

        // // These two number multiplied should be > nbins_
        // L_pt_plot.n_columns_      = 5  ;
        // L_pt_plot.n_rows_         = 5  ;
        
        // L_pt_plot.ymin_           = global_ymin ;
        // L_pt_plot.ymax_           = global_ymax ;

        // L_pt_plot.MakeSlicePlot();


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


    else if (highpt) {

        // // ======================================
        // // pt slices

        // BinPlot H_pt_plot = BinPlot();

        // const Int_t H_pt_nbounds = 14;
        // Double_t H_pt_xbins[H_pt_nbounds] = {
        //     0.,    20.,   40.,   60.,   80.,
        //     100.,  150.,  200.,  300.,  500.,
        //     700.,  1000., 1500., 2000.
        //     };
        // H_pt_plot.Set_xbins( H_pt_nbounds , H_pt_xbins );

        // // Load variables into class
        // H_pt_plot.hdata_          = hdata ;
        // H_pt_plot.rawvar_         = rawvar ;
        // H_pt_plot.ecorvar_        = ecorvar ;
        // H_pt_plot.tgtvar_         = tgtvar ;

        // if (isElectron){
        //     H_pt_plot.draw_old_regression_ = true ;
        //     H_pt_plot.ecor74var_ = ecor74var ;
        //     }

        // H_pt_plot.slicevar_       = pt ;
        // H_pt_plot.slicevartitle_  = "pt"  ;
        // H_pt_plot.slicevarname_   = "pt"  ;

        // // These two number multiplied should be > nbins_
        // H_pt_plot.n_columns_      = 5  ;
        // H_pt_plot.n_rows_         = 5  ;
        
        // H_pt_plot.ymin_           = global_ymin ;
        // H_pt_plot.ymax_           = global_ymax ;

        // H_pt_plot.MakeSlicePlot();


        // ======================================
        // genpt slices

        BinPlot H_genpt_plot = BinPlot();

        // const Int_t H_genpt_nbounds = 14;
        // Double_t H_genpt_xbins[H_genpt_nbounds] = {
        //     0.,    20.,   40.,   60.,   80.,
        //     100.,  150.,  200.,  300.,  500.,
        //     700.,  1000., 1500., 2000.
        //     };

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
        // genE slices

        // BinPlot H_genE_plot = BinPlot();

        // const Int_t H_genE_nbounds = 14;
        // Double_t H_genE_xbins[H_genE_nbounds] = {
        //     0.,    50.,   100.,  150.,  200.,
        //     250.,  300.,  500.,  750.,  1000.,
        //     1250., 1500., 2000., 3000.
        //     };
        // H_genE_plot.Set_xbins( H_genE_nbounds , H_genE_xbins );

        // // Load variables into class
        // H_genE_plot.hdata_          = hdata ;
        // H_genE_plot.rawvar_         = rawvar ;
        // H_genE_plot.ecorvar_        = ecorvar ;
        // H_genE_plot.tgtvar_         = tgtvar ;

        // // if (isElectron){
        //     H_genE_plot.draw_old_regression_ = true ;
        //     H_genE_plot.ecor74var_ = ecor74var ;
        //     // }

        // H_genE_plot.slicevar_       =  genE ;
        // H_genE_plot.slicevartitle_  = "genEnergy"  ;
        // H_genE_plot.slicevarname_   = "genEnergy"  ;

        // // These two number multiplied should be > nbins_
        // H_genE_plot.n_columns_      = 4  ;
        // H_genE_plot.n_rows_         = 4  ;
        
        // H_genE_plot.ymin_           = global_ymin ;
        // H_genE_plot.ymax_           = global_ymax ;

        // H_genE_plot.MakeSlicePlot();


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





    }