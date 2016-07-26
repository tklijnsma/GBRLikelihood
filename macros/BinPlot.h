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

    bool add_extension_to_filename_ = false;
    TString extension_ = "";

    Double_t fitxmin_ = 0.9, fitxmax_ = 1.1;

    // For plotting purposes
    Double_t ymin_ = 0.9, ymax_ = 1.1 ;
    Double_t ymin_sigma_ = 0.0, ymax_sigma_ = 0.1 ;
    Double_t sliceplot_hmargin_ = 0.07;
    Double_t sliceplot_legheight_ = 0.12;

    // Dimensions of main slice canvas
    Int_t c_width_ = 1000, c_height_ = 800;

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
    void Append_filename( TString extension );
    void SaveCanvas( TCanvas *c, TString filename );

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
        }
    }

void BinPlot::Append_filename( TString extension ){
    add_extension_to_filename_ = true;
    extension_ = extension;
    }

void BinPlot::SaveCanvas( TCanvas *c, TString filename ){

    cout << "    Saving " << filename << ".png" ;
    if (draw_pdf_) cout << " (and .pdf)";
    cout << " from canvas" << c->GetName() << endl;

    c->SaveAs( filename + ".png" );
    if (draw_pdf_) c->SaveAs( filename + ".pdf" );
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

    pdfCB.fitTo( *hdata_var, Range( fitxmin_, fitxmax_ ), PrintEvalErrors(-1), PrintLevel(-1) );

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

    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.1);

    RooPlot *var_datapoints = var->frame(0.,2, 250);
    
    // var_datapoints->SetTitle( "Distribution of " + varname + "/true " + sel_str );
    var_datapoints->SetTitle( varname + "/true " + sel_str );
    var_datapoints->SetTitleSize(0.06);

    var_datapoints->GetXaxis()->SetTitle( varname + "/true" );
    var_datapoints->GetXaxis()->SetLimits( 0., 1.3 );

    var_datapoints->GetXaxis()->SetLabelSize(0.05);
    var_datapoints->GetXaxis()->SetTitleSize(0.06);

    var_datapoints->GetYaxis()->SetLabelSize(0.05);
    var_datapoints->GetYaxis()->SetTitleSize(0.06);

    
    hdata_reduced->plotOn( var_datapoints, Name( "datapoints" + varname ), MarkerSize(0.02));
    pdfCB.plotOn( var_datapoints, Name( "fit" + varname ), LineColor(kRed), PrintEvalErrors(-1) );
    var_datapoints->Draw();


    // --------------------------------
    // Labels

    TLatex *l = new TLatex();
    l->SetTextAlign(13);
    l->SetNDC();
    l->SetTextSize(0.05);

    Double_t text_x     = 0.18;
    Double_t textheight = 0.65;
    Double_t nextline   = 0.07;
    Double_t column_shift = 0.1;

    TString TString_effsigma        ;
    TString TString_mean            ;
    TString TString_error_on_mean   ;
    TString TString_sigma           ;
    TString TString_error_on_sigma  ;

    TString_effsigma        .Form( "%.4f", effsigma        );   
    TString_mean            .Form( "%.4f #pm %.4f", mean.getVal(), mean.getError() );   
    TString_sigma           .Form( "%.4f #pm %.4f", sig.getVal(),  sig.getError() );   

    l->DrawLatex( text_x                , textheight, "#sigma_{eff}"         );
    l->DrawLatex( text_x + column_shift , textheight, TString_effsigma       ); textheight -= nextline;
    l->DrawLatex( text_x                , textheight, "#mu_{CB}"             );
    l->DrawLatex( text_x + column_shift , textheight, TString_mean           ); textheight -= nextline;
    l->DrawLatex( text_x                , textheight, "#sigma_{CB}"          );
    l->DrawLatex( text_x + column_shift , textheight, TString_sigma          ); textheight -= nextline;
    
    TLegend *legend = new TLegend( 0.18, 0.75, 0.6, 0.9 );
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

    // craw->SetLeftMargin(0.14);
    // craw->SetBottomMargin(0.14);
    // craw->SetTopMargin(0.02);
    // craw->SetRightMargin(0.02);


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

    // craw->SaveAs( extension_ + slicevarname_ + "_PerBinFits_raw.png" );
    // ccor->SaveAs( extension_ + slicevarname_ + "_PerBinFits_cor.png" );
    // if (draw_old_regression_) ccor74->SaveAs( extension_ + slicevarname_ + "_PerBinFits_cor74.png" );

    // if (draw_pdf_){
    //     craw->SaveAs( extension_ + slicevarname_ + "_PerBinFits_raw.pdf" );
    //     ccor->SaveAs( extension_ + slicevarname_ + "_PerBinFits_cor.pdf" );
    //     if (draw_old_regression_) ccor74->SaveAs( extension_ + slicevarname_ + "_PerBinFits_cor74.pdf" );
    //     }

    SaveCanvas( craw, extension_ + slicevarname_ + "_PerBinFits_raw" );
    SaveCanvas( ccor, extension_ + slicevarname_ + "_PerBinFits_cor" );
    if (draw_old_regression_) SaveCanvas( ccor74, extension_ + slicevarname_ + "_PerBinFits_cor74" );


    //#######################################
    // Make main slice plot
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* c = new TCanvas( "c_" + slicevarname_, "c_" + slicevarname_, c_width_, c_height_ );
    c->cd();

    c->SetLeftMargin(   0.14 );
    c->SetRightMargin(  0.02 );
    c->SetBottomMargin( 0.14 );
    c->SetTopMargin(    0.01 + sliceplot_legheight_ );

    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);


    // ======================================
    // Draw the means of DSCB

    h_raw->SetMarkerStyle(8);
    h_raw->SetMarkerColor(kRed);
    // h_raw->SetMarkerSize(0);
    h_raw->SetLineWidth(2);
    h_raw->SetLineStyle(1);
    h_raw->SetLineColor(kRed);
    h_raw->SetFillColorAlpha(kRed,0.15);
    h_raw->SetName("h_raw");
    h_raw->Draw("HISTE2L");

    // h_cor->SetMarkerSize(0);
    h_cor->SetMarkerStyle(22);
    h_cor->SetMarkerColor(kBlue);
    h_cor->SetLineWidth(2);
    h_cor->SetLineStyle(1);
    h_cor->SetLineColor(kBlue);
    h_cor->SetFillColorAlpha(kBlue,0.15);
    h_cor->SetName("h_cor");
    h_cor->Draw("HISTSAMELE2");

    // h_raw->Draw("e p same");
    // h_cor->Draw("e p same");
    if (draw_old_regression_){
        h_cor74->SetMarkerStyle(22);
        h_cor74->SetMarkerColor(kGreen);
        // h_cor74->SetMarkerSize(0);
        h_cor74->SetLineWidth(2);
        h_cor74->SetLineStyle(1);
        h_cor74->SetLineColor(kGreen);
        h_cor74->SetFillColorAlpha(kGreen,0.15);
        h_cor74->SetName("h_cor74");
        h_cor74->Draw("HISTSAMELE2");
        // h_cor74->Draw("e p same");
        }

    // Set axis properties here
    h_raw->GetYaxis()->SetRangeUser( ymin_, ymax_ );
    h_raw->GetYaxis()->SetTitle("#mu_{CB}");

    h_raw->GetXaxis()->SetRangeUser( xbins_[0], xbins_[nbins_] );
    h_raw->GetXaxis()->SetTitle( slicevartitle_ );

    h_raw->GetYaxis()->SetTitleOffset(1.1);

    h_raw->GetXaxis()->SetLabelSize(0.05);
    h_raw->GetXaxis()->SetTitleSize(0.06);
    h_raw->GetYaxis()->SetLabelSize(0.05);
    h_raw->GetYaxis()->SetTitleSize(0.06);


    // ======================================
    // Draw both profiles ('averages')

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
    // sliceplot_legend->SetBorderSize(0);
    
    sliceplot_legend->AddEntry( "h_raw",      "E_{raw}/E_{true}:  #mu_{CB} #pm #sigma_{eff}", "pfl" );
    sliceplot_legend->AddEntry( "h_cor",      "E_{corr}/E_{true}:  #mu_{CB} #pm #sigma_{eff}", "pfl" );
    
    if (draw_old_regression_) sliceplot_legend->AddEntry( "h_cor74", "E_{corr}/E_{true}:  #mu_{CB} #pm #sigma_{eff} (74X)", "pfl" ) ;

    // sliceplot_legend->AddEntry( "ProfX_raw",  "E_{raw}/E_{true}:  mean #pm RMS", "pf" );
    // sliceplot_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean #pm RMS", "pf" );
    if (draw_meanRMS_) sliceplot_legend->AddEntry(  "ProfX_raw",  "E_{raw}/E_{true}:  mean", "p" );
    if (draw_meanRMS_) sliceplot_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean", "p" );
    
    sliceplot_legend->Draw("same");

    // Output to file

    c->SaveAs( extension_ + slicevarname_ + "_Sliced.png" );
    if (draw_pdf_) c->SaveAs( extension_ + slicevarname_ + "_Sliced.pdf" );


    //#######################################
    // Same type of plot, but now with errors on the mean
    //#######################################

    // Canvas for the main sliced plot
    TCanvas* c2 = new TCanvas( "c2_" + slicevarname_, "c2_" + slicevarname_, c_width_, c_height_ );
    c2->cd();

    c2->SetLeftMargin(   0.14 );
    c2->SetRightMargin(  0.02 );
    c2->SetBottomMargin( 0.14 );
    c2->SetTopMargin(    0.01 + sliceplot_legheight_ );

    gPad->SetGridx();
    gPad->SetGridy();
    // gStyle->SetOptStat(0);

    // if (draw_meanRMS_) ProfX_raw->Draw("HISTSAMEP");
    // if (draw_meanRMS_) ProfX_cor->Draw("HISTSAMEP");

    // Overwrite bin errors with error on mean
    for(Int_t ibin=0; ibin<nbins_; ibin++) {
        h_raw->SetBinError( ibin+1, errors_on_mean_raw[ibin] );
        h_cor->SetBinError( ibin+1, errors_on_mean_cor[ibin] );
        if (draw_old_regression_) h_cor74->SetBinError( ibin+1, errors_on_mean_cor74[ibin] );
        }

    h_raw->Draw("E2L");
    h_cor->Draw("HISTSAMELE2");

    if (draw_old_regression_) h_cor74->Draw("HISTSAMELE2");


    // ======================================
    // Legend
    
    TLegend *sliceplot_errorsonmean_legend = new TLegend( 0.14, 0.99 - sliceplot_legheight_,  1.0-0.02 , 0.99 );
    sliceplot_errorsonmean_legend->SetNColumns(3);
    sliceplot_errorsonmean_legend->SetFillStyle(0);
    // sliceplot_errorsonmean_legend->SetBorderSize(0);
    sliceplot_errorsonmean_legend->AddEntry( "h_raw",      "E_{raw}/E_{true}:  #mu_{CB} #pm #Delta#mu_{CB}", "plf" );
    sliceplot_errorsonmean_legend->AddEntry( "h_cor",      "E_{corr}/E_{true}:  #mu_{CB} #pm #Delta#mu_{CB}", "plf" );
    if (draw_old_regression_) sliceplot_errorsonmean_legend->AddEntry( "h_cor74", "E_{corr}/E_{true}:  #mu_{CB} #pm #Delta#mu_{CB} (74X)", "plf" ) ;

    // sliceplot_errorsonmean_legend->AddEntry( "ProfX_raw",  "E_{raw}/E_{true}:  mean #pm RMS", "pf" );
    // sliceplot_errorsonmean_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean #pm RMS", "pf" );
    if (draw_meanRMS_) sliceplot_errorsonmean_legend->AddEntry(  "ProfX_raw",  "E_{raw}/E_{true}:  mean", "p" );
    if (draw_meanRMS_) sliceplot_errorsonmean_legend->AddEntry( "ProfX_cor",  "E_{corr}/E_{true}:  mean", "p" );
    
    sliceplot_errorsonmean_legend->Draw("same");

    // Output to file
    c2->SaveAs( extension_ + slicevarname_ + "_Sliced_errorsonmean.png" );
    if (draw_pdf_) c2->SaveAs( extension_ + slicevarname_ + "_Sliced_errorsonmean.pdf" );


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
    csigma->SaveAs( extension_ + slicevarname_ + "_Sliced_sigma.png" );
    if (draw_pdf_) csigma->SaveAs( extension_ + slicevarname_ + "_Sliced_sigma.pdf" );

    }