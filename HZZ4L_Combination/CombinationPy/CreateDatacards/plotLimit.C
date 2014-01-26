#include <iostream>
#include <fstream>
using namespace std;



void plotLimit(TString file, TString POI,TString title=""){
/*
  title = gSystem->Getenv("PLOT_TAG");
  std::cout <<"PLOT_TAG = " <<title << std::endl;*/
//   cout << Form("The 3.001 and 3.001 are the same numbers = %d", int(AreSame(3.001,3.001))) << endl;
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  TString parameter = POI;
//   TString parameter="fa3";
  
//   TFile *f=new TFile(Form("higgsCombine1D.MultiDimFit.mH125.7.%s.root",parameter.Data()),"READ");
  
  
  TFile *f=new TFile(file,"READ"); 
//   TFile *f=new TFile("higgsCombine1D.MultiDimFit.mH125.7.k31.10000points.root","READ");
  TTree *t=(TTree*)f->Get("limit");
//   t->Draw("2*deltaNLL:x", "2*deltaNLL <= 3.84");
  t->Draw(Form("2*deltaNLL:%s",POI.Data()));
  TGraph *gr_xnll=(TGraph*) gROOT->FindObject("Graph")->Clone();
//   gr_xnll->SetMarkerStyle(34); 
  gr_xnll->SetMarkerSize(2.0);


  

  TCanvas *c1=new TCanvas("cc1","CANVAS1",600,600);
  c1->cd();
  gPad->SetRightMargin(0.0085);
  TString ytitle;
  if (parameter=="k3k1_ratio") ytitle = "k_{3}/k_{1}";
  if (parameter=="k2k1_ratio") ytitle = "k_{2}/k_{1}";
  
//   TString ytitle = "|k_{3}/k_{1}|" ;
  gr_xnll->GetXaxis()->SetTitle(ytitle);
  gr_xnll->GetYaxis()->SetTitle("-2 #Delta ln L");
  gr_xnll->GetXaxis()->SetLabelSize(0.04);
  gr_xnll->SetMarkerStyle(20);
  gr_xnll->SetMarkerSize(.5);
  
  gr_xnll->GetYaxis()->SetTitleOffset(1.3);
  gr_xnll->GetXaxis()->SetTitleOffset(1.3);
  gr_xnll->SetTitle("");
  
  gr_xnll->GetYaxis()->SetRangeUser(0,5);
  
  FitResult_t fitres = getFitResults(t,POI);
//   FitResult_t fitres_v0 = getFitResults_v0(t,POI);
  
  Double_t from_fit[] = {fitres.bf_global,fitres.ul68,fitres.ul95};
  
  
//   cout << Form("@@@@@ FitResult: bf_glob=%3.2f | bf_scan=%3.2f | wf=%3.2f | ll68=%3.2f | ll95=%3.2f | ul68=%3.2f | ul95=%3.2f \n", 
//              fitres.bf_global,fitres.bf,fitres.wf,fitres.ll68,fitres.ll95,fitres.ul68,fitres.ul95); 
  
  
  
  ofstream myfile;
  myfile.open (Form("%s.limits",file.Data()));
  myfile << Form("BF:%f\n",fitres.bf_global);
  myfile << Form("WF:%f\n",fitres.wf);
  myfile << Form("LL68:%f\n",fitres.ll68);
  myfile << Form("LL95:%f\n",fitres.ll95);
  myfile << Form("UL68:%f\n",fitres.ul68);
  myfile << Form("UL95:%f\n",fitres.ul95);
  myfile.close();
  
//   if (POI != "k3k1_ratio") {
//     gr_xnll->Draw("AP");
//     c1->SaveAs(Form("%s.png",file.Data()));
//   }
  
//
  
  
  if (POI == "k3k1_ratio") {
    
    gr_xnll->Draw("AP");
    putInformationOnPlot(fitres,c1, gr_xnll->GetXaxis(), gr_xnll->GetYaxis());
    putLimitLinesOnPlot(&fitres, c1);
    c1->SaveAs(Form("%s.AxisK3K1.png",file.Data()));
    c1->SaveAs(Form("%s.AxisK3K1.root",file.Data()));
    
    TCanvas *c2=new TCanvas("cc2","CANVAS2",600,600);
    c2->cd();
    gPad->SetRightMargin(0.0085);
    
    TGraph * gr_xnll_arctan = getRescaledGraph(gr_xnll,0.03676929746947648);
    gr_xnll_arctan->GetXaxis()->SetTitle("1/#pi arctan(#sqrt{#gamma_{33}} k_{3}/k_{1})");
    gr_xnll_arctan->GetYaxis()->SetTitle("-2 #Delta ln L");
    gr_xnll_arctan->GetYaxis()->SetTitleOffset(1.3);
    gr_xnll_arctan->GetXaxis()->SetTitleOffset(1.3);
    gr_xnll_arctan->SetTitle("");
    gr_xnll_arctan->GetXaxis()->SetLabelSize(0.04);
    gr_xnll_arctan->SetMarkerStyle(20);
    gr_xnll_arctan->SetMarkerSize(.5);
    gr_xnll_arctan->GetYaxis()->SetRangeUser(0,5);
    gr_xnll_arctan->Draw("AP");
//     max_y = gr_xnll_arctan->GetYaxis()->GetXmax();
//     min_x = gr_xnll_arctan->GetXaxis()->GetXmin();
//     latex.DrawLatex(min_x, max_y,fit_info);
//     putInformationOnPlot(t, POI, c2, gr_xnll_arctan->GetXaxis(), gr_xnll_arctan->GetYaxis());
    putInformationOnPlot(fitres, c2, gr_xnll_arctan->GetXaxis(), gr_xnll_arctan->GetYaxis());
    putLimitLinesOnPlot(&getRescaledFitResult(fitres,0.03676929746947648), c2);
    
    c2->SaveAs(Form("%s.AxisArctanGamK3K1.png",file.Data()));

    
    
    
    TCanvas *c3=new TCanvas("cc3","CANVAS3",600,600);
    c3->cd();
    gPad->SetRightMargin(0.0085);
    
    TGraph * gr_xnll_arctan_2 = getRescaledGraph(gr_xnll,1);
    gr_xnll_arctan_2->GetXaxis()->SetTitle("1/#pi arctan(k_{3}/k_{1})");
    gr_xnll_arctan_2->GetYaxis()->SetTitle("-2 #Delta ln L");
    gr_xnll_arctan_2->GetYaxis()->SetTitleOffset(1.3);
    gr_xnll_arctan_2->GetXaxis()->SetTitleOffset(1.3);
    gr_xnll_arctan_2->SetTitle("");
    gr_xnll_arctan_2->GetXaxis()->SetLabelSize(0.04);
    gr_xnll_arctan_2->SetMarkerStyle(20);
    gr_xnll_arctan_2->SetMarkerSize(.5);
    gr_xnll_arctan_2->GetYaxis()->SetRangeUser(0,5);
    gr_xnll_arctan_2->Draw("AP");
    putInformationOnPlot(fitres, c3, gr_xnll_arctan_2->GetXaxis(), gr_xnll_arctan_2->GetYaxis());
    putLimitLinesOnPlot(&getRescaledFitResult(fitres,1), c3);
    c3->SaveAs(Form("%s.AxisArctanK3K1.png",file.Data()));
    
    
    //fa3 result
    TCanvas *c4=new TCanvas("cc4","CANVAS4",600,600);
    c4->cd();
    gPad->SetRightMargin(0.0085);
    
    TGraph * gr_xnll_fa3 = getfa3RescaledGraph(gr_xnll);
    gr_xnll_fa3->GetXaxis()->SetTitle("f_{a3}");
    gr_xnll_fa3->GetYaxis()->SetTitle("-2 #Delta ln L");
    gr_xnll_fa3->GetYaxis()->SetTitleOffset(1.3);
    gr_xnll_fa3->GetXaxis()->SetTitleOffset(1.3);
    gr_xnll_fa3->SetTitle("");
    gr_xnll_fa3->GetXaxis()->SetLabelSize(0.04);
    gr_xnll_fa3->SetMarkerStyle(20);
    gr_xnll_fa3->SetMarkerSize(.5);
    gr_xnll_fa3->GetYaxis()->SetRangeUser(0,5);
    gr_xnll_fa3->Draw("AP");
//     max_y = gr_xnll_fa3->GetYaxis()->GetXmax();
// //     min_x = gr_xnll_fa3->GetXaxis()->GetXmin();
//     latex.DrawLatex(min_x, max_y,fit_info);
//     putInformationOnPlot(t, POI, c4, gr_xnll_fa3->GetXaxis(), gr_xnll_fa3->GetYaxis());
    putInformationOnPlot(getfa3RescaledFitResult(fitres), c4, gr_xnll_fa3->GetXaxis(), gr_xnll_fa3->GetYaxis(),"f_{a3}");
    putLimitLinesOnPlot(&getfa3RescaledFitResult(fitres), c4);
    c4->SaveAs(Form("%s.AxisFa3.png",file.Data()));
  }
  
  
  if (POI == "k2k1_ratio") {
    
    gr_xnll->Draw("AP");
    putInformationOnPlot(fitres,c1, gr_xnll->GetXaxis(), gr_xnll->GetYaxis());
    putLimitLinesOnPlot(&fitres, c1);
    c1->SaveAs(Form("%s.AxisK2K1.png",file.Data()));
    c1->SaveAs(Form("%s.AxisK2K1.root",file.Data()));
  }
  
  
  
}




Bool_t AreSame(Double_t a, Double_t b) {
//     return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
//     return fabs(a - b) < EPSILON;
    return fabs(a - b) < 1.0E-8;
};

struct FitResult_t{
  Double_t bf_global;
  Double_t bf;
  Double_t wf;
  Double_t ll68;
  Double_t ll95;
  Double_t ul68;
  Double_t ul95;
} ;

FitResult_t default_fr = {-99999,-99999,-99999,99999,99999,-99999,-99999 };

FitResult_t getFitResults(TTree *t, TString POI){
//   1) get best/worst fit
//   2) starting to left/right from best fit value get lower/upeer limits at 68(95)%
//       - check if the value is not found (there should be one value greater than NLL of limit or the value is equal to NLL limit) and set to infinity
   FitResult_t fitres = {-99999,-99999,-99999,99999,99999,-99999,-99999 };
   FitResult_t fitres_entry = {-99,-99,-99,-99,-99,-99,-99 };  //here we will have integer entries to fitres
//    Double_t ret=-9999;
   Long64_t nentries = t->GetEntriesFast();
   Float_t x, deltaNLL,quantileExpected;
   
   t->SetBranchAddress("quantileExpected", &quantileExpected);
   t->SetBranchAddress(POI.Data(), &x);
   t->SetBranchAddress("deltaNLL", &deltaNLL);
//    std::cout <<"Getting limits"<< std::endl;
   
   
   //find two closest x points to the quantileExpected
   Long64_t nbytes = 0, nb = 0;

   //find best/worst fit
   Double_t lowestDeltaNLL=99999;
   Double_t highestDeltaNLL=-99999;
   nb = t->GetEntry(0);   //nbytes += nb;
   fitres.bf_global  = x;
   for (Long64_t jentry=1; jentry<nentries;jentry++) {
      nb = t->GetEntry(jentry);   //nbytes += nb;
      //best fit
      if ( deltaNLL < lowestDeltaNLL) {
        lowestDeltaNLL = deltaNLL;
        fitres.bf  = x;
        fitres_entry.bf  = jentry;
      }
      //worst fit
      if ( deltaNLL > highestDeltaNLL) {
        highestDeltaNLL = deltaNLL;
        fitres.wf  = x;
        fitres_entry.wf  = jentry;
      }
   }
 
   //find upper limits
   Float_t highestDeltaNLL_right=-99999;
   Float_t min_diff_deltaNLL = 0.1;
   
   TGraph *ul_inverse = new TGraph();
   Int_t inverse_i=0; 
   for (Long64_t jentry=int(fitres_entry.bf); jentry<nentries;jentry++) {
       
      nb = t->GetEntry(jentry);   
      
      if (fabs(2*deltaNLL - 3.84)<0.5) {
          ul_inverse->SetPoint(inverse_i, 2*deltaNLL, x);    
//           cout << Form("2*dNLL=%4.3f x=%4.3f Eval=%4.3f minx=%4.3f maxx=%4.3f", 2*deltaNLL, x, ul_inverse->Eval(3.84),ul_inverse->GetXaxis()->GetXmin(),ul_inverse->GetXaxis()->GetXmax() )<< endl;
          inverse_i++;
      }
      
      if (fabs(2*deltaNLL - 1)<0.5) {
          ul_inverse->SetPoint(inverse_i, 2*deltaNLL, x);    
//           cout << Form("2*dNLL=%4.3f x=%4.3f Eval=%4.3f minx=%4.3f maxx=%4.3f", 2*deltaNLL, x, ul_inverse->Eval(1),ul_inverse->GetXaxis()->GetXmin(),ul_inverse->GetXaxis()->GetXmax() )<< endl;
          inverse_i++;
      }
   }
   
//    Double_t ul_inverse_val = ul_inverse->Eval(3.84);
   if (ul_inverse->GetN() > 1 && 
       3.84 < ul_inverse->GetXaxis()->GetXmax() && 
       3.84 > ul_inverse->GetXaxis()->GetXmin()) fitres.ul95 = ul_inverse->Eval(3.84);

   if (ul_inverse->GetN() > 1 && 
       1 < ul_inverse->GetXaxis()->GetXmax() && 
       1 > ul_inverse->GetXaxis()->GetXmin()) fitres.ul68 = ul_inverse->Eval(1);

   
   
   //find lower limits
   Float_t highestDeltaNLL_left=-99999;
   inverse_i=0;
   TGraph *ll_inverse = new TGraph();
   for (Long64_t jentry=int(fitres_entry.bf); jentry>0;jentry--) {
      nb = t->GetEntry(jentry);   
      if (fabs(2*deltaNLL - 3.84)<0.5) {
          ll_inverse->SetPoint(inverse_i, 2*deltaNLL, x);    
//           cout << Form("2*dNLL=%4.3f x=%4.3f Eval=%4.3f minx=%4.3f maxx=%4.3f", 2*deltaNLL, x, ll_inverse->Eval(3.84),ll_inverse->GetXaxis()->GetXmin(),ll_inverse->GetXaxis()->GetXmax() )<< endl;
          inverse_i++;
      }
      
      if (fabs(2*deltaNLL - 1)<0.5) {
          ll_inverse->SetPoint(inverse_i, 2*deltaNLL, x);    
//           cout << Form("2*dNLL=%4.3f x=%4.3f Eval=%4.3f minx=%4.3f maxx=%4.3f", 2*deltaNLL, x, ll_inverse->Eval(1),ll_inverse->GetXaxis()->GetXmin(),ll_inverse->GetXaxis()->GetXmax() )<< endl;
          inverse_i++;
      }
   }
   
//    Double_t ll_inverse_val = ll_inverse->Eval(3.84);
   if (ll_inverse->GetN() > 1 && 
       3.84 < ll_inverse->GetXaxis()->GetXmax() && 
       3.84 > ll_inverse->GetXaxis()->GetXmin()) fitres.ll95 = ll_inverse->Eval(3.84);

   if (ll_inverse->GetN() > 1 && 
       1 < ll_inverse->GetXaxis()->GetXmax() && 
       1 > ll_inverse->GetXaxis()->GetXmin()) fitres.ll68 = ll_inverse->Eval(1);

//     set to beginning of interval
   nb = t->GetEntry(1);   
//    
   cout << Form("@@@@@ NEW FitResult in getFitResults: bf_glob=%3.2f | bf_scan=%3.2f | wf=%3.2f | ll68=%3.2f | ll95=%3.2f | ul68=%3.2f | ul95=%3.2f \n", 
                fitres.bf_global,fitres.bf,fitres.wf,fitres.ll68,fitres.ll95,fitres.ul68,fitres.ul95);
   
   return fitres;
  
}


FitResult_t getFitResults_v0(TTree *t, TString POI){
//   1) get best/worst fit
//   2) starting to left/right from best fit value get lower/upeer limits at 68(95)%
//       - check if the value is not found (there should be one value greater than NLL of limit or the value is equal to NLL limit) and set to infinity
   FitResult_t fitres = {-99999,-99999,-99999,99999,99999,-99999,-99999 };
   FitResult_t fitres_entry = {-99,-99,-99,-99,-99,-99,-99 };  //here we will have integer entries to fitres
//    Double_t ret=-9999;
   Long64_t nentries = t->GetEntriesFast();
   Float_t x, deltaNLL,quantileExpected;
   
   t->SetBranchAddress("quantileExpected", &quantileExpected);
   t->SetBranchAddress(POI.Data(), &x);
   t->SetBranchAddress("deltaNLL", &deltaNLL);
//    std::cout <<"Getting limits"<< std::endl;
   
   
   //find two closest x points to the quantileExpected
   Long64_t nbytes = 0, nb = 0;

   //find best/worst fit
   Double_t lowestDeltaNLL=99999;
   Double_t highestDeltaNLL=-99999;
   nb = t->GetEntry(0);   //nbytes += nb;
   fitres.bf_global  = x;
   for (Long64_t jentry=1; jentry<nentries;jentry++) {
      nb = t->GetEntry(jentry);   //nbytes += nb;
      //best fit
      if ( deltaNLL < lowestDeltaNLL) {
	lowestDeltaNLL = deltaNLL;
	fitres.bf  = x;
	fitres_entry.bf  = jentry;
      }
      //worst fit
      if ( deltaNLL > highestDeltaNLL) {
	highestDeltaNLL = deltaNLL;
	fitres.wf  = x;
	fitres_entry.wf  = jentry;
      }
   }
 
   //find upper limits
   Float_t distQ68=99999.;
   Float_t distQ95=99999.;
   Float_t highestDeltaNLL_right=-99999;
   Float_t min_diff_deltaNLL = 0.1;
   
   for (Long64_t jentry=int(fitres_entry.bf); jentry<nentries;jentry++) {
       
      nb = t->GetEntry(jentry);   
      //find entry with Q closest to one provided
      if ( fabs(quantileExpected-0.32) > distQ68 && 2*deltaNLL > 1) distQ68=0; //stop condition
      if ( fabs(quantileExpected-0.32) < distQ68 ) {
// 	closestEntry = jentry;
	distQ68 = fabs(quantileExpected-0.32);
        if (fabs(2*deltaNLL - 1)<min_diff_deltaNLL){
            fitres.ul68=x;
            fitres_entry.ul68=jentry;
        }
	
      }
      
      if ( fabs(quantileExpected-0.05) > distQ95 && 2*deltaNLL > 3.84) distQ95=0; //stop condition
      if ( fabs(quantileExpected-0.05) < distQ95 ) {
// 	closestEntry = jentry;
	distQ95 = fabs(quantileExpected-0.05);
        if (fabs(2*deltaNLL - 3.84)<min_diff_deltaNLL){
            fitres.ul95=x;
            fitres_entry.ul95=jentry;
        }
      }
      
      if ( deltaNLL > highestDeltaNLL_right) {
	highestDeltaNLL_right = deltaNLL;
      }
   }
   
   
   
   //find lower limits
   distQ68=99999;
   distQ95=99999;
   Float_t highestDeltaNLL_left=-99999;
   for (Long64_t jentry=int(fitres_entry.bf); jentry>0;jentry--) {
      nb = t->GetEntry(jentry);   
      //find entry with Q closest to one provided
      if ( fabs(quantileExpected-0.32) > distQ68 && 2*deltaNLL > 1) distQ68=0; //stop condition
      if ( fabs(quantileExpected-0.32) < distQ68 ) {
	distQ68 = fabs(quantileExpected-0.32);
        if (fabs(2*deltaNLL - 1)<min_diff_deltaNLL){
            fitres.ll68=x;
            fitres_entry.ll68=jentry;
        }
      }
      
      if ( fabs(quantileExpected-0.05) > distQ95 && 2*deltaNLL > 3.84) distQ95=0; //stop condition
      if ( fabs(quantileExpected-0.05) < distQ95 ) {
// 	closestEntry = jentry;
	distQ95 = fabs(quantileExpected-0.05);
        if (fabs(2*deltaNLL - 3.84)<min_diff_deltaNLL){
            fitres.ll95=x;
            fitres_entry.ll95=jentry;
        }
      }
      
      if ( deltaNLL > highestDeltaNLL_left) {
	highestDeltaNLL_left = deltaNLL;
      }
   }
//     set to beginning of interval
   nb = t->GetEntry(1);   
   if (2*highestDeltaNLL_left < 3.84 && fabs(2*highestDeltaNLL_left - 3.84) < min_diff_deltaNLL) {
       fitres.ll95=x;
//        if (fabs(2*highestDeltaNLL_left - 3.84) < min_diff_deltaNLL) 
   }
   if (2*highestDeltaNLL_left < 1    && fabs(2*highestDeltaNLL_left - 1) < min_diff_deltaNLL ) fitres.ll68=x;
   
   //set to end of interval
   nb = t->GetEntry(nentries-1);     
   if (2*highestDeltaNLL_right < 3.84  && fabs(2*highestDeltaNLL_right - 3.84) < min_diff_deltaNLL) fitres.ul95=x;
   if (2*highestDeltaNLL_right < 1     && fabs(2*highestDeltaNLL_right - 1) < min_diff_deltaNLL) fitres.ul68=x;
   cout << Form("Highest NLL_L=%3.2f NLL_R=%3.2f\n",2*highestDeltaNLL_left,2*highestDeltaNLL_right);
   
   cout << Form("@@@@@ OLD FitResult in getFitResults: bf_glob=%3.2f | bf_scan=%3.2f | wf=%3.2f | ll68=%3.2f | ll95=%3.2f | ul68=%3.2f | ul95=%3.2f \n", 
		fitres.bf_global,fitres.bf,fitres.wf,fitres.ll68,fitres.ll95,fitres.ul68,fitres.ul95);
   
   //check left and right ... to get the region
   return fitres;
  
}


void putLimitLinesOnPlot(FitResult_t * fitresult, TCanvas * c){
    
//   //1 sigma and 2 sigma lines
  c->cd();
  TLine line;
  FitResult_t fitres = *fitresult;
  Double_t bestFit = fitres.bf_global;
  TArrow bf_arr;

  bf_arr.SetLineColor(kRed);
  bf_arr.SetLineWidth(2);
  bf_arr.DrawArrow(bestFit,-0.5,bestFit,-0.25,0.02,"|->");
  
  
  Double_t ulimit68 = fitres.ul68;
  Double_t ulimit95 = fitres.ul95;
  Double_t llimit68 = fitres.ll68;
  Double_t llimit95 = fitres.ll95;
  
  Bool_t doLL95=0, doLL68=0, doUL68=0, doUL95=0;
  doLL95=!(AreSame(llimit95,default_fr.ll95));
  doLL68=!(AreSame(llimit68,default_fr.ll68));
  doUL68=!(AreSame(ulimit68,default_fr.ul68));
  doUL95=!(AreSame(ulimit95,default_fr.ul95));
  
  
  cout << Form("ul68=%f, ul95=%f,ll68=%f, ll95=%f", ulimit68, ulimit95, llimit68,llimit95)<< endl;
//   cout << Form("@@@@ Do lines: ul68=%f, ul95=%f,ll68=%f, ll95=%f", doLL68, doUL95, doLL68,doLL95)<< endl;
//   cout << Form("@@@@ Do lines: ul68=%s, ul95=%s,ll68=%s, ll95=%s", doLL68, doUL95, doLL68,doLL95)<< endl;
  cout << Form("@@@@ Do lines: ul68=%d, ul95=%d,ll68=%d, ll95=%d", int(doLL68), int(doUL95), int(doLL68),int(doLL95))<< endl;
   TAxis * gra = ((TGraph*)c->GetPrimitive("Graph"))->GetXaxis();
   Double_t x_min = gra->GetXmin();
   Double_t x_max = gra->GetXmax();
//    Double_t x_min = -10;
//    Double_t x_max = 10;
   
//  //68% C.L
    line.SetLineStyle(kDashed);
    line.SetLineColor(kRed);
    line.SetLineWidth(2);
    if (!doLL68) llimit68=x_min;
    if (!doUL68) ulimit68=x_max;
    if (doLL68 || doUL68) line.DrawLine(llimit68,1,ulimit68,1);
    if (doLL68) line.DrawLine(llimit68,0,llimit68,1);
    if (doLL68) line.DrawLine(ulimit68,0,ulimit68,1);
    
    //95 %CL
    line.SetLineColor(kBlue);
    if (!doLL95) llimit95=x_min;
    if (!doUL95) ulimit95=x_max;
    if (doLL95 || doUL95) line.DrawLine(llimit95,3.84,ulimit95,3.84);
    if (doLL95) line.DrawLine(llimit95,0,llimit95,3.84);
    if (doUL95) line.DrawLine(ulimit95,0,ulimit95,3.84);

    return;  
}


void putInformationOnPlot(FitResult_t val, TCanvas * c, TAxis *x, TAxis *y, TString POI="|k_{3}/k_{1}|"){
  
  TLatex latex;
  latex.SetTextSize(0.025);
  latex.SetTextAlign(10);  //align at special bottom
//   Double_t max_y = y->GetXmax();
//   Double_t max_y_2 = y->GetXmax() - 6*(y->GetXmax() - y->GetXmin())/100;
//   Double_t max_y_3 = y->GetXmax() - 10*(y->GetXmax() - y->GetXmin())/100;
  
  Double_t max_y = 5;
  Double_t max_y_2 = max_y - 6*(max_y)/100;
  Double_t max_y_3 = max_y - 10*(max_y)/100;
  Double_t max_y_4 = max_y - 14*(max_y)/100;
  
//   Double_t max_y_2 = 10*(max_y)/100;
//   Double_t max_y_3 = 6*(max_y)/100;
//   Double_t max_y_4 = 2*(max_y)/100;
  
  
  Double_t min_x = x->GetXmin() + 2*(x->GetXmax() - x->GetXmin())/100;
  std::cout <<"min_x = " <<min_x <<" max_y = " <<max_y << std::endl;
//   Double_t limit68 = getLimitFromTree(t,0.68, POI);
//   Double_t limit95 = getLimitFromTree(t,0.95, POI);
//   Float_t best_fit = findBestFitValue(t,POI);
//   TString fit_info = Form("%s | Best fit %s = %4.3f | U.L. 68(95)\% = %4.3f(%4.3f)",title.Data(), ytitle.Data(),best_fit,limit68, limit95);
  TString title = gSystem->Getenv("PLOT_TAG");
  std::cout <<"PLOT_TAG = " <<title << std::endl;
  TString fit_info = Form("%s",title.Data());
//   TString fit_info_2 = Form("Best fit |k_{3}/k_{1}| = %4.3f | U.L. 68(95)\% = %4.3f(%4.3f)", best_fit,limit68, limit95);
//   TString fit_info_2 = Form("Best fit |k_{3}/k_{1}| = %4.3f", best_fit);
//   TString fit_info_3 = Form("U.L. 68(95)\% = %4.3f(%4.3f)", limit68, limit95);
  
        
    Double_t ulimit68 = val.ul68;
    Double_t ulimit95 = val.ul95;
    Double_t llimit68 = val.ll68;
    Double_t llimit95 = val.ll95;
    TString ll68 = (AreSame(llimit68,default_fr.ll68)) ? "N.A." :  Form("%4.3f",val.ll68);
    TString ll95 = (AreSame(llimit95,default_fr.ll95)) ? "N.A." :  Form("%4.3f",val.ll95);
    TString ul68 = (AreSame(ulimit68,default_fr.ul68)) ? "N.A." :  Form("%4.3f",val.ul68);
    TString ul95 = (AreSame(ulimit95,default_fr.ul95)) ? "N.A." :  Form("%4.3f",val.ul95);
  
  
  TString fit_info_2 = Form("Best fit %s = %4.3f", POI.Data(), val.bf_global);
  TString fit_info_3 = Form("U.L. 68(95)\% = %s(%s)", ul68.Data(), ul95.Data());
  TString fit_info_4 = Form("L.L. 68(95)\% = %s(%s)", ll68.Data(), ll95.Data());
  
  std::cout << fit_info << std::endl;
  
  c->cd();
  latex.DrawLatex(min_x, max_y,fit_info);
//   min_x = 0;
  latex.DrawLatex(min_x, max_y_2,fit_info_2);
  latex.DrawLatex(min_x, max_y_3,fit_info_3);
  latex.DrawLatex(min_x, max_y_4,fit_info_4);
  
  return;
  
}



Double_t getArcTanParam(Double_t k3k1,  Double_t gamma=1){
  return (1/TMath::Pi())*TMath::ATan(TMath::Sqrt(gamma)*k3k1);
}
FitResult_t getRescaledFitResult(FitResult_t fitres, Double_t gamma){
//     rescale all values in fit result
    FitResult_t nfr;
    nfr.bf_global = getArcTanParam(fitres.bf_global,gamma);
    nfr.bf = getArcTanParam(fitres.bf,gamma);
    nfr.wf = getArcTanParam(fitres.wf,gamma);
    
    Double_t ulimit68 = fitres.ul68;
    Double_t ulimit95 = fitres.ul95;
    Double_t llimit68 = fitres.ll68;
    Double_t llimit95 = fitres.ll95;
    nfr.ll68 = (AreSame(llimit68,default_fr.ll68)) ? default_fr.ll68 :  getArcTanParam(fitres.ll68,gamma);
    nfr.ll95 = (AreSame(llimit95,default_fr.ll95)) ? default_fr.ll95 :  getArcTanParam(fitres.ll95,gamma);
    nfr.ul68 = (AreSame(ulimit68,default_fr.ul68)) ? default_fr.ul68 :  getArcTanParam(fitres.ul68,gamma);
    nfr.ul95 = (AreSame(ulimit95,default_fr.ul95)) ? default_fr.ul95 :  getArcTanParam(fitres.ul95,gamma);
return nfr;    
}

FitResult_t getfa3RescaledFitResult(FitResult_t fitres){
//     rescale all values in fit result
    FitResult_t nfr= {-99999,-99999,-99999,99999,99999,-99999,-99999 };
    nfr.bf_global = fa3(fitres.bf_global,0);
    nfr.bf = fa3(fitres.bf,0);
    nfr.wf = fa3(fitres.wf,0);
    
    Double_t ulimit68 = fitres.ul68;
    Double_t ulimit95 = fitres.ul95;
    Double_t llimit68 = fitres.ll68;
    Double_t llimit95 = fitres.ll95;
    nfr.ll68 = (AreSame(llimit68,default_fr.ll68)) ? default_fr.ll68 :  fa3(fitres.ll68,0);
    nfr.ll95 = (AreSame(llimit95,default_fr.ll95)) ? default_fr.ll95 :  fa3(fitres.ll95,0);
    nfr.ul68 = (AreSame(ulimit68,default_fr.ul68)) ? default_fr.ul68 :  fa3(fitres.ul68,0);
    nfr.ul95 = (AreSame(ulimit95,default_fr.ul95)) ? default_fr.ul95 :  fa3(fitres.ul95,0);
return nfr;    
}


TGraph * getfa3RescaledGraph(TGraph *gr){
  Int_t N = gr->GetN();
  std::cout <<"N = " <<N << std::endl;
  Double_t *x;
  std::vector<Double_t> new_x;
  Double_t *y;
  
  x = gr->GetX();
  y = gr->GetY();
//   std::cout <<"N = " <<N << std::endl;
  for (Int_t i=0; i<N; i++){
    new_x.push_back(fa3(x[i],0));
//     std::cout <<"x, new_x = " << x[i] << "; "<< new_x[i]<< std::endl;
  }
  
  TGraph * gr_ret = new TGraph(N,&(new_x.front()),y);
//   TGraph * gr_ret = new TGraph(N,x,y);
  return gr_ret;
}


TGraph * getRescaledGraph(TGraph *gr, Double_t gamma){
  Int_t N = gr->GetN();
  std::cout <<"N = " <<N << std::endl;
  Double_t *x;
  std::vector<Double_t> new_x;
  Double_t *y;
  
  x = gr->GetX();
  y = gr->GetY();
//   std::cout <<"N = " <<N << std::endl;
  for (Int_t i=0; i<N; i++){
    new_x.push_back(getArcTanParam(x[i],gamma));
//     std::cout <<"x, new_x = " << x[i] << "; "<< new_x[i]<< std::endl;
  }
  
  TGraph * gr_ret = new TGraph(N,&(new_x.front()),y);
//   TGraph * gr_ret = new TGraph(N,x,y);
  return gr_ret;
}


Double_t getGamma33(Int_t channel=0){
    Double_t g=0.038;
    if (channel==0)  
        g = 0.040;
//         g = (0.46154957824608*0.040+0.3543010974318054*0.034+0.18414932432211456*0.034 ); //weighted mean in case of calculation for all channels
    else if (channel==3) g = 0.040;
    else g=0.034 ;
    
//     cout<< Form("@@@@@ gamma33 (ch=%d)=%f",channel,g)<< endl;
    return g;
    }

Double_t fa3(Double_t k3k1, Int_t channel){
    Short_t sgn = 1;
    if (k3k1<0) sgn = -1; 
//     cout <<  "@@@@ fa3Calculator: " << endl;
    Double_t g=getGamma33(channel);
    return sgn*g*k3k1*k3k1/(1+g*k3k1*k3k1);
}

    
