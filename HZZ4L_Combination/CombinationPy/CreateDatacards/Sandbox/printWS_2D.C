void printWS_2D(){

        gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so");
//         gROOT->ProcessLine(".L /afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/tdrstyle.cc");
//         setTDRStyle(true);
        gStyle->SetOptStat(0);
//         TFile *fin = new TFile("/afs/cern.ch/user/r/roko/public/4Mingshui/workspaceWithAsimov_k2k1_ratio_0_lumi_25.root"); //WORKSPACE
        TFile *fin = new TFile("/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.PedjaIn.IntOn.D0Ph.Dint12.8TeV.M4l121to131.RECO/HCG/126/workspaceWithAsimov_k2k1_ratio_0_lumi_19.7.root");
        RooWorkspace *w = (RooWorkspace*) fin->Get("w");
//         w->Print();
//         return;

//TH2D* hh_data = data->createHistogram("hh_data",x,Binning(20),YVar(y,Binning(20))) ;
//   TH1* hh_data = data->createHistogram("x,y",20,20)
  
        RooRealVar x = *(w->var("CMS_zz4l_Djcp"));
        RooRealVar y = *(w->var("CMS_zz4l_Djcp_int"));

        TH1 *hInt = (TH1*) w->data("sigTempDataHist_ggInt_12N_3_8")->createHistogram("ggInt_12N", x, RooFit::YVar(y));
        TH1 *hggH = (TH1*) w->data("sigTempDataHist_ggH_3_8")->createHistogram("ggH", x, RooFit::YVar(y));
        TH1 *hgg0Ph = (TH1*) w->data("sigTempDataHist_gg0Ph_3_8")->createHistogram("gg0Ph", x, RooFit::YVar(y));
        hInt->Print("V");
        cout<<" hInt x:  "<<hInt->GetXaxis()->GetXmin() <<" to "<< hInt->GetXaxis()->GetXmax() << " with "<< hInt->GetXaxis()->GetNbins() << " bins " <<endl; 
        cout<<" hggH x:  "<<hggH->GetXaxis()->GetXmin() <<" to "<< hggH->GetXaxis()->GetXmax() << " with "<< hggH->GetXaxis()->GetNbins() << " bins " <<endl; 
        cout<<" hgg0Ph x:  "<<hgg0Ph->GetXaxis()->GetXmin() <<" to "<< hgg0Ph->GetXaxis()->GetXmax() << " with "<< hgg0Ph->GetXaxis()->GetNbins()<<  " bins " <<endl; 

        cout<<" hInt y:  "<<hInt->GetYaxis()->GetXmin() <<" to "<< hInt->GetYaxis()->GetXmax() << " with "<< hInt->GetYaxis()->GetNbins() << " bins " <<endl; 
        cout<<" hggH y:  "<<hggH->GetYaxis()->GetXmin() <<" to "<< hggH->GetYaxis()->GetXmax() << " with "<< hggH->GetYaxis()->GetNbins() << " bins " <<endl; 
        cout<<" hgg0Ph y:  "<<hgg0Ph->GetYaxis()->GetXmin() <<" to "<< hgg0Ph->GetYaxis()->GetXmax() << " with "<< hgg0Ph->GetYaxis()->GetNbins()<<  " bins " <<endl; 

        
        cout<<"hInt integral = "<<   hInt -> Integral()<<endl;
        cout<<"hggH integral = "<<   hggH -> Integral()<<endl;
        cout<<"hgg0Ph integral = "<< hgg0Ph -> Integral()<<endl;
        TH2D * h_total_neg = (TH2D*)hggH->Clone("h_total_neg");
        h_total_neg->GetXaxis()->SetTitle("X");
        h_total_neg->GetYaxis()->SetTitle("Y");
        h_total_neg->GetZaxis()->SetTitle("");
        h_total_neg->Reset();
        
        
        TH2D * h_total_pos = (TH2D*)hggH->Clone("h_total_pos");
        h_total_pos->GetXaxis()->SetTitle("X");
        h_total_pos->GetYaxis()->SetTitle("Y");
        h_total_pos->GetZaxis()->SetTitle("");
        h_total_pos->Reset();
        
        cout << "Integral of h_total_neg = " << h_total_neg->Integral()<< endl;
//         cout << "h_total_neg title= " << h_total_neg->GetTitle() << endl;
//         cout << "h_total_neg name= " << h_total_neg->GetName() << endl;
        // check normalization of each hists 
        double tmp = 0;
        TCanvas * c1 = new TCanvas("c1","c1",1400,1650);
        gPad->SetRightMargin(0.2);
        gStyle->SetOptTitle(1);
        Double_t the_most_negative=999;
        Double_t the_most_negative_integral=999;
        Double_t the_most_negative_integral_k2k1;
        
        Double_t the_most_positive=-999;
        Double_t the_most_positive_k2k1,the_most_negative_k2k1 ;
//         exit();
        for(double k2k1=1; k2k1<=5; k2k1+=0.1){
                w->var("k2k1_ratio")->setVal(k2k1);
                double ggInt = w->function("n_exp_final_binch3_proc_ggInt_12N")->getVal();
                double gg0Ph = w->function("n_exp_final_binch3_proc_gg0Ph")->getVal();
                double ggH = w->function("n_exp_final_binch3_proc_ggH")->getVal();
                int nnegbins = 0;

                for(int i=0; i<=50; i++){
                    for(int j=0; j<=50; j++){
                            tmp = ggInt * hInt->GetBinContent(i,j) + ggH * hggH->GetBinContent(i,j) + gg0Ph * hgg0Ph->GetBinContent(i,j) ;
                            if(tmp<0) { 
                                    cout<< "bin "<<i<<","<< j <<": sum is negative "<<tmp<<".  ggInt="<<ggInt<<" x "<<hInt->GetBinContent(i,j)<<" ggH="<<ggH<<" x "<<hggH->GetBinContent(i,j)<<" gg0Ph="<<gg0Ph<<" x "<<hgg0Ph->GetBinContent(i,j)<<endl; 
                                    nnegbins++;
                                    h_total_neg->SetBinContent(i,j,tmp);
                                    if (tmp < the_most_negative) {
                                        the_most_negative=tmp;
                                        the_most_negative_k2k1 = k2k1;
                                    }
                                    
                                    
                                    
                            }
                            else {
                                h_total_pos->SetBinContent(i,j,tmp);
                                if (tmp > the_most_positive) {
                                    the_most_positive=tmp;
                                    the_most_positive_k2k1 = k2k1;
                                }
                                    
                            }
                    }
                }
                Double_t tmp_integral = h_total_neg->Integral();
                if (tmp_integral < the_most_negative_integral){
                    the_most_negative_integral = tmp_integral;
                    the_most_negative_integral_k2k1 = k2k1;
                }
                cout<<"--------> In total, there are "<<nnegbins<<" bins with 0 or negative sum"<<endl;                
                c1->Clear();
                c1->Divide(3,3);
                c1->cd(1);
                TString title = Form("PDF < 0 for k2k1 = %3.2f",k2k1);
                if (fabs(k2k1-3.1)<0.0001) title = Form("PDF < 0 for k2k1 = %3.2f  !!! The most negative!!!",k2k1);
                h_total_neg->SetTitle(title);
//                 h_total_neg->SetMaximum(-0.001);
//                 h_total_neg->SetMinimum(-0.126);
                h_total_neg->Scale(-1);
                h_total_neg->SetMaximum(0.126);
                h_total_neg->SetMinimum(0.);
//                 c1->SetTitle(Form("k_{2}/k_{1} = %3.2",k2k1));
                h_total_neg->Draw("colz");
              
                c1->cd(2);
                title = Form("PDF > 0 for k2k1 = %3.2f",k2k1);
                h_total_pos->SetTitle(title);
                h_total_pos->SetMaximum(0.092);
                h_total_pos->SetMinimum(0);
//                 c1->SetTitle(Form("k_{2}/k_{1} = %3.2",k2k1));
                h_total_pos->Draw("colz");
                
                c1->cd(3);
                TBox k2k1_box;
//                 k2k1_box.SetNDC(1);
                k2k1_box.SetFillColor(kRed);
                k2k1_box.DrawBox(0.1,0.1,0.1+0.7*(k2k1)/5, 0.2 );
                TLatex latex;
                latex.SetTextSize(0.1);
                latex.SetTextAlign(10);
                TString label = Form("k_{2}/k_{1} = %3.2f",k2k1);
                latex.DrawLatex(0.1, 0.25,label.Data());
                
                TH2D *hggH_2 = ((TH2D*)hggH->Clone("hggH_2"));
                TH2D *hgg0Ph_2 = ((TH2D*)hgg0Ph->Clone("hgg0Ph_2"));
                TH2D *hInt_2 = ((TH2D*)hInt->Clone("hInt_2"));
                c1->cd(4);
                hggH_2->Scale(ggH);
                hggH_2->SetMaximum(0.025*95);
                TString template_title = Form("%4.2f x ggH template",ggH);
                hggH_2->SetTitle(template_title);
                hggH_2->Draw("colz");
                
                c1->cd(5);
                hgg0Ph_2->Scale(gg0Ph);
                hgg0Ph_2->SetMaximum(0.025*95);
                template_title = Form("%4.2f x gg0Ph template",gg0Ph);
                hgg0Ph_2->SetTitle(template_title);
                hgg0Ph_2->Draw("colz");

                c1->cd(6);
                hInt_2->Scale(-1*ggInt);
                hInt_2->SetMaximum(0.025*95);
                template_title = Form("(-1) x %4.2f x ggInt_12N template",ggInt);
                hInt_2->SetTitle(template_title);
                hInt_2->Draw("colz");

                
// // // // // // // // // // // // // // // // // // //                 
                
                c1->cd(7);
                hggH->SetMaximum(0.025);
                template_title = Form("ggH template");
                hggH->SetTitle(template_title);
                hggH->Draw("colz");
                
                c1->cd(8);
                hgg0Ph->SetMaximum(0.025);
                template_title = Form("gg0Ph template");
                hgg0Ph->SetTitle(template_title);
                hgg0Ph->Draw("colz");

                c1->cd(9);
                hInt->SetMaximum(0.025);
                template_title = Form("ggInt_12N template");
                hInt->SetTitle(template_title);
                hInt->Draw("colz");
                
                
                
                
                TString fig_name = Form("totalPDF_is_neg_k2k1_%3.2f",k2k1);
                c1->SaveAs(fig_name + ".png");
                c1->SaveAs(fig_name + ".root");
                c1->SaveAs(fig_name + ".gif");
                h_total_neg->Reset();
                h_total_pos->Reset();
//                 break;
                
                delete hggH_2;
                delete hgg0Ph_2;
                delete hInt_2;
                
                
        }
        
        cout <<"The most negative value : " << the_most_negative << " for k2k1 = " << the_most_negative_k2k1 << endl;
        cout <<"The most positive value : " << the_most_positive << " for k2k1 = " << the_most_positive_k2k1 << endl;
        cout <<"The most negative integral : " << the_most_negative_integral << " for k2k1 = " << the_most_negative_integral_k2k1 << endl;
        
        
       

}
