from ROOT import *
import ROOT



def norm(k2k1_ratio, k3k1_ratio=0):
        norm = (0.581984238389*k2k1_ratio*(-1))/(0.354301097432*(1+0.09*k2k1_ratio*k2k1_ratio+2*((-0.264118324659*k2k1_ratio)))\
                                        +0.184149324322*(1+0.09*k2k1_ratio*k2k1_ratio+2*((-0.269303399267*k2k1_ratio)))\
                                        +0.461549578246*(1+0.09*k2k1_ratio*k2k1_ratio+2*((-0.290992119195*k2k1_ratio))))
        return norm
        

        
if __name__ == '__main__':        
    #k2k1_ratio = 0.4
    #the_list = []
    #for k2k1_multiplier in range(1000):
        #k2k1_ratio = -20 + 0.1*k2k1_multiplier
        #the_list.append((k2k1_ratio,norm(k2k1_ratio)))
        
        
    ##print "Norm({0}) = {1}".format(k2k1_ratio, norm(k2k1_ratio))
    #print 'The list = ', the_list
    #print "Max item = {0}".format(max(the_list, key=lambda p: p[1]))
    #print "Min item = {0}".format(min(the_list, key=lambda p: p[1]))

    gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    gROOT.SetBatch(1)

    
    #fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.0.PedjaIn.IntOn.D0Ph.Dint12.8TeV.M4l121to131.RECO/HCG/126/workspaceWithAsimov_k2k1_ratio_0_lumi_25.root'
    #fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.PedjaIn.IntOn.D0Ph.Dint12.8TeV.M4l121to131.RECO/HCG/126/combine.ws.4l.v1.root'
    #fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.0.PedjaIn.IntOn.D0Ph.Dint12.7TeV.M4l121to131.RECO.NEW/HCG/126/combine.ws.4l.v1.root'
    #fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.0.PedjaIn.IntOn.D0Ph.Dint12.7TeV.M4l106to141.RECO.NEW/HCG/126/combine.ws.4l.v1.root'
    fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.PedjaIn.IntOn.D0Ph.Dint12.8TeV.M4l106to141.RECO/HCG/126/combine.ws.4l.v1.root'
    #fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_test_2Dtemplates/HCG/126/workspaceWithAsimov_k3k1_ratio_0_lumi_25.root'
    #f = TFile.Open('~/public/4Mingshui/k2k1_ws.root','READ')
    f = TFile.Open(fname,'READ')
    
    w.var('cmshzz4l_lumi_8').setVal(19.79)
    #w.var('cmshzz4l_lumi_7').setVal(5.051)
    w.var('k2k1_ratio').setRange(-20,20)
    
    total_norms = { 1 : {'norm': None, 'color':kBlue, 'name': '4mu'},
                    2 : {'norm': None, 'color':kGreen, 'name': '4e'},
                    3 : {'norm': None, 'color':kRed, 'name': '2e2mu'}
                  }
                    
                    
    for channel in total_norms.keys():
        plot = w.var('k2k1_ratio').frame()
        w.function("n_exp_final_binch{0}_proc_ggInt_12N".format(channel)).plotOn(plot, RooFit.LineColor(kGreen),RooFit.Name('ggInt_12N'))
        w.function("n_exp_final_binch{0}_proc_ggH".format(channel)).plotOn(plot, RooFit.LineColor(kRed),RooFit.Name('ggH'))
        w.function("n_exp_final_binch{0}_proc_gg0Ph".format(channel)).plotOn(plot, RooFit.LineColor(kBlue),RooFit.Name('gg0Ph'))
        total = RooFormulaVar('total','@0+@1+@2',RooArgList(w.function("n_exp_final_binch{0}_proc_gg0Ph".format(channel)), w.function("n_exp_final_binch{0}_proc_ggH".format(channel)),w.function("n_exp_final_binch{0}_proc_ggInt_12N".format(channel))) )
        
        total.plotOn(plot, RooFit.LineColor(kBlack), RooFit.Name('total'))
        leg = TLegend(0.2,0.3,0.5,0.45)
        leg.SetFillColor(kWhite)
        leg.SetLineColor(kWhite)
        leg.AddEntry(plot.getObject(0), 'ggInt_12N term', 'L')
        leg.AddEntry(plot.getObject(1), 'ggH term', 'L')
        leg.AddEntry(plot.getObject(2), 'gg0Ph term', 'L')
        leg.AddEntry(plot.getObject(3), 'Sum of terms, total', 'L')
        canv = TCanvas('Norms','Normalization for k2/k1',800,800)
        canv.cd()
        plot.SetTitle('Term normalizations for {0} channel'.format(total_norms[channel]['name']))
        plot.SetYTitle('Expected yield (@25 fb^{-1})')
        plot.SetXTitle('k_{2}/k_{1}')
        plot.Draw()
        leg.Draw()
        canv.SaveAs('Normalizations_k2k1_{0}.SeparatedDen.png'.format(total_norms[channel]['name']))
        canv.SaveAs('Normalizations_k2k1_{0}.SeparatedDen.root'.format(total_norms[channel]['name']))
    
    
    
    #Now plot the channel norms
    
    leg_chan = TLegend(0.65,0.6,0.86,0.75)
    plot_chan = w.var('k2k1_ratio').frame()
    for channel in total_norms.keys():
        total_norms[channel]['norm'] = RooFormulaVar('total_norm_{0}'.format(channel),'@0+@1+@2',RooArgList(w.function("n_exp_final_binch{0}_proc_gg0Ph".format(channel)), w.function("n_exp_final_binch{0}_proc_ggH".format(channel)),w.function("n_exp_final_binch{0}_proc_ggInt_12N".format(channel))) )
        total_norms[channel]['norm'].plotOn(plot_chan, RooFit.LineColor(total_norms[channel]['color']), RooFit.Name('total_{0}'.format(channel)))
        leg_chan.AddEntry(plot_chan.findObject('total_{0}'.format(channel)), total_norms[channel]['name'], 'l')
        
        
    sum_norms = RooFormulaVar('sum_norms','@0+@1+@2',RooArgList(total_norms[1]['norm'], total_norms[2]['norm'],total_norms[3]['norm']) )
    sum_norms.plotOn(plot_chan, RooFit.LineColor(kBlack), RooFit.Name('sum_norms'))    
    leg_chan.AddEntry(plot_chan.findObject('sum_norms'), 'Total yield','l')

    
    k2k1 = 0
    denominator = TFormula('denominator',"(0.354301097432*(1.0+0.09*x*x+0.034*0.0*0.0+2*((-0.264118324659*x)+(0.0*0.0)+(0.0*x*0.0)))+0.184149324322*(1.0+0.09*x*x+0.034*0.0*0.0+2*((-0.269303399267*x)+(0.0*0.0)+(0.0*x*0.0)))+0.461549578246*(1.0+0.09*x*x+0.04*0.0*0.0+2*((-0.290992119195*x)+(0.0*0.0)+(0.0*x*0.0))))")
    while True:
        w.var('k2k1_ratio').setVal(k2k1)
        print 'k2k1: {4} total:{0} 2e2mu: {1} 4e: {2} 4mu: {3} denom = {5}'.format(sum_norms.getVal(), 
                    total_norms[3]['norm'].getVal(),
                    total_norms[2]['norm'].getVal(),
                    total_norms[1]['norm'].getVal(),
                    w.var('k2k1_ratio').getVal(),
                    denominator.Eval(k2k1)
        )
        k2k1 +=0.1
        if k2k1 > 0.7: break    
        
    canv_chan = TCanvas('Norms_channel','Normalization for k2/k1 in subchannels',800,800)
    canv_chan.cd()
    plot_chan.SetTitle('Channel normalizations')
    plot_chan.SetYTitle('Yield contribution (@25 fb^{-1})')
    plot_chan.SetXTitle('k_{2}/k_{1}')
    plot_chan.Draw()
    leg_chan.SetFillColor(kWhite)
    leg_chan.SetLineColor(kWhite)
    leg_chan.Draw()
    canv_chan.SaveAs('ChannelNormalizations_k2k1.SeparatedDen.png')
    canv_chan.SaveAs('ChannelNormalizations_k2k1.SeparatedDen.root')
    
    
    bkg=[w.function('n_exp_binch3_proc_bkg2d_qqzz'),w.function('n_exp_binch2_proc_bkg2d_qqzz'),w.function('n_exp_binch1_proc_bkg2d_qqzz')]
    print 'qqZZ: total:{0} 2e2mu: {1} 4e: {2} 4mu: {3}'.format(bkg[0].getVal()+bkg[1].getVal()+bkg[2].getVal(), bkg[0].getVal(),bkg[1].getVal(),bkg[2].getVal())                  
    bkg=[w.function('n_exp_binch3_proc_bkg2d_ggzz'),w.function('n_exp_binch2_proc_bkg2d_ggzz'),w.function('n_exp_binch1_proc_bkg2d_ggzz')]
    print 'ggZZ: total:{0} 2e2mu: {1} 4e: {2} 4mu: {3}'.format(bkg[0].getVal()+bkg[1].getVal()+bkg[2].getVal(), bkg[0].getVal(),bkg[1].getVal(),bkg[2].getVal())                  
    bkg=[w.function('n_exp_binch3_proc_bkg2d_zjets'),w.function('n_exp_binch2_proc_bkg2d_zjets'),w.function('n_exp_binch1_proc_bkg2d_zjets')]
    print 'zjets: total:{0} 2e2mu: {1} 4e: {2} 4mu: {3}'.format(bkg[0].getVal()+bkg[1].getVal()+bkg[2].getVal(), bkg[0].getVal(),bkg[1].getVal(),bkg[2].getVal())                  
    
    
    

    #### k3k1_ratio
    
    #f.Close()
    #fname = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_test_2Dtemplates/HCG/126/workspaceWithAsimov_k3k1_ratio_0_lumi_25.root'
    ##f = TFile.Open('~/public/4Mingshui/k2k1_ws.root','READ')
    #fk3k1 = TFile.Open(fname,'READ')
    #fk3k1.cd()
    
    #w = fk3k1.Get('w')
    ##w.Print()
    #w.var('cmshzz4l_lumi').setVal(25)
    #w.var('k3k1_ratio').setRange(-30,30)
    
    #total_norms = {1 : {'norm': None, 'color':kBlue, 'name': '4mu'},
                    #2 : {'norm': None, 'color':kGreen, 'name': '4e'},
                    #3:{'norm': None, 'color':kRed, 'name': '2e2mu'}
                    #}
                    
                    
    #for channel in total_norms.keys():
        #plot = w.var('k3k1_ratio').frame()
        #w.function("n_exp_final_binch{0}_proc_ggInt_13N".format(channel)).plotOn(plot, RooFit.LineColor(kGreen),RooFit.Name('ggInt_13N'))
        #w.function("n_exp_final_binch{0}_proc_ggInt_13P".format(channel)).plotOn(plot, RooFit.LineColor(kGreen),RooFit.LineStyle(kDashed), RooFit.Name('ggInt_13P'))
        #w.function("n_exp_final_binch{0}_proc_ggH".format(channel)).plotOn(plot, RooFit.LineColor(kRed),RooFit.Name('ggH'))
        #w.function("n_exp_final_binch{0}_proc_gg0M".format(channel)).plotOn(plot, RooFit.LineColor(kBlue),RooFit.Name('gg0M'))
        #total = RooFormulaVar('total','@0+@1+@2+@3',RooArgList(w.function("n_exp_final_binch{0}_proc_gg0M".format(channel)), w.function("n_exp_final_binch{0}_proc_ggH".format(channel)),w.function("n_exp_final_binch{0}_proc_ggInt_13N".format(channel)),w.function("n_exp_final_binch{0}_proc_ggInt_13P".format(channel))) )
        #total.plotOn(plot, RooFit.LineColor(kBlack), RooFit.Name('total'))
        #leg = TLegend(0.15,0.45,0.45,0.6)
        #leg.SetFillColor(0)
        #leg.SetFillStyle(0)
        #leg.SetLineColor(kWhite)
        #leg.AddEntry(plot.getObject(0), 'ggInt_13N term', 'L')
        #leg.AddEntry(plot.getObject(1), 'ggInt_13P term', 'L')
        #leg.AddEntry(plot.getObject(2), 'ggH term', 'L')
        #leg.AddEntry(plot.getObject(3), 'gg0M term', 'L')
        #leg.AddEntry(plot.getObject(4), 'Sum of terms, total', 'L')
        #canv = TCanvas('Norms','Normalization for k2/k1',800,800)
        #canv.cd()
        #plot.SetTitle('Term normalizations for {0} channel'.format(total_norms[channel]['name']))
        #plot.SetYTitle('Expected yield (@25 fb^{-1})')
        #plot.SetXTitle('k_{3}/k_{1}')
        #plot.Draw()
        #leg.Draw()
        #canv.SaveAs('Normalizations_k3k1_{0}.png'.format(total_norms[channel]['name']))
        #canv.SaveAs('Normalizations_k3k1_{0}.root'.format(total_norms[channel]['name']))
    
    
    
    ##Now plot the channel norms
    
    #leg_chan = TLegend(0.65,0.6,0.86,0.75)
    #plot_chan = w.var('k3k1_ratio').frame()
    #for channel in total_norms.keys():
        #total_norms[channel]['norm'] = RooFormulaVar('total_norm_{0}'.format(channel),'@0+@1+@2+@3',RooArgList(w.function("n_exp_final_binch{0}_proc_gg0M".format(channel)), w.function("n_exp_final_binch{0}_proc_ggH".format(channel)),w.function("n_exp_final_binch{0}_proc_ggInt_13N".format(channel)),w.function("n_exp_final_binch{0}_proc_ggInt_13P".format(channel))) )
        #total_norms[channel]['norm'].plotOn(plot_chan, RooFit.LineColor(total_norms[channel]['color']), RooFit.Name('total_{0}'.format(channel)))
        #leg_chan.AddEntry(plot_chan.findObject('total_{0}'.format(channel)), total_norms[channel]['name'], 'l')
    #sum_norms = RooFormulaVar('sum_norms','@0+@1+@2',RooArgList(total_norms[1]['norm'], total_norms[2]['norm'],total_norms[3]['norm']) )
    #sum_norms.plotOn(plot_chan, RooFit.LineColor(kBlack), RooFit.Name('sum_norms'))    
    #leg_chan.AddEntry(plot_chan.findObject('sum_norms'), 'Total yield','l')
    
    #canv_chan = TCanvas('Norms_channel','Normalization for k2/k1 in subchannels',800,800)
    #canv_chan.cd()
    #plot_chan.SetTitle('Channel normalizations')
    #plot_chan.SetYTitle('Yield contribution (@25 fb^{-1})')
    #plot_chan.SetXTitle('k_{3}/k_{1}')
    #plot_chan.Draw()
    #leg_chan.SetFillColor(kWhite)
    #leg_chan.SetLineColor(kWhite)
    #leg_chan.Draw()
    #canv_chan.SaveAs('ChannelNormalizations_k3k1.png')
    #canv_chan.SaveAs('ChannelNormalizations_k3k1.root')
    

    
    
    
    
    
        












