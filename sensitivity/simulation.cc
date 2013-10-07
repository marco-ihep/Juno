#include <map>
#include <utility>
#include "simulation.h"

using namespace std;


Simulation::Simulation(char *outputname)
{
    cout << "Simulation::Simulation()" << endl;
    output = new TFile((string(outputname) + "/output.root").c_str(), "recreate");
    path_ = string(outputname);

    //RooMsgService::instance().getStream(1).addTopic(Workspace) ;
    // RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveIntegratorND").setRealValue("maxWarn",10) ; 

    // RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30) ;
    // RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveIntegratorND").setRealValue("maxEval2D", 2000000) ; 
    //RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D") ;
    //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg",100) ;

    RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D") ;
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","21Points") ;
    RooAbsReal::defaultIntegratorConfig()->Print("v") ;

    /// DEBUG
    // RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-2);
    // RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-2);


    cout << "Integrator Configured" << endl;

    spectrum_index = 0;
    smeared_spectrum_basename_ = "smeared_spectrum";
    multiple_spectra_basename_ = "multiple_spectra_pdf";
    
    map_baseline_reactorPower_[string("50500.0")] = 17.4; // m - GW
    map_baseline_reactorPower_[string("51000.0")] = 18.4; // m - GW

    cout << "constructor::smeared_spectrum_basename : " << smeared_spectrum_basename_ << endl;

    w_nh = 0;
    w_ih = 0;
}



objectsToPrint
Simulation::generate_multiple_spectra( Hierarchy thisHierarchy ) 
{
    char thisHierarchy_char[2];
    char ws_name[20];
    char ws_title[50];
    double tot_power = 0.0;
    int nSpectra = 0;
    string factory_cmd;
    string multiple_spectra_name;
    map<string, string> map_baseline_spectrumName;


    sprintf(thisHierarchy_char, "%s", thisHierarchy==NH ? "nh": "ih");
    sprintf(ws_name, "w_%s", thisHierarchy_char);
    sprintf(ws_title, "%s oscillation models", thisHierarchy_char);

    /// initialize workspace pointer
    RooWorkspace *ws = new RooWorkspace(ws_name, ws_title) ;
    if(thisHierarchy == NH)
        w_nh = ws;
    else
        w_ih = ws;


    /// determine tot power to compute relative fractions among spectra
    for(map<string, double>::iterator it = map_baseline_reactorPower_.begin(); it != map_baseline_reactorPower_.end(); ++it)
        tot_power += it->second;

    /// initialize common vars
    init_common_vars(ws);
    cout << ws->GetName() << ":: common variables initialized" << endl;

    /// add spectra to ws
    for(map<string, double>::iterator it = map_baseline_reactorPower_.begin(); it != map_baseline_reactorPower_.end(); ++it, ++nSpectra)
    {
        map_baseline_spectrumName[it->first] = addSpectrum(ws, thisHierarchy, (it->first).c_str());

        cout << ws->GetName() << " :: added " 
             << smeared_spectrum_basename_
             << map_baseline_spectrumName[it->first] 
             << " with baseline "
             << it->first 
             << " and power "
             << map_baseline_reactorPower_[it->first] << endl;

        /// create RooAddPdf string to be processed by RooWorkspace::factory()
        char fraction[20];
        sprintf(fraction, (nSpectra < int(map_baseline_reactorPower_.size() - 1))? "frac%d[%f]*" : "", nSpectra, it->second / tot_power );

        factory_cmd += string(fraction);
        factory_cmd += smeared_spectrum_basename_ + map_baseline_spectrumName[it->first]; 
        factory_cmd += nSpectra < int(map_baseline_reactorPower_.size() - 1) ? ", " : " ";

    }

    /// process RooAddPdf only if there is more than one component
    if(map_baseline_reactorPower_.size()>1)
    {
        multiple_spectra_name = multiple_spectra_basename_ + "_" + string(thisHierarchy_char);
        factory_cmd = "SUM::" + multiple_spectra_name + "("+ factory_cmd + ")"; 
        ws->factory(factory_cmd.c_str());

        cout << factory_cmd << endl;
    }


    /*******************
     *                 *
     *     caching     *
     *                 *
     *******************/

    //ws->var("baseline_m")->setBins(100, "cache");
    //ws->var("e_nu")->setBins(500, "cache");
    //ws->var("e_obs")->setBins(500, "cache");
    // ws->pdf("raw_spectrum_NH")->setStringAttribute("CACHEPARAMINT","e_nu:baseline_m");
    // ws->pdf("survival_NH_baseline")->setStringAttribute("CACHEPARAMINT","baseline_m");
    // ws->pdf("smeared_spectrum_2D_NH")->setStringAttribute("CACHEPARAMINT","e_nu:e_obs");
    //ws->pdf(multiple_spectra_name.c_str())->setStringAttribute("CACHEPARAMINT","e_obs");


    /*************************
     *                       *
     *      asimov data      *
     *                       *
     *************************/

    ws->var("e_obs")->setBins(150, "data_binning");
    RooDataHist *binnedData;

    /// dependinf on the number of components, generate asimov data
    /// according to a single-source pdf or to a multi-sources pdf
    if(map_baseline_reactorPower_.size()>1)
        binnedData = ws->pdf(multiple_spectra_name.c_str())->generateBinned(*ws->var("e_obs"), NumEvents(21214), Asimov());
    else
        binnedData =  ws->pdf((smeared_spectrum_basename_ + map_baseline_spectrumName.begin()->second).c_str() \
                                            )->generateBinned(*ws->var("e_obs"), NumEvents(21214), Asimov());

    ws->import(*binnedData, Rename((string("AsimovData") + (thisHierarchy==NH ? "NH": "IH")).c_str()));


    /// absolute normalization NH
    // ws->var("e_nu")->setRange("int_range_nh", 0.1, 9.9);
    // RooAbsReal *myInt_NH =  ws->function("spectrum_NH_norm")->createIntegral(RooArgSet(*ws->var("e_nu")), "int_range_nh");
    // double nEvents_NH =  myInt_NH->getVal(); 
    // cout << "nevents NH = " << nEvents_NH << endl;


    output->cd();
    ws->Write();

    objectsToPrint outputnames(map_baseline_spectrumName, multiple_spectra_name, thisHierarchy);
    return outputnames;
}


void
Simulation::drawBaseline(objectsToPrint & names)
{
    output->cd();

    RooWorkspace *ws = names.thisHierarchy_==NH ? w_nh : w_ih;

    char thisHierarchy_char[2];
    sprintf(thisHierarchy_char, "%s", names.thisHierarchy_==NH ? "nh": "ih");

    RooPlot *b_frame = ws->var((string("baseline_m") + names.map_baseline_spectrumName_.begin()->second).c_str())->frame();
    b_frame->SetName("baseline_frame");
    b_frame->SetTitle("Survival Prob (NH) vs Baseline");
    ws->pdf( (string("survival_prob_pdf") + names.map_baseline_spectrumName_.begin()->second + "_baseline").c_str())->plotOn(b_frame);

    TCanvas *baseline_c = new TCanvas("baseline_c", "baseline_c", 800, 600);
    b_frame->Draw();
    b_frame->Write();
    baseline_c->Write();
    baseline_c->SaveAs((path_ + "/baseline_" + thisHierarchy_char + ".png").c_str());
}




void 
Simulation::drawJointHierachies(objectsToPrint & names_nh, objectsToPrint & names_ih)
{
    if(!w_nh || !w_ih)
    {
        cout << "ERROR :: drawJointHierachies :: at least one among w_nh and w_ih is null" << endl;
        return;
    }

    map<string, string> map_nh(names_nh.map_baseline_spectrumName_);
    map<string, string> map_ih(names_ih.map_baseline_spectrumName_);

    string spectrum_nh(names_nh.multiple_spectra_name_);
    string spectrum_ih(names_ih.multiple_spectra_name_);

    bool isMultipleSpectra = (map_nh.size() > 1) && (map_ih.size() > 1);


    /// fix baselines
    for(map<string, string>::iterator it = map_nh.begin(); it != map_nh.end(); ++it)
        w_nh->var( (string("baseline_m") + it->second).c_str() )->setConstant();

    for(map<string, string>::iterator it = map_ih.begin(); it != map_ih.end(); ++it)
        w_ih->var( (string("baseline_m") + it->second).c_str() )->setConstant();


    output->cd();
    
    RooPlot * nh_data_frame = w_nh->var("e_obs")->frame();
    RooPlot * ih_data_frame = w_ih->var("e_obs")->frame();
    RooPlot * nh_spectrum_only_frame = w_nh->var("e_obs")->frame();
    RooPlot * ih_spectrum_only_frame = w_ih->var("e_obs")->frame();

    nh_data_frame->SetName("nh_data_frame");
    ih_data_frame->SetName("ih_data_frame");
    nh_spectrum_only_frame->SetName("nh_spectrum_only_frame");
    ih_spectrum_only_frame->SetName("ih_spectrum_only_frame");

    nh_data_frame->SetTitle("NH Asimov Data");
    ih_data_frame->SetTitle("IH Asimov Data");
    nh_spectrum_only_frame->SetTitle("NH/IH Spectra");
    ih_spectrum_only_frame->SetTitle("NH/IH Spectra");

    // RooPlot *nh_spectrum_nh_data_frame = ws->var("e_obs")->frame();
    //          nh_spectrum_nh_data_frame->SetTitle((string("Spectrum") + thisHierarchy_char).c_str() );
    //          nh_spectrum_nh_data_frame->SetName((string(thisHierarchy_char) + "_spectrum_" + thisHierarchy_char + "_data_frame").c_str());


    w_nh->data("AsimovDataNH")->plotOn(nh_data_frame, DataError(RooAbsData::Poisson));
    w_ih->data("AsimovDataIH")->plotOn(ih_data_frame, DataError(RooAbsData::Poisson));

    if(isMultipleSpectra)
    {
        w_nh->pdf(spectrum_nh.c_str())->plotOn(nh_data_frame,          LineColor(kBlue), Normalization(21214, RooAbsReal::NumEvent));
        w_nh->pdf(spectrum_nh.c_str())->plotOn(nh_spectrum_only_frame, LineColor(kBlue), Normalization(21214, RooAbsReal::NumEvent));

        w_ih->pdf(spectrum_ih.c_str())->plotOn(ih_data_frame,          LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent));
        w_ih->pdf(spectrum_ih.c_str())->plotOn(ih_spectrum_only_frame, LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent));
    }
    else
    {
        w_nh->pdf((smeared_spectrum_basename_ + map_nh.begin()->second).c_str())->plotOn( nh_data_frame, \
                LineColor(kBlue), Normalization(21214, RooAbsReal::NumEvent));
        w_nh->pdf((smeared_spectrum_basename_ + map_nh.begin()->second).c_str())->plotOn( nh_spectrum_only_frame, \
                LineColor(kBlue), Normalization(21214, RooAbsReal::NumEvent));

        w_ih->pdf((smeared_spectrum_basename_ + map_ih.begin()->second).c_str())->plotOn( ih_data_frame, \
                LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent));
        w_ih->pdf((smeared_spectrum_basename_ + map_ih.begin()->second).c_str())->plotOn( ih_spectrum_only_frame, \
                LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent));
    }

    nh_data_frame->Write();
    ih_data_frame->Write();
    nh_spectrum_only_frame->Write();
    ih_spectrum_only_frame->Write();


    TCanvas *joint_spectra_nh_data_canvas = new TCanvas("joint_spectra_nh_data_canvas", "joint_spectra_nh_data_canvas", 800, 700);
    nh_data_frame->Draw();
    ih_spectrum_only_frame->Draw("same");
    joint_spectra_nh_data_canvas->Write();
    joint_spectra_nh_data_canvas->SaveAs((path_ + "/joint_spectra_w_nh_data.png").c_str());


    TCanvas *joint_spectra_ih_data_canvas = new TCanvas("joint_spectra_ih_data_canvas", "joint_spectra_ih_data_canvas", 800, 700);
    ih_data_frame->Draw();
    nh_spectrum_only_frame->Draw("same");
    joint_spectra_ih_data_canvas->Write();
    joint_spectra_ih_data_canvas->SaveAs((path_ + "/joint_spectra_w_ih_data.png").c_str());


    TCanvas *both_spectra_canvas = new TCanvas("both_spectra_canvas", "both_spectra_canvas", 800, 600);
    nh_spectrum_only_frame->Draw();
    ih_spectrum_only_frame->Draw("same");
    both_spectra_canvas->Write();
    both_spectra_canvas->SaveAs((path_ + "/joint_spectra.png").c_str());

}





void 
Simulation::drawSpectra(objectsToPrint & names)
{
    int iComponent = 0;
    int colors[] = {kRed-7, kRed-3};
    char thisHierarchy_char[2];

    RooWorkspace *ws = names.thisHierarchy_==NH ? w_nh : w_ih;

    sprintf(thisHierarchy_char, "%s", names.thisHierarchy_==NH ? "nh": "ih");

    map<string, string> map_baseline_spectrumName(names.map_baseline_spectrumName_);
    string multiple_spectra_name(names.multiple_spectra_name_);

    bool isMultipleSpectra = multiple_spectra_name.size() > 1;


    /// fix baselines
    for(map<string, string>::iterator it = map_baseline_spectrumName.begin(); it != map_baseline_spectrumName.end(); ++it)
        ws->var( (string("baseline_m") + it->second).c_str() )->setConstant();



    output->cd();
    
    RooPlot* e_obs_f[10]; 
    e_obs_f[0] = ws->var("e_obs")->frame();
    e_obs_f[1] = ws->var("e_obs")->frame();
    e_obs_f[2] = ws->var("e_obs")->frame();

    RooPlot *nh_spectrum_nh_data_frame = ws->var("e_obs")->frame();
             nh_spectrum_nh_data_frame->SetTitle((string("Spectrum") + thisHierarchy_char).c_str() );
             nh_spectrum_nh_data_frame->SetName((string(thisHierarchy_char) + "_spectrum_" + thisHierarchy_char + "_data_frame").c_str());

    for(map<string, string>::iterator it = map_baseline_spectrumName.begin(); it != map_baseline_spectrumName.end(); ++it, iComponent++)
    {
        ws->pdf( (smeared_spectrum_basename_ + it->second).c_str() \
                 )->plotOn(e_obs_f[iComponent], LineColor(colors[iComponent%2]), Normalization(21214, RooAbsReal::NumEvent)); 

        if(isMultipleSpectra)
            ws->pdf( multiple_spectra_name.c_str() \
                     )->plotOn(nh_spectrum_nh_data_frame, LineColor(colors[iComponent%2]), Normalization(21214, RooAbsReal::NumEvent), \
                               Components(*ws->pdf( (smeared_spectrum_basename_ + it->second).c_str() ))); 
    }

    //binnedData->plotOn(nh_spectrum_nh_data_frame);
    ws->data((string("AsimovData") + (names.thisHierarchy_==NH ? "NH": "IH")).c_str())->plotOn(nh_spectrum_nh_data_frame);

    if(isMultipleSpectra)
        ws->pdf(multiple_spectra_name.c_str())->plotOn(nh_spectrum_nh_data_frame, LineColor(kBlue), LineStyle(kDashed), Normalization(21214, RooAbsReal::NumEvent));
    else
        ws->pdf((smeared_spectrum_basename_ + map_baseline_spectrumName.begin()->second).c_str())->plotOn(nh_spectrum_nh_data_frame, LineColor(kBlue), Normalization(21214, RooAbsReal::NumEvent));

    e_obs_f[0]->SetName("component0_f");
    e_obs_f[1]->SetName("component1_f");

    nh_spectrum_nh_data_frame->Write();
    e_obs_f[0]->Write();
    e_obs_f[1]->Write();

    TCanvas *nh_spectrum_nh_data_canvas = new TCanvas((string(thisHierarchy_char) + "_spectrum_" + thisHierarchy_char + "_data_canvas").c_str(), \
                                           (string(thisHierarchy_char) + "_spectrum_" + thisHierarchy_char + "_data_canvas").c_str() , 800, 700);
    nh_spectrum_nh_data_frame->Draw();
    nh_spectrum_nh_data_canvas->Write();
    nh_spectrum_nh_data_canvas->SaveAs((path_ + "/spectrum_w_data_" + thisHierarchy_char + ".png").c_str());

    TCanvas *component0_c = new TCanvas((string("component0_") + thisHierarchy_char + "_c").c_str(), "component0_c", 800, 600);
    e_obs_f[0]->Draw();
    component0_c->Write();
    component0_c->SaveAs((path_ + "/component0_" + thisHierarchy_char + ".png").c_str());

    TCanvas *component1_c = new TCanvas((string("component1_") + thisHierarchy_char + "_c").c_str(), "component1_c", 800, 600);
    e_obs_f[1]->Draw();
    component1_c->Write();
    component1_c->SaveAs((path_ + "/component1_" + thisHierarchy_char + ".png").c_str());

    // TCanvas *nh_baseline_canvas = new TCanvas("nh_baseline_canvas", "nh_baseline_canvas", 800, 700);
    // b_frame->Draw();
    // b_frame->Write();
    // nh_baseline_canvas->Write();


    //output->Close();

}


string 
Simulation::addSpectrum(RooWorkspace *ws, Hierarchy hierarchy, const char *baseline )
{
    char suffix_s[50];
    sprintf(suffix_s, "_spect%d_%s", spectrum_index, (hierarchy==NH ? "nh" : "ih"));
    string suffix(suffix_s);
    
    vector<string> factory_cmds(15, "");

    factory_cmds.at(0) = "baseline_m" + suffix + "[" + baseline + ", 0., 100000.]";
    factory_cmds.at(1) = "expr::delta21" + suffix + "('1.27*dm221*baseline_m" + suffix + "/e_nu', {dm221, e_nu, baseline_m" + suffix + "})";
    factory_cmds.at(2) = "expr::delta31" + suffix + "('1.27*dm231*baseline_m" + suffix + "/e_nu', {dm231, e_nu, baseline_m" + suffix + "})";
    factory_cmds.at(3) = "expr::survival_1" + suffix + "('(e_nu>1.8) * cos2theta13^2 * sin22teta12 * ";
    factory_cmds.at(3) +="(sin(delta21" + suffix + "))^2', {e_nu, cos2theta13, sin22teta12, delta21" + suffix + "})";

    factory_cmds.at(4) = "expr::survival_2" + suffix + "('(e_nu>1.8) * sin22teta13* (sin(abs(delta31" + suffix + ")))^2', ";
    factory_cmds.at(4) +="{e_nu, sin22teta13, delta31" + suffix + "})";

    factory_cmds.at(5) = "expr::survival_3" + suffix + "('(e_nu>1.8) * sin2theta12 * sin22teta13 * (sin(delta21" + suffix + "))^2 * ";
    factory_cmds.at(5) +="cos(2.0*abs(delta31" + suffix + "))', {e_nu, sin2theta12, sin22teta13, delta21" + suffix + ", delta31" + suffix + "})";

    factory_cmds.at(6) = "expr::survival_4" + suffix + "('(e_nu>1.8) * 0.5 * sin2theta12 * sin22teta13 * sin(2.0*delta21" + suffix + ") * ";
    factory_cmds.at(6) += "sin(2.0*abs(delta31" + suffix + "))', {e_nu, sin2theta12, sin22teta13, delta21" + suffix + ", delta31" + suffix + "})";

    factory_cmds.at(7) = "expr::spectrum_norm" + suffix + "('detector_mass *1000000. * proton_fraction * exposure_time_seconds / (1.6726231 * ";
    factory_cmds.at(7) += "4.0 * myPi * baseline_m" + suffix + "^2)', {detector_mass, proton_fraction, exposure_time_seconds, myPi, baseline_m";
    factory_cmds.at(7) += suffix + "})";

    factory_cmds.at(8) =  "expr::survival_prob_formula" + suffix + "('1-survival_1" + suffix + "-survival_2" + suffix + "-survival_3" + suffix;  
    factory_cmds.at(8) += (hierarchy==NH ? "+" : "-");
    factory_cmds.at(8) += "survival_4" + suffix + "',{survival_1" + suffix + ", survival_2" + suffix + ", survival_3" + suffix;  
    factory_cmds.at(8) += ", survival_4" + suffix + "})";

    factory_cmds.at(9) =  "expr::spectrum_formula" + suffix + "('flux_expr * ibd_xsec_expr * survival_prob_formula" + suffix;  
    factory_cmds.at(9) += " * spectrum_norm" + suffix + " ', {flux_expr, ibd_xsec_expr, survival_prob_formula" + suffix;
    factory_cmds.at(9) += ", spectrum_norm" + suffix + "})";

    factory_cmds.at(10) = "EXPR::survival_prob_pdf" + suffix + "('1-survival_1" + suffix + "-survival_2" + suffix + "-survival_3" + suffix;
    factory_cmds.at(10) += (hierarchy==NH ? "+" : "-");
    factory_cmds.at(10) += "survival_4" + suffix + "',{survival_1" + suffix + ", survival_2" + suffix + ", survival_3"; 
    factory_cmds.at(10) += suffix + ", survival_4" + suffix + "})";

    factory_cmds.at(11) = "PROD::true_spectrum" + suffix + "(survival_prob_pdf" + suffix + ", reactorSpectrum)";

    factory_cmds.at(12) = "PROD::smeared_spectrum_2D" + suffix + "(response|e_nu, true_spectrum" + suffix + ")";

    factory_cmds.at(13) = "PROJ::survival_prob_pdf" + suffix + "_baseline(true_spectrum" + suffix + ", e_nu)";

    factory_cmds.at(14) = "PROJ::smeared_spectrum" + suffix + "(smeared_spectrum_2D" + suffix + ", e_nu)";


    /// process commands 
    //for(unsigned int i=0; i < factory_cmds.size(); ++i)
    for(vector<string>::iterator it=factory_cmds.begin(); it != factory_cmds.end(); ++it)
    {
        //cout  << *it << endl;
        ws->factory(it->c_str());
        //cout << endl << endl;
    }

    spectrum_index++;
    return string(suffix);
}


void 
Simulation::init_common_vars(RooWorkspace *ws)
{
    ws->factory("myPi[3.1415926535897]");
	ws->factory("proton_fraction[0.12]");
	ws->factory("reactor_power[20.]"); // GW
	ws->factory("detector_mass[5.]"); // Kt
	ws->factory("exposure_time[5.]"); // years
    ws->factory("sin22teta12[0.857, 0., 1.]");
    ws->factory("sin22teta12_fix[0.857]");
    ws->factory("sin22teta12_err[0.024]");
    ws->factory("sin22teta13[0.089, 0., 1.]");
    ws->factory("sin22teta13_fix[0.089]");
    ws->factory("sin22teta13_err[0.005]");
    ws->factory("dm221[0.000075, 0., 0.0005]");
    ws->factory("dm221_fix[0.000075]");
    ws->factory("dm221_err[0.000002]");
    ws->factory("dm231[0.00232, 0., 0.005]");
    ws->factory("dm231_fix[0.00232]");
    ws->factory("dm231_err[0.0001]");
    ws->factory("sys_fit[0.03, 0.0, 2.0]");
    ws->factory("sys_err[0.03]");
    ws->factory("e_nu[1.0, 10.0]");
    ws->factory("e_obs[0.0, 10.0]");
    ws->factory("weighted_fission_energy[0.58*201.7+0.3*210.+0.07*205+0.05*212.4]");

    ws->factory("RooPolyVar::e_vis(e_nu, {e_shift[-0.8], e_ratio[1.0]})");

    ws->factory("expr::exposure_time_seconds('exposure_time*220.0*86400.0',{exposure_time})");
    ws->factory("expr::flux_a('exp(0.87 -0.16 *e_nu-0.091 *(e_nu^2))',{e_nu})");
    ws->factory("expr::flux_b('exp(0.896-0.239*e_nu-0.0981*(e_nu^2))',{e_nu})");
    ws->factory("expr::flux_c('exp(0.976-0.162*e_nu-0.079 *(e_nu^2))',{e_nu})");
    ws->factory("expr::flux_d('exp(0.793-0.08 *e_nu-0.1085*(e_nu^2))',{e_nu})");
    ws->factory("expr::sin2theta12('0.5*(1.-sqrt(1.-sin22teta12))',{sin22teta12})"); 
    ws->factory("expr::cos2theta13('0.5*(1.+sqrt(1.-sin22teta13))',{sin22teta13})");
    ws->factory("expr::resolution('0.03 * sqrt(e_vis)', {e_vis})");

    /// needed for normalization
    ws->factory("expr::ibd_xsec_expr('(e_nu>1.8) * 9.52*(e_nu-1.29)*sqrt((e_nu-1.29)^2 - 0.511^2)', {e_nu})");
	ws->factory("expr::flux_expr('0.58*flux_a + 0.3*flux_b + 0.07*flux_c + 0.05*flux_d',{flux_a, flux_b, flux_c, flux_d})");

    /// smearing
    ws->factory("Gaussian::response(e_obs, e_vis, resolution)");

    /// pdfs merging all the previous formulas
    ws->factory("EXPR::ibd_xsec('(e_nu>1.8) * 9.52*(e_nu-1.29)*sqrt((e_nu-1.29)^2 - 0.511^2)', {e_nu})");
	ws->factory("EXPR::flux('0.58*flux_a + 0.3*flux_b + 0.07*flux_c + 0.05*flux_d',{flux_a, flux_b, flux_c, flux_d})");

    ws->factory("PROD::reactorSpectrum(flux, ibd_xsec)");
}




RooFormulaVar*
Simulation::getPenalizedChi2(RooWorkspace *ws, RooChi2Var *chi2_data)
{
    string penalizedChi2(chi2_data->GetName());
           penalizedChi2 += " + ";
           penalizedChi2 += "((sin22teta12-sin22teta12_fix)/sin22teta12_err)^2 + ";
           penalizedChi2 += "((sin22teta13-sin22teta13_fix)/sin22teta13_err)^2 + ";
           penalizedChi2 += "((dm221-dm221_fix)/dm221_err)^2 + ";
           penalizedChi2 += "((dm231-dm231_fix)/dm231_err)^2 + ";
           penalizedChi2 += "((sys_fit-1.0)/sys_err)^2";

    cout << endl << "******* PENALIZED CHI2 ********"  << endl;
    cout << penalizedChi2 << endl << endl;


    /****************************************************
     *   Standalone constructor for PENALIZED CHI2
     *   RooArgList has to be built separately
     *   The chi2 variables have to be shared
     *   with the data chi2 ->  RecycleConflictNodes()
     ****************************************************/

    RooArgList list(*chi2_data);
    list.add(*ws->var("sin22teta12"));
    list.add(*ws->var("sin22teta12_fix"));
    list.add(*ws->var("sin22teta12_err"));
    list.add(*ws->var("sin22teta13"));
    list.add(*ws->var("sin22teta13_fix"));
    list.add(*ws->var("sin22teta13_err"));
    list.add(*ws->var("dm221"));
    list.add(*ws->var("dm221_fix"));
    list.add(*ws->var("dm221_err"));
    list.add(*ws->var("dm231"));
    list.add(*ws->var("dm231_fix"));
    list.add(*ws->var("dm231_err"));
    list.add(*ws->var("sys_fit"));
    list.add(*ws->var("sys_err"));

    const char *cc = penalizedChi2.c_str();
    RooFormulaVar *chi2_penalized = new RooFormulaVar("chi2_penalized", "penalized chi2", cc , list);

    return chi2_penalized;
}



//      void
//      Simulation::fit(const char *nh_filename, const char *ih_filename)
//      {
//      
//          TFile *nh_file = TFile::Open(nh_filename);
//          TFile *ih_file = TFile::Open(ih_filename);
//      
//          if(!nh_file) {cout << "ERROR: couldn't open " << nh_file << endl; exit (-1);} 
//          if(!ih_file) {cout << "ERROR: couldn't open " << ih_file << endl; exit (-1);}
//      
//          w_nh = (RooWorkspace*)nh_file->Get("w_nh");
//          w_ih = (RooWorkspace*)ih_file->Get("w_ih");
//      
//          if(!w_nh) { cout << "ERROR: couldn't find w_nh in " << nh_filename << endl; exit(-1); }
//          if(!w_ih) { cout << "ERROR: couldn't find w_ih in " << ih_filename << endl; exit(-1); }
//      
//      
//          w_nh->var("baseline_m")->setConstant();
//          w_ih->var("baseline_m")->setConstant();
//      
//      
//      
//          /******************************************
//           *                                        *
//           *    chi2 definition and minimization    *
//           *                                        *
//           ******************************************/
//      
//          w_nh->var("e_obs")->setRange("fit_range", 1.1, 8.0);
//      
//          //RooChi2Var chi2_data_NH("chi2_data_NH","chi2_data_NH", *w_nh->pdf("smeared_spectrum_NH"), *((RooDataHist*)w_nh->data("AsimovDataNH")), DataError(RooAbsData::Poisson)) ;
//          RooChi2Var chi2_data_NH("chi2_data_NH","chi2_data_NH", *w_nh->pdf("smeared_spectrum_NH"), *((RooDataHist*)w_nh->data("AsimovDataNH")), DataError(RooAbsData::Poisson), Range("fit_range"), Verbose()) ;
//          RooFormulaVar *chi2_penalized_NH = getPenalizedChi2(w_nh, &chi2_data_NH); 
//          chi2_penalized_NH->SetName("chi2_penalized_NH");
//      
//          RooChi2Var chi2_data_IH("chi2_data_IH","chi2_data_IH", *w_ih->pdf("smeared_spectrum_IH"), *((RooDataHist*)w_nh->data("AsimovDataNH")), DataError(RooAbsData::Poisson), Range("fit_range"), Verbose()) ;
//          RooFormulaVar *chi2_penalized_IH = getPenalizedChi2(w_ih, &chi2_data_IH); 
//          chi2_penalized_IH->SetName("chi2_penalized_IH");
//          
//          cout << endl << "************* MINUIT ON CHI2 NH *************" << endl;
//          RooMinuit minimizer_nh(*chi2_penalized_NH) ;
//                minimizer_nh.migrad() ;
//                minimizer_nh.hesse() ;
//      
//                RooFitResult* r_chi2_nh = minimizer_nh.save() ;
//                cout << endl  << endl << endl << endl << endl << endl;
//                cout << "************* FIT RESULT NH *************" << endl;
//                r_chi2_nh->Print("v");
//      
//      
//          /// IH chi2 built with NH asimov dataset
//      
//          cout << endl << "************* MINUIT ON CHI2 IH *************" << endl;
//          RooMinuit minimizer_ih(*chi2_penalized_IH) ;
//                minimizer_ih.migrad() ;
//                minimizer_ih.hesse() ;
//        
//                RooFitResult* r_chi2_ih = minimizer_ih.save() ;
//                cout << endl  << endl << endl << endl << endl << endl;
//                cout << "************* FIT RESULT IH *************" << endl;
//                r_chi2_ih->Print("v");
//         
//      
//          /***********************
//           *                     *
//           *    chi2 plotting    *
//           *                     *
//           ***********************/
//      
//      
//          RooPlot *nh_chi2_frame = w_nh->var("dm231")->frame(); 
//                   nh_chi2_frame->SetName("nh_chi2_frame");
//                   nh_chi2_frame->SetTitle("#Chi^{2}(NH, AsimovNH) vs #Deltam^{2}_{31}");
//                   cout << "DEBUG :: w_nh->var(dm231)->GetName()  " << w_nh->var("dm231")->GetName() << endl;
//                   cout << "DEBUG :: nh_chi2_frame " << nh_chi2_frame << "  " << nh_chi2_frame->GetName() << endl;
//                   chi2_penalized_NH->plotOn(nh_chi2_frame, LineColor(kBlue), LineStyle(kDashed));
//                   output->cd();
//                   //nh_chi2_frame->Write();
//      
//          TCanvas *nh_chi2_canvas = new TCanvas("nh_chi2_canvas", "nh_chi2_canvas", 800, 600 );
//                   nh_chi2_frame->Draw();
//                   nh_chi2_canvas->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/nh_chi2_canvas.png");
//                   output->cd();
//                   nh_chi2_canvas->Write();
//      
//      
//          
//          RooPlot *ih_chi2_frame = w_ih->var("dm231")->frame(); 
//                   ih_chi2_frame->SetName("ih_chi2_frame");
//                   ih_chi2_frame->SetTitle("#Chi^{2}(IH, AsimovNH) vs #Deltam^{2}_{31}");
//                   chi2_penalized_IH->plotOn(ih_chi2_frame, LineColor(kRed));
//                   output->cd();
//                   //ih_chi2_frame->Write();
//      
//          TCanvas *ih_chi2_canvas = new TCanvas("ih_chi2_canvas", "ih_chi2_canvas", 800, 600 );
//                   ih_chi2_frame->Draw();
//                   ih_chi2_canvas->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/ih_chi2_canvas.png");
//                   output->cd();
//                   ih_chi2_canvas->Write();
//          
//      
//          TCanvas *both_chi2_canvas = new TCanvas("both_chi2_canvas", "both_chi2_canvas", 800, 600);
//                   ih_chi2_frame->Draw();
//                   nh_chi2_frame->Draw("same");
//                   both_chi2_canvas->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/both_chi2_canvas.png");
//                   output->cd();
//                   both_chi2_canvas->Write();
//      
//      
//          output->cd();
//          /// combined spectra
//          TCanvas *data_c = new TCanvas("data_c", "data_c", 800, 600);
//          RooPlot *data_f = w_nh->var("e_obs")->frame();
//          data_f->SetName("data_f");
//          // w_ih->pdf("smeared_spectrum_IH")->plotOn(data_f, Normalization(21214, RooAbsReal::NumEvent) , LineColor(8), Range("fit_range"));  
//          // w_nh->pdf("smeared_spectrum_NH")->plotOn(data_f, Normalization(21214, RooAbsReal::NumEvent) , LineColor(kRed), Range("fit_range")); 
//          w_ih->pdf("smeared_spectrum_IH")->plotOn(data_f, LineColor(8), Range("fit_range"));  
//          w_nh->pdf("smeared_spectrum_NH")->plotOn(data_f, LineColor(kRed), Range("fit_range")); 
//          w_nh->data("AsimovDataNH")->plotOn(data_f);
//          data_f->Draw();
//          data_c->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/data_c.png");
//          output->cd();
//          //data_f->Write();
//          data_c->Write();
//      
//          // TCanvas *both_spectra_nh_data_canvas = new TCanvas("both_spectra_nh_data_canvas", "both_spectra_nh_data_canvas", 800, 700);
//          // ((RooPlot*)nh_file->Get("nh_spectrum_nh_data_frame"))->Draw();
//          // ((RooPlot*)nh_file->Get("ih_spectrum_frame"))->Draw("same");
//          // both_spectra_nh_data_canvas->Write();
//      
//      
//          /// standalone fitted spectrum to overcome normalization problems
//          output->cd();
//          TCanvas *fitted_spectrum_canvas = new TCanvas("fitted_spectrum_canvas","fitted_spectrum_canvas", 800,600); 
//          RooPlot *fitted_spectrum_frame = w_nh->var("e_obs")->frame();
//                   fitted_spectrum_frame->SetName("fitted_spectrum_frame");
//          w_nh->pdf("smeared_spectrum_NH")->plotOn(fitted_spectrum_frame, LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent), Range("fit_range")); 
//          fitted_spectrum_frame->Draw();
//          fitted_spectrum_canvas->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/NHFittedSpectrum.png");
//          output->cd();
//          fitted_spectrum_canvas->Write();
//          //fitted_spectrum_frame->Write();
//      
//      
//      
//          nh_file->Close();
//          ih_file->Close();
//          output->Close();
//      }



void
Simulation::fit(objectsToPrint & names_nh, objectsToPrint & names_ih, Hierarchy dataHierarchy, Hierarchy fcnHierarchy)
{

    if(!w_nh || !w_ih)
    {
        cout << "ERROR :: fit_nh :: at least one among w_nh and w_ih is null" << endl;
        return;
    }


    map<string, string> map_nh(names_nh.map_baseline_spectrumName_);
    map<string, string> map_ih(names_ih.map_baseline_spectrumName_);

    string spectrum_nh(names_nh.multiple_spectra_name_);
    string spectrum_ih(names_ih.multiple_spectra_name_);

    bool isMultipleSpectra = (map_nh.size() > 1) && (map_ih.size() > 1);


    /// fix baselines
    for(map<string, string>::iterator it = map_nh.begin(); it != map_nh.end(); ++it)
        w_nh->var( (string("baseline_m") + it->second).c_str() )->setConstant();

    for(map<string, string>::iterator it = map_ih.begin(); it != map_ih.end(); ++it)
        w_ih->var( (string("baseline_m") + it->second).c_str() )->setConstant();



    /******************************************
     *                                        *
     *    chi2 definition and minimization    *
     *                                        *
     ******************************************/

    // w_nh->var("e_obs")->setRange("fit_range", 1.1, 8.0);
    // w_nh->var("e_obs")->setRange("full_range", 0.0, 10.0);
    // w_ih->var("e_obs")->setRange("fit_range", 1.1, 8.0);
    // w_ih->var("e_obs")->setRange("full_range", 0.0, 10.0);

    string dataHierarchy_name  = dataHierarchy == NH ? "NH" : "IH"; 
    string fcnHierarchy_name   = fcnHierarchy == NH ? "NH" : "IH"; 
    string chi2_name           = "chi2_data_" + dataHierarchy_name + "_fcn_" + fcnHierarchy_name;
    string chi2_penalized_name = "chi2_penalized_data_" + dataHierarchy_name + "_fcn_" + fcnHierarchy_name;
    string data_name           = "AsimovData" + dataHierarchy_name;

    RooAbsPdf *fit_fcn;
    if(isMultipleSpectra)
    {
        fit_fcn = fcnHierarchy == NH ? w_nh->pdf(spectrum_nh.c_str()) : w_ih->pdf(spectrum_ih.c_str());
    }
    else
    {
        fit_fcn = fcnHierarchy == NH ? w_nh->pdf((smeared_spectrum_basename_ + map_nh.begin()->second).c_str()) : \
                                       w_ih->pdf((smeared_spectrum_basename_ + map_ih.begin()->second).c_str());
    }

    RooWorkspace *actual_ws = fcnHierarchy == NH ? w_nh : w_ih;
    RooWorkspace *data_ws = dataHierarchy == NH ? w_nh : w_ih;

    actual_ws->var("e_obs")->setRange("fit_range", 1.1, 8.0);
    actual_ws->var("e_obs")->setRange("full_range", 0.0, 10.0);

    cout << "DEBUG: (*fit_fcn).GetName() = " << (*fit_fcn).GetName() << endl;
    cout << "DEBUG: (*((RooDataHist*)data_ws->data(data_name.c_str()))).GetName() = " 
         << (*((RooDataHist*)data_ws->data(data_name.c_str()))).GetName() << endl;

    RooChi2Var chi2_data(chi2_name.c_str(), chi2_name.c_str(), \
            *fit_fcn, \
            *((RooDataHist*)data_ws->data(data_name.c_str())), \
            DataError(RooAbsData::Poisson), \
            //Range("fit_range"), Verbose()) ;
             Verbose()) ;


    RooFormulaVar *chi2_penalized = getPenalizedChi2(actual_ws, &chi2_data); 

    chi2_penalized->SetName(chi2_penalized_name.c_str());

    
    cout << endl << "************* MINUIT ON CHI2 *************" << endl;
    RooMinuit minimizer(*chi2_penalized) ;
          minimizer.setVerbose();
          minimizer.migrad() ;
          minimizer.hesse() ;

          RooFitResult* r_chi2 = minimizer.save() ;
          cout << endl  << endl << endl << endl << endl << endl;
          cout << "************* FIT RESULT  *************" << endl;
          r_chi2->Print("v");


   

    /***********************
     *                     *
     *    chi2 plotting    *
     *                     *
     ***********************/

    string plot_title =  "#chi^{2}(" + dataHierarchy == NH ? "NH" : "IH";
           plot_title += ", Asimov"  +  fcnHierarchy == NH ? "NH" : "IH";  
           plot_title += ") vs #Deltam^{2}_{31}";

    RooPlot *chi2_frame = actual_ws->var("dm231")->frame(); 
             chi2_frame->SetName((chi2_name + "_frame").c_str());
             chi2_frame->SetTitle(plot_title.c_str());
             chi2_penalized->plotOn(chi2_frame, LineColor(kBlue));
             //nh_chi2_frame->Write();

             output->cd();

    TCanvas *chi2_canvas = new TCanvas((chi2_name + "_canvas").c_str(), (chi2_name + "_canvas").c_str(), 800, 600 );
             chi2_frame->Draw();
             chi2_canvas->SaveAs((path_ + "/" + chi2_name + ".png").c_str());
             chi2_canvas->Write();

}



void
Simulation::fit_ih(const char *nh_filename, const char *ih_filename)
{

    TFile *nh_file = TFile::Open(nh_filename);
    TFile *ih_file = TFile::Open(ih_filename);

    if(!nh_file) {cout << "ERROR: couldn't open " << nh_file << endl; exit (-1);} 
    if(!ih_file) {cout << "ERROR: couldn't open " << ih_file << endl; exit (-1);}

    w_nh = (RooWorkspace*)nh_file->Get("w_nh");
    w_ih = (RooWorkspace*)ih_file->Get("w_ih");

    if(!w_nh) { cout << "ERROR: couldn't find w_nh in " << nh_filename << endl; exit(-1); }
    if(!w_ih) { cout << "ERROR: couldn't find w_ih in " << ih_filename << endl; exit(-1); }


    w_nh->var("baseline_m")->setConstant();
    w_ih->var("baseline_m")->setConstant();



    /******************************************
     *                                        *
     *    chi2 definition and minimization    *
     *                                        *
     ******************************************/

    w_ih->var("e_obs")->setRange("fit_range", 1.1, 8.0);
    w_ih->var("e_obs")->setRange("full_range", 0.0, 10.0);

    RooChi2Var chi2_data_IH("chi2_data_IH","chi2_data_IH", *w_ih->pdf("smeared_spectrum_IH"), *((RooDataHist*)w_nh->data("AsimovDataNH")), DataError(RooAbsData::Poisson), Range("fit_range"), Verbose()) ;
    RooFormulaVar *chi2_penalized_IH = getPenalizedChi2(w_ih, &chi2_data_IH); 
    chi2_penalized_IH->SetName("chi2_penalized_IH");
    
    cout << endl << "************* MINUIT ON CHI2 IH *************" << endl;
    RooMinuit minimizer_ih(*chi2_penalized_IH) ;
          minimizer_ih.migrad() ;
          minimizer_ih.hesse() ;
  
          RooFitResult* r_chi2_ih = minimizer_ih.save() ;
          cout << endl  << endl << endl << endl << endl << endl;
          cout << "************* FIT RESULT IH *************" << endl;
          r_chi2_ih->Print("v");
   

    /***********************
     *                     *
     *    chi2 plotting    *
     *                     *
     ***********************/

    output->cd();

    RooPlot *ih_chi2_frame = w_ih->var("dm231")->frame(); 
             ih_chi2_frame->SetName("ih_chi2_frame");
             ih_chi2_frame->SetTitle("#chi^{2}(IH, AsimovNH) vs #Deltam^{2}_{31}");
             chi2_penalized_IH->plotOn(ih_chi2_frame, LineColor(kRed));
             //ih_chi2_frame->Write();

    TCanvas *ih_chi2_canvas = new TCanvas("ih_chi2_canvas", "ih_chi2_canvas", 800, 600 );
             ih_chi2_frame->Draw();
             ih_chi2_canvas->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/ih_chi2_canvas.png");
             ih_chi2_canvas->Write();
    

    output->cd();

    RooPlot *fitted_IH_spectrum_frame = w_ih->var("e_obs")->frame();
             fitted_IH_spectrum_frame->SetName("fitted_IH_spectrum_frame");

    TCanvas *fitted_IH_spectrum_canvas = new TCanvas("fitted_IH_spectrum_canvas","fitted_IH_spectrum_canvas", 800,600); 
             w_ih->pdf("smeared_spectrum_IH")->plotOn(fitted_IH_spectrum_frame, LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent), Range("fit_range"), NormRange("full_range")); 
             fitted_IH_spectrum_frame->Draw();
             fitted_IH_spectrum_canvas->SaveAs("/scratchfs/dyw/marco/generate_and_fit/fit_v2/NHFittedSpectrum.png");
             fitted_IH_spectrum_canvas->Write();
             //fitted_IH_spectrum_frame->Write();



    nh_file->Close();
    ih_file->Close();
    output->Close();
}




    /*************************************
     *                                   *
     *             OLD STUFF             *
     *                                   *
     *************************************/

void
Simulation::generate_ih_ws() 
{

    /*********************************
     *                               *
     *    create & fill workspace    *
     *                               *
     *********************************/

    w_ih = new RooWorkspace("w_ih","IH oscillation models") ;

    this->fill_common_ws(w_ih);

    /// avoid baseline to be fitted as a floating parameter
    w_ih->var("baseline_m")->setConstant();
    
    w_ih->factory("expr::survival_IH_expr('1-survival_1-survival_2-survival_3-survival_4',{survival_1, survival_2, survival_3, survival_4})");
    w_ih->factory("expr::spectrum_IH_norm('flux_expr * ibd_xsec_expr * survival_IH_expr * spectrum_norm ', {flux_expr, ibd_xsec_expr, survival_IH_expr, spectrum_norm})");
    w_ih->factory("EXPR::survival_IH('1-survival_1-survival_2-survival_3-survival_4',{survival_1, survival_2, survival_3, survival_4})");
    w_ih->factory("PROD::raw_spectrum_IH(survival_IH, reactorSpectrum)");
    w_ih->factory("PROD::smeared_spectrum_2D_IH(response|e_nu, raw_spectrum_IH)");
    w_ih->factory("PROJ::survival_IH_baseline(raw_spectrum_IH, e_nu)");
    w_ih->factory("PROJ::smeared_spectrum_IH(smeared_spectrum_2D_IH, e_nu)");


    /*********************
     *                   *
     *      caching      *
     *                   *
     *********************/

    w_ih->var("baseline_m")->setBins(100, "cache");
    w_ih->var("e_nu")->setBins(500, "cache");
    w_ih->var("e_obs")->setBins(500, "cache");
    w_ih->pdf("raw_spectrum_IH")->setStringAttribute("CACHEPARAMINT","e_nu:baseline_m");
    w_ih->pdf("survival_IH_baseline")->setStringAttribute("CACHEPARAMINT","baseline_m");
    w_ih->pdf("smeared_spectrum_2D_IH")->setStringAttribute("CACHEPARAMINT","e_nu:e_obs");
    w_ih->pdf("smeared_spectrum_IH")->setStringAttribute("CACHEPARAMINT","e_obs");



    /// absolute normalization IH
    w_ih->var("e_nu")->setRange("int_range_ih", 0.1, 9.9);
    RooAbsReal *myInt_IH =  w_ih->function("spectrum_IH_norm")->createIntegral(RooArgSet(*w_ih->var("e_nu")), "int_range_ih");
    double nEvents_IH =  myInt_IH->getVal(); 
    cout << "nevents IH = " << nEvents_IH << endl;


   /***********************
    *                     *
    *  plotting & saving  *
    *                     *
    ***********************/

    output->cd();

    RooPlot *ih_spectrum_frame = w_ih->var("e_obs")->frame();
    w_ih->pdf("smeared_spectrum_IH")->plotOn(ih_spectrum_frame, LineColor(kRed), Normalization(21214, RooAbsReal::NumEvent));
    ih_spectrum_frame->SetName("ih_spectrum_frame");
    ih_spectrum_frame->SetTitle("IH Spectrum");
    ih_spectrum_frame->Write();

    TCanvas *ih_spectrum_canvas = new TCanvas("ih_spectrum_canvas", "ih_spectrum_canvas", 800, 700);
    ih_spectrum_frame->Draw();
    ih_spectrum_canvas->Write();


    w_ih->Write();
    output->Close();

}


void
Simulation::generate_nh_ws() 
{

    /*********************************
     *                               *
     *    create & fill workspace    *
     *                               *
     *********************************/

    w_nh = new RooWorkspace("w_nh","NH oscillation models") ;

    this->fill_common_ws(w_nh);

    /// normal hierarchy
    w_nh->factory("expr::survival_NH_expr('1-survival_1-survival_2-survival_3+survival_4',{survival_1, survival_2, survival_3, survival_4})");
    w_nh->factory("expr::spectrum_NH_norm('flux_expr * ibd_xsec_expr * survival_NH_expr * spectrum_norm ', {flux_expr, ibd_xsec_expr, survival_NH_expr, spectrum_norm})");
    w_nh->factory("EXPR::survival_NH('1-survival_1-survival_2-survival_3+survival_4',{survival_1, survival_2, survival_3, survival_4})");
    w_nh->factory("PROD::raw_spectrum_NH(survival_NH, reactorSpectrum)");
    w_nh->factory("PROD::smeared_spectrum_2D_NH(response|e_nu, raw_spectrum_NH)");
    w_nh->factory("PROJ::survival_NH_baseline(raw_spectrum_NH, e_nu)");
    w_nh->factory("PROJ::smeared_spectrum_NH(smeared_spectrum_2D_NH, e_nu)");


    /*************************
     *                       *
     *      asimov data      *
     *                       *
     *************************/

     w_nh->var("e_obs")->setBins(150, "data_binning");
     RooDataHist *binnedData = w_nh->pdf("smeared_spectrum_NH")->generateBinned(*w_nh->var("e_obs"), NumEvents(21214), Asimov());
     w_nh->import(*binnedData, Rename("AsimovDataNH"));



    /*******************
     *                 *
     *     caching     *
     *                 *
     *******************/

    w_nh->var("baseline_m")->setBins(100, "cache");
    w_nh->var("e_nu")->setBins(500, "cache");
    w_nh->var("e_obs")->setBins(500, "cache");
    w_nh->pdf("raw_spectrum_NH")->setStringAttribute("CACHEPARAMINT","e_nu:baseline_m");
    w_nh->pdf("survival_NH_baseline")->setStringAttribute("CACHEPARAMINT","baseline_m");
    w_nh->pdf("smeared_spectrum_2D_NH")->setStringAttribute("CACHEPARAMINT","e_nu:e_obs");
    w_nh->pdf("smeared_spectrum_NH")->setStringAttribute("CACHEPARAMINT","e_obs");


   //   w->pdf("model")->setNormValueCaching(3) ;
   //   // Evaluate p.d.f. once to trigger filling of cache
   //   RooArgSet normSet(*w->var("x"),*w->var("y"),*w->var("z")) ;
   //   w->pdf("model")->getVal(&normSet) ;
   //   w->writeToFile("rf903_numintcache.root") ;


    /*************************
     *                       *
     *     baseline plot     *
     *                       *
     *************************/

    RooPlot *b_frame = w_nh->var("baseline_m")->frame();
    b_frame->SetName("baseline_frame");
    b_frame->SetTitle("Survival Prob (NH) vs Baseline");
    w_nh->pdf("survival_NH_baseline")->plotOn(b_frame);
    ///  after making this plot, fix the baseline 
    w_nh->var("baseline_m")->setConstant();


    /// absolute normalization NH
    w_nh->var("e_nu")->setRange("int_range_nh", 0.1, 9.9);
    RooAbsReal *myInt_NH =  w_nh->function("spectrum_NH_norm")->createIntegral(RooArgSet(*w_nh->var("e_nu")), "int_range_nh");
    double nEvents_NH =  myInt_NH->getVal(); 
    cout << "nevents NH = " << nEvents_NH << endl;



   /***********************
    *                     *
    *  plotting & saving  *
    *                     *
    ***********************/

    output->cd();

    RooPlot *nh_spectrum_nh_data_frame = w_nh->var("e_obs")->frame();
             nh_spectrum_nh_data_frame->SetTitle("Spectrum NH");
             nh_spectrum_nh_data_frame->SetName("nh_spectrum_nh_data_frame");

    binnedData->plotOn(nh_spectrum_nh_data_frame);
    w_nh->pdf("smeared_spectrum_NH")->plotOn(nh_spectrum_nh_data_frame, LineColor(kBlue), Normalization(21214, RooAbsReal::NumEvent));

    nh_spectrum_nh_data_frame->Write();

    TCanvas *nh_spectrum_nh_data_canvas = new TCanvas("nh_spectrum_nh_data_canvas", "nh_spectrum_nh_data_canvas", 800, 700);
    nh_spectrum_nh_data_frame->Draw();
    nh_spectrum_nh_data_canvas->Write();

    TCanvas *nh_baseline_canvas = new TCanvas("nh_baseline_canvas", "nh_baseline_canvas", 800, 700);
    b_frame->Draw();
    b_frame->Write();
    nh_baseline_canvas->Write();

    w_nh->Write();

    output->Close();



    /***********************
     *                     *
     *    old canvases     *
     *                     *
     ***********************/

/*
    RooPlot *b_frame  = baseline_m.frame();
    RooPlot *b_frame2 = baseline_m.frame();

    RooPlot *e_frame  = e_nu.frame();
    RooPlot *e_frame2 = e_nu.frame();
    RooPlot *e_frame3 = e_nu.frame();
    RooPlot *e_frame4 = e_nu.frame();
    RooPlot *e_frame5 = e_nu.frame();
    RooPlot *e_frame6 = e_nu.frame();

    RooPlot *obs_frame = e_obs.frame();
    RooPlot *obs_frame2 = e_obs.frame();
    RooPlot *obs_frame3 = e_obs.frame();
    RooPlot *obs_frame4 = e_obs.frame();


    pdf_NH_2.plotOn(e_frame);
    flux.plotOn(e_frame2);
    ibd_xsec.plotOn(e_frame3);
    reactorSpectrum.plotOn(e_frame4);
    survival_NH.plotOn(e_frame5);
    raw_spectrum_NH.plotOn(e_frame6, LineColor(kRed));

    survival_NH.plotOn(b_frame);
    survival_NH_baseline->plotOn(b_frame2);
    raw_spectrum_NH.plotOn(b_frame2, LineColor(kRed), LineStyle(kDashed), ProjectionRange("A"));

    // pdf_NH_smear_cont->plotOn(obs_frame, LineColor(kRed));
    // response_proj->plotOn(obs_frame2);
    smeared_spectrum_NH->plotOn(obs_frame3);
    smeared_spectrum_IH->plotOn(obs_frame4, LineColor(kRed));

    //pdf_NH_int->plotOn(b_frame);

    TCanvas *c8 = new TCanvas("c8","c8",1024,700);
    c8->Divide(2,2);
    c8->cd(1);
    e_frame->Draw();
    c8->cd(2);
    e_frame2->Draw();
    c8->cd(3);
    e_frame3->Draw();
    c8->cd(4);
    e_frame4->Draw();


    TCanvas *c9 = new TCanvas("c9", "c9", 1024, 700);
    c9->Divide(2,2);
    c9->cd(1);
    e_frame5->Draw();
    c9->cd(2);
    e_frame6->Draw();
    c9->cd(3);
    //b_frame->Draw();
    reactor_b_h->Draw("lego");
    c9->cd(4);
    b_frame2->Draw();
    //b_frame->Draw();



    // TCanvas *c10 = new TCanvas("c10", "c10", 1024, 700);
    // c10->Divide(2,2);
    // c10->cd(1);
    // response_h->Draw("lego");
    // c10->cd(2);
    // /// response projection
    // obs_frame2->Draw();
    // c10->cd(3);
    // /// smeared reactor spectrum
    // e_frame4->Draw();
    // obs_frame->Draw("same");
    // c10->cd(4);
    // reactor_h->Draw("lego");

    TCanvas *c11 = new TCanvas("c11", "c11", 1024, 700);
    obs_frame3->Draw();
    e_frame6->Draw("same");

    TCanvas *c12 = new TCanvas("c12", "c12", 1024, 700);
    obs_frame4->Draw();
    obs_frame3->Draw("same");


    output->cd();
    c8->Write();
    c9->Write();
    c11->Write();
    c12->Write();


    output->Close();

    return 0;
    */
}



void 
Simulation::fill_common_ws(RooWorkspace *ws)
{
    ws->factory("myPi[3.1415926535897]");
	ws->factory("proton_fraction[0.12]");
	ws->factory("reactor_power[20.]"); // GW
	ws->factory("baseline_m[50000., 0., 100000.]");
	ws->factory("detector_mass[5.]"); // Kt
	ws->factory("exposure_time[5.]"); // years
    ws->factory("sin22teta12[0.857, 0., 1.]");
    ws->factory("sin22teta12_fix[0.857]");
    ws->factory("sin22teta12_err[0.024]");
    ws->factory("sin22teta13[0.089, 0., 1.]");
    ws->factory("sin22teta13_fix[0.089]");
    ws->factory("sin22teta13_err[0.005]");
    ws->factory("dm221[0.000075, 0., 0.0005]");
    ws->factory("dm221_fix[0.000075]");
    ws->factory("dm221_err[0.000002]");
    ws->factory("dm231[0.00232, 0., 0.005]");
    ws->factory("dm231_fix[0.00232]");
    ws->factory("dm231_err[0.0001]");
    ws->factory("sys_fit[0.03, 0.0, 2.0]");
    ws->factory("sys_err[0.03]");
    ws->factory("e_nu[1.0, 10.0]");
    ws->factory("e_obs[0.0, 10.0]");
    ws->factory("weighted_fission_energy[0.58*201.7+0.3*210.+0.07*205+0.05*212.4]");

    ws->factory("RooPolyVar::e_vis(e_nu, {e_shift[-0.8], e_ratio[1.0]})");

    ws->factory("expr::exposure_time_seconds('exposure_time*220.0*86400.0',{exposure_time})");
    ws->factory("expr::flux_a('exp(0.87 -0.16 *e_nu-0.091 *(e_nu^2))',{e_nu})");
    ws->factory("expr::flux_b('exp(0.896-0.239*e_nu-0.0981*(e_nu^2))',{e_nu})");
    ws->factory("expr::flux_c('exp(0.976-0.162*e_nu-0.079 *(e_nu^2))',{e_nu})");
    ws->factory("expr::flux_d('exp(0.793-0.08 *e_nu-0.1085*(e_nu^2))',{e_nu})");
    ws->factory("expr::sin2theta12('0.5*(1.-sqrt(1.-sin22teta12))',{sin22teta12})"); 
    ws->factory("expr::cos2theta13('0.5*(1.+sqrt(1.-sin22teta13))',{sin22teta13})");
    ws->factory("expr::delta21('1.27*dm221*baseline_m/e_nu', {dm221, e_nu, baseline_m})");
    ws->factory("expr::delta31('1.27*dm231*baseline_m/e_nu', {dm231, e_nu, baseline_m})");
    ws->factory("expr::survival_1('(e_nu>1.8) * cos2theta13^2 * sin22teta12 * (sin(delta21))^2', {e_nu, cos2theta13, sin22teta12, delta21})");
    ws->factory("expr::survival_2('(e_nu>1.8) * sin22teta13* (sin(abs(delta31)))^2', {e_nu, sin22teta13, delta31})");
    ws->factory("expr::survival_3('(e_nu>1.8) * sin2theta12 * sin22teta13 * (sin(delta21))^2 * cos(2.0*abs(delta31))', {e_nu, sin2theta12, sin22teta13, delta21, delta31})");
    ws->factory("expr::survival_4('(e_nu>1.8) * 0.5 * sin2theta12 * sin22teta13 * sin(2.0*delta21) * sin(2.0*abs(delta31))', {e_nu, sin2theta12, sin22teta13, delta21, delta31})");
    ws->factory("expr::spectrum_norm('detector_mass *1000000. * proton_fraction * exposure_time_seconds / (1.6726231 * 4.0 * myPi * baseline_m^2)', {detector_mass, proton_fraction, exposure_time_seconds, myPi, baseline_m})");
    ws->factory("expr::resolution('0.03 * sqrt(e_vis)', {e_vis})");

    /// needed for normalization
    ws->factory("expr::ibd_xsec_expr('(e_nu>1.8) * 9.52*(e_nu-1.29)*sqrt((e_nu-1.29)^2 - 0.511^2)', {e_nu})");
	ws->factory("expr::flux_expr('0.58*flux_a + 0.3*flux_b + 0.07*flux_c + 0.05*flux_d',{flux_a, flux_b, flux_c, flux_d})");

    /// smearing
    ws->factory("Gaussian::response(e_obs, e_vis, resolution)");

    /// pdfs merging all the previous formulas
    ws->factory("EXPR::ibd_xsec('(e_nu>1.8) * 9.52*(e_nu-1.29)*sqrt((e_nu-1.29)^2 - 0.511^2)', {e_nu})");
	ws->factory("EXPR::flux('0.58*flux_a + 0.3*flux_b + 0.07*flux_c + 0.05*flux_d',{flux_a, flux_b, flux_c, flux_d})");

    ws->factory("PROD::reactorSpectrum(flux, ibd_xsec)");
}








    /*************************************
     *                                   *
     *     PENALIZED CHI2 DEFINITION     *
     *                                   *
     *************************************/

 //    string penalizedChi2("RooChi2Var::chi2p('chi2_data + ");
 //           penalizedChi2 += "((sin22teta12-sin22teta12_fix)/sin22teta12_err)^2 + ";
 //           penalizedChi2 += "((sin22teta13-sin22teta13_fix)/sin22teta13_err)^2 + ";
 //           penalizedChi2 += "((dm221-dm221_fix)/dm221_err)^2 + ";
 //           penalizedChi2 += "((dm231-dm231_fix)/dm231_err)^2 + ";
 //           penalizedChi2 += "((sys_fit-1.0)/sys_err)^2', ";
 //           penalizedChi2 += "{sin22teta12, sin22teta12_fix, sin22teta12_err, ";
 //           penalizedChi2 += "sin22teta13,sin22teta13_fix, sin22teta13_err, ";
 //           penalizedChi2 += "dm221, dm221_fix, dm221_err, ";
 //           penalizedChi2 += "dm231, dm231_fix, dm231_err, ";
 //           penalizedChi2 += "sys_fit, sys_err})";
 // 
 //
 // w_nh->factory(penalizedChi2.c_str());





 



    /************************************
     *                                  *
     *     CACHING in ROOT v5.34/07     *
     *                                  *
     ************************************/

    /*
     
    w_nh->var("baseline_m")->setBins(100, "cache");
    w_nh->var("e_nu")->setBins(500, "cache");
    w_nh->var("e_obs")->setBins(500, "cache");
    w_nh->pdf("raw_spectrum_NH")->setStringAttribute("CACHEPARAMINT","e_nu:baseline_m");
    w_nh->pdf("survival_NH_baseline")->setStringAttribute("CACHEPARAMINT","baseline_m");
    w_nh->pdf("smeared_spectrum_2D_NH")->setStringAttribute("CACHEPARAMINT","e_nu:e_obs");
    w_nh->pdf("smeared_spectrum_NH")->setStringAttribute("CACHEPARAMINT","e_obs");

    w_ih->var("baseline_m")->setBins(100, "cache");
    w_ih->var("e_nu")->setBins(500, "cache");
    w_ih->var("e_obs")->setBins(500, "cache");
    w_ih->pdf("raw_spectrum_IH")->setStringAttribute("CACHEPARAMINT","e_nu:baseline_m");
    w_ih->pdf("survival_IH_baseline")->setStringAttribute("CACHEPARAMINT","baseline_m");
    w_ih->pdf("smeared_spectrum_2D_IH")->setStringAttribute("CACHEPARAMINT","e_nu:e_obs");
    w_ih->pdf("smeared_spectrum_IH")->setStringAttribute("CACHEPARAMINT","e_obs");


    RooArgSet normSet(*w_nh->var("baseline_m"));
    cout << "baseline norm: " << w_nh->pdf("survival_NH_baseline")->getVal(&normSet) << endl;

    RooArgSet normSet2(*w_nh->var("e_obs"));
    cout << "spectrum norm: " << w_nh->pdf("smeared_spectrum_NH")->getVal(&normSet2) << endl;

    */




    /************************************
     *                                  *
     *    STORE CHI2 IN WS : TRIALS     *
     *                                  *
     ************************************/


    /***************************************
     *            test 1:                  *
     *   chi2 built within the workspace   *
     ***************************************/

    //w_nh->factory("Gaussian::gg(e_obs, mm[0], ss[1])");
    //RooDataHist *binnedData = w_nh->pdf("gg")->generateBinned(*w_nh->var("e_obs"), NumEvents(100000), Asimov());
    //w_nh->import(*binnedData, Rename("AsimovDataNH"));
    //w_nh->factory("chi2::chi2_data(gg, AsimovDataNH)");


    /***************************************
     *            test 2:                  *
     *   standalone chi2 imported in ws    *
     ***************************************/

    // RooGaussian gg("gg", "gg", *w_nh->var("e_obs"), RooConst(0), RooConst(1));
    // RooDataHist *binnedData = gg.generateBinned(*w_nh->var("e_obs"), NumEvents(100000), Asimov());
    // w_nh->import(*binnedData, Rename("AsimovDataNH"));
    // RooChi2Var chi2_data("chi2_data","chi2_data", gg, *binnedData, DataError(RooAbsData::Poisson)) ;
    // w_nh->import(chi2_data);


    /***************************************
     *            test 3:                  *
     *       no storage in workspace       *
     ***************************************/

     // RooGaussian gg("gg", "gg", *w_nh->var("e_obs"), RooConst(0), RooConst(1));
     // RooDataHist *binnedData = gg.generateBinned(*w_nh->var("e_obs"), NumEvents(100000), Asimov());
     // w_nh->import(*binnedData, Rename("AsimovDataNH"));
     // RooChi2Var chi2_data("chi2_data","chi2_data", gg, *binnedData, DataError(RooAbsData::Poisson)) ;



    /***************************************
     *            test 4:                  *
     *    data from smeared NH spectrum    *
     ***************************************/

    //RooDataHist* binnedData = w_nh->pdf("smeared_spectrum_NH")->generateBinned(*w_nh->var("e_obs"), NumEvents(100000), Asimov()); //data->binnedClone() ;


    /*****************************************
     *  not working because roochi2var       *
     *  cannot be built using RooArgs        *
     *  within the factory method            *
     *****************************************/

    //w_nh->factory("RooChi2Var::chi2(smeared_spectrum_NH, AsimovDataNH, DataError(RooAbsData::Poisson))");
    //w_nh->factory("chi2_data(gg, AsimovDataNH, DataError(RooAbsData::Poisson))");
