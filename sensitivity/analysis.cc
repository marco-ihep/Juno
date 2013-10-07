#include <iostream>
#include <string>
#include "simulation.h"

using namespace std;

int main(int argc, char *argv[] )
{
    if(argc != 3) {
        cout << "usage: " << argv[0] << " mode path" << endl;
        cout << "   mode is: " << endl;
        cout << "   - gen_nh: " << endl;
        cout << "   - gen_ih: " << endl;
        cout << "   - fit_nh: " << endl;
        cout << "   - fit_ih: " << endl;
        exit(-1);
    }

    // string nh_file("/scratchfs/dyw/marco/generate_and_fit/class_generate_nh/analysis_submission_cache_500_nEvents_21214/output_nh.root");
    // string ih_file("/scratchfs/dyw/marco/generate_and_fit/class_generate_ih/analysis_submission_cache_500_nEvents_21214/output_ih.root");
    string nh_file("/scratchfs/dyw/marco/generate_and_fit/fit_v2/gen_nh_submission/output_nh.root");
    string ih_file("/scratchfs/dyw/marco/generate_and_fit/fit_v2/gen_ih_submission/output_ih.root");

    Simulation sim(argv[2]);
    string mode(argv[1]);

    objectsToPrint objs_nh = sim.generate_multiple_spectra(NH);
    cout << endl << "********** DONE WITH GENERATION NH **********" << endl << endl ;
    // sim.drawBaseline(objs_nh);
    sim.drawSpectra(objs_nh);


    //objectsToPrint objs_ih = sim.generate_multiple_spectra(IH);
    //cout << endl << "********** DONE WITH GENERATION IH **********" << endl << endl ;
    // sim.drawBaseline(objs_ih);
    // sim.drawSpectra(objs_ih);

    
    // cout << endl << endl << "********** FITTING NH DATA with NH FCN **********" << endl;
    // sim.fit(objs_nh, objs_ih, NH, NH);
    // 

    // cout << endl << endl << "********** FITTING NH DATA with IH FCN**********" << endl;
    // sim.fit(objs_nh, objs_ih, NH, IH);


    // cout << "********** PLOTTING JOINT SPECTRA **********" << endl;
    // sim.drawJointHierachies(objs_nh, objs_ih);

    return 0;
}
