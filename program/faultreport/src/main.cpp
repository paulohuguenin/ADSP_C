

#include "datasignal.h"
#include "pll.h"

using namespace std;

int main(int argc, char *argv[])
{
    // --------------------------------------------
    // program - faultreport
    // argv[1] - file name
    // --------------------------------------------

    // Define the ".cfg" file path
    // _MAX_PATH defined in stdlib.h
    char InputFile[_MAX_PATH];
    double rateTarget;
    char auxStr[_MAX_PATH];
    if (argc > 1)
    {
        strcpy( InputFile, argv[1]);
        
    }
    else
    {
        //strcpy( InputFile, "../../osc/x001.cfg");
    }
    
    //////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;

    dataSignal = new CComtradeSignal;
    dataSignal->setFileName(InputFile);

    //dataSignal->setBlockSize(genData->getBlockSize());
    //dataSignal->setBlockHop(genData->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // -----------------------------------------------------------
    // run PLL
    //        
    double Fs;
    
    Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
    
    
    cout << "Fs:" << Fs << endl;
    
    int numSignals = dataSignal->getNumSignal();
    int sigSize = dataSignal->getSignalSize();
    
    double** signal;
    signal = dataSignal->getSignal();
    
    int i,j;
        
    
        
    fstream file_report("pllreport.out",ios::out);
    file_report <<  setw (10) << setfill(' ') << "Signal"    << " "
                <<  setw (10) << setfill(' ') << "Sample"    << " "
                <<  setw (15) << setfill(' ') << "Value"   << " "
                <<  setw (15) << setfill(' ') << "VEstim"   << " "
                <<  setw (15) << setfill(' ') << "FreqEstim"<< " "
                <<  setw (15) << setfill(' ') << "AmpEstim" << " "
                <<  setw (15) << setfill(' ') << "WtEstim"  << " "
                << endl;    

    for (i=0; i<numSignals; i++)
    {
        double* V_estim = new double[sigSize];
        double* freq_estim = new double[sigSize];
        double* amp_estim = new double[sigSize]; 
        double* wt_estim =new double[sigSize];
        cout << "Run PLL channel: " <<  i+1 << endl;
        
        pll_lovo__(signal[i], V_estim, sigSize, freq_estim, 
                    amp_estim, wt_estim,Fs);
        for (j=0; j<sigSize; j++)
        {   
            file_report <<  setw (10) << setfill(' ') << i+1    << " "
                        <<  setw (10) << setfill(' ') << j+1    << " "
                        <<  setw (15) << setfill(' ') << signal[i][j]   << " "
                        <<  setw (15) << setfill(' ') << V_estim[j]   << " "
                        <<  setw (15) << setfill(' ') << freq_estim[j]<< " "
                        <<  setw (15) << setfill(' ') << amp_estim[j] << " "
                        <<  setw (15) << setfill(' ') << wt_estim[j]  << " "    
                        << endl;
        }
                    
        delete [] V_estim;
        delete [] freq_estim;
        delete [] amp_estim; 
        delete [] wt_estim;
    }
    file_report.close();
        
    delete dataSignal;
    
    return 0;
}