#include "mpursuit.h"

void matchingPursuit(   cgMatrix<double>& residue,
                        strtParameter* chosenParm,
                        int dicSize,
                        CDataSignal* dataSignal,
                        CFileDictionary* dicData,
                        CFileDecomp* genData, int step)
{
    int decincAsymmFlag;
    int i,j;
    int k;
    int s;
    int chosenTau = 0;
    double chosenXi = 0.0;
    int delta_tau;
    double chosenOptPhase = 0.0;
    double Fs;
    // double Ffund;
    // double delta_f = 0;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    int dicType;
    // int signalSize;
    double  innerProd;
    double* decRealAtom;
    double* incRealAtom;
    double* realAtomWin = new double[dicSize];
    CComplex* complexAtomXi = new CComplex[dicSize];
    CComplex* complexAtom2Xi = new CComplex[dicSize];
    CDictionary* dic;

    strtParameter* parm = new strtParameter;
    // strtParameter* chosenParm = new strtParameter;
    ((strtParameter*)chosenParm)->innerProduct=0.0;
    ((strtParameter*)chosenParm)->rho = 0.0;
    ((strtParameter*)chosenParm)->xi = 0.0;
    ((strtParameter*)chosenParm)->phase = 0.0;
    ((strtParameter*)chosenParm)->a = 0;
    ((strtParameter*)chosenParm)->b = 0;
    double xi;
    int tau;
    double maxInnerProd = 0.0;

    int MPType;
    MPType = genData->getMPType();


    FILE* stream;   
    char* fileName;

    // if (MPType == 1)
    // {
    //     fileName = "MPTradicional.out";
        
    // }

    // if (MPType == 2)
    // {
    //     fileName = "fastMPKolasa.out";
        
    // }

    // else if (MPType == 3)
    // {
    //     fileName = "fastMPKolasaModified.out";
    // }

    for(j=0;j<dicData->getNumDicBlock();j++)
    {   
        // stream = fopen( fileName, "a" );
        // fprintf(stream,"%d --------------------------------------------------------------\n", j);
        // fprintf(stream,"Coef.           Rho             Xi              Phase           tau             a               b               innerProd_xp    innerProd_xq    innerProd_pp    innerProd_qq    innerProd_pq    \n");
        // fflush(stream);
        // fclose(stream);
        // Using the dictionary data
        s = (dicData->getDicBlock())[j].scale;
        delta_tau = (dicData->getDicBlock())[j].delta_tau;
        fdiscrtype = (dicData->getDicBlock())[j].fdiscrtype;
        freqi = (dicData->getDicBlock())[j].freqi;
        freqf = (dicData->getDicBlock())[j].freqf;
        dicType = (dicData->getDicBlock())[j].dicType;

        int N = dicSize;

        if ((s > dicSize)&&(s!=88888))
        {
            cout << "Scale greater than block size!!!" << endl;
            printf("Scale %d - block size %d\n",s,dicSize);
            exit(1);
        }
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////

        if (dataSignal->getType() == 1)
        {
            Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
        }
        else if (dataSignal->getType() == 2) // Audio
        {
            Fs = ((CAudioSignal*)dataSignal)->getSamplingRate();
        }
        else if (dataSignal->getType() == 3) // Noise for Audio
        {
            Fs = ((CNoiseSignal*)dataSignal)->getSamplingRate();
        }

        if ( (freqf==9999999999) && (freqi==9999999999) )
        {
            if (dataSignal->getType() == 1)
            {
                freqi = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
                freqf = Fs;
            }
            else if (dataSignal->getType() == 2) // Audio
            {
                freqi =  27.5; // A0 = 27.5 Hz
                freqf = Fs;
                //double A0 = 27.5; // Hz
                //double C8 = 4186.01;
                //double freqi= A0 * pow(2.0,-23/12);
                //double freqf= C8 * pow(2.0,23/12); // B9
            }
            else if (dataSignal->getType() == 3) // Noise for Audio
            {
                freqi =  Fs/100;
                freqf = Fs;
            }

        }
        else if (freqf==9999999999)
        {
            freqf = Fs;
        }
        else if (freqi==9999999999)
        {
            freqi = freqf/s;
        }

        if (freqf>Fs)
        {
            cout << "Final frequency greater than sampling frequency !!!" << endl;
            printf("final freq. %f - Fs %f\n",freqf,Fs);
            exit(1);
        }

        if (fdiscrtype==1) // linear
        {
            // delta_f = (2*pi/freqf)* freqi;
            Nfreq = (int)(freqf/(2*freqi));
            //Nfreq = (int)ceil(freqf/freqi);
            xi_vec = new double[Nfreq];
            for (i=0;i<Nfreq;i++)
            {
                // xi_vec[i] = delta_f * i;
                xi_vec[i] = (2*pi/Fs) * (freqi * i );
            }
        }
        else if (fdiscrtype==2) // geometric with quarter-tone discretization
        {
            Nfreq = (int)ceil(24 * ( log10(freqf/freqi)/log10(2.0) ) )+1;

            xi_vec = new double[Nfreq];
            xi_vec[0] = 0.0;
            for (i=1;i<Nfreq;i++)
            {
                xi_vec[i] = (2*pi/Fs) * (freqi * pow (2.0, (double)(i-1)/24) );
            }
        }

        if (s == 1) // impulse
        {
            for(i=0;i<N;i++)
            {
                innerProd = residue[0][i];
                if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                {
                    maxInnerProd = innerProd;
                    setParameters ( chosenParm,
                                    dicType,
                                    maxInnerProd,
                                    s,
                                    0.0,
                                    0.0,
                                    i,
                                    i);
                    // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                    // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                    // ((strtParameter*)chosenParm)->xi = 0.0;
                    // ((strtParameter*)chosenParm)->phase = 0.0;
                    // ((strtParameter*)chosenParm)->a = i;
                    // ((strtParameter*)chosenParm)->b = i; 
                }
            }
        }

        else
        {
            if (MPType == 1)
            {
                ((strtParameter*)chosenParm)->dicType = dicType;
                dic = new CExpDictionary;
                ((CExpDictionary*)dic)->setSignalSize(dicSize);

                decincAsymmFlag = 1;
                MPTradicional(  residue,
                                maxInnerProd,
                                chosenOptPhase,
                                chosenTau,
                                chosenXi,
                                dic,
                                N, //int dicSize
                                decincAsymmFlag,
                                s,
                                delta_tau,
                                xi_vec,
                                Nfreq,
                                fileName);
                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                {
                    setParameters ( chosenParm,
                                    dicType,
                                    maxInnerProd,
                                    s,
                                    chosenXi,
                                    chosenOptPhase,
                                    chosenTau,
                                    N-1);
                    // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                    // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                    // ((strtParameter*)chosenParm)->xi = chosenXi;
                    // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                    // ((strtParameter*)chosenParm)->a = chosenTau;
                    // ((strtParameter*)chosenParm)->b = N-1;
                }

                // decincAsymmFlag = -1;
                // MPTradicional(  residue,
                //                 maxInnerProd,
                //                 chosenOptPhase,
                //                 chosenTau,
                //                 chosenXi,
                //                 dic,
                //                 N, //int dicSize
                //                 decincAsymmFlag,
                //                 s,
                //                 delta_tau,
                //                 xi_vec,
                //                 Nfreq,
                //                 fileName);

                // if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                // {
                //     ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                //     ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                //     ((strtParameter*)chosenParm)->xi = chosenXi;
                //     ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                //     ((strtParameter*)chosenParm)->a = 0.0;
                //     ((strtParameter*)chosenParm)->b = chosenTau;
                // }
            }

            else if (MPType==2)
            {
                ((strtParameter*)chosenParm)->dicType = dicType;
                dic = new CExpDictionary;
                ((CExpDictionary*)dic)->setSignalSize(dicSize);

                int delta_tau = s;
                for (tau=0; tau<(double)N; tau=tau+delta_tau)
                {
                    //decreasing
                    ((strtParameter*)parm)->rho = (double)(1.0/(double)s);
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->a = tau;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)dic)->setRealAtom(parm);
                    decRealAtom = dic->getRealAtom();
                    // cout << (((strtParameter*)parm)->rho) << endl;
                    // cout << *decRealAtom << endl;

                    fastMPKolasa(   residue,
                                    maxInnerProd,
                                    chosenOptPhase,
                                    chosenXi,
                                    dicSize,
                                    tau,
                                    s,
                                    decRealAtom,
                                    fileName);

                    if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    {
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        chosenXi,
                                        chosenOptPhase,
                                        tau,
                                        N-1);
                        // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                        // ((strtParameter*)chosenParm)->xi = chosenXi;
                        // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                        // ((strtParameter*)chosenParm)->a = tau;
                        // ((strtParameter*)chosenParm)->b = N-1; 
                        // cout << ((strtParameter*)chosenParm)->rho << endl;
                        // cout << "AQUI" << endl;
                    }

                    //increasing
                    // ((strtParameter*)parm)->rho = (double)(1/s);
                    // ((strtParameter*)parm)->xi = 0.0;
                    // ((strtParameter*)parm)->phase = 0.0;
                    // ((strtParameter*)parm)->a = 0.0;
                    // ((strtParameter*)parm)->b = tau;
                    // ((CExpDictionary*)dic)->setRealAtom(parm);
                    // incRealAtom = dic->getRealAtom();
                    
                    // fastMPKolasa(   residue,
                    //                 maxInnerProd,
                    //                 chosenOptPhase,
                    //                 chosenXi,
                    //                 dicSize,
                    //                 tau,
                    //                 s,
                    //                 incRealAtom,
                    //                 fileName);

                    // if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    // {
                    //     ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                    //     ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                    //     ((strtParameter*)chosenParm)->xi = chosenXi;
                    //     ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                    //     ((strtParameter*)chosenParm)->a = 0.0;
                    //     ((strtParameter*)chosenParm)->b = tau; 
                    // }
                }


            }

            else if (MPType==3)
            {
                ((strtParameter*)chosenParm)->dicType = dicType;
                dic = new CExpDictionary;
                ((CExpDictionary*)dic)->setSignalSize(dicSize);

                ((strtParameter*)parm)->rho = 1.0/(double)s;
                ((strtParameter*)parm)->xi = 0.0;
                ((strtParameter*)parm)->phase = 0.0;
                ((strtParameter*)parm)->a = 0;
                ((strtParameter*)parm)->b = N-1; 
                ((CExpDictionary*)dic)->setRealAtom(parm);
                // cout << "AQUI" << endl;
                // realAtomWin = dic->getRealAtom();
                memcpy(realAtomWin, (dic->getRealAtom()), sizeof(double)*dicSize);

                // for (int k1=0; k1<dicSize;k1++)
                // {
                //     realAtomWin[k1] = *(dic->getRealAtom()+k1);
                // }

                for(k=0;k<Nfreq;k++)
                {
                    //xi = k * ((2*pi)/s);
                    xi = xi_vec[k];
                    // cout << xi_vec[0] << endl << endl;
                    if ( (xi==0.0) || (xi>=((2*pi)/s) ))
                    {
                        // z1
                        ((strtParameter*)parm)->rho = 0.0;
                        ((strtParameter*)parm)->xi = xi;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1; 
                        ((CExpDictionary*)dic)->setComplexAtom(parm);

                        memcpy(complexAtomXi, (dic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                        // for (int k1=0; k1<dicSize;k1++)
                        // {
                        //     complexAtomXi[k1] = *(dic->getComplexAtom()+k1);
                        // }

                        // z2
                        ((strtParameter*)parm)->rho = 0.0;
                        ((strtParameter*)parm)->xi = 2*xi;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1; 
                        ((CExpDictionary*)dic)->setComplexAtom(parm);

                        memcpy(complexAtom2Xi, (dic->getComplexAtom()), sizeof(CComplex) * (dicSize));


                        // for (int k1=0; k1<dicSize;k1++)
                        // {
                        //     complexAtom2Xi[k1] = *(dic->getComplexAtom()+k1);
                        // }


                        //Decreasing exponencial
                        decincAsymmFlag = 1;
                        fastMPKolasaModified(   residue,
                                                maxInnerProd,
                                                chosenOptPhase,
                                                chosenTau,
                                                dicSize,
                                                decincAsymmFlag,
                                                delta_tau,
                                                xi,
                                                realAtomWin,
                                                complexAtomXi,
                                                complexAtom2Xi,
                                                fileName,
                                                s);
                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {   
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            xi,
                                            chosenOptPhase,
                                            chosenTau,
                                            N-1);
                            // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                            // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                            // ((strtParameter*)chosenParm)->xi = xi;
                            // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                            // ((strtParameter*)chosenParm)->a = chosenTau;
                            // ((strtParameter*)chosenParm)->b = N-1; 
                        }
                        // cout << *realAtomWin << endl;
                    
                        
                        //Increasing exponencial
                        // decincAsymmFlag = -1;
                        // fastMPKolasaModified(   residue,
                        //                         maxInnerProd,
                        //                         chosenOptPhase,
                        //                         chosenTau,
                        //                         dicSize,
                        //                         decincAsymmFlag,
                        //                         delta_tau,
                        //                         xi,
                        //                         realAtomWin,
                        //                         complexAtomXi,
                        //                         complexAtom2Xi,
                        //                         fileName);
                        // if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        // {
                        //     ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        //     ((strtParameter*)chosenParm)->rho = -1.0/(double)s;
                        //     ((strtParameter*)chosenParm)->xi = xi;
                        //     ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                        //     ((strtParameter*)chosenParm)->a = 0;
                        //     ((strtParameter*)chosenParm)->b = chosenTau; 
                        // }
                    }
                }   
            }
        }

        if (xi_vec!=NULL) 
        {
            delete [] xi_vec;
        }
        if (dic!= NULL) 
        {
            delete dic;
        }
    }
    if (parm!=NULL) 
    {
        delete [] parm; 
    }
    delete complexAtomXi;
    delete complexAtom2Xi;

    return;
}

void MPTradicional( cgMatrix<double>& residue,
                    double& maxInnerProd,
                    double& chosenOptPhase,
                    int& chosenTau,
                    double& chosenXi,
                    CDictionary* dic,
                    int N, //int dicSize
                    int decincAsymmFlag,
                    int s,
                    int delta_tau,
                    double* xi_vec,
                    int Nfreq,
                    char* fileName)
{
    int k;
    double opt_phase;
    CComplex* complexAtom;
    int tau;
    double xi;
    strtParameter* parm;
    parm = new strtParameter;
    double innerProd = 0;

    for(k=0;k<Nfreq;k++)
    {
        // cout << "k - " << k << endl;
        for(tau=0;tau<N;tau+=delta_tau)
        {   
            if (decincAsymmFlag == 1)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];
                ((strtParameter*)parm)->rho = (1.0/(double)s);
                ((strtParameter*)parm)->xi = xi;
                ((strtParameter*)parm)->phase = 0.0;
                ((strtParameter*)parm)->a = tau;
                ((strtParameter*)parm)->b = N-1; 
                ((CExpDictionary*)dic)->setComplexAtom(parm);
                complexAtom = dic->getComplexAtom();
                computeOptimumPhase(residue,opt_phase,innerProd,N,xi,complexAtom, fileName,s,tau,N);


                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    chosenXi = xi;
                    chosenOptPhase = opt_phase;
                    chosenTau = tau;
                }
            }
            if (decincAsymmFlag == -1)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];
                ((strtParameter*)parm)->rho = -(1.0/(double)s);
                ((strtParameter*)parm)->xi = xi;
                ((strtParameter*)parm)->phase = 0.0;
                ((strtParameter*)parm)->a = 0.0;
                ((strtParameter*)parm)->b = tau; 
                ((CExpDictionary*)dic)->setComplexAtom(parm);
                complexAtom = dic->getComplexAtom();

                computeOptimumPhase(residue,opt_phase,innerProd,N,xi,complexAtom, fileName,s,tau,N);


                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    chosenXi = xi;
                    chosenOptPhase = opt_phase;
                    chosenTau = tau;
                }

            }
        }
    }
    delete parm;
    return;
}

void fastMPKolasa(  cgMatrix<double>& residue,
                    double& maxInnerProd,
                    double& chosenOptPhase,
                    double& chosenXi,
                    int N, //dicSize
                    int tau,
                    int s,
                    double* realAtom,
                    char* fileName)
{
    // FILE* stream1;
    // FILE* stream2;
    int i,k;
    double opt_phase=0;
    int NC;

    // double *in1,*in2;
    fftw_complex *in1, *in2;
    fftw_complex *out1, *out2;
    fftw_plan plan_forward1, plan_forward2;

    NC = ( N/2 ) +1;

    in1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * N );
    in2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * N );
    // in1 = (double*) fftw_malloc ( sizeof ( double ) * N);
    out1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * N );
    // in2 = (double*) fftw_malloc ( sizeof ( double ) * N );
    out2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * N );

    // plan_forward1 = fftw_plan_dft_r2c_1d ( N, in1, out1, FFTW_ESTIMATE );
    // plan_forward2 = fftw_plan_dft_r2c_1d ( N, in2, out2, FFTW_ESTIMATE );
    plan_forward1 = fftw_plan_dft_1d ( N, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE );
    plan_forward2 = fftw_plan_dft_1d ( N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE );
    

    double xi;
    int k_pi;
    double C;
    double innerProd = 0;
    double innerProd_xp = 0;
    double innerProd_xq = 0;
    double innerProd_pp = 0;
    double innerProd_qq = 0;
    double innerProd_pq = 0;
    double a1,b1;


    for (i=0;i<N;i++)
    {
        in1[i][0] = 0.0;
        in1[i][1] = 0.0;
        in2[i][0] = 0.0;
        in2[i][1] = 0.0;
    }

    for (i=0;i<N;i++)
    {
        in1[i][0] = residue[0][i] * realAtom[i];
        in2[i][0] = realAtom[i] * realAtom[i];

    }

    fftw_execute ( plan_forward1 );

    fftw_execute ( plan_forward2 );
    C = out2[0][0]; 
    k_pi=NC-1;
    int delta_k = (int)(ceil((double)N/(double)s));

    for (k=0;k<NC-1;k=k+delta_k)
    {
        // cout << k << "AQUI" << endl;
        innerProd_xp = out1[k][0];  // /sqrt((double)N);

        innerProd_xq = -out1[k][1]; // /sqrt((double)N);

        if ( (2*k)<NC )
        {
            innerProd_pp = 0.5 * (C + out2[2*k][0]);// /sqrt((double)N) );
            innerProd_qq = 0.5 * (C - out2[2*k][0]);// /sqrt((double)N) );
            innerProd_pq = -0.5 * (out2[2*k][1]);// /sqrt((double)N) );
        }
        else
        {
            innerProd_pp = 0.5 * (C + out2[N-(2*k)][0]);// /sqrt((double)N) );
            innerProd_qq = 0.5 * (C - out2[N-(2*k)][0]);// /sqrt((double)N) );
            innerProd_pq = -0.5 * ( - out2[N-(2*k)][1]);// /sqrt((double)N) );
        }

        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

        if ( (k == 0) || (k == k_pi) )
        {
            opt_phase = 0;
            innerProd = innerProd_xp  / sqrt(innerProd_pp);
        }
        else if (a1 == 0)
        {
            opt_phase = (double)(pi/2);
            innerProd = -innerProd_xq / sqrt(innerProd_qq);
        }
        else if ( (a1!=0) && (k != 0) )
        {
            opt_phase = atan( -(b1/a1) );
            innerProd = (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                        sqrt(a1*a1*innerProd_pp +
                             b1*b1*innerProd_qq +
                             2*a1*b1*innerProd_pq);
        }

        xi = (k*2*pi)/N;
        if (fabs(innerProd)>fabs(maxInnerProd))
        {
            maxInnerProd = innerProd;
            chosenXi = xi,
            chosenOptPhase = opt_phase;
        }

        // For debug
        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
        //                innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
        //                innerProd_pq, s, (double)tau, xi,opt_phase);

        // FILE* stream;
        // stream = fopen(fileName,"a");
        // fprintf (   stream,"%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n",
        //             innerProd,
        //             1.0/(double)s,
        //             xi,
        //             opt_phase,
        //             (double)tau,
        //             (double)tau,
        //             (double)(N-1),
        //             innerProd_xp,
        //             innerProd_xq,
        //             innerProd_pp,
        //             innerProd_qq,
        //             innerProd_pq);
        // fflush(stream);
        // fclose(stream);
    }

    fftw_destroy_plan ( plan_forward1 );
    fftw_destroy_plan ( plan_forward2 );


    fftw_free ( in1 );
    fftw_free ( out1 );
    fftw_free ( in2 );
    fftw_free ( out2 );
    return;
}




void fastMPKolasaModified(  cgMatrix<double>& residue,
                            double& maxInnerProd,
                            double& chosenOptPhase,
                            int& chosenTau,
                            int N, //int dicSize
                            int decincAsymmFlag,
                            int delta_tau,
                            double xi,
                            double* realAtomWin,
                            CComplex* complexAtomXi,
                            CComplex* complexAtom2Xi,
                            char* fileName,
                            int s)

{   
    int i;
    int tau;
    double opt_phase = 0;
    double innerProd = 0;
    double innerProd_xp = 0;
    double innerProd_xq = 0;
    double innerProd_pp = 0;
    double innerProd_qq = 0;
    double innerProd_pq = 0;
    double a1, b1;

    fftw_complex    *z1, *z2, *z2_xi0;
    fftw_complex    *w1, *w2/*, *w2_zxi0*/;
    fftw_complex    *z1_fft, *z2_fft, *z2_xi0_fft, *w1_fft, *w2_fft/*, *w2_zxi0_fft*/;
    fftw_complex    *conv_zw1, *conv_zw2, *conv_zxi0_w2/*, *conv_zw2_xi0_expinc*/;
    fftw_complex    *prod_cvec_zw1, *prod_cvec_zw2, *prod_cvec_zw2_xi0;
    
    fftw_plan       plan_forward1a, plan_forward1b, plan_backward1,
                    plan_forward2a, plan_forward2a_xi0, plan_forward2b, plan_backward2, plan_backward2_xi0;
    //FILE* stream;
    //  ---------------------------------------------------
    //stream = fopen( "results_modified_kolasa.out", "w" );
    
    // Memory allocation
    z1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_xi0 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    // w2_zxi0 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_xi0_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    // w2_zxi0_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    // conv_zw2_xi0_expinc = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    prod_cvec_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    prod_cvec_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    prod_cvec_zw2_xi0 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    // PLAN FFT
    plan_forward1a = fftw_plan_dft_1d ( 2*N, z1, z1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward1b = fftw_plan_dft_1d ( 2*N, w1, w1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward1b = fftw_plan_dft_r2c_1d ( 2*N, w1,  w1_fft, FFTW_ESTIMATE);
    plan_backward1 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw1, conv_zw1, FFTW_BACKWARD, FFTW_ESTIMATE); 

    plan_forward2a = fftw_plan_dft_1d ( 2*N, z2, z2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2a_xi0 = fftw_plan_dft_1d ( 2*N, z2_xi0, z2_xi0_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2b = fftw_plan_dft_1d ( 2*N, w2, w2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    // plan_forward2b = fftw_plan_dft_1d ( 2*N, w2_zxi0, w2_zxi0_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward2b = fftw_plan_dft_r2c_1d ( 2*N, w2, w2_fft,FFTW_ESTIMATE );
    plan_backward2 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw2, conv_zw2, FFTW_BACKWARD, FFTW_ESTIMATE); 
    plan_backward2_xi0 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw2_xi0, conv_zxi0_w2, FFTW_BACKWARD, FFTW_ESTIMATE);
    // cout << &(prod_cvec_zw2) << endl;
    for (i=0;i<2*N;i++)
    {
        z1[i][0] = 0.0;
        z1[i][1] = 0.0;
        z2[i][0] = 0.0;
        z2[i][1] = 0.0;
        z2_xi0[i][0] = 1.0;
        z2_xi0[i][1] = 0.0;
        w1[i][0] = 0.0;
        w1[i][1] = 0.0;
        w2[i][0] = 0.0;
        w2[i][1] = 0.0;
        // w2_zxi0[i][0] = 0.0;
        // w2_zxi0[i][1] = 1.0;
    }
    // cout << "--------------------------" << endl;
    double tol = 0;//1e-10;

    for (i=0;i<N;i++)
    {
        w1[i][0] = realAtomWin[N-1-i];
        w2[i][0] = realAtomWin[N-1-i] * realAtomWin[N-1-i];
        // w2_zxi0[i][0] = realAtomWin[N-1-i] * realAtomWin[N-1-i];
        // cout << w2_zxi0[i][0] << " " << w2_zxi0[i][1] << endl;
        // w2_zxi0[i][1] = 1.0;
        // w2[i][0] = realAtomWin[N-1-i] * realAtomWin[N-1-i];
        // if (step == 4)
        // {
        //    cout << i << " - w2: " << w2[i][0] << endl; 
        // }
        // 

        // cout << w1[i][0] << endl;
    }
    // cout << endl << endl;
    for (i=0;i<N;i++)
    {
        z1[i][0] = complexAtomXi[i].Real() * residue[0][i];
        z1[i][1] = complexAtomXi[i].Imag() * residue[0][i];
    }
    // cout << "--------------------------" << endl;
    for (i=0;i<N;i++)
    {
        z2[i][0] = complexAtom2Xi[i].Real();
        z2[i][1] = complexAtom2Xi[i].Imag();
        // if (step == 4)
        // {
        //    cout << i << " - z2: " << z2[i][0] << " & " << z2[i][1] << endl;
        // }
    }

    // zw1
    fftw_execute (plan_forward1a);
    fftw_execute (plan_forward1b);

    // product between complex vectors
    for (i=0;i<2*N;i++)
    {
        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*decincAsymmFlag*(w1_fft[i][1]));
        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*decincAsymmFlag*(w1_fft[i][1]));
        // if (step == 4)
        // {
        //    cout << prod_cvec_zw1[i][0]<< " + j* " << prod_cvec_zw1[i][1] << endl;
        // }
    }

    fftw_execute ( plan_backward1);

    // zw2
    fftw_execute (plan_forward2a);
    fftw_execute (plan_forward2a_xi0);
    fftw_execute (plan_forward2b);

    // product between complex vectors
    // cout << "---------------------------" << endl;
    // for (i = 0; i<2*N; i++)
    // {
    //     // cout << N;
    //     w2_fft[i][0] =      0.5 * (w2_zxi0_fft[i][0] + w2_zxi0_fft[((2*N)-i)%(2*N)][0]);
    //     w2_fft[i][1] =      0.5 * (w2_zxi0_fft[i][1] - w2_zxi0_fft[((2*N)-i)%(2*N)][1]);                
    //     // if (i == 0)
    //     // {
    //         z2_xi0_fft[i][0] =  0.5 * (w2_zxi0_fft[i][1] + w2_zxi0_fft[((2*N)-i)%(2*N)][1]);
    //         z2_xi0_fft[i][1] =  0.5 * (w2_zxi0_fft[((2*N)-i)%(2*N)][0] - w2_zxi0_fft[i][0]);
        // }
        // else
        // {
        //     z2_xi0_fft[i][0] = 0.0;
        //     z2_xi0_fft[i][1] = 0.0;
        // }

        // z2_xi0_fft[i][0] = 0.0;
        // z2_xi0_fft[i][1] = 0.0;

        // z2_xi0_fft[0][0] =  0.5 * (w2_zxi0_fft[i][1] + w2_zxi0_fft[((2*N)-i)%(2*N)][1]);
        // z2_xi0_fft[0][1] =  0.5 * (w2_zxi0_fft[((2*N)-i)%(2*N)][0] - w2_zxi0_fft[i][0]);

        // cout << "w2_zxi0_fft[" << i << "][1]: " << w2_zxi0_fft[i][0] << " - w2_zxi0_fft[" << (2*N)%((2*N)-i) << "][1]: " << w2_zxi0_fft[i][1] << endl;
        // cout << "w2_zxi0_fft[" << i << "][0]: " << w2_zxi0_fft[i][0] << " - w2_zxi0_fft[" << i << "][1]: " << w2_zxi0_fft[i][1] << endl;
        // cout << "W2["<<i<<"][0]:   " << w2_fft[i][0] << " - W2["<<i<<"][1]:   " << w2_fft[i][1] << endl;
        // cout << "Z2_xi0["<<i<<"][0]:   " << z2_xi0_fft[i][0] << " - Z2_xi0["<<i<<"][1]:   " << z2_xi0_fft[i][1] << endl;

    // }
    for (i=0;i<=2*N;i++)
    {
        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*decincAsymmFlag*(w2_fft[i][1]));
        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*decincAsymmFlag*(w2_fft[i][1]));
        // if (step == 4)
        // {
           // cout << i << " - prod_cvec_zw2[i][0]: " << prod_cvec_zw2[i][0] << " - prod_cvec_zw2[i][1]: " << prod_cvec_zw2[i][1] << endl;
        // }
    }

    for (i=0;i<2*N;i++)
    {
        prod_cvec_zw2_xi0[i][0]= (z2_xi0_fft[i][0]*w2_fft[i][0] - z2_xi0_fft[i][1]*decincAsymmFlag*(w2_fft[i][1]));
        prod_cvec_zw2_xi0[i][1]= (z2_xi0_fft[i][1]*w2_fft[i][0] + z2_xi0_fft[i][0]*decincAsymmFlag*(w2_fft[i][1]));
        // if (xi == 0)
        // {
        //     cout << "conv_zxi0_w2["<<i<<"][0]:   " << conv_zxi0_w2[i][0] << " - conv_zxi0_w2["<<i<<"][1]:   " << conv_zxi0_w2[i][1] << endl;
        // }
    }
    fftw_execute ( plan_backward2);
    fftw_execute ( plan_backward2_xi0);

    // if (xi==0.0) //(k==0)
    // {
    //     for (i=0;i<2*N;i++)
    //     {
    //         // cout << "Zxio["<<i<<"][0]: " << z2_xi0_fft[i][1] << " - Zxio["<<i<<"][1]: " << z2_xi0_fft[i][1] << endl;
    //     }

    //     // fftw_execute(plan_backward2_xi0);
    //     // memcpy(conv_zxi0_w2, conv_zw2, sizeof(fftw_complex) * (2*N));
    //     // memcpy(conv_zxi0_w2_expinc, conv_zw2, sizeof(fftw_complex) * (2*N));
    // }

    if (decincAsymmFlag==1) //Decreasing
    {
        // Compute inner product and optimum phase
        for (tau=0; tau<(double)N; tau=tau+delta_tau)
        {
            innerProd_xp = conv_zw1[tau+N-1][0]/(double)(2*N);
            innerProd_xq = conv_zw1[tau+N-1][1]/(double)(2*N);
            innerProd_pp = 0.5*(conv_zxi0_w2[tau+N-1][0] + conv_zw2[tau+N-1][0])/(double)(2*N);
            innerProd_qq = 0.5*(conv_zxi0_w2[tau+N-1][0] - conv_zw2[tau+N-1][0])/(double)(2*N);
            innerProd_pq = 0.5*(conv_zw2[tau+N-1][1]/(double)(2*N));
            // if (step == 4)
            // {
            //    cout << "s: " << s << " - xi: " << xi << " - tau: " << tau << "conv_zxi0_w2: " << conv_zw2_xi0[tau+N-1][0] << " - conv_zw2:" << conv_zw2[tau+N-1][0] << endl; 
            // }

            a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
            b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;
            
            if (( xi == 0.0 )||( (int)xi*10000 == (int)pi*10000)) 
            {
                opt_phase = 0;
                if (fabs(innerProd_pp) > tol)
                {
                    innerProd = innerProd_xp / sqrt(innerProd_pp);
                }
                else
                {
                    innerProd = 0.0;
                    cout << " 1 - " << innerProd_pp << endl;
                }
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
                if (fabs(innerProd_qq) > tol)
                {
                    innerProd = -innerProd_xq / sqrt(innerProd_qq);
                }
                else
                {
                    innerProd = 0.0;
                    cout << " 2 - " << innerProd_qq  << endl;
                }
            }
            else //if ( (a1 != 0) && (xi != 0) )
            {
                opt_phase = atan( -(b1/a1) );
                if ( fabs(a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) > tol)
                {
                    innerProd =  (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) / 
                                 sqrt(a1*a1*innerProd_pp + 
                                 b1*b1*innerProd_qq +
                                 2*a1*b1*innerProd_pq);
                }
                else
                {
                    innerProd = 0.0;
                    cout << " 3 - " << (a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) << endl;
                }

            }
            if (fabs(innerProd)>fabs(maxInnerProd))
            {
                maxInnerProd = innerProd;
                chosenTau = tau;
                chosenOptPhase = opt_phase;
            }
            // FILE* stream;
            // stream = fopen(fileName,"a");
            // fprintf (   stream,"%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n",
            //         innerProd,
            //         1/(double)s,
            //         xi,
            //         opt_phase,
            //         (double)tau,
            //         (double)tau,
            //         (double)(N-1),
            //         innerProd_xp,
            //         innerProd_xq,
            //         innerProd_pp,
            //         innerProd_qq,
            //         innerProd_pq);
            // fflush(stream);
            // fclose(stream);
            // if (step == 4)
            // {
            //    cout << "innerProd: " << innerProd << endl;
            // }
        }   
    }

    // else if (decincAsymmFlag==-1) //Increasing
    // {
    //     // Compute inner product and optimum phase
    //     for (int tau=N-1; tau>0; tau=tau-delta_tau)
    //     {
    //         innerProd_xp = conv_zw1[(tau+N+1)%(2*N)][0]/(double)(2*N);
    //         innerProd_xq = conv_zw1[(tau+N+1)%(2*N)][1]/(double)(2*N);
    //         innerProd_pp = 0.5*(conv_zw2_xi0_expinc[(tau+N+1)%(2*N)][0] + conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
    //         innerProd_qq = 0.5*(conv_zw2_xi0_expinc[(tau+N+1)%(2*N)][0] - conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
    //         innerProd_pq = 0.5*(conv_zw2[(tau+N+1)%(2*N)][1]/(double)(2*N));

    //         a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
    //         b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

    //         if (( xi == 0.0 )||( (int)xi*10000 == (int)pi*10000)) 
    //         {
    //             opt_phase = 0;
    //             if (fabs(innerProd_pp) > tol)
    //             {
    //                 innerProd = innerProd_xp / sqrt(innerProd_pp);
    //             }
    //             else
    //             {
    //                 innerProd = 0.0;
    //                 cout << " 1 - " << innerProd_pp << endl;
    //             }
    //         }
    //         else if (a1 == 0)
    //         {
    //             opt_phase = (double)(pi/2);
    //             if (fabs(innerProd_qq) > tol)
    //             {
    //                 innerProd = -innerProd_xq / sqrt(innerProd_qq);
    //             }
    //             else
    //             {
    //                 innerProd = 0.0;
    //                 cout << " 2 - " << innerProd_qq  << endl;
    //             }
    //         }
    //         else //if ( (a1 != 0) && (xi != 0) )
    //         {
    //             opt_phase = atan( -(b1/a1) );
    //             if ( fabs(a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) > tol)
    //             {
    //                 innerProd =  (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) / 
    //                              sqrt(a1*a1*innerProd_pp + 
    //                              b1*b1*innerProd_qq +
    //                              2*a1*b1*innerProd_pq);
    //             }
    //             else
    //             {
    //                 innerProd = 0.0;
    //                 cout << " 3 - " << (a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) << endl;
    //             }

    //         }
    //         if (fabs(innerProd)>fabs(maxInnerProd))
    //         {
    //             maxInnerProd = innerProd;
    //             chosenTau = tau;
    //             chosenOptPhase = opt_phase;
    //         }

            // FILE* stream;
            // stream = fopen(fileName,"a");
            // fprintf (   stream,"%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n",
            //         innerProd,
            //         1/(double)s,
            //         xi,
            //         opt_phase,
            //         (double)tau,
            //         (double)tau,
            //         (double)(N-1),
            //         innerProd_xp,
            //         innerProd_xq,
            //         innerProd_pp,
            //         innerProd_qq,
            //         innerProd_pq);
            // fflush(stream);
            // fclose(stream);
    //     }
    // }
    // cout << "AQUI 1" << endl;
    fftw_destroy_plan ( plan_forward1a );
    fftw_destroy_plan ( plan_forward1b );
    fftw_destroy_plan ( plan_backward1 );
    fftw_destroy_plan ( plan_forward2a );
    fftw_destroy_plan ( plan_forward2a_xi0);
    fftw_destroy_plan ( plan_forward2b );
    fftw_destroy_plan ( plan_backward2 );
    fftw_destroy_plan ( plan_backward2_xi0);
    
    fftw_free(z1);
    fftw_free(z2);
    fftw_free(z2_xi0);
    fftw_free(w1);
    fftw_free(w2);
    // fftw_free(w2_zxi0);
    fftw_free(z1_fft);
    fftw_free(z2_fft);
    fftw_free(z2_xi0_fft);
    fftw_free(w1_fft);
    fftw_free(w2_fft);
    // fftw_free(w2_zxi0_fft);
    fftw_free(conv_zw1);
    fftw_free(conv_zw2);
    fftw_free(conv_zxi0_w2);
    // fftw_free(conv_zw2_xi0_expinc);
    fftw_free(prod_cvec_zw1);
    // cout << "AQUI 2" << endl;
    fftw_free(prod_cvec_zw2);
    // cout << "AQUI 3" << endl;
    fftw_free(prod_cvec_zw2_xi0);


    return;
}

void setParameters (strtParameter* chosenParm,
                    int dicType,
                    double innerProd,
                    int s,
                    double xi,
                    double opt_phase,
                    int a,
                    int b)
{
    if (dicType == 1)
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
    } 

    return;    
}




// double computeOptimumPhase( cgMatrix<double>& residue,
//                             double xi,
//                             double& innerProd, 
//                             int signalSize,
//                             CComplex* complexAtom)
// {
//     double opt_phase = 0;
//     int i;
//     cgMatrix<double> realPartComplexDic(1,signalSize,0.0);
//     cgMatrix<double> imagPartComplexDic(1,signalSize,0.0);

//     for ( i=0;i<signalSize;i++)
//     {
//         realPartComplexDic[0][i] = complexAtom[i].Real();
//         imagPartComplexDic[0][i] = complexAtom[i].Imag();
//     }

//     double innerProdReal=0;
//     double innerProdImag=0;
//     double innerProdRealImag=0;
//     for (i=0;i< signalSize;i++)
//     {
//         innerProdReal += residue.getData(0,i) * complexAtom[i].Real();
//         innerProdImag += residue.getData(0,i) * complexAtom[i].Imag();
//         innerProdRealImag += complexAtom[i].Real() * complexAtom[i].Imag();
//     }

//     cgMatrix<double> innerProductReal(1,1,innerProdReal);
//     cgMatrix<double> innerProductImag(1,1,innerProdImag);
//     //innerProductReal = residue * (realPartComplexDic.transpose());
//     //innerProductImag = residue * (imagPartComplexDic.transpose());

//     double p,q ;
//     p = realPartComplexDic.norm();
//     q = imagPartComplexDic.norm();

//     cgMatrix<double> innerProductRealImag(1,1,innerProdRealImag);
//     //innerProductRealImag = realPartComplexDic * (imagPartComplexDic.transpose());

//     cgMatrix<double> a1;
//     cgMatrix<double> b1;
//     a1 = innerProductReal*(q*q) - innerProductImag * innerProductRealImag;
//     b1 = innerProductImag*(p*p) - innerProductReal * innerProductRealImag;

//     if ( (xi == 0) ||
//         ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
//     {
//         opt_phase = 0;
//         innerProd = innerProductReal[0][0]/p;
//     }
//     else if (a1.getData(0,0) == 0)
//     {
//         opt_phase = (double)(pi/2);
//         innerProd = -innerProductImag[0][0]/q;
//     }
//     else if (   (a1.getData(0,0)!=0) && (xi!=0) )
//     {
//         opt_phase = atan( -(b1.getData(0,0)/a1.getData(0,0)) );
//         innerProd = (a1[0][0]/fabs(a1[0][0]))*(innerProductReal[0][0]*a1[0][0] + innerProductImag[0][0]*b1[0][0])/
//             sqrt(a1[0][0]*a1[0][0]*p*p+b1[0][0]*b1[0][0]*q*q+2*a1[0][0]*b1[0][0]*innerProductRealImag[0][0]);
//     }
//     return opt_phase;
// }

void computeOptimumPhase(   cgMatrix<double>& residue,
                            double& opt_phase,
                            double& innerProd, 
                            int signalSize,
                            double xi,
                            CComplex* complexAtom,
                            char* fileName,
                            int s,
                            double tau, 
                            int N)
{
    // double opt_phase = 0;
    int i;
    cgMatrix<double> realPartComplexDic(1,signalSize,0.0);
    cgMatrix<double> imagPartComplexDic(1,signalSize,0.0);

    for ( i=0;i<signalSize;i++)
    {
        realPartComplexDic[0][i] = complexAtom[i].Real();
        imagPartComplexDic[0][i] = complexAtom[i].Imag();
    }

    double innerProdReal=0;
    double innerProdImag=0;
    double innerProdRealImag=0;
    for (i=0;i< signalSize;i++)
    {
        innerProdReal += residue.getData(0,i) * complexAtom[i].Real();
        innerProdImag += residue.getData(0,i) * complexAtom[i].Imag();
        innerProdRealImag += complexAtom[i].Real() * complexAtom[i].Imag();
    }

    cgMatrix<double> innerProductReal(1,1,innerProdReal);
    cgMatrix<double> innerProductImag(1,1,innerProdImag);
    //innerProductReal = residue * (realPartComplexDic.transpose());
    //innerProductImag = residue * (imagPartComplexDic.transpose());

    double p,q ;
    p = realPartComplexDic.norm();
    q = imagPartComplexDic.norm();

    cgMatrix<double> innerProductRealImag(1,1,innerProdRealImag);
    //innerProductRealImag = realPartComplexDic * (imagPartComplexDic.transpose());

    cgMatrix<double> a1;
    cgMatrix<double> b1;
    a1 = innerProductReal*(q*q) - innerProductImag * innerProductRealImag;
    b1 = innerProductImag*(p*p) - innerProductReal * innerProductRealImag;

    if ( (xi == 0) ||
        ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
    {
        opt_phase = 0;
        innerProd = innerProductReal[0][0]/p;
    }
    else if (a1.getData(0,0) == 0)
    {
        opt_phase = (double)(pi/2);
        innerProd = -innerProductImag[0][0]/q;
    }
    else if (   (a1.getData(0,0)!=0) && (xi!=0) )
    {
        opt_phase = atan( -(b1.getData(0,0)/a1.getData(0,0)) );
        innerProd = (a1[0][0]/fabs(a1[0][0]))*(innerProductReal[0][0]*a1[0][0] + innerProductImag[0][0]*b1[0][0])/
            sqrt(a1[0][0]*a1[0][0]*p*p+b1[0][0]*b1[0][0]*q*q+2*a1[0][0]*b1[0][0]*innerProductRealImag[0][0]);
    }
    // FILE* stream;
    // stream = fopen(fileName,"a");
    // fprintf (   stream,"%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n",
    //             innerProd,
    //             1.0/(double)s,
    //             xi,
    //             opt_phase,
    //             (double)tau,
    //             (double)tau,
    //             (double)(N-1),
    //             innerProdReal,
    //             innerProdImag,
    //             (p*p),
    //             (q*q),
    //             innerProdRealImag);
    // fflush(stream);
    // fclose(stream);
    return;
}

void adjustParameters ( cgMatrix<double>& residue,
                        strtParameter* parm)
{
    CDictionary* dic;
    
    // if ( parm->dicType==1 )
    // {
        // cout << "AQUI" << endl;

        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*)dic)->adjustParameters(residue,parm);
        // setDictionary(dic);
    // }
    
    if (dic!=NULL) delete dic;   
}

void updateResidue (cgMatrix<double>& residue,
                    double& norm,
                    int dicSize,
                    strtParameter* parm)
{
    CDictionary* dic;
    double* realAtom;
    cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);
    // if ( parm->dicType==1 )
    // {   
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CExpDictionary*) dic)->getRealAtom();
    // }
    cgRealAtom.fillVector(realAtom);
    cgRealAtomAux = cgRealAtom*(((strtParameter*)parm)->innerProduct);
    residue = residue - cgRealAtom*(((strtParameter*)parm)->innerProduct);
    norm = residue.norm();
    if(dic!=NULL) delete dic;
    return;
}