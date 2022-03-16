#include "mpursuit.h"

strtParameter* matchingPursuit( cgMatrix<double>& residue,
                                int dicSize,
                                CDataSignal* dataSignal,
                                CFileDictionary* dicData)
{

    int i,j;
    int k;
    int s;
    // double xi;
    // int tau;
    int delta_tau;
    double Fs;
    double Ffund;
    double delta_f;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    int dicType;
    int signalSize;
    //double* realAtom;
    CComplex* complexAtom;
    //realAtom = dic->getRealAtom();
    //complexAtom = dic->getComplexAtom();

    CDictionary* dic;

    // CParameter* parm;
    //CParameter* chosenParm;
    strtParameter* parm = new strtParameter;
    strtParameter* chosenParm = new strtParameter;
    
    double xi;

    double innerProd = 0;
    double maxInnerProd = 0;
    double opt_phase;


    for(j=0;j<dicData->getNumDicBlock();j++)
    {   
        // Using the dictionary data
        s = (dicData->getDicBlock())[j].scale;
        delta_tau = (dicData->getDicBlock())[j].delta_tau;
        fdiscrtype = (dicData->getDicBlock())[j].fdiscrtype;
        freqi = (dicData->getDicBlock())[j].freqi;
        freqf = (dicData->getDicBlock())[j].freqf;
        dicType = (dicData->getDicBlock())[j].dicType;


        
        int N = signalSize;

         if ((s > dataSignal->getBlockSize())&&(s!=88888))
        {
            cout << "Scale greater than block size!!!" << endl;
            printf("Scale %d - block size %d\n",s,signalSize);
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


        if (freqf>Fs)
        {
            cout << "Final frequency greater than sampling frequency !!!" << endl;
            printf("final freq. %f - Fs %f\n",freqf,Fs);
            exit(1);
        }

        if (fdiscrtype==1) // linear
        {
            delta_f = (2*pi/freqf)* freqi;
            Nfreq = (int)(freqf/(2*freqi));
            //Nfreq = (int)ceil(freqf/freqi);
            xi_vec = new double[Nfreq];
            for (i=0;i<Nfreq;i++)
            {
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

        // if (dicType==1)
        // {
        //     ((strtParameter*)chosenParm)->innerProduct = 0;
        //     ((strtParameter*)chosenParm)->rho = 0;
        //     ((strtParameter*)chosenParm)->xi = 0;
        //     ((strtParameter*)chosenParm)->phase = 0;
        //     ((strtParameter*)chosenParm)->a = 0;
        //     ((strtParameter*)chosenParm)->b = 0;
        // }

        // if (s==1)
        // {
        //     // s=1 Impulse
        //     for(i=0;i<N;i+=delta_tau)
        //     {
        //         innerProd = residue[0][i];
        //         if (fabs(innerProd)>fabs(maxInnerProd))
        //         {
        //             maxInnerProd = innerProd;
        //             ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
        //             ((strtParameter*)chosenParm)->rho = 0.0;
        //             ((strtParameter*)chosenParm)->xi = 0.0;
        //             ((strtParameter*)chosenParm)->phase = 0.0;
        //             ((strtParameter*)chosenParm)->a = i;
        //             ((strtParameter*)chosenParm)->b = i;
        //             // cout << N << " " << delta_tau <<" "<< i << " " << maxInnerProd << endl;
        //         }
        //         //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
        //         //                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
        //         //                  0.0, 1.0, (double)i, 0.0,0.0);
        //     }
        // }
        // else if (s==88888)
        // {
        //     // s = N Aprox. as Pure Sinusoid
        //     cgMatrix<double> innerProduct;
        //     //k_pi = (int) s/2.0;
        //     //for(k=0;k<(int)s;k++)
        //     for(k=0;k<Nfreq;k++)
        //     {
        //         //xi = k * ((2*pi)/s);
        //         xi = xi_vec[k];

        //         if (dicType==1)
        //         {
        //             ((strtParameter*)parm)->rho =     0.0;
        //             ((strtParameter*)parm)->xi =      xi;
        //             ((strtParameter*)parm)->phase =   0.0;
        //             ((strtParameter*)parm)->a =       0;
        //             ((strtParameter*)parm)->b =       N-1;
        //             ((CExpDictionary*)dic)->setComplexAtom(parm);
        //         }


        //         complexAtom = dic->getComplexAtom();

        //         // Computing optimum phase
        //         opt_phase = computeOptimumPhase(    residue,
        //                                             xi,
        //                                             innerProd,
        //                                             signalSize,
        //                                             complexAtom);

        //         if (fabs(innerProd)>fabs(maxInnerProd))
        //         {        
        //             maxInnerProd = innerProd;
        //             {
        //                 if (dicType==1)
        //                 {
        //                     ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
        //                     ((strtParameter*)chosenParm)->rho = 0.0;
        //                     ((strtParameter*)chosenParm)->xi = xi;
        //                     ((strtParameter*)chosenParm)->phase = opt_phase;
        //                     ((strtParameter*)chosenParm)->a = 0;
        //                     ((strtParameter*)chosenParm)->b = N-1;    
        //                 }
        //             }
        //         }
        //     }
        // }

        // else // s =[2,4,...,N/2] Damped Sinusoids
        // {

        if (dicType == 1)
        {
            dic = new CExpDictionary;
            ((CExpDictionary*)dic)->setSignalSize(dicSize);
            signalSize = dic->getSignalSize();
        }  

        if (dicType==1)
        {
            ((strtParameter*)parm)->rho = 1.0/(double)s;
            ((strtParameter*)parm)->xi = 0.0;
            ((strtParameter*)parm)->phase = 0.0;
            ((strtParameter*)parm)->a = 0;
            ((strtParameter*)parm)->b = N-1; 
            ((CExpDictionary*)dic)->setRealAtom(parm);
        }   

        realAtom = dic->getRealAtom();

        // for (i=0;i<N;i++)
        // {
        //     w1[i][0] = realAtom[N-1-i];
        //     w2[i][0] = realAtom[N-1-i] * realAtom[N-1-i];
        // }
        // //k_pi = (int) s/2.0;
        // //for(k=0;k<(int)s;k++)
        for(k=0;k<Nfreq;k++)
        {
            //xi = k * ((2*pi)/s);
            xi = xi_vec[k];

            if ( (xi==0.0) || (xi>=((2*pi)/s) ))
            {
                
                // z1
                if (dicType==1)
                {
                    ((strtParameter*)parm)->rho = 0.0;
                    ((strtParameter*)parm)->xi = xi;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = N-1; 
                    ((CExpDictionary*)dic)->setComplexAtom(parm);
                }

                complexAtom = dic->getComplexAtom();



                cout << "s: " << s <<  
                fastMPKolasaModified(   residue, dic, i, s, delta_tau, 
                                        Nfreq, dicType, N, k, xi, xi_vec, 
                                        innerProd, maxInnerProd, opt_phase, 
                                        complexAtom, parm, chosenParm);  
            }
        if (xi_vec!=NULL) 
        {
            delete [] xi_vec;
        }
        if (dic!= NULL) delete dic;
    }
    
    if ( dicType==1 )
    {
        ((strtParameter*)chosenParm)->dicType = 1;
    }

    delete parm;
    return chosenParm;
}
}


//===============================================================
// Function: fastMPKolasa
// Goal: 
// Return:
//==============================================================


/*CParameter* fastMPKolasa(   cgMatrix<double>& residue,
                            CDataSignal* dataSignal,
                            CFileDictionary* dicData,
                            CDictionary* dic)
{
    CParameter* parm;
    CParameter* chosenParm;
    if (dicData->getDicType()==1)
    {
        parm = new CExpParm;
        chosenParm = new CExpParm;
        ((CExpParm*)chosenParm)->innerProd = 0;
        ((CExpParm*)chosenParm)->rho = 0;
        ((CExpParm*)chosenParm)->xi = 0;
        ((CExpParm*)chosenParm)->phase = 0;
        ((CExpParm*)chosenParm)->a = 0;
        ((CExpParm*)chosenParm)->b = 0;
    }

    int i,k;
    double opt_phase=0;
    int NC;
    int N;

    double *in1,*in2;
    fftw_complex *out1, *out2;
    fftw_plan plan_forward1, plan_forward2;

    FILE* stream;
//  ---------------------------------------------------
    stream = fopen( "results_kolasa.out", "w" );


    N = m_signalSize;
    NC = ( N/2 ) +1;

    in1 = (double*) fftw_malloc ( sizeof ( double ) * N);
    out1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * NC );
    in2 = (double*) fftw_malloc ( sizeof ( double ) * N );
    out2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * NC );

    plan_forward1 = fftw_plan_dft_r2c_1d ( N, in1, out1, FFTW_ESTIMATE );
    plan_forward2 = fftw_plan_dft_r2c_1d ( N, in2, out2, FFTW_ESTIMATE );

    double s, tau, xi;
    int k_pi;
    double C;
    double innerProd_xp;
    double innerProd_xq;
    double innerProd_pp;
    double innerProd_qq;
    double innerProd_pq;
    double a1,b1;
    double innerProd;
    double maxInnerProd = 0;


    // s=1 Impulse
    for(i=0;i<N;i++)
    {
        innerProd = residue[0][i];
        if (fabs(innerProd)>fabs(maxInnerProd))
        {
            maxInnerProd = innerProd;
            if (dicData->getDicType()==1)
            {
                ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                ((CExpParm*)chosenParm)->rho = 1.0;
                ((CExpParm*)chosenParm)->xi = 0.0;
                ((CExpParm*)chosenParm)->phase = 0.0 ;
                ((CExpParm*)chosenParm)->a = i;
                ((CExpParm*)chosenParm)->b = i;
            }
        }
        fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                                  0.0, 1.0, (double)i, 0.0,0.0);
    }
    // Damped sinusoids
    s=2;
    while (s<N)
    {
        int delta_tau = s;
        // Decreasing damped sinusoids
        for (tau=0; tau<(double)N; tau=tau+delta_tau)
        {
            if (dicData->getDicType()==1)
            {
                ((CExpParm*)parm)->rho = (double)(1/s);
                ((CExpParm*)parm)->xi = 0.0;
                ((CExpParm*)parm)->phase = 0.0;
                ((CExpParm*)parm)->a = tau;
                ((CExpParm*)parm)->b = N-1;
                ((CExpDictionary*)dic)->setRealAtom(parm);
            }

            for (i=0;i<N;i++)
            {
                in1[i] = residue[0][i] * m_realAtom[i];
                in2[i] = m_realAtom[i] * m_realAtom[i];
            }

            fftw_execute ( plan_forward1 );

            fftw_execute ( plan_forward2 );

            C = out2[0][0];
            k_pi=NC-1;
            int delta_k = N/s;
            for (k=0;k<NC;k=k+delta_k)
            {
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
                    if (dicData->getDicType()==1)
                    {
                        ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                        ((CExpParm*)chosenParm)->rho = (double)1/s;
                        ((CExpParm*)chosenParm)->xi = xi;
                        ((CExpParm*)chosenParm)->phase = opt_phase ;
                        ((CExpParm*)chosenParm)->a = tau;
                        ((CExpParm*)chosenParm)->b = N-1;    
                    }  
                }

                // For debug
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                //                innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                innerProd_pq, s, (double)tau, xi,opt_phase);
            }
        }

        // Increasing
        for (tau=1; tau<(double)N; tau=tau+delta_tau)
        {
            if (dicData->getDicType()==1)
            {
                ((CExpParm*)parm)->rho = -(double)(1/s);
                ((CExpParm*)parm)->xi = 0.0;
                ((CExpParm*)parm)->phase = 0.0;
                ((CExpParm*)parm)->a = 0.0;
                ((CExpParm*)parm)->b = tau;
                ((CExpDictionary*)dic)->setRealAtom(parm);    
            }
            
            for (i=0;i<N;i++)
            {
                in1[i] = residue[0][i] * m_realAtom[i];
                in2[i] = m_realAtom[i] * m_realAtom[i];
            }

            fftw_execute ( plan_forward1 );

            fftw_execute ( plan_forward2 );

            C = out2[0][0];
            k_pi=NC-1;
            int delta_k = N/s;
            for (k=0;k<NC;k=k+delta_k)
            {
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
                    if (dicDaa->getDicType()==1)
                    {
                        ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                        ((CExpParm*)chosenParm)->rho = -(double)1/s;
                        ((CExpParm*)chosenParm)->xi = xi;
                        ((CExpParm*)chosenParm)->phase = opt_phase ;
                        ((CExpParm*)chosenParm)->a = 0.0;
                        ((CExpParm*)chosenParm)->b = tau;    
                    }
                    
                }

                // For debug
                // fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                //                   innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                   innerProd_pq, s, (double)tau, xi,opt_phase);
            }
        }

        s=s*2;
    }

    // s = N Aprox. as Pure Sinusoid
    tau =0;
    if (dicData->getDicType()==1)
    {
        ((CExpParm*)parm)->rho = 0.0;
        ((CExpParm*)parm)->xi = 0.0;
        ((CExpParm*)parm)->phase = 0.0;
        ((CExpParm*)parm)->a = 0;
        ((CExpParm*)parm)->b = N-1;
        ((CExpDictionary*)dic)->setRealAtom(parm);
    }

    for (i=0;i<N;i++)
    {
        in1[i] = residue[0][i] * m_realAtom[i];
        in2[i] = m_realAtom[i] * m_realAtom[i];
    }

    fftw_execute ( plan_forward1 );

    fftw_execute ( plan_forward2 );

    C = out2[0][0];
    k_pi=NC-1;
    int delta_k = N/s;
    for (k=0;k<NC;k=k+delta_k)
    {
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
            innerProd = innerProd_xp / sqrt(innerProd_pp);
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
            if (dicData->getDicType()==1)
            {
                ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                ((CExpParm*)chosenParm)->rho = -1.0/(double)s;
                ((CExpParm*)chosenParm)->xi = xi;
                ((CExpParm*)chosenParm)->phase = opt_phase ;
                ((CExpParm*)chosenParm)->a = 0;
                ((CExpParm*)chosenParm)->b = tau;
            }
        }
        // For debug
        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
        //                  innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
        //                  innerProd_pq, s, (double)tau, xi,opt_phase);
    }

    fclose(stream);

    fftw_destroy_plan ( plan_forward1 );
    fftw_destroy_plan ( plan_forward2 );

    fftw_free ( in1 );
    fftw_free ( out1 );
    fftw_free ( in2 );
    fftw_free ( out2 );

    delete parm;
    return chosenParm;
}*/


//===============================================================
// Function: fastMPKolasaModified
// Goal: 
// Return:
//==============================================================

/*CParameter* fastMPKolasaModified( cgMatrix<double>& residue,
                                    CDataSignal* dataSignal,
                                    CFileDictionary* dicData,
                                    CDictionary* dic)
{
    CParameter* parm;
    CParameter* chosenParm;
    if (dicData->getDicType()==1)
    {
        parm = new CExpParm;
        chosenParm = new CExpParm;
        ((CExpParm*)chosenParm)->innerProd = 0;
        ((CExpParm*)chosenParm)->rho = 0;
        ((CExpParm*)chosenParm)->xi = 0;
        ((CExpParm*)chosenParm)->phase = 0;
        ((CExpParm*)chosenParm)->a = 0;
        ((CExpParm*)chosenParm)->b = 0;
    }
    //double            *in_vec1, *in_vec2, *out_vec3;
    fftw_complex    *z1, *z2;
    fftw_complex    *w1, *w2;
    fftw_complex    *z1_fft, *z2_fft, *w1_fft, *w2_fft ;
    fftw_complex    *conv_zw1, *conv_zw2, *conv_zxi0_w2, *conv_zxi0_w2_expinc;
    fftw_complex    *prod_cvec_zw1, *prod_cvec_zw2;
    
    fftw_plan       plan_forward1a, plan_forward1b, plan_backward1,
                    plan_forward2a, plan_forward2b, plan_backward2; 
    
    int N = signalSize;

    //FILE* stream;
    //  ---------------------------------------------------
    //stream = fopen( "results_modified_kolasa.out", "w" );
    
    // Memory allocation
    z1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    //w1 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    //w2 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    w1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2_expinc = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    prod_cvec_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    prod_cvec_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    // PLAN FFT
    plan_forward1a = fftw_plan_dft_1d ( 2*N, z1, z1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward1b = fftw_plan_dft_1d ( 2*N, w1, w1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward1b = fftw_plan_dft_r2c_1d ( 2*N, w1,  w1_fft, FFTW_ESTIMATE);
    plan_backward1 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw1, conv_zw1, FFTW_BACKWARD, FFTW_ESTIMATE); 

    plan_forward2a = fftw_plan_dft_1d ( 2*N, z2, z2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2b = fftw_plan_dft_1d ( 2*N, w2, w2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward2b = fftw_plan_dft_r2c_1d ( 2*N, w2, w2_fft,FFTW_ESTIMATE );
    plan_backward2 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw2, conv_zw2, FFTW_BACKWARD, FFTW_ESTIMATE); 

    double innerProd_xp;
    double innerProd_xq;
    double innerProd_pp;
    double innerProd_qq;
    double innerProd_pq;
    double a1, b1;
    double innerProd = 0;
    double maxInnerProd = 0;
    double opt_phase;
    
    int i,j;
    int k;
    int s;
    double xi;
    int delta_tau;
    int tau;
    double Fs;
    double Ffund;
    double delta_f;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    
    //dicData->printToScreen();

    for(j=0;j<dicData->getNumDicBlock();j++)
    {
        // Using the dictionary data
        s = (dicData->getDicBlock())[j].scale;
        delta_tau = (dicData->getDicBlock())[j].delta_tau;
        fdiscrtype = (dicData->getDicBlock())[j].fdiscrtype;
        freqi = (dicData->getDicBlock())[j].freqi;
        freqf = (dicData->getDicBlock())[j].freqf;
        
        //printf(" %d %d %d %f %f \n",s,delta_tau,fdiscrtype,freqi,freqf);
        
        if ((s > dic->getSignalSize())&&(s!=88888))
        {
            cout << "Scale greater than block size!!!" << endl;
            printf("Scale %d - block size %d\n",s,dic->getSignalSize());
            exit(1);
        }

        // Set linear frequencies for COMTRADE files
        if (dataSignal->getType() == 1) // Comtrade
        {
            if ( (freqf==9999999999) && (freqi==9999999999) )
            {
                Ffund = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
            } 
            else if (freqf==9999999999)
            {
                Ffund = freqi;
                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
            }
            delta_f = (2*pi/Fs)* Ffund;
            Nfreq = (int)(Fs/(2*Ffund));
            xi_vec = new double[Nfreq];
            for (i=0;i<Nfreq;i++)
            {
                xi_vec[i] = (double)i * delta_f;
            }
        }
        // Set linear/geometric frequencies space for audio or noise
        if ( (dataSignal->getType() == 2) || // Audio
             (dataSignal->getType() == 3) )   // Noise for Audio
        {
            if (dataSignal->getType() == 2) Fs = ((CAudioSignal*)dataSignal)->getSamplingRate();
            if (dataSignal->getType() == 3) Fs = ((CNoiseSignal*)dataSignal)->getSamplingRate();
            //double A0 = 27.5; // Hz
            //double C8 = 4186.01;
            //double freqi= A0 * pow(2.0,-23/12);
            //double freqf= C8 * pow(2.0,23/12); // B9

            if (fdiscrtype==1) // linear
            {
                Nfreq = (int)ceil(freqf/freqi);
                xi_vec = new double[Nfreq];
                for (i=0;i<Nfreq;i++)
                {
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
        }

        if (s==1)
        {

            // s=1 Impulse
            for(i=0;i<N;i+=delta_tau)
            {
                innerProd = residue[0][i];
                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                    ((CExpParm*)chosenParm)->rho = 0.0;
                    ((CExpParm*)chosenParm)->xi = 0.0;
                    ((CExpParm*)chosenParm)->phase = 0.0 ;
                    ((CExpParm*)chosenParm)->a = i;
                    ((CExpParm*)chosenParm)->b = i;

                }
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                //                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                  0.0, 1.0, (double)i, 0.0,0.0);
            }
        }
        else if (s==88888)
        {
            // s = N Aprox. as Pure Sinusoid
            cgMatrix<double> innerProduct;
            //k_pi = (int) s/2.0;
            //for(k=0;k<(int)s;k++)
            for(k=0;k<Nfreq;k++)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];

                if (dicData->getDicType()==1)
                {
                    ((CExpParm*)parm)->rho      = 0.0;
                    ((CExpParm*)parm)->xi       = xi;
                    ((CExpParm*)parm)->phase    = 0.0;
                    ((CExpParm*)parm)->a        = 0;
                    ((CExpParm*)parm)->b        = N-1;
                    ((CExpDictionary*)dic)->setComplexAtom(parm);
                }                
                                    
                // Computing optimum phase
                opt_phase = computeOptimumPhase(    residue,
                                                    xi,
                                                    innerProd);

                if (fabs(innerProd)>fabs(maxInnerProd))
                {        
                    maxInnerProd = innerProd;
                    if (dicData->getDicType()==1)
                    {
                        ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                        ((CExpParm*)chosenParm)->rho = 0.0;
                        ((CExpParm*)chosenParm)->xi = xi;
                        ((CExpParm*)chosenParm)->phase = opt_phase ;
                        ((CExpParm*)chosenParm)->a = 0;
                        ((CExpParm*)chosenParm)->b = N-1;
                    }
                }

            }
        
        }
        else // s =[2,4,...,N/2] Damped Sinusoids
        {
            for (i=0;i<2*N;i++)
            {
                z1[i][0] = 0.0;
                z1[i][1] = 0.0;
                z2[i][0] = 0.0;
                z2[i][1] = 0.0;
                w1[i][0] = 0.0;
                w1[i][1] = 0.0;
                w2[i][0] = 0.0;
                w2[i][1] = 0.0;
            }
            double tol = 0;//1e-10;
            //  int k_pi;

            if (dicData->getDicType()==1)
            {
                ((CExpParm*)parm)->rho      = 1.0/(double)s;
                ((CExpParm*)parm)->xi       = 0.0;
                ((CExpParm*)parm)->phase    = 0.0;
                ((CExpParm*)parm)->a        = 0;
                ((CExpParm*)parm)->b        = N-1;
                ((CExpDictionary*)dic)->setRealAtom(parm);
            }   

            for (i=0;i<N;i++)
            {
                w1[i][0] = realAtom[N-1-i];
                w2[i][0] = realAtom[N-1-i] * realAtom[N-1-i];
            }
            //k_pi = (int) s/2.0;
            //for(k=0;k<(int)s;k++)
            for(k=0;k<Nfreq;k++)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];

                if ( (xi==0.0) || (xi>=((2*pi)/s) ))
                {
                    
                    // z1
                    if (dicData->getDicType()==1)
                    {
                        ((CExpParm*)parm)->rho      = 0.0;
                        ((CExpParm*)parm)->xi       = xi;
                        ((CExpParm*)parm)->phase    = 0.0;
                        ((CExpParm*)parm)->a        = 0;
                        ((CExpParm*)parm)->b        = N-1;
                        ((CExpDictionary*)dic)->setComplexAtom(parm);
                    }
                    

                    for (i=0;i<N;i++)
                    {
                        z1[i][0] = complexAtom[i].Real() * residue[0][i];
                        z1[i][1] = complexAtom[i].Imag() * residue[0][i];
                    }

                    // z2
                    if (dicData->getDicType()==1)
                    {
                        ((CExpParm*)parm)->rho      = 0.0;
                        ((CExpParm*)parm)->xi       = 2*xi;
                        ((CExpParm*)parm)->phase    = 0.0;
                        ((CExpParm*)parm)->a        = 0;
                        ((CExpParm*)parm)->b        = N-1;
                        ((CExpDictionary*)dic)->setComplexAtom(parm);
                    }

                    for (i=0;i<N;i++)
                    {
                        z2[i][0] = complexAtom[i].Real();
                        z2[i][1] = complexAtom[i].Imag();
                    }

                    // zw1
                    fftw_execute (plan_forward1a);
                    fftw_execute (plan_forward1b);

                    // product between complex vectors
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(w1_fft[i][1]));
                        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(w1_fft[i][1]));
                        //cout << prod_cvec_zw1[i][0]<< " + j* " << prod_cvec_zw1[i][1] << endl;
                    }

                    fftw_execute ( plan_backward1);

                    // zw2
                    fftw_execute (plan_forward2a);
                    fftw_execute (plan_forward2b);

                    // product between complex vectors
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(w2_fft[i][1]));
                        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(w2_fft[i][1]));
                    }

                    fftw_execute ( plan_backward2);

                    if (k==0) // (xi==0)
                    {
                        memcpy(conv_zxi0_w2, conv_zw2, sizeof(fftw_complex) * (2*N));
                    }

                    // Compute inner product and optimum phase
                    //for (int tau=-N+1; tau<N+1; tau+=delta_tau)
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        innerProd_xp = conv_zw1[tau+N-1][0]/(double)(2*N);
                        innerProd_xq = conv_zw1[tau+N-1][1]/(double)(2*N);
                        innerProd_pp = 0.5*(conv_zxi0_w2[tau+N-1][0] + conv_zw2[tau+N-1][0])/(double)(2*N);
                        innerProd_qq = 0.5*(conv_zxi0_w2[tau+N-1][0] - conv_zw2[tau+N-1][0])/(double)(2*N);
                        innerProd_pq = 0.5*(conv_zw2[tau+N-1][1]/(double)(2*N));

                        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                        if (( k == 0 )||( (int)xi*10000 == (int)pi*10000)) 
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
                            //printf("Decreasing\n");
                            //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                            maxInnerProd = innerProd;
                            if (dicData->getDicType()==1)
                            {
                                ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                                ((CExpParm*)chosenParm)->rho = 1.0/(double)s;
                                ((CExpParm*)chosenParm)->xi = xi;
                                ((CExpParm*)chosenParm)->phase = opt_phase ;
                                ((CExpParm*)chosenParm)->a = tau;
                                ((CExpParm*)chosenParm)->b = N-1;
                            }
                        }
                        // For debug
                        //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                        //                  innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                        //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N), 
                                          //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0], 
                                          //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                    }

                    // Increasing exponential

                    // product between complex vectors (z1_fft e  conjugate(w1_fft))  
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(-w1_fft[i][1]));
                        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(-w1_fft[i][1]));
                    }

                    fftw_execute ( plan_backward1);

                    // product between complex vectors (z2_fft e  conjugate(w2_fft))  
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(-w2_fft[i][1]));
                        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(-w2_fft[i][1]));
                    }

                    fftw_execute ( plan_backward2);

                    if (k==0) // (xi==0)
                    {
                        memcpy(conv_zxi0_w2_expinc, conv_zw2, sizeof(fftw_complex) * (2*N));
                    }

                    // Compute inner product and optimum phase
                    for (tau=N-1; tau>0; tau=tau-delta_tau)
                    {
                        innerProd_xp = conv_zw1[(tau+N+1)%(2*N)][0]/(double)(2*N);
                        innerProd_xq = conv_zw1[(tau+N+1)%(2*N)][1]/(double)(2*N);
                        innerProd_pp = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] + conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                        innerProd_qq = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] - conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                        innerProd_pq = 0.5*(conv_zw2[(tau+N+1)%(2*N)][1]/(double)(2*N));

                        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                        if (( k == 0 )||( (int)xi*10000 == (int)pi*10000)) 
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
                            //printf("Increasing\n");
                            //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                            maxInnerProd = innerProd;
                            if (dicData->getDicType()==1)
                            {
                                ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                                ((CExpParm*)chosenParm)->rho = -1.0/(double)s;
                                ((CExpParm*)chosenParm)->xi = xi;
                                ((CExpParm*)chosenParm)->phase = opt_phase ;
                                ((CExpParm*)chosenParm)->a = 0;
                                ((CExpParm*)chosenParm)->b = tau;
                            }
                        }
                        // For debug
                        //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                        //                 innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                        //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N), 
                                          //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0], 
                                          //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                    }
                }
            }
            //s=s*2;
        }
        if (xi_vec!=NULL) 
        {
            delete [] xi_vec;
        }
    }
    
    //fclose(stream);
    

    fftw_destroy_plan ( plan_forward1a );
    fftw_destroy_plan ( plan_forward1b);
    fftw_destroy_plan ( plan_backward1 );
    fftw_destroy_plan ( plan_forward2a );
    fftw_destroy_plan ( plan_forward2b);
    fftw_destroy_plan ( plan_backward2 );
    
    fftw_free(z1);
    fftw_free(z2);
    fftw_free(w1);
    fftw_free(w2);
    fftw_free(z1_fft);
    fftw_free(z2_fft);
    fftw_free(w1_fft);
    fftw_free(w2_fft);
    fftw_free(conv_zw1);
    fftw_free(conv_zw2);
    fftw_free(conv_zxi0_w2);
    fftw_free(conv_zxi0_w2_expinc);
    fftw_free(prod_cvec_zw1);
    fftw_free(prod_cvec_zw2);

    delete parm;
    return chosenParm;
}*/

void fastMPKolasaModified(  cgMatrix<double>& residue,
                            CDictionary* dic,
                            int i,
                            int s,
                            int delta_tau,
                            int Nfreq,
                            int dicType,
                            int N,
                            int k,
                            double xi,
                            double* xi_vec,
                            double innerProd,
                            double maxInnerProd,
                            double opt_phase,
                            CComplex* complexAtom,
                            strtParameter* parm,
                            strtParameter* chosenParm)
{
    // int k;
    int tau;
    // double xi;
    // strtParameter* parm = new strtParameter;

    double* realAtom;
    //CComplex* complexAtom;

    /*if (dicType==1)
    {
        // strtParameter* chosenParm;
        ((strtParameter*)chosenParm)->innerProduct = 0;
        ((strtParameter*)chosenParm)->rho = 0;
        ((strtParameter*)chosenParm)->xi = 0;
        ((strtParameter*)chosenParm)->phase = 0;
        ((strtParameter*)chosenParm)->a = 0;
        ((strtParameter*)chosenParm)->b = 0;
    }*/
    //double            *in_vec1, *in_vec2, *out_vec3;
    fftw_complex    *z1, *z2;
    fftw_complex    *w1, *w2;
    fftw_complex    *z1_fft, *z2_fft, *w1_fft, *w2_fft ;
    fftw_complex    *conv_zw1, *conv_zw2, *conv_zxi0_w2, *conv_zxi0_w2_expinc;
    fftw_complex    *prod_cvec_zw1, *prod_cvec_zw2;
    
    fftw_plan       plan_forward1a, plan_forward1b, plan_backward1,
                    plan_forward2a, plan_forward2b, plan_backward2; 

    //FILE* stream;
    //  ---------------------------------------------------
    //stream = fopen( "results_modified_kolasa.out", "w" );
    
    // Memory allocation
    z1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    //w1 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    //w2 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    w1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2_expinc = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    prod_cvec_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    prod_cvec_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    // PLAN FFT
    plan_forward1a = fftw_plan_dft_1d ( 2*N, z1, z1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward1b = fftw_plan_dft_1d ( 2*N, w1, w1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward1b = fftw_plan_dft_r2c_1d ( 2*N, w1,  w1_fft, FFTW_ESTIMATE);
    plan_backward1 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw1, conv_zw1, FFTW_BACKWARD, FFTW_ESTIMATE); 

    plan_forward2a = fftw_plan_dft_1d ( 2*N, z2, z2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2b = fftw_plan_dft_1d ( 2*N, w2, w2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward2b = fftw_plan_dft_r2c_1d ( 2*N, w2, w2_fft,FFTW_ESTIMATE );
    plan_backward2 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw2, conv_zw2, FFTW_BACKWARD, FFTW_ESTIMATE); 

    double innerProd_xp;
    double innerProd_xq;
    double innerProd_pp;
    double innerProd_qq;
    double innerProd_pq;
    double a1, b1;

    /*double innerProd = 0;
    double maxInnerProd = 0;
    double opt_phase;

    if (s==1)
        {
            // s=1 Impulse
            for(i=0;i<N;i+=delta_tau)
            {
                innerProd = residue[0][i];
                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                    ((CExpParm*)chosenParm)->rho = 0.0;
                    ((CExpParm*)chosenParm)->xi = 0.0;
                    ((CExpParm*)chosenParm)->phase = 0.0 ;
                    ((CExpParm*)chosenParm)->a = i;
                    ((CExpParm*)chosenParm)->b = i;
                    ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                    ((strtParameter*)chosenParm)->rho = 0.0;
                    ((strtParameter*)chosenParm)->xi = 0.0;
                    ((strtParameter*)chosenParm)->phase = 0.0;
                    ((strtParameter*)chosenParm)->a = i;
                    ((strtParameter*)chosenParm)->b = i;
                }
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                //                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                  0.0, 1.0, (double)i, 0.0,0.0);
            }
        }
    else if (s==88888)
    {
        // s = N Aprox. as Pure Sinusoid
        cgMatrix<double> innerProduct;
        //k_pi = (int) s/2.0;
        //for(k=0;k<(int)s;k++)
        for(k=0;k<Nfreq;k++)
        {
            //xi = k * ((2*pi)/s);
            xi = xi_vec[k];

            if (dicType==1)
            {
                ((CExpParm*)parm)->rho      = 0.0;
                ((CExpParm*)parm)->xi       = xi;
                ((CExpParm*)parm)->phase    = 0.0;
                ((CExpParm*)parm)->a        = 0;
                ((CExpParm*)parm)->b        = N-1;
                ((strtParameter*)parm)->rho =     0.0;
                ((strtParameter*)parm)->xi =      xi;
                ((strtParameter*)parm)->phase =   0.0;
                ((strtParameter*)parm)->a =       0;
                ((strtParameter*)parm)->b =       N-1;
                ((CExpDictionary*)dic)->setComplexAtom(parm);
            }
            // Computing optimum phase
            opt_phase = computeOptimumPhase(    residue,
                                                xi,
                                                innerProd,
                                                signalSize,
                                                complexAtom);

            if (fabs(innerProd)>fabs(maxInnerProd))
            {        
                maxInnerProd = innerProd;
                if (dicType==1)
                {
                    ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                    ((CExpParm*)chosenParm)->rho = 0.0;
                    ((CExpParm*)chosenParm)->xi = xi;
                    ((CExpParm*)chosenParm)->phase = opt_phase ;
                    ((CExpParm*)chosenParm)->a = 0;
                    ((CExpParm*)chosenParm)->b = N-1;
                    ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                    ((strtParameter*)chosenParm)->rho = 0.0;
                    ((strtParameter*)chosenParm)->xi = xi;
                    ((strtParameter*)chosenParm)->phase = opt_phase;
                    ((strtParameter*)chosenParm)->a = 0;
                    ((strtParameter*)chosenParm)->b = N-1; 
                }
            }

        }
    
    }
    else // s =[2,4,...,N/2] Damped Sinusoids
    {
    */
    for (i=0;i<2*N;i++)
    {
        z1[i][0] = 0.0;
        z1[i][1] = 0.0;
        z2[i][0] = 0.0;
        z2[i][1] = 0.0;
        w1[i][0] = 0.0;
        w1[i][1] = 0.0;
        w2[i][0] = 0.0;
        w2[i][1] = 0.0;
    }
    double tol = 0;//1e-10;
    //  int k_pi;

    // if (dicType==1)
    // {
    //     ((strtParameter*)parm)->rho = 1.0/(double)s;
    //     ((strtParameter*)parm)->xi = 0.0;
    //     ((strtParameter*)parm)->phase = 0.0;
    //     ((strtParameter*)parm)->a = 0;
    //     ((strtParameter*)parm)->b = N-1; 
    //     ((CExpDictionary*)dic)->setRealAtom(parm);
    // }   

    // realAtom = dic->getRealAtom();

    for (i=0;i<N;i++)
    {
        w1[i][0] = realAtomWin[N-1-i];
        w2[i][0] = realAtomWin[N-1-i] * realAtomWin[N-1-i];
    }
    //k_pi = (int) s/2.0;
    //for(k=0;k<(int)s;k++)
    for(k=0;k<Nfreq;k++)
    {
        //xi = k * ((2*pi)/s);
        xi = xi_vec[k];

        if ( (xi==0.0) || (xi>=((2*pi)/s) ))
        {
            
            // z1
            if (dicType==1)
            {
                ((strtParameter*)parm)->rho = 0.0;
                ((strtParameter*)parm)->xi = xi;
                ((strtParameter*)parm)->phase = 0.0;
                ((strtParameter*)parm)->a = 0;
                ((strtParameter*)parm)->b = N-1; 
                ((CExpDictionary*)dic)->setComplexAtom(parm);
            }

            complexAtom = dic->getComplexAtom();
            

            for (i=0;i<N;i++)
            {
                z1[i][0] = complexAtomxi[i].Real() * residue[0][i];
                z1[i][1] = complexAtomxi[i].Imag() * residue[0][i];
            }

            // z2
            if (dicType==1)
            {
                ((strtParameter*)parm)->rho = 0.0;
                ((strtParameter*)parm)->xi = 2*xi;
                ((strtParameter*)parm)->phase = 0.0;
                ((strtParameter*)parm)->a = 0;
                ((strtParameter*)parm)->b = N-1; 
                ((CExpDictionary*)dic)->setComplexAtom(parm);
            }

            complexAtom2xi = dic->getComplexAtom();

            for (i=0;i<N;i++)
            {
                z2[i][0] = complexAtom2xi[i].Real();
                z2[i][1] = complexAtom2xi[i].Imag();
            }

            // zw1
            fftw_execute (plan_forward1a);
            fftw_execute (plan_forward1b);

            // product between complex vectors
            for (i=0;i<2*N;i++)
            {
                prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(w1_fft[i][1]));
                prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(w1_fft[i][1]));
                //cout << prod_cvec_zw1[i][0]<< " + j* " << prod_cvec_zw1[i][1] << endl;
            }

            fftw_execute ( plan_backward1);

            // zw2
            fftw_execute (plan_forward2a);
            fftw_execute (plan_forward2b);

            // product between complex vectors
            for (i=0;i<2*N;i++)
            {
                prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(w2_fft[i][1]));
                prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(w2_fft[i][1]));
            }

            fftw_execute ( plan_backward2);

            if (k==0) // (xi==0)
            {
                memcpy(conv_zxi0_w2, conv_zw2, sizeof(fftw_complex) * (2*N));
            }

            // Compute inner product and optimum phase
            //for (int tau=-N+1; tau<N+1; tau+=delta_tau)
            for (tau=0; tau<(double)N; tau=tau+delta_tau)
            {
                innerProd_xp = conv_zw1[tau+N-1][0]/(double)(2*N);
                innerProd_xq = conv_zw1[tau+N-1][1]/(double)(2*N);
                innerProd_pp = 0.5*(conv_zxi0_w2[tau+N-1][0] + conv_zw2[tau+N-1][0])/(double)(2*N);
                innerProd_qq = 0.5*(conv_zxi0_w2[tau+N-1][0] - conv_zw2[tau+N-1][0])/(double)(2*N);
                innerProd_pq = 0.5*(conv_zw2[tau+N-1][1]/(double)(2*N));

                a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                if (( k == 0 )||( (int)xi*10000 == (int)pi*10000)) 
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
                    //printf("Decreasing\n");
                    //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                    maxInnerProd = innerProd;
                    if (dicType==1)
                    {
                        ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                        ((strtParameter*)chosenParm)->xi = xi;
                        ((strtParameter*)chosenParm)->phase = opt_phase;
                        ((strtParameter*)chosenParm)->a = tau;
                        ((strtParameter*)chosenParm)->b = N-1; 
                    }
                }
                // For debug
                //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                //                  innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N), 
                                  //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0], 
                                  //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
            }

            // Increasing exponential

            // product between complex vectors (z1_fft e  conjugate(w1_fft))  
            for (i=0;i<2*N;i++)
            {
                prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(-w1_fft[i][1]));
                prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(-w1_fft[i][1]));
            }

            fftw_execute ( plan_backward1);

            // product between complex vectors (z2_fft e  conjugate(w2_fft))  
            for (i=0;i<2*N;i++)
            {
                prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(-w2_fft[i][1]));
                prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(-w2_fft[i][1]));
            }

            fftw_execute ( plan_backward2);

            if (k==0) // (xi==0)
            {
                memcpy(conv_zxi0_w2_expinc, conv_zw2, sizeof(fftw_complex) * (2*N));
            }

            // Compute inner product and optimum phase
            for (tau=N-1; tau>0; tau=tau-delta_tau)
            {
                innerProd_xp = conv_zw1[(tau+N+1)%(2*N)][0]/(double)(2*N);
                innerProd_xq = conv_zw1[(tau+N+1)%(2*N)][1]/(double)(2*N);
                innerProd_pp = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] + conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                innerProd_qq = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] - conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                innerProd_pq = 0.5*(conv_zw2[(tau+N+1)%(2*N)][1]/(double)(2*N));

                a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                if (( k == 0 )||( (int)xi*10000 == (int)pi*10000)) 
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
                    //printf("Increasing\n");
                    //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                    maxInnerProd = innerProd;
                    if (dicType==1)
                    {
                        ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                        ((strtParameter*)chosenParm)->xi = xi;
                        ((strtParameter*)chosenParm)->phase = opt_phase;
                        ((strtParameter*)chosenParm)->a = 0;
                        ((strtParameter*)chosenParm)->b = tau; 
                    }
                }
                    // For debug
                    //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                    //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                    //                 innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                    //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N), 
                                      //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0], 
                                      //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                //}
            }
        }
        //s=s*2;
    }

/*    if (xi_vec!=NULL) 
    {
        delete [] xi_vec;
    }
    */
    //fclose(stream);
    

    fftw_destroy_plan ( plan_forward1a );
    fftw_destroy_plan ( plan_forward1b);
    fftw_destroy_plan ( plan_backward1 );
    fftw_destroy_plan ( plan_forward2a );
    fftw_destroy_plan ( plan_forward2b);
    fftw_destroy_plan ( plan_backward2 );
    
    fftw_free(z1);
    fftw_free(z2);
    fftw_free(w1);
    fftw_free(w2);
    fftw_free(z1_fft);
    fftw_free(z2_fft);
    fftw_free(w1_fft);
    fftw_free(w2_fft);
    fftw_free(conv_zw1);
    fftw_free(conv_zw2);
    fftw_free(conv_zxi0_w2);
    fftw_free(conv_zxi0_w2_expinc);
    fftw_free(prod_cvec_zw1);
    fftw_free(prod_cvec_zw2);

    // delete parm;
}


double computeOptimumPhase( cgMatrix<double>& residue,
                            double xi,
                            double& innerProd, 
                            int signalSize,
                            CComplex* complexAtom)
{
    double opt_phase = 0;
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

//  FILE* stream;
//  ---------------------------------------------------
/*  stream = fopen( "results_modified_kolasa.out", "a" );
    fprintf(stream,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                                    innerProd,
                                                    innerProductReal[0][0],
                                                    innerProductImag[0][0],
                                                    p*p,
                                                    q*q,
                                                    innerProductRealImag[0][0],
                                                    (double)m_signalSize,//pow(2.0, (double)(expDic.getGammaExp())[atomIndex].j),
                                                    0.0,//((expDic.getGammaExp())[atomIndex].p)*pow(2.0, (double)(expDic.getGammaExp())[atomIndex].j-1),
                                                    xi,//(expDic.getGammaExp())[atomIndex].k*2*pi*pow(2.0, -(double)(expDic.getGammaExp())[atomIndex].j));
                                                    opt_phase);
    fclose(stream);
*/

    return opt_phase;
}


void adjustParameters ( cgMatrix<double>& residue,
                        strtParameter* parm)
{
    CDictionary* dic;
    
    if ( parm->dicType==1 )
    {
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*)dic)->adjustParameters(residue,parm);
        // setDictionary(dic);
    }
    
    if (dic!=NULL) delete dic;   
}

double* updateResidue ( cgMatrix<double>& residue,
                        strtParameter* parm)
{
    CDictionary* dic;
    double* realAtom;
    if ( parm->dicType==1 )
    {   
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CExpDictionary*) dic)->getRealAtom();
    }
    return realAtom;
}



////////////////////////////////////////////////////////////////////////

/*CParameter* fastMPKolasaModified(   cgMatrix<double>& residue,
                                    CDataSignal* dataSignal,
                                    CFileDictionary* dicData,
                                    CDictionary* dic)
{
    CParameter* parm;
    CParameter* chosenParm;
    if ((dicData->getDicBlock())[j].dicType==1)
    {
        parm = new CExpParm;
        chosenParm = new CExpParm;
        ((CExpParm*)chosenParm)->innerProd = 0;
        ((CExpParm*)chosenParm)->rho = 0;
        ((CExpParm*)chosenParm)->xi = 0;
        ((CExpParm*)chosenParm)->phase = 0;
        ((CExpParm*)chosenParm)->a = 0;
        ((CExpParm*)chosenParm)->b = 0;
    }
    //double            *in_vec1, *in_vec2, *out_vec3;
    fftw_complex    *z1, *z2;
    fftw_complex    *w1, *w2;
    fftw_complex    *z1_fft, *z2_fft, *w1_fft, *w2_fft ;
    fftw_complex    *conv_zw1, *conv_zw2, *conv_zxi0_w2, *conv_zxi0_w2_expinc;
    fftw_complex    *prod_cvec_zw1, *prod_cvec_zw2;
    
    fftw_plan       plan_forward1a, plan_forward1b, plan_backward1,
                    plan_forward2a, plan_forward2b, plan_backward2; 
    
    int N = signalSize;

    //FILE* stream;
    //  ---------------------------------------------------
    //stream = fopen( "results_modified_kolasa.out", "w" );
    
    // Memory allocation
    z1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    //w1 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    //w2 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    w1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2_expinc = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    prod_cvec_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    prod_cvec_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    // PLAN FFT
    plan_forward1a = fftw_plan_dft_1d ( 2*N, z1, z1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward1b = fftw_plan_dft_1d ( 2*N, w1, w1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward1b = fftw_plan_dft_r2c_1d ( 2*N, w1,  w1_fft, FFTW_ESTIMATE);
    plan_backward1 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw1, conv_zw1, FFTW_BACKWARD, FFTW_ESTIMATE); 

    plan_forward2a = fftw_plan_dft_1d ( 2*N, z2, z2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2b = fftw_plan_dft_1d ( 2*N, w2, w2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward2b = fftw_plan_dft_r2c_1d ( 2*N, w2, w2_fft,FFTW_ESTIMATE );
    plan_backward2 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw2, conv_zw2, FFTW_BACKWARD, FFTW_ESTIMATE); 

    double innerProd_xp;
    double innerProd_xq;
    double innerProd_pp;
    double innerProd_qq;
    double innerProd_pq;
    double a1, b1;
    double innerProd = 0;
    double maxInnerProd = 0;
    double opt_phase;
    
    int i,j;
    int k;
    int s;
    double xi;
    int delta_tau;
    int tau;
    double Fs;
    double Ffund;
    double delta_f;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    
    //dicData->printToScreen();


    for(j=0;j<dicData->getNumDicBlock();j++)
    {
        // Using the dictionary data
        s = (dicData->getDicBlock())[j].scale;
        delta_tau = (dicData->getDicBlock())[j].delta_tau;
        fdiscrtype = (dicData->getDicBlock())[j].fdiscrtype;
        freqi = (dicData->getDicBlock())[j].freqi;
        freqf = (dicData->getDicBlock())[j].freqf;
        dicType = (dicData->getDicBlock())[j].dicType;

         if ((s > dic->getSignalSize())&&(s!=88888))
        {
            cout << "Scale greater than block size!!!" << endl;
            printf("Scale %d - block size %d\n",s,signalSize);
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


        if (freqf>Fs)
        {
            cout << "Final frequency greater than sampling frequency !!!" << endl;
            printf("final freq. %f - Fs %f\n",freqf,Fs);
            exit(1);
        }

        if (fdiscrtype==1) // linear
        {
            delta_f = (2*pi/freqf)* freqi;
            Nfreq = (int)(freqf/(2*freqi));
            //Nfreq = (int)ceil(freqf/freqi);
            xi_vec = new double[Nfreq];
            for (i=0;i<Nfreq;i++)
            {
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

        if (s==1)
        {

            // s=1 Impulse
            for(i=0;i<N;i+=delta_tau)
            {
                innerProd = residue[0][i];
                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                    ((CExpParm*)chosenParm)->rho = 0.0;
                    ((CExpParm*)chosenParm)->xi = 0.0;
                    ((CExpParm*)chosenParm)->phase = 0.0 ;
                    ((CExpParm*)chosenParm)->a = i;
                    ((CExpParm*)chosenParm)->b = i;

                }
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                //                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                  0.0, 1.0, (double)i, 0.0,0.0);
            }
        }
        else if (s==88888)
        {
            // s = N Aprox. as Pure Sinusoid
            cgMatrix<double> innerProduct;
            //k_pi = (int) s/2.0;
            //for(k=0;k<(int)s;k++)
            for(k=0;k<Nfreq;k++)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];

                if ((dicData->getDicBlock())[j].dicType==1)
                {
                    ((CExpParm*)parm)->rho      = 0.0;
                    ((CExpParm*)parm)->xi       = xi;
                    ((CExpParm*)parm)->phase    = 0.0;
                    ((CExpParm*)parm)->a        = 0;
                    ((CExpParm*)parm)->b        = N-1;
                    ((CExpDictionary*)dic)->setComplexAtom(parm);
                }                
                                    
                // Computing optimum phase
                opt_phase = computeOptimumPhase(    residue,
                                                    xi,
                                                    innerProd);

                if (fabs(innerProd)>fabs(maxInnerProd))
                {        
                    maxInnerProd = innerProd;
                    if ((dicData->getDicBlock())[j].dicType==1)
                    {
                        ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                        ((CExpParm*)chosenParm)->rho = 0.0;
                        ((CExpParm*)chosenParm)->xi = xi;
                        ((CExpParm*)chosenParm)->phase = opt_phase ;
                        ((CExpParm*)chosenParm)->a = 0;
                        ((CExpParm*)chosenParm)->b = N-1;
                    }
                }

            }
        
        }
        else // s =[2,4,...,N/2] Damped Sinusoids
        {
            for (i=0;i<2*N;i++)
            {
                z1[i][0] = 0.0;
                z1[i][1] = 0.0;
                z2[i][0] = 0.0;
                z2[i][1] = 0.0;
                w1[i][0] = 0.0;
                w1[i][1] = 0.0;
                w2[i][0] = 0.0;
                w2[i][1] = 0.0;
            }
            double tol = 0;//1e-10;
            //  int k_pi;

            if ((dicData->getDicBlock())[j].dicType==1)
            {
                ((CExpParm*)parm)->rho      = 1.0/(double)s;
                ((CExpParm*)parm)->xi       = 0.0;
                ((CExpParm*)parm)->phase    = 0.0;
                ((CExpParm*)parm)->a        = 0;
                ((CExpParm*)parm)->b        = N-1;
                ((CExpDictionary*)dic)->setRealAtom(parm);
            }   

            for (i=0;i<N;i++)
            {
                w1[i][0] = realAtom[N-1-i];
                w2[i][0] = realAtom[N-1-i] * realAtom[N-1-i];
            }
            //k_pi = (int) s/2.0;
            //for(k=0;k<(int)s;k++)
            for(k=0;k<Nfreq;k++)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];

                if ( (xi==0.0) || (xi>=((2*pi)/s) ))
                {
                    
                    // z1
                    if ((dicData->getDicBlock())[j].dicType==1)
                    {
                        ((CExpParm*)parm)->rho      = 0.0;
                        ((CExpParm*)parm)->xi       = xi;
                        ((CExpParm*)parm)->phase    = 0.0;
                        ((CExpParm*)parm)->a        = 0;
                        ((CExpParm*)parm)->b        = N-1;
                        ((CExpDictionary*)dic)->setComplexAtom(parm);
                    }
                    

                    for (i=0;i<N;i++)
                    {
                        z1[i][0] = complexAtom[i].Real() * residue[0][i];
                        z1[i][1] = complexAtom[i].Imag() * residue[0][i];
                    }

                    // z2
                    if ((dicData->getDicBlock())[j].dicType==1)
                    {
                        ((CExpParm*)parm)->rho      = 0.0;
                        ((CExpParm*)parm)->xi       = 2*xi;
                        ((CExpParm*)parm)->phase    = 0.0;
                        ((CExpParm*)parm)->a        = 0;
                        ((CExpParm*)parm)->b        = N-1;
                        ((CExpDictionary*)dic)->setComplexAtom(parm);
                    }

                    for (i=0;i<N;i++)
                    {
                        z2[i][0] = complexAtom[i].Real();
                        z2[i][1] = complexAtom[i].Imag();
                    }

                    // zw1
                    fftw_execute (plan_forward1a);
                    fftw_execute (plan_forward1b);

                    // product between complex vectors
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(w1_fft[i][1]));
                        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(w1_fft[i][1]));
                        //cout << prod_cvec_zw1[i][0]<< " + j* " << prod_cvec_zw1[i][1] << endl;
                    }

                    fftw_execute ( plan_backward1);

                    // zw2
                    fftw_execute (plan_forward2a);
                    fftw_execute (plan_forward2b);

                    // product between complex vectors
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(w2_fft[i][1]));
                        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(w2_fft[i][1]));
                    }

                    fftw_execute ( plan_backward2);

                    if (k==0) // (xi==0)
                    {
                        memcpy(conv_zxi0_w2, conv_zw2, sizeof(fftw_complex) * (2*N));
                    }

                    // Compute inner product and optimum phase
                    //for (int tau=-N+1; tau<N+1; tau+=delta_tau)
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        innerProd_xp = conv_zw1[tau+N-1][0]/(double)(2*N);
                        innerProd_xq = conv_zw1[tau+N-1][1]/(double)(2*N);
                        innerProd_pp = 0.5*(conv_zxi0_w2[tau+N-1][0] + conv_zw2[tau+N-1][0])/(double)(2*N);
                        innerProd_qq = 0.5*(conv_zxi0_w2[tau+N-1][0] - conv_zw2[tau+N-1][0])/(double)(2*N);
                        innerProd_pq = 0.5*(conv_zw2[tau+N-1][1]/(double)(2*N));

                        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                        if (( k == 0 )||( (int)xi*10000 == (int)pi*10000)) 
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
                            //printf("Decreasing\n");
                            //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                            maxInnerProd = innerProd;
                            if ((dicData->getDicBlock())[j].dicType==1)
                            {
                                ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                                ((CExpParm*)chosenParm)->rho = 1.0/(double)s;
                                ((CExpParm*)chosenParm)->xi = xi;
                                ((CExpParm*)chosenParm)->phase = opt_phase ;
                                ((CExpParm*)chosenParm)->a = tau;
                                ((CExpParm*)chosenParm)->b = N-1;
                            }
                        }
                        // For debug
                        //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                        //                  innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                        //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N), 
                                          //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0], 
                                          //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                    }

                    // Increasing exponential

                    // product between complex vectors (z1_fft e  conjugate(w1_fft))  
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(-w1_fft[i][1]));
                        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(-w1_fft[i][1]));
                    }

                    fftw_execute ( plan_backward1);

                    // product between complex vectors (z2_fft e  conjugate(w2_fft))  
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(-w2_fft[i][1]));
                        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(-w2_fft[i][1]));
                    }

                    fftw_execute ( plan_backward2);

                    if (k==0) // (xi==0)
                    {
                        memcpy(conv_zxi0_w2_expinc, conv_zw2, sizeof(fftw_complex) * (2*N));
                    }

                    // Compute inner product and optimum phase
                    for (tau=N-1; tau>0; tau=tau-delta_tau)
                    {
                        innerProd_xp = conv_zw1[(tau+N+1)%(2*N)][0]/(double)(2*N);
                        innerProd_xq = conv_zw1[(tau+N+1)%(2*N)][1]/(double)(2*N);
                        innerProd_pp = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] + conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                        innerProd_qq = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] - conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                        innerProd_pq = 0.5*(conv_zw2[(tau+N+1)%(2*N)][1]/(double)(2*N));

                        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                        if (( k == 0 )||( (int)xi*10000 == (int)pi*10000)) 
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
                            //printf("Increasing\n");
                            //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                            maxInnerProd = innerProd;
                            if ((dicData->getDicBlock())[j].dicType==1)
                            {
                                ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                                ((CExpParm*)chosenParm)->rho = -1.0/(double)s;
                                ((CExpParm*)chosenParm)->xi = xi;
                                ((CExpParm*)chosenParm)->phase = opt_phase ;
                                ((CExpParm*)chosenParm)->a = 0;
                                ((CExpParm*)chosenParm)->b = tau;
                            }
                        }
                        // For debug
                        //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                        //                 innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                        //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N), 
                                          //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0], 
                                          //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                    }
                }
            }
            //s=s*2;
        }
        if (xi_vec!=NULL) 
        {
            delete [] xi_vec;
        }
    }
    
    //fclose(stream);
    

    fftw_destroy_plan ( plan_forward1a );
    fftw_destroy_plan ( plan_forward1b);
    fftw_destroy_plan ( plan_backward1 );
    fftw_destroy_plan ( plan_forward2a );
    fftw_destroy_plan ( plan_forward2b);
    fftw_destroy_plan ( plan_backward2 );
    
    fftw_free(z1);
    fftw_free(z2);
    fftw_free(w1);
    fftw_free(w2);
    fftw_free(z1_fft);
    fftw_free(z2_fft);
    fftw_free(w1_fft);
    fftw_free(w2_fft);
    fftw_free(conv_zw1);
    fftw_free(conv_zw2);
    fftw_free(conv_zxi0_w2);
    fftw_free(conv_zxi0_w2_expinc);
    fftw_free(prod_cvec_zw1);
    fftw_free(prod_cvec_zw2);

    delete parm;
    return chosenParm;
}*/