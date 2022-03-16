#include "mpursuit.h"



void matchingPursuit(   cgMatrix<double>& residue,
                        strtParameter* chosenParm,
                        int dicSize,
                        CDataSignal* dataSignal,
                        CFileDictionary* dicData,
                        CFileDecomp* genData, int step, int chosenDic,
                        strtParameter** dicAtoms,
                        int& a0,
                        int& b0,
                        int flagOMP)
{
    int decincAsymmFlag;
    int i,j;
    int k;
    // int s;
    double s;
    int a;
    int b;
    int e = 8;

    // int a0 = -1;
    // int b0 = dicSize;

    int chosenAtom;

    int chosenTau = 0;
    double chosenXi = 0.0;
    int delta_tau;
    double chosenOptPhase = 0.0;
    double Fs;
    double beta = 0.0;
    double rho = 0.0;
    double decay=0.0;
    double rise= 0.0;
    // double Ffund;
    // double delta_f = 0;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    int dicType;
    // int signalSize;
    double innerProd;
    double  maxInnerProd2 = 0.0;

    double* realAtom = new double[dicSize];
    double* realAtomDec = new double [dicSize];
    double* realAtomInc = new double [dicSize];
    double* realAtomWin = new double[2*dicSize];

    CComplex* complexAtom = new CComplex[dicSize];
    CComplex* complexAtomDec = new CComplex[dicSize];
    CComplex* complexAtomInc = new CComplex[dicSize];
    CComplex* complexAtomXi = new CComplex[dicSize];
    CComplex* complexAtom2Xi = new CComplex[dicSize];


    double* conv_zxi0_w2 = new double[2*dicSize];
    double* conv_zxi0_w2_expinc = new double[2*dicSize];

    CDictionary *expDic = new CExpDictionary;
    CDictionary * gaborDic = new CGaborDictionary;
    CDictionary * triangDic = new CTriangDictionary;
    CDictionary * triangDic2 = new CTriangDictionary;
    CDictionary * bateDic = new CBatemanDictionary;

    ((CExpDictionary*)expDic)->setSignalSize(dicSize);
    ((CGaborDictionary*)gaborDic)->setSignalSize(dicSize);
    ((CTriangDictionary*)triangDic)->setSignalSize(dicSize);
    ((CTriangDictionary*)triangDic2)->setSignalSize(2*dicSize);
    ((CBatemanDictionary*)bateDic)->setSignalSize(dicSize);

    strtParameter* parm = new strtParameter;
    ((strtParameter*)chosenParm)->innerProduct=0.0;
    ((strtParameter*)chosenParm)->rho = 0.0;
    ((strtParameter*)chosenParm)->xi = 0.0;
    ((strtParameter*)chosenParm)->phase = 0.0;
    ((strtParameter*)chosenParm)->a = 0;
    ((strtParameter*)chosenParm)->b = 0;
    double xi = 0.0;
    int tau = 0;
    double maxInnerProd = 0.0;


    int MPType;
    MPType = genData->getMPType();

    FILE* stream;
    char* fileName;

    int d = 0;


    for(j=0;j<dicData->getNumDicBlock();j++)
    {
        // stream = fopen( fileName, "a" );
        // fprintf(stream,"%d --------------------------------------------------------------\n", step);
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
        rise = (dicData->getDicBlock())[j].rise;
        decay = (dicData->getDicBlock())[j].decay;
        int N = dicSize;
        for (int n=0; n<dicSize; n++)
        {
            realAtom[n] = 0.0;
            realAtomDec[n] = 0.0;
            realAtomInc[n] = 0.0;
        }

        for (int n=0; n<2*dicSize; n++)
        {
            realAtomWin[n] = 0.0;
            conv_zxi0_w2[n] = 0.0;
            conv_zxi0_w2_expinc[n] = 0.0;
        }


        if (chosenDic == dicType || chosenDic == 0)
        {
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
            //else if (dataSignal->getType() == 4) // ECG
            //{
            //    Fs = ((CECGSignal*)dataSignal)->getSamplingRate();
            //}
            //else if (dataSignal->getType()==5) //EDA
            //Fs=((CEDASignal*)dataSignal)->getSamplingRate();


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
                else if (dataSignal->getType() == 4) // ECG
                {
                    freqi = 0;
                    freqf = Fs;
                }
                /*
                else if (dataSignal->getType()==5)
                {
                    freqi= ;
                    freqf= ;

                }

                */

            }
            else if (freqf==9999999999)
            {
                freqf = Fs;
            }
            else if (freqi==9999999999)
            {
                freqi = freqf/s;
            }
            else if (freqi==0000000000)
            {
                freqi = freqf/2;
            }

            if (freqf>Fs)
            {
                cout << "Final frequency greater than sampling frequency !!!" << endl;
                printf("final freq. %f - Fs %f\n",freqf,Fs);
                exit(1);
            }

            if (fdiscrtype==1) // linear
            {
                // delta_f = (2*PI/freqf)* freqi;
                Nfreq = (int)(freqf/(2*freqi));
                // cout << Nfreq << endl;
                //Nfreq = (int)ceil(freqf/freqi);
                xi_vec = new double[Nfreq];
                for (i=0;i<Nfreq;i++)
                {
                    // xi_vec[i] = delta_f * i;
                    xi_vec[i] = (2*PI/Fs) * (freqi * i );
                    // cout << xi_vec[i] << endl;
                }
            }
            else if (fdiscrtype==2) // geometric with quarter-tone discretization
            {
                Nfreq = (int)ceil(24 * ( log10(freqf/freqi)/log10(2.0) ) )+1;

                xi_vec = new double[Nfreq];
                xi_vec[0] = 0.0;
                for (i=1;i<Nfreq;i++)
                {
                    xi_vec[i] = (2*PI/Fs) * (freqi * pow (2.0, (double)(i-1)/24) );
                }
            }

            if (dicType == 1) //exponential
            {
                if (MPType == 1)
                {
                    for(k=0;k<Nfreq;k++)
                    {
                        for(tau=0;tau<(double)N;tau+=delta_tau)
                        {
                            // cout << tau << endl;

                            //Decreasing
                            //xi = k * ((2*PI)/s);
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = tau;
                            ((strtParameter*)parm)->b = N-1;

                            a = tau;
                            b = tau + e*s*log(10);

                            if (b>N-1) b = N-1;

                            //save atom parameters

                            if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                            {
                                ((CExpDictionary*)expDic)->setComplexAtom(parm);
                                complexAtomDec = expDic->getComplexAtom();
                                // cout << s  << ' ' << xi << ' ' << tau << endl;
                                // MPTradicional(  residue,
                                //                 complexAtomDec,
                                //                 maxInnerProd,
                                //                 chosenOptPhase,
                                //                 tau,
                                //                 xi,
                                //                 N, //int dicSize
                                //                 s,
                                //                 fileName);

                                computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtomDec, fileName,s,tau,N);

                                dicAtoms[d]->innerProduct = maxInnerProd;
                                dicAtoms[d]->s = s;
                                dicAtoms[d]->rho = 1.0/s;
                                dicAtoms[d]->xi = chosenXi;
                                dicAtoms[d]->phase = chosenOptPhase;
                                dicAtoms[d]->u = tau;
                                dicAtoms[d]->a = tau;
                                dicAtoms[d]->b = N-1;
                                dicAtoms[d]->beta = 0.0;
                                dicAtoms[d]->dicType = dicType;

                            /*if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                tau,
                                                N-1,
                                                0.0);
                                // cout << "AQUI" << endl;
                                // cout << ((strtParameter*)chosenParm)->s  << ' ' << ((strtParameter*)chosenParm)->xi << ' ' << ((strtParameter*)chosenParm)->u << endl;
                            }*/
                            }
                            d++;
                        }
                        for (tau=N-1;tau>0; tau=tau-delta_tau)
                        {
                            // Increasing
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = -1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = tau;

                            a = tau - e*s*log(10);
                            b = tau;

                            if (a<0) a=0;

                            // cout << a << endl;
                            // cout << b << endl;
                            // cout << a0 << endl;
                            // cout << b0 << endl;
                            // cout << "AQUI1" << endl;
                            // cout << d << endl;
                            // cout << a << " " << b << endl;
                            // cout << "AQUI2" << endl;
                            if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                            {
                                ((CExpDictionary*)expDic)->setComplexAtom(parm);
                                complexAtomInc = expDic->getComplexAtom();

                                // for (int k1=0; k1<dicSize;k1++)
                                // {
                                //     complexAtomInc[k1] = *(expDic->getComplexAtom()+k1);
                                //     cout << complexAtomInc[k1].Real() << endl;
                                // }

                                // MPTradicional(  residue,
                                //                 complexAtomInc,
                                //                 maxInnerProd,
                                //                 chosenOptPhase,
                                //                 chosenTau,
                                //                 chosenXi,
                                //                 dic,
                                //                 N, //int dicSize
                                //                 // decincAsymmFlag,
                                //                 s,
                                //                 fileName);

                                computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtomInc,fileName,s,tau,N);

                                dicAtoms[d]->innerProduct = maxInnerProd;
                                dicAtoms[d]->s = -s;
                                dicAtoms[d]->rho = -1.0/s;
                                dicAtoms[d]->xi = chosenXi;
                                dicAtoms[d]->phase = chosenOptPhase;
                                dicAtoms[d]->u = tau;
                                dicAtoms[d]->a = 0.0;
                                dicAtoms[d]->b = tau;
                                dicAtoms[d]->beta = 0.0;
                                dicAtoms[d]->dicType = dicType;


                            /*if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                -s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                0,
                                                tau,
                                                0.0);
                            }*/
                            }
                            d++;
                        }
                    }
                }

                else if (MPType==2)
                {
                    // int delta_tau = s;
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        //decreasing

                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = tau;
                        ((strtParameter*)parm)->a = tau;
                        ((strtParameter*)parm)->b = N-1;


                        a = tau;
                        b = tau + e*s*log(10);

                        if (b>N-1) b = N-1;

                        //save atom parameters

                        if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                        {
                            // cout << d << endl;
                            ((CExpDictionary*)expDic)->setRealAtom(parm);
                            realAtomDec = expDic->getRealAtom();
                            fastMPKolasa(   residue,
                                            maxInnerProd,
                                            chosenOptPhase,
                                            chosenXi,
                                            dicSize,
                                            tau,
                                            s,
                                            realAtomDec,
                                            fileName);

                            dicAtoms[d]->innerProduct = maxInnerProd;
                            dicAtoms[d]->s = s;
                            dicAtoms[d]->rho = 1.0/s;
                            dicAtoms[d]->xi = chosenXi;
                            dicAtoms[d]->phase = chosenOptPhase;
                            dicAtoms[d]->u = tau;
                            dicAtoms[d]->a = tau;
                            dicAtoms[d]->b = N-1;
                            dicAtoms[d]->beta = 0.0;
                            dicAtoms[d]->dicType = dicType;


                            // if (d == 150) cout << maxInnerProd << endl;

                            // memcpy(dicAtoms[d],parm,sizeof(parm));
                            // if (d == 150) cout << dicAtoms[d]->xi << endl;
                        }

                        d++;

                            /*if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                // cout << maxInnerProd << endl;
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                chosenXi,
                                                chosenOptPhase,
                                                tau,
                                                tau,
                                                N-1,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = chosenXi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = tau;
                                // ((strtParameter*)chosenParm)->b = N-1;
                            }*/
                    }

                    for (tau=1; tau<=(double)N; tau=tau+delta_tau)
                    {
                        // increainsg
                        ((strtParameter*)parm)->s = -s;
                        ((strtParameter*)parm)->rho = -1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->a = 0.0;
                        ((strtParameter*)parm)->b = tau;
                        ((strtParameter*)parm)->u = tau;

                        a = tau - e*s*log(10);
                        b = tau;

                        if (a<0) a=0;

                        // cout << a << endl;
                        // cout << b << endl;
                        // cout << a0 << endl;
                        // cout << b0 << endl;
                        // cout << "AQUI1" << endl;
                        // cout << d << endl;
                        // cout << a << " " << b << endl;
                        // cout << "AQUI2" << endl;
                        if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                        {
                            // cout << d << endl;
                            ((CExpDictionary*)expDic)->setRealAtom(parm);
                            realAtomInc = expDic->getRealAtom();

                            fastMPKolasa(   residue,
                                            maxInnerProd,
                                            chosenOptPhase,
                                            chosenXi,
                                            dicSize,
                                            tau,
                                            s,
                                            realAtomInc,
                                            fileName);

                            dicAtoms[d]->innerProduct = maxInnerProd;
                            dicAtoms[d]->s = -s;
                            dicAtoms[d]->rho = -1.0/s;
                            dicAtoms[d]->xi = chosenXi;
                            dicAtoms[d]->phase = chosenOptPhase;
                            dicAtoms[d]->u = tau;
                            dicAtoms[d]->a = 0.0;
                            dicAtoms[d]->b = tau;
                            dicAtoms[d]->beta = 0.0;
                            dicAtoms[d]->dicType = dicType;

                            // memcpy(dicAtoms[d],parm,sizeof(parm));
                        }

                        d++;

                        // cout << chosenXi << endl;
                        // cout << tau << endl;
                        /*if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << s << endl;
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            -s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            0,
                                            tau,
                                            0.0);
                            // cout << endl << "ATOM" << endl;
                            // cout << -1/double(s) << " " << tau << endl;

                            // for (int p = 0; p<512; p++)
                            // {
                            //     cout << realAtomInc[p] << endl;
                            // }

                            // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                            // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                            // ((strtParameter*)chosenParm)->xi = chosenXi;
                            // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                            // ((strtParameter*)chosenParm)->a = 0.0;
                            // ((strtParameter*)chosenParm)->b = tau;
                        }*/
                    }
                }

                else if (MPType==3)
                {
                    // cout << "AQUI" << endl;
                    ((strtParameter*)parm)->s = s;
                    ((strtParameter*)parm)->rho = 1.0/s;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = 0;
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)expDic)->setRealAtom(parm);
                    // realAtomWin = expDic->getRealAtom();
                    // memcpy(realAtomWin, (expDic->getRealAtom()), sizeof(double)*dicSize);
                    for (int k1=0; k1<dicSize;k1++)
                    {
                        realAtomWin[k1+N] = *(expDic->getRealAtom()+k1);
                    }

                    for(k=0;k<Nfreq;k++)
                    {
                        //xi = k * ((2*PI)/s);
                        xi = xi_vec[k];
                        if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                        {
                            // z1
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            // complexAtomXi = expDic->getComplexAtom();
                            // memcpy(complexAtomXi, (expDic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                            }

                            // z2
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = 2*xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            // complexAtom2Xi = expDic->getComplexAtom();
                            // memcpy(complexAtom2Xi, (expDic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                            }


                            //Decreasing exponential
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
                                                    conv_zxi0_w2,
                                                    fileName,
                                                    s);
                            // if (k == 0){
                            // }

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                chosenTau,
                                                N-1,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = xi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = chosenTau;
                                // ((strtParameter*)chosenParm)->b = N-1;
                            }

                            //Increasing exponential
                            decincAsymmFlag = -1;
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
                                                    conv_zxi0_w2_expinc,
                                                    fileName,
                                                    s);
                            // // if (k == 0){
                            // // }

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                -s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                0,
                                                chosenTau,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = -1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = xi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = 0;
                                // ((strtParameter*)chosenParm)->b = chosenTau;
                            }
                        }
                    }
                }
            }

            else if (dicType == 2) // pure sine
            {
                // ((CExpDictionary*)expDic)->setSignalSize(dicSize);

                for (tau=0; tau<(double)N; tau=tau+delta_tau)
                {
                    ((strtParameter*)parm)->rho = 0.0;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = tau;
                    ((strtParameter*)parm)->a = tau;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)expDic)->setRealAtom(parm);
                    realAtom = expDic->getRealAtom();
                    fastMPKolasa(   residue,
                                    maxInnerProd,
                                    chosenOptPhase,
                                    chosenXi,
                                    dicSize,
                                    tau,
                                    s,
                                    realAtom,
                                    fileName);
                    if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    {
                        // cout << maxInnerProd << endl;
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        chosenXi,
                                        chosenOptPhase,
                                        tau,
                                        tau,
                                        N-1,
                                        0.0);
                    }
                }
            }

            else if (dicType == 3) // impulse
            {
                for(i=0;i<N;i++)
                {
                    a = i;
                    b = i;
                    if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                    {
                        innerProd = residue[0][i];
/*                    if (fabs(innerProd)>fabs(maxInnerProd))
                    {
                        // cout << maxInnerProd << endl;
                        maxInnerProd = innerProd;
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        0.0,
                                        0.0,
                                        i,
                                        i,
                                        i,
                                        0.0);
                        // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                        // ((strtParameter*)chosenParm)->xi = 0.0;
                        // ((strtParameter*)chosenParm)->phase = 0.0;
                        // ((strtParameter*)chosenParm)->a = i;
                        // ((strtParameter*)chosenParm)->b = i;
                    }*/
                    }
                    d++;
                }
            }

            else if (dicType == 4) // Gabor
            {
                if (MPType == 1)
                {
                    for(k=0;k<Nfreq;k++)
                    {
                        for(tau=0;tau<N;tau+=delta_tau)
                        {
                            //Decreasing
                            //xi = k * ((2*PI)/s);
                            // cout << s << endl;
                            // cout << tau << endl;
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;

                            a = tau - s*sqrt(e*log(10)/PI);
                            b = tau + s*sqrt(e*log(10)/PI);

                            if (a<0) a=0;
                            if (b>N-1) b=N-1;

                            if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                            {

                                ((CGaborDictionary*)gaborDic)->setComplexAtom(parm);
                                complexAtom = gaborDic->getComplexAtom();


                                // MPTradicional(  residue,
                                //                 complexAtomDec,
                                //                 maxInnerProd,
                                //                 chosenOptPhase,
                                //                 tau,
                                //                 xi,
                                //                 N, //int dicSize
                                //                 s,
                                //                 fileName);

                                computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                                // cout << maxInnerProd << endl;

                                dicAtoms[d]->innerProduct = maxInnerProd;
                                dicAtoms[d]->s = s;
                                dicAtoms[d]->rho = 1.0/s;
                                dicAtoms[d]->xi = chosenXi;
                                dicAtoms[d]->phase = chosenOptPhase;
                                dicAtoms[d]->u = tau;
                                dicAtoms[d]->a = 0.0;
                                dicAtoms[d]->b = N-1;
                                dicAtoms[d]->beta = 0.0;
                                dicAtoms[d]->dicType = dicType;

                                /*if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                                {
                                    setParameters ( chosenParm,
                                                    dicType,
                                                    maxInnerProd,
                                                    s,
                                                    xi,
                                                    chosenOptPhase,
                                                    tau,
                                                    0,
                                                    N-1,
                                                    0.0);
                                }*/
                            }
                            d++;
                        }
                    }

                    // for(k=0;k<Nfreq;k++)
                    // {
                    //     for(tau=0;tau<N;tau+=delta_tau)
                    //     {
                    //         // cout << xi /*<< " " << xi << " " << dicType */<< endl;
                    //         //xi = k * ((2*PI)/s);
                    //         xi = xi_vec[k];
                    //         ((strtParameter*)parm)->s = s;
                    //         ((strtParameter*)parm)->rho = (1.0/(double)s);
                    //         ((strtParameter*)parm)->xi = xi;
                    //         ((strtParameter*)parm)->phase = 0.0;
                    //         ((strtParameter*)parm)->a = tau;
                    //         ((strtParameter*)parm)->b = N-1;
                    //         ((CExpDictionary*)dic)->setComplexAtom(parm);
                    //         complexAtom = dic->getComplexAtom();
                    //         cout << xi << " " << tau << " " << &complexAtom[0] << endl;

                    //         computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                    //         // MPTradicional(  residue,
                    //         //                 complexAtom,
                    //         //                 maxInnerProd,
                    //         //                 chosenOptPhase,
                    //         //                 tau,
                    //         //                 xi,
                    //         //                 // dic,
                    //         //                 N, //int dicSize
                    //         //                 s,
                    //         //                 fileName);
                    //         // cout << maxInnerProd << endl;
                    //         if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    //         {
                    //             setParameters ( chosenParm,
                    //                             dicType,
                    //                             maxInnerProd,
                    //                             s,
                    //                             xi,
                    //                             chosenOptPhase,
                    //                             tau,
                    //                             0,
                    //                             N-1);
                    //         }
                    //     }
                    // }
                }

                else if (MPType==2)
                {
                    int delta_tau = (int)s;
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = tau;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;

                        a = tau - s*sqrt(e*log(10)/PI);
                        b = tau + s*sqrt(e*log(10)/PI);

                        if (a<0) a=0;
                        if (b>N-1) b=N-1;

                        if ((flagOMP == 1) || ((a>=a0 && a<=b0)||(b>=a0 && b<=b0)||(a<=a0 && b>=b0)))
                        {
                            // cout << d << endl;

                            ((CGaborDictionary*)gaborDic)->setRealAtom(parm);
                            realAtom = gaborDic->getRealAtom();
                            // cout << "AQUI" << endl;
                            fastMPKolasa(   residue,
                                            maxInnerProd,
                                            chosenOptPhase,
                                            chosenXi,
                                            dicSize,
                                            tau,
                                            s,
                                            realAtom,
                                            fileName);

                            dicAtoms[d]->innerProduct = maxInnerProd;
                            dicAtoms[d]->s = s;
                            dicAtoms[d]->rho = 1.0/s;
                            dicAtoms[d]->xi = chosenXi;
                            dicAtoms[d]->phase = chosenOptPhase;
                            dicAtoms[d]->u = tau;
                            dicAtoms[d]->a = 0.0;
                            dicAtoms[d]->b = N-1;
                            dicAtoms[d]->beta = 0.0;
                            dicAtoms[d]->dicType = dicType;


                            // memcpy(dicAtoms[d],parm,sizeof(parm));
                        }
                        // cout << "AQUI2" << endl;

                        /*if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            0,
                                            N-1,
                                            0.0);
                        }*/
                        d++;
                    }
                }


                else if (MPType==3)
                {
                    ((strtParameter*)parm)->s = s;
                    ((strtParameter*)parm)->rho = 1.0/s;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = N-1;
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = N-1;
                    ((CGaborDictionary*)gaborDic)->setRealAtom(parm);
                    // memcpy(realAtomWin, (gaborDic->getRealAtom()), sizeof(double)*dicSize);
                    // cout << "realAtomWin" << endl;
                    for (int k1=0; k1<dicSize;k1++)
                    {
                        realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                        realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);
                    }

                    // cout << "AQUI" << endl;
                    // for (int k1=0; k1<2*dicSize; k1++)
                    // {
                    //     cout << realAtomWin[k1] << endl;
                    // }
                    // cout << "ATOM" << endl;
                    // for (int k1=0;k1<2*dicSize;k1++)
                    // {
                    //     cout << realAtomWin[k1] << endl;
                    // }

                    for(k=0;k<Nfreq;k++)
                    {
                        xi = xi_vec[k];
                        if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                        {
                            // cout << s << endl;
                            // cout << xi << endl;
                            // z1
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0; //generating complex exponential with rectangular window
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CGaborDictionary*)expDic)->setComplexAtom(parm);
                            // cout << "complexAtomXi" << endl;
                            // memcpy(complexAtomXi, (gaborDic->getComplexAtom()), sizeof(CComplex)*(dicSize));
                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                                // cout << complexAtomXi[K1] << endl;
                            }
                            // cout << "AQUI" << endl;
                            // for (int k1=0; k1<dicSize; k1++)
                            // {
                            //     cout << complexAtomXi[k1].Real() << " " << complexAtomXi[k1].Imag() << endl;
                            // }


                            // z2
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;  //generating complex exponential with rectangular window
                            ((strtParameter*)parm)->xi = 2*xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CGaborDictionary*)expDic)->setComplexAtom(parm);
                            // memcpy(complexAtom2Xi, (gaborDic->getComplexAtom()), sizeof(CComplex)*(dicSize));
                            // cout << "complexAtom2Xi" << endl;
                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                                // cout << complexAtom2Xi[K1] << endl;
                            }

                            fastMPKolasaModified(   residue,
                                                    maxInnerProd,
                                                    chosenOptPhase,
                                                    chosenTau,
                                                    dicSize,
                                                    1,
                                                    delta_tau,
                                                    xi,
                                                    realAtomWin,
                                                    complexAtomXi,
                                                    complexAtom2Xi,
                                                    conv_zxi0_w2,
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
                                                0,
                                                N-1,
                                                0.0);
                            }
                        }
                    }
                }
            }

            else if (dicType == 5) // triangular
            {
                if (MPType == 1 || MPType == 2)
                {
                    for(tau=0;tau<N;tau+=delta_tau)
                    {
                        // cout << "a: " << a << endl;
                        for(a=0;a<s;a+=delta_tau)
                        {
                            // cout << "tau: " << tau << endl;
                            //Decreasing
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = a;
                            ((strtParameter*)parm)->b = (int)s-1-a;
                            ((CTriangDictionary*)triangDic)->setComplexAtom(parm);
                            complexAtom = triangDic->getComplexAtom();

                            // cout << "AQUI" << endl;
                            // cout <<" s- " << s << "tau- " << tau << "a- " << a << endl;

                            // for (int k1=0;k1<N;k1++)
                            // {
                            //     // complexAtomXi[k1] = *(triangDic->getComplexAtom()+k1);
                            //     // cout << residue[0][k1] << endl;
                            //     cout << complexAtom[k1].Real() << endl;
                            // }

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                            // cout << maxInnerProd << endl;
                            // cout << "maxInnerProd " << maxInnerProd << endl;
                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                // cout << "s- " << s << "tau- " << tau << "a- " << a << endl;
                                // cout << maxInnerProd << endl;
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                a,
                                                (int)s-a-1,
                                                0.0);
                            }
                            d++;
                        }
                    }
                }

                else if (MPType==3)
                {
                    // cout << s << endl;
                    for(a=1;a<(int)s;a+=delta_tau)
                    {
                        // cout << a << endl;
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = N;
                        ((strtParameter*)parm)->a = a;
                        ((strtParameter*)parm)->b = (int)s-a-1;
                        ((CTriangDictionary*)triangDic2)->setRealAtom(parm);

                        // cout << "s " << s << " - a " << a << endl;
                        for (int k1=0; k1<2*dicSize;k1++)
                        {
                            realAtomWin[k1] = *(triangDic2->getRealAtom()+k1);
                            // cout << realAtomWin[k1] << endl;
                            // realAtomWin[(2*N-1)-k1] = *(triangDic2->getRealAtom()+k1);
                        }
                        // z1
                        ((strtParameter*)parm)->rho = 0.0;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CTriangDictionary*)expDic)->setComplexAtom(parm);

                        for (int k1=0; k1<dicSize;k1++)
                        {
                            complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                        }

                        // z2
                        ((strtParameter*)parm)->rho = 0.0;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CTriangDictionary*)expDic)->setComplexAtom(parm);

                        // cout << "AQUI1" << endl;
                        for (int k1=0; k1<dicSize;k1++)
                        {
                            complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                            // cout << residue[1][k1] << endl;
                        }
                        // cout << "AQUI2" << endl;

                        fastMPKolasaModified(   residue,
                                                maxInnerProd,
                                                chosenOptPhase,
                                                chosenTau,
                                                dicSize,
                                                1,
                                                delta_tau,
                                                xi,
                                                realAtomWin,
                                                complexAtomXi,
                                                complexAtom2Xi,
                                                conv_zxi0_w2,
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
                                        a,
                                        (int)s-a-1,
                                        0.0);
                        }
                    }
                }
            }

            if (dicType == 6) //Bateman
            {
                if (MPType == 1)
                {
                    for (beta = 2.0; beta > 1.0/s; beta/=2.0)
                    {
                        for(k=0;k<Nfreq;k++)
                        {
                            for(tau=0;tau<(double)N;tau+=delta_tau)
                            {
                                xi = xi_vec[k];
                                ((strtParameter*)parm)->rho = 1.0/s;
                                ((strtParameter*)parm)->xi = xi;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = tau;
                                ((strtParameter*)parm)->a = tau;
                                ((strtParameter*)parm)->b = N-1;
                                ((strtParameter*)parm)->beta = beta;
                                ((CBatemanDictionary*)bateDic)->setComplexAtom(parm);
                                complexAtom = bateDic->getComplexAtom();

                                computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);

                                if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                                {
                                    setParameters ( chosenParm,
                                                    dicType,
                                                    maxInnerProd,
                                                    s,
                                                    xi,
                                                    chosenOptPhase,
                                                    tau,
                                                    tau,
                                                    N-1,
                                                    beta);
                                }
                                d++;
                            }
                        }
                    }
                }

                else if (MPType==2)
                {
                    // int delta_tau = s;
                    for (beta=1.0; beta > 1.0/s; beta/=(2.0))
                    {
                        for (tau=0; tau<(double)N; tau=tau+delta_tau)
                        {
                            ((strtParameter*)parm)->rho = 1/s;
                            ((strtParameter*)parm)->beta = beta;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = tau;
                            ((strtParameter*)parm)->b = N-1;
                            ((CBatemanDictionary*)bateDic)->setRealAtom(parm);
                            realAtom = bateDic->getRealAtom();
                            fastMPKolasa(   residue,
                                            maxInnerProd,
                                            chosenOptPhase,
                                            chosenXi,
                                            dicSize,
                                            tau,
                                            s,
                                            realAtom,
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
                                                tau,
                                                N-1,
                                                beta);
                            }
                        }
                    }
                }

                else if (MPType==3)
                {
                    for (beta=1.0; beta > 1.0/s; beta/=2.0)
                    {
                        ((strtParameter*)parm)->rho = 1/s;
                        ((strtParameter*)parm)->beta = beta;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CBatemanDictionary*)bateDic)->setRealAtom(parm);

                        for (int k1=0; k1<dicSize;k1++)
                        {
                            realAtomWin[k1+N] = *(bateDic->getRealAtom()+k1);
                        }

                        for(k=0;k<Nfreq;k++)
                        {
                            //xi = k * ((2*PI)/s);
                            xi = xi_vec[k];
                            if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                            {
                                // z1
                                ((strtParameter*)parm)->rho = 0.0;
                                ((strtParameter*)parm)->beta = 0.0;
                                ((strtParameter*)parm)->xi = 0.0;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = 0;
                                ((strtParameter*)parm)->a = 0;
                                ((strtParameter*)parm)->b = N-1;
                                ((CExpDictionary*)expDic)->setComplexAtom(parm);

                                for (int k1=0; k1<dicSize;k1++)
                                {
                                    complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                                }

                                // z2
                                ((strtParameter*)parm)->rho = 0.0;
                                ((strtParameter*)parm)->beta = 0.0;
                                ((strtParameter*)parm)->xi = 0.0;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = 0;
                                ((strtParameter*)parm)->a = 0;
                                ((strtParameter*)parm)->b = N-1;
                                ((CExpDictionary*)expDic)->setComplexAtom(parm);

                                for (int k1=0; k1<dicSize;k1++)
                                {
                                    complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                                }


                                fastMPKolasaModified(   residue,
                                                        maxInnerProd,
                                                        chosenOptPhase,
                                                        chosenTau,
                                                        dicSize,
                                                        1,
                                                        delta_tau,
                                                        xi,
                                                        realAtomWin,
                                                        complexAtomXi,
                                                        complexAtom2Xi,
                                                        conv_zxi0_w2,
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
                                                    chosenTau,
                                                    N-1,
                                                    beta);
                                }
                            }
                        }
                    }
                }
            }

            if (xi_vec) delete [] xi_vec;
        }
    }


    // cout << "AQUI" << endl;

    for (int atom=0; atom<d; atom++)
    {
        if (fabs(dicAtoms[atom]->innerProduct)>fabs(maxInnerProd2))
        {
            chosenAtom = atom;
            maxInnerProd2 = dicAtoms[atom]->innerProduct;


            // cout << dicAtoms[atom]->innerProduct << endl;
        }
    }

    // cout << "AQUI2" << endl;
    // cout << dicAtoms[chosenAtom]->innerProduct << endl;

    // cout << dicAtoms[chosenAtom]->dicType << endl;

    if (dicAtoms[chosenAtom]->dicType == 1)
    {
        if (dicAtoms[chosenAtom]->s > 0)
        {
            a0 = dicAtoms[chosenAtom]->a;
            b0 = dicAtoms[chosenAtom]->u + e*dicAtoms[chosenAtom]->s*log(10);

            if (b0>dicSize-1) b0 = dicSize-1;
        }
        else if (dicAtoms[chosenAtom]->s < 0)
        {
            a0 = dicAtoms[chosenAtom]->u + e*dicAtoms[chosenAtom]->s*log(10);
            b0 = dicAtoms[chosenAtom]->b;
            if (a0<0) a0=0;
        }
    }

    else if (dicAtoms[chosenAtom]->dicType == 3)
    {
        a0 = dicAtoms[chosenAtom]->u;
        b0 = dicAtoms[chosenAtom]->u;
    }

    else if (dicAtoms[chosenAtom]->dicType == 4)
    {
        a0 = dicAtoms[chosenAtom]->u - dicAtoms[chosenAtom]->s*sqrt(e*log(10)/PI);
        b0 = dicAtoms[chosenAtom]->u + dicAtoms[chosenAtom]->s*sqrt(e*log(10)/PI);
        // cout << s << " " << a0 << " " << b0 << endl;
        if (a0<0) a0=0;
        if (b0>dicSize-1) b0=dicSize-1;

    }

    // cout << a0 << " " << b0 << endl;
    // cout << chosenAtom << endl;


    setParameters(  chosenParm,
                    dicAtoms[chosenAtom]->dicType,
                    dicAtoms[chosenAtom]->innerProduct,
                    dicAtoms[chosenAtom]->s,
                    dicAtoms[chosenAtom]->xi,
                    dicAtoms[chosenAtom]->phase,
                    dicAtoms[chosenAtom]->u,
                    dicAtoms[chosenAtom]->a,
                    dicAtoms[chosenAtom]->b,
                    dicAtoms[chosenAtom]->beta);

    if (parm) delete parm;

    if (expDic) delete expDic;
    if (gaborDic) delete gaborDic;
    if (triangDic) delete triangDic;
    if (triangDic) delete triangDic2;
    if (bateDic) delete bateDic;

    if (conv_zxi0_w2_expinc) delete[] conv_zxi0_w2_expinc;
    if (conv_zxi0_w2) delete[] conv_zxi0_w2;

    return;
}

void MPTradicional( cgMatrix<double>& residue,
                    CComplex* complexAtom,
                    double& maxInnerProd,
                    double& chosenOptPhase,
                    int& chosenTau,
                    double& chosenXi,
                    // CDictionary* dic,
                    int N, //int dicSize
                    // int decincAsymmFlag,
                    double s,
                    // int delta_tau,
                    // double* xi_vec,
                    // int Nfreq,
                    char* fileName)
{
    int k;
    double opt_phase = 0;
    double innerProd = 0.0;

    computeOptimumPhase(residue,chosenOptPhase,innerProd,N,chosenXi,complexAtom, fileName,s,chosenTau,N);

    if (fabs(innerProd)>fabs(maxInnerProd))
    {
        maxInnerProd = innerProd;
        chosenOptPhase = opt_phase;
    }
    return;
}

// void MPTradicional( cgMatrix<double>& residue,
//                     CComplex* complexAtom,
//                     double& maxInnerProd,
//                     double& chosenOptPhase,
//                     int& chosenTau,
//                     double& chosenXi,
//                     // CDictionary* dic,
//                     int N, //int dicSize
//                     // int decincAsymmFlag,
//                     int s,
//                     // int delta_tau,
//                     // double* xi_vec,
//                     // int Nfreq,
//                     char* fileName)
// {
//     int k;
//     double opt_phase = 0;
//     // CComplex* complexAtom;
//     int tau = 0;
//     double xi = 0;
//     // strtParameter* parm;
//     // parm = new strtParameter;
//     double innerProd = 0.0;
//     // for(k=0;k<Nfreq;k++)
//     // {
//     //     // cout << "k - " << k << endl;
//     //     for(tau=0;tau<N;tau+=delta_tau)
//     //     {
//     //         if (decincAsymmFlag == 1)
//     //         {
//     //             //xi = k * ((2*PI)/s);
//     //             xi = xi_vec[k];
//     //             ((strtParameter*)parm)->rho = (1.0/(double)s);
//     //             ((strtParameter*)parm)->xi = xi;
//     //             ((strtParameter*)parm)->phase = 0.0;
//     //             ((strtParameter*)parm)->a = tau;
//     //             ((strtParameter*)parm)->b = N-1;
//     //             ((CExpDictionary*)dic)->setComplexAtom(parm);
//     //             complexAtom = dic->getComplexAtom();
//     computeOptimumPhase(residue,opt_phase,innerProd,N,xi,complexAtom, fileName,s,tau,N);
//     // cout << innerProd << endl;

//     if (fabs(innerProd)>fabs(maxInnerProd))
//     {
//         // cout << innerProd << endl;
//         maxInnerProd = innerProd;
//         chosenXi = xi;
//         chosenOptPhase = opt_phase;
//         chosenTau = tau;
//     }
//             // }
//             // if (decincAsymmFlag == -1)
//             // {
//             //     //xi = k * ((2*PI)/s);
//             //     xi = xi_vec[k];
//             //     ((strtParameter*)parm)->rho = -(1.0/(double)s);
//             //     ((strtParameter*)parm)->xi = xi;
//             //     ((strtParameter*)parm)->phase = 0.0;
//             //     ((strtParameter*)parm)->a = 0.0;
//             //     ((strtParameter*)parm)->b = tau;
//             //     ((CExpDictionary*)dic)->setComplexAtom(parm);
//             //     complexAtom = dic->getComplexAtom();

//                 // computeOptimumPhase(residue,opt_phase,innerProd,N,xi,complexAtom, fileName,s,tau,N);


//                 // if (fabs(innerProd)>fabs(maxInnerProd))
//                 // {
//                 //     maxInnerProd = innerProd;
//                 //     chosenXi = xi;
//                 //     chosenOptPhase = opt_phase;
//                 //     chosenTau = tau;
//                 // }

//     //         }
//     //     }
//     // }
//     // delete parm;
//     return;
// }

void fastMPKolasa(  cgMatrix<double>& residue,
                    double& maxInnerProd,
                    double& chosenOptPhase,
                    double& chosenXi,
                    int N, //dicSize
                    int tau,
                    double s,
                    double* realAtom,
                    char* fileName)
{
     //cout << "AQUI3" << endl;
     //FILE* stream;
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

    maxInnerProd = 0.0;


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
    int delta_k = (int)(ceil((double)N/s));

    for (k=0;k<NC-1;k=k+delta_k)
    {
        // cout << NC << endl;
        // cout << k << "AQUI" << endl;
        innerProd_xp = out1[k][0];  // /sqrt((double)N);

        innerProd_xq = -out1[k][1]; // /sqrt((double)N);
        //cout<<k<<endl<<out2[2*k][0]<<endl;
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

        if (fabs(a1) < 1e-160)
        {
            a1 = 0.0;
        }

        if ( (k == 0) || (k == k_pi) )
        {
            opt_phase = 0;
            innerProd = innerProd_xp  / sqrt(innerProd_pp);
            // cout << innerProd << endl;
        }
        else if (a1 == 0)
        {
            opt_phase = (double)(PI/2);
            innerProd = -innerProd_xq / sqrt(innerProd_qq);
            // cout << innerProd << endl;
        }
        else if ( (a1!=0) && (k != 0) )
        {
            opt_phase = atan( -(b1/a1) );
            innerProd = (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                        sqrt(a1*a1*innerProd_pp +
                             b1*b1*innerProd_qq +
                             2*a1*b1*innerProd_pq);
            // if (innerProd > 1e6)
            // {
                // cout << innerProd << endl;
                // cout << a1 << endl;
                // cout << a1*a1*innerProd_pp << endl;
                // cout << b1*b1*innerProd_qq << endl;
                // cout << 2*a1*b1*innerProd_pq << endl;
            // }

        }

        xi = (k*2*PI)/N;
        // cout << xi << endl;
        if (fabs(innerProd)>fabs(maxInnerProd))
        {
            maxInnerProd = innerProd;
            chosenXi = xi,
            chosenOptPhase = opt_phase;
        }


        // For debug
        // fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
        //                innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
        //                innerProd_pq, s, (double)tau, xi,opt_phase);

        /* FILE* stream;
         stream = fopen("kolasa.out","a");
         fprintf (   stream,"IP- %15.8f rho - %15.8f xi - %15.8f optph - %15.8f tau - %15.8f tau - %15.8f N-1 %15.8f xp - %15.8f xq - %15.8f pp - %15.8f qq - %15.8f pq - %15.8f \n",
                     innerProd,
                     1.0/(double)s,
                     xi,
                     opt_phase,
                     (double)tau,
                     (double)tau,
                     (double)(N-1),
                     innerProd_xp,
                     innerProd_xq,
                     innerProd_pp,
                     innerProd_qq,
                     innerProd_pq);
         fflush(stream);
         fclose(stream);
        */
    }

    fftw_destroy_plan ( plan_forward1 );
    fftw_destroy_plan ( plan_forward2 );


    fftw_free ( in1 );
    fftw_free ( out1 );
    fftw_free ( in2 );
    fftw_free ( out2 );
    // cout << "AQUI4" << endl;
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
                            double* conv_zxi0_w2,
                            char* fileName,
                            double s)

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

    fftw_complex    *z1, *z2;
    fftw_complex    *w1, *w2;
    fftw_complex    *z1_fft, *z2_fft, *w1_fft, *w2_fft ;
    fftw_complex    *conv_zw1, *conv_zw2;
    fftw_complex    *prod_cvec_zw1, *prod_cvec_zw2;

    fftw_plan       plan_forward1a, plan_forward1b, plan_backward1,
                    plan_forward2a, plan_forward2b, plan_backward2;


 
    // Memory allocation
         /*
         z1=fftw_complex *fftw_alloc_complex(size_t (2*N));
         z2 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         w1 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         w2 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         z1_fft = fftw_complex *fftw_alloc_complex(size_t (2*N));
         z2_fft = fftw_complex *fftw_alloc_complex(size_t (2*N));
         w1_fft = fftw_complex *fftw_alloc_complex(size_t (2*N));
         w2_fft = fftw_complex *fftw_alloc_complex(size_t (2*N));
         conv_zw1 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         conv_zw2 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         prod_cvec_zw1 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         prod_cvec_zw2 = fftw_complex *fftw_alloc_complex(size_t (2*N));
         */
         //(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n)
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

    maxInnerProd = 0.0;

    for (i=0;i<2*N;i++) //for (i=0;i<2*N;i++)
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

        double tol = 1e-10;
        //  int k_pi;

    // if (s==2 && xi==0.0 && chosenTau==0)
    // {
         //cout << endl << endl << endl << "AQUIk" << endl;

    //     cout << endl << "realAtomWin1: " << " - realAtomWin2: " << endl;
    // // }
    //w1[0][0]=realAtomWin[0];
    //w2[0][0]=realAtomWin[0];
    for (i=0;i<N;i++) //for (i=0;i<2*N;i++)
    {
        w1[i][0] = realAtomWin[N-1-i]; //w1[i][0] = realAtomWin[2*N-1-i];
        w2[i][0] = realAtomWin[N-1-i] * realAtomWin[N-1-i];  //w2[i][0] = realAtomWin[2*N-1-i] * realAtomWin[2*N-1-i];
         //cout << "1- "<<w1[i][0] << "2- "<<w2[i][0]<< endl;
         //if (s==5 && xi==0.0 && chosenTau==0)
         //{
           //  cout << w1[i][0] << "----- " << w2[i][0] << endl;
         //}
    }
    // for (i=0;i<N;i++)
    // {
    //     w1[i][0] = realAtomWin[N-1-i];
    //     w2[i][0] = realAtomWin[N-1-i] * realAtomWin[N-1-i];
    //     // if (s==2 && xi==0.0 && chosenTau==0)
    //     // {
    //     //     cout << w1[i][0] << " " << w2[i][0] << endl;
    //     // }
    // }


    // if (s==2 && xi==0.0 && chosenTau==0)
    // {
    //     cout << endl << "complexAtom1: " << " - complexAtom2: " << endl;
    // }
    for (i=0;i<N;i++) //for (i=0;i<2*N;i++)
    {

        z1[i][0] = residue[0][i];//complexAtomXi[i].Real() * residue[0][i];
        z1[i][1] = complexAtomXi[i].Imag() * residue[0][i];
        z2[i][0] = complexAtom2Xi[i].Real();
        z2[i][1] = complexAtom2Xi[i].Imag();
         //if (xi==0.0 && chosenTau==0)
         //{
         //   cout << "CAX"<< endl;
         //cout <<"z1Re---"<< z1[i][0] << "z1Im--- " << z1[i][1] << "z2Re---- " << z2[i][0] << "z2Im----- " << z2[i][1] <<"res--"<<residue[0][i]<<"CA"<<complexAtom2Xi[i].Real() << endl;
         //}
    }
   /*FILE* stream;
             stream = fopen("BATMP.out","a");
             for(int k1=0;k1<N;k1++)
             {
             fprintf (   stream," CAR- %15.8f CAR2 - %15.8f  Z2R - %15.8f \n",
                     complexAtomXi[k1].Real(),
                     complexAtom2Xi[k1].Real(),
                     z2[k1][0]);
                }
             fflush(stream);
             fclose(stream);
            */
    // if (s==2 && xi==0.0 && chosenTau==0)
    // {
    //     cout << endl << endl << endl;
    // }

    fftw_execute (plan_forward1a);
    fftw_execute (plan_forward1b);

    for (i=0;i<2*N;i++)  //for (i=0;i<2*N;i++)
    {
        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*decincAsymmFlag*(w1_fft[i][1]));
        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*decincAsymmFlag*(w1_fft[i][1]));
        // if (step == 4)
        // if (s==2 && xi==0.0 && chosenTau==0)
        // {
        //    cout << prod_cvec_zw1[i][0]<< " + j* " << prod_cvec_zw1[i][1] << endl;
        // }
    }

    fftw_execute ( plan_backward1);

    // zw2
    fftw_execute (plan_forward2a);
    fftw_execute (plan_forward2b);

    // product between complex vectors
    // cout << "---------------------------" << endl;
    for (i=0;i<2*N;i++)  //for (i=0;i<2*N;i++)
    {
        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*decincAsymmFlag*(w2_fft[i][1]));
        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*decincAsymmFlag*(w2_fft[i][1]));
        // if (step == 4)
        // {
           // cout << i << " - prod_cvec_zw2[i][0]: " << prod_cvec_zw2[i][0] << " - prod_cvec_zw2[i][1]: " << prod_cvec_zw2[i][1] << endl;
        // // }
        // cout << "W2["<<i<<"][0]:   " << w2_fft[i][0] << " - W2["<<i<<"][1]:   " << w2_fft[i][1] << endl;
    }

    fftw_execute ( plan_backward2);

    if (xi==0.0) //(k==0)
    {
        for (i = 0;i<2*N;i++)  //for (i=0;i<2*N;i++)
        {
            conv_zxi0_w2[i] = conv_zw2[i][0];
        }
    }

    if (decincAsymmFlag==1) //Decreasing
    {
        // Compute inner product and optimum phase
        for (tau=0; tau<(double)N; tau=tau+delta_tau)
        {
            //cout << "Tau - "<< tau << "N - " << (double)N << endl;
            innerProd_xp = conv_zw1[tau+N-1][0]/(double)(2*N);
            innerProd_xq = conv_zw1[tau+N-1][1]/(double)(2*N);
            innerProd_pp = 0.5*(conv_zxi0_w2[tau+N-1] + conv_zw2[tau+N-1][0])/(double)(2*N);
            innerProd_qq = 0.5*(conv_zxi0_w2[tau+N-1] - conv_zw2[tau+N-1][0])/(double)(2*N);
            innerProd_pq = 0.5*(conv_zw2[tau+N-1][1]/(double)(2*N));
             //if (s == 5 && xi == 0.0 && tau == 0)
             //{
                //cout << endl << "s: " << conv_zxi0_w2[tau] << " - xi: " << conv_zw2[tau][0] << " - tau: " << (double)(2*N) << " - xp: " << innerProd_xp << " - xq:" << innerProd_xq << " - pp:" << innerProd_pp << " - qq:" << innerProd_qq << " - pq:" << innerProd_pq << endl << endl;
             //}

            a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
            b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

            //if (s == 5 && xi == 0.0 && tau == 0)
            //{
             //   cout << " 1 - "<< a1 <<"  "<<"2 - " <<b1<<endl;//<< innerProd_qq << endl;
            //}
             //cout << " 3 - " << innerProd_xp << "4- " << innerProd_pp << endl;
             //cout << " 5 - " << innerProd_xq << "6- " << innerProd_pq << endl;

            if (( xi == 0.0 )||( (int)(xi*10000) == (int)(PI*10000)))
            {
                opt_phase = 0;
                if (fabs(innerProd_pp) > tol)
                {
                    innerProd = innerProd_xp / sqrt(innerProd_pp);
                }
                else
                {
                    innerProd = 0.0;
                    // cout << " 1 - " << innerProd_pp << endl;
                }
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(PI/2);
                if (fabs(innerProd_qq) > tol)
                {
                    innerProd = -innerProd_xq / sqrt(innerProd_qq);
                }
                else
                {
                    innerProd = 0.0;
                    // cout << " 2 - " << innerProd_qq  << endl;
                }
            }
            else //if ( (a1 != 0) && (xi != 0) )
            {
                opt_phase = atan( -(b1/a1) );
                if ( fabs(a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) > tol)
                {
                    innerProd =  (((a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1)) /
                                 (sqrt(a1*a1*innerProd_pp +
                                 b1*b1*innerProd_qq +
                                 2*a1*b1*innerProd_pq)));
                                 //cout<< "2 - " << tol << endl;
                }
                else
                {
                    innerProd = 0.0;
                     //cout << " 3 - " << (a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) << endl;
                }
            }

            if (fabs(innerProd)>fabs(maxInnerProd))
            {
                maxInnerProd = innerProd;
                chosenTau = tau;
                chosenOptPhase = opt_phase;
                //cout << "IP - "<< innerProd<<endl;
            }

             //cout << a1   << endl;
             /*
             FILE* stream;
             stream = fopen("fastkolasa.out","a");
             fprintf (   stream," IP- %15.8f rho - %15.8f x1- %15.8f optph- %15.8f tau - %15.8f tau - %15.8f N-1 - %15.8f xp -  %15.8f xq - %15.8f pp -  %15.8f qq - %15.8f pq -%15.8f \n",
                     innerProd,
                     1/(double)s,
                     xi,
                     opt_phase,
                     (double)tau,
                     (double)tau,
                     (double)(N-1),
                     innerProd_xp,
                     innerProd_xq,
                     innerProd_pp,
                     innerProd_qq,
                     innerProd_pq);
             fflush(stream);
             fclose(stream);
             */
             //if (step == 4)
             //{
            //  cout << "innerProd: " << innerProd << endl;
            // }

        }
    }

    else if (decincAsymmFlag==-1) //Increasing
    {
        // Compute inner product and optimum phase
        for (tau=N-1; tau>0; tau=tau-delta_tau)
        {
            innerProd_xp = conv_zw1[(tau+N+1)%(2*N)][0]/(double)(2*N);
            innerProd_xq = conv_zw1[(tau+N+1)%(2*N)][1]/(double)(2*N);
            innerProd_pp = 0.5*(conv_zxi0_w2[(tau+N+1)%(2*N)] + conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
            innerProd_qq = 0.5*(conv_zxi0_w2[(tau+N+1)%(2*N)] - conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
            innerProd_pq = 0.5*(conv_zw2[(tau+N+1)%(2*N)][1]/(double)(2*N));

            a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
            b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

            if (( xi == 0.0 )||( (int)(xi*10000) == (int)(PI*10000)))
            {
                opt_phase = 0;
                if (fabs(innerProd_pp) > tol)
                {
                    innerProd = innerProd_xp / sqrt(innerProd_pp);
                }
                else
                {
                    innerProd = 0.0;
                    // cout << " 1 - " << innerProd_pp << endl;
                }
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(PI/2);
                if (fabs(innerProd_qq) > tol)
                {
                    innerProd = -innerProd_xq / sqrt(innerProd_qq);
                }
                else
                {
                    innerProd = 0.0;
                    // cout << " 2 - " << innerProd_qq  << endl;
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
                    // cout << " 3 - " << (a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) << endl;
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
        }
    }
    // cout << "AQUI 1" << endl;
    fftw_destroy_plan ( plan_forward1a );
    fftw_destroy_plan ( plan_forward1b );
    fftw_destroy_plan ( plan_backward1 );
    fftw_destroy_plan ( plan_forward2a );
    fftw_destroy_plan ( plan_forward2b );
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
    fftw_free(prod_cvec_zw1);
    fftw_free(prod_cvec_zw2);

    return;
}

void setParameters (strtParameter* chosenParm,
                    int dicType,
                    double innerProd,
                    double s,
                    double xi,
                    double opt_phase,
                    int tau,
                    int a,
                    int b,
                    double beta)
{

    // cout << "AQUI" << endl;
    if (dicType == 1) //exp
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->s = (int)s;
        ((strtParameter*)chosenParm)->rho = 1.0/s;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->u = tau;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
        ((strtParameter*)chosenParm)->beta = 0.0;
        ((strtParameter*)chosenParm)->dicType = dicType;
    }
    else if (dicType == 2) //sine
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->s = 512;
        ((strtParameter*)chosenParm)->rho = 0.0;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->u = tau;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
        ((strtParameter*)chosenParm)->beta = 0.0;
        ((strtParameter*)chosenParm)->dicType = dicType;
    }

    else if (dicType == 3) //impulse
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->s = (int)s;
        ((strtParameter*)chosenParm)->rho = 1.0/s;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->u = tau;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
        ((strtParameter*)chosenParm)->beta = 0.0;
        ((strtParameter*)chosenParm)->dicType = dicType;
    }

    else if (dicType == 4) //Gabor
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->s = (int)s;
        ((strtParameter*)chosenParm)->rho = 1.0/s;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->u = tau;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
        ((strtParameter*)chosenParm)->beta = 0.0;
        ((strtParameter*)chosenParm)->dicType = dicType;
    }

    else if (dicType == 5) //Triangular
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->s = (int)s;
        ((strtParameter*)chosenParm)->rho = 1.0/s;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->u = tau;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
        ((strtParameter*)chosenParm)->beta = 0.0;
        ((strtParameter*)chosenParm)->dicType = dicType;
    }

    else if (dicType == 6) //Bateman
    {
        ((strtParameter*)chosenParm)->innerProduct = innerProd;
        ((strtParameter*)chosenParm)->s = s;
        ((strtParameter*)chosenParm)->rho = 1/s;
        ((strtParameter*)chosenParm)->beta = beta;
        ((strtParameter*)chosenParm)->xi = xi;
        ((strtParameter*)chosenParm)->phase = opt_phase;
        ((strtParameter*)chosenParm)->u = tau;
        ((strtParameter*)chosenParm)->a = a;
        ((strtParameter*)chosenParm)->b = b;
        ((strtParameter*)chosenParm)->dicType = dicType;
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
//         ((int)(10000*xi) == (int)(10000*PI))) //caso n�o haja sen�ide
//     {
//         opt_phase = 0;
//         innerProd = innerProductReal[0][0]/p;
//     }
//     else if (a1.getData(0,0) == 0)
//     {
//         opt_phase = (double)(PI/2);
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
                            double s,
                            double tau,
                            int N)
{
    // double opt_phase = 0;
    int i;
    cgMatrix<double> realPartComplexDic(1,signalSize,0.0);
    cgMatrix<double> imagPartComplexDic(1,signalSize,0.0);

    innerProd = 0.0;

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
    // cout << p << " " << q << endl;
    cgMatrix<double> innerProductRealImag(1,1,innerProdRealImag);
    //innerProductRealImag = realPartComplexDic * (imagPartComplexDic.transpose());

    cgMatrix<double> a1;
    cgMatrix<double> b1;
    a1 = innerProductReal*(q*q) - innerProductImag * innerProductRealImag;
    b1 = innerProductImag*(p*p) - innerProductReal * innerProductRealImag;

    if ( (xi == 0) ||
        ((int)(10000*xi) == (int)(10000*PI))) //caso n�o haja sen�ide
    {
        opt_phase = 0;
        innerProd = innerProductReal[0][0]/p;
    }
    else if (a1.getData(0,0) == 0)
    {
        // cout << "AQUI" << endl;
        opt_phase = (double)(PI/2);
        innerProd = -innerProductImag[0][0]/q;
        // cout << "innerProd: " << innerProd << endl;
    }
    else if (   (a1.getData(0,0)!=0) && (xi!=0) )
    {
        // cout << "AQUI3" << endl;
        opt_phase = atan( -(b1.getData(0,0)/a1.getData(0,0)) );
        innerProd = (a1[0][0]/fabs(a1[0][0]))*(innerProductReal[0][0]*a1[0][0] + innerProductImag[0][0]*b1[0][0])/
        sqrt(a1[0][0]*a1[0][0]*p*p+b1[0][0]*b1[0][0]*q*q+2*a1[0][0]*b1[0][0]*innerProductRealImag[0][0]);
    }
    /*
     cout << innerProd << endl;
     cout << innerProd << endl;
     FILE* stream;
     stream = fopen("MPTrad.out","a");
      fprintf (   stream," IP- %15.8f rho - %15.8f xi- %15.8f optph- %15.8f tau - %15.8f tau - %15.8f N-1 - %15.8f xp -  %15.8f xq - %15.8f pp -  %15.8f qq - %15.8f pq -%15.8f \n",
                 innerProd,
                 1.0/(double)s,
                 xi,
                 opt_phase,
                 (double)tau,
                 (double)tau,
                 (double)(N-1),
                 innerProdReal,
                 innerProdImag,
                 (p*p),
                 (q*q),
                 innerProdRealImag);
     fflush(stream);
     fclose(stream);
     */
    return;
}

void adjustParameters ( cgMatrix<double>& residue,
                        strtParameter* parm)
{
    CDictionary* dic;

    if ( parm->dicType==1 || parm->dicType==2 || parm->dicType==3 )
    {
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*)dic)-> adjustParameters(residue,parm);
    }

    else if ( parm->dicType==4 )
    {
        dic = new CGaborDictionary;
        ((CGaborDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CGaborDictionary*)dic)-> adjustParameters(residue,parm);
    }

    else if ( parm->dicType==5 )
    {
        dic = new CTriangDictionary;
        ((CTriangDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CTriangDictionary*)dic)-> adjustParameters(residue,parm);
    }

    else if ( parm->dicType==6 )
    {
        dic = new CBatemanDictionary;
        ((CBatemanDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CBatemanDictionary*)dic)-> adjustParameters(residue,parm);
    }


    if (dic) delete dic;
}

void optimizeContinuousParms(   cgMatrix<double>& residue,
                                                strtParameter* parm)

{
    CDictionary* dic;

    if ( parm->dicType==1 || parm->dicType==2 || parm->dicType==3 )
    {
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
       // ((CExpDictionary*)dic)-> optimizeContinuousParms(residue,parm);
    }

    else if ( parm->dicType==4 )
    {
        dic = new CGaborDictionary;
        ((CGaborDictionary*)dic)-> setSignalSize(residue.getColumns());
        //((CGaborDictionary*)dic)-> optimizeContinuousParms(residue,parm);
    }

    else if ( parm->dicType==5 )
    {
        dic = new CTriangDictionary;
        ((CTriangDictionary*)dic)-> setSignalSize(residue.getColumns());
        //((CTriangDictionary*)dic)-> optimizeContinuousParms(residue,parm);
    }

    else if ( parm->dicType==6 )
    {
        dic = new CBatemanDictionary;
        ((CBatemanDictionary*)dic)->setSignalSize(residue.getColumns());
        ((CBatemanDictionary*)dic)-> optimizeContinuousParms(residue,parm);
    }


    if (dic) delete dic;

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

    if ( parm->dicType==1 || parm->dicType==2 || parm->dicType==3 )
    {
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CExpDictionary*) dic)->getRealAtom();
    }

    else if ( parm->dicType==4)
    {
        dic = new CGaborDictionary;
        ((CGaborDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CGaborDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CGaborDictionary*) dic)->getRealAtom();
    }


    else if ( parm->dicType==5)
    {
        dic = new CTriangDictionary;
        ((CTriangDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CTriangDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CTriangDictionary*) dic)->getRealAtom();
    }

    else if ( parm->dicType==6)
    {
        dic = new CBatemanDictionary;
        ((CBatemanDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CBatemanDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CBatemanDictionary*) dic)->getRealAtom();
    }

    cgRealAtom.fillVector(realAtom);
    cgRealAtomAux = cgRealAtom*(((strtParameter*)parm)->innerProduct);
        residue = residue - cgRealAtom*(((strtParameter*)parm)->innerProduct);
    norm = residue.norm();
    // cout << "AQUI" << norm << endl;
    if(dic) delete dic;
    return;
}

void updateResidue (    cgMatrix<double>& residue,
                        cgMatrix<double>& signal,
                        cgMatrix<double>& b,
                        cgMatrix<double>& v,
                        cgMatrix<double>& Ai,
                        cgMatrix<double>& a,
                        double& norm,
                        int dicSize,
                        int step,
                        cgMatrix<double>& prevAtoms,
                        strtParameter* parm)
{
    double* u = new double[dicSize];
    double* approx = new double[dicSize];
    double vb = 0.0;
    // double bb = 0.0;
    double beta = 0.0;
    double norm2u = 0.0;
    double* bg = new double[dicSize];
    CDictionary* dic;
    CDictionary* dicAux;
    double* realAtom;
    double* realAtomAux;
    cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    cgMatrix<double> cgApprox(1,dicSize,0.0);
    // cgMatrix<double> bb(step,step,0.0);
    strtParameter* gGammaAux = new strtParameter;

    // for(int i=0;i<dicSize;i++)
    // {
    //     cout << residue[0][i] << endl;
    // }


    for (int i=0; i<step; i++)
    {
        v[i][0] = 0.0;
        b[i][0] = 0.0;
    }

    for (int i=0; i<dicSize; i++)
    {
        u[i] = 0.0;
        approx[i] = 0.0;
        bg[i] = 0.0;
    }

    // cgMatrix<strtParameter*> gGamma(step,1,parm);

    prevAtoms[step][0]  = parm->innerProduct ;
    prevAtoms[step][1]  = parm->s            ;
    prevAtoms[step][2]  = parm->rho          ;
    prevAtoms[step][3]  = parm->xi           ;
    prevAtoms[step][4]  = parm->phase        ;
    prevAtoms[step][5]  = parm->u            ;
    prevAtoms[step][6]  = parm->a            ;
    prevAtoms[step][7]  = parm->b            ;
    prevAtoms[step][8]  = parm->nextAtom     ;
    prevAtoms[step][9]  = parm->prevAtom     ;
    prevAtoms[step][10] = parm->origAtomIndex;
    prevAtoms[step][11] = parm->beta         ;
    prevAtoms[step][12] = parm->dicType      ;

    // prevAtoms[step] = parm;
    // cout << step << endl;
    // cout << prevAtoms[0][1] << endl;
    // cout << prevAtoms[step][1] << endl;
    // for (int i=0;i<=step;i++)
    // {
    //     cout << prevAtoms[i]->rho << endl;
    // }


    // strtParameter* parmAux;

    if ( parm->dicType==1 || parm->dicType==2 || parm->dicType==3 )
    {
        dic = new CExpDictionary;
        ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CExpDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CExpDictionary*) dic)->getRealAtom();
    }

    else if ( parm->dicType==4)
    {
        dic = new CGaborDictionary;
        ((CGaborDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CGaborDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CGaborDictionary*) dic)->getRealAtom();
    }


    else if ( parm->dicType==5)
    {
        dic = new CTriangDictionary;
        ((CTriangDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CTriangDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CTriangDictionary*) dic)->getRealAtom();
    }

    else if ( parm->dicType==6)
    {
        dic = new CBatemanDictionary;
        ((CBatemanDictionary*)dic)-> setSignalSize(residue.getColumns());
        ((CBatemanDictionary*) dic)->setRealAtom(parm);
        realAtom = ((CBatemanDictionary*) dic)->getRealAtom();
    }

    cgRealAtom.fillVector(realAtom);





    // cgRealAtomAux = cgRealAtom*(((strtParameter*)parm)->innerProduct);
    if (step == 0)
    {
        Ai[0][0] = 1;
        a[0][step] = parm->innerProduct;
        cgApprox = cgRealAtom*(((strtParameter*)parm)->innerProduct);
    }
    else
    {
        for (int i=0; i<step; i++)
        {
            gGammaAux = new strtParameter;
            gGammaAux->innerProduct = prevAtoms[i][0];
            gGammaAux->s            = prevAtoms[i][1];
            gGammaAux->rho          = prevAtoms[i][2];
            gGammaAux->xi           = prevAtoms[i][3];
            gGammaAux->phase        = prevAtoms[i][4];
            gGammaAux->u            = prevAtoms[i][5];
            gGammaAux->a            = prevAtoms[i][6];
            gGammaAux->b            = prevAtoms[i][7];
            gGammaAux->nextAtom     = prevAtoms[i][8];
            gGammaAux->prevAtom     = prevAtoms[i][9];
            gGammaAux->origAtomIndex= prevAtoms[i][10];
            gGammaAux->beta         = prevAtoms[i][11];
            gGammaAux->dicType      = prevAtoms[i][12];

            if ( gGammaAux->dicType==1 || gGammaAux->dicType==2 || gGammaAux->dicType==3 )
            {
                dicAux = new CExpDictionary;
                ((CExpDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CExpDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CExpDictionary*) dicAux)->getRealAtom();
            }

            else if ( gGammaAux->dicType==4)
            {
                dicAux = new CGaborDictionary;
                ((CGaborDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CGaborDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CGaborDictionary*) dicAux)->getRealAtom();
            }


            else if ( gGammaAux->dicType==5)
            {
                dicAux = new CTriangDictionary;
                ((CTriangDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CTriangDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CTriangDictionary*) dicAux)->getRealAtom();
            }

            else if ( gGammaAux->dicType==6)
            {
                dicAux = new CBatemanDictionary;
                ((CBatemanDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CBatemanDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CBatemanDictionary*) dicAux)->getRealAtom();
            }

            for (int j=0; j<dicSize; j++)
            {
                v[i][0] += cgRealAtom[0][j]*realAtomAux[j];
            }
            // }

            if (dicAux) delete dicAux;

            // v[i][1] = cgRealAtom.fast_dprod(realAtomAux,0,dicSize);
        }
        for (int i=0; i<step; i++)
        {
            for (int j=0; j<step; j++)
            {
                b[i][0] += Ai[i][j]*v[j][0];
            }
            vb += v[i][0]*b[i][0];
            // bb += b[i][0]*b[i][0];
        }
        // cout << vb << endl;
        // cout << bbt << endl;

        beta = 1/(1-vb);
        // cout << beta << endl;
        // calculate Ainv
        // for (int i=0; i<step; i++)
        // {
        //     for (int j=0; j<step; j++)
        //     {
        //         Ai[i][j] += beta*bb;
        //     }
        // }


        for (int i=0; i<step; i++)
        {
            for (int j=0; j<step; j++)
            {
                Ai[i][j] += beta*b[i][0]*b[j][0];
            }
        }

        for (int i=0; i<step; i++)
        {
            Ai[i][step] = -beta*b[i][0];
        }

        for (int j=0; j<step; j++)
        {
            Ai[step][j] = -beta*b[j][0];
        }

        Ai[step][step] = beta;
        // cout << "AQUI" << endl;

        for (int i=0; i<step; i++)
        {
            gGammaAux->innerProduct = prevAtoms[i][0];
            gGammaAux->s            = prevAtoms[i][1];
            gGammaAux->rho          = prevAtoms[i][2];
            gGammaAux->xi           = prevAtoms[i][3];
            gGammaAux->phase        = prevAtoms[i][4];
            gGammaAux->u            = prevAtoms[i][5];
            gGammaAux->a            = prevAtoms[i][6];
            gGammaAux->b            = prevAtoms[i][7];
            gGammaAux->nextAtom     = prevAtoms[i][8];
            gGammaAux->prevAtom     = prevAtoms[i][9];
            gGammaAux->origAtomIndex= prevAtoms[i][10];
            gGammaAux->beta         = prevAtoms[i][11];
            gGammaAux->dicType      = prevAtoms[i][12];

            if ( gGammaAux->dicType==1 || gGammaAux->dicType==2 || gGammaAux->dicType==3 )
            {
                dicAux = new CExpDictionary;
                ((CExpDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CExpDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CExpDictionary*) dicAux)->getRealAtom();
            }

            else if ( gGammaAux->dicType==4)
            {
                dicAux = new CGaborDictionary;
                ((CGaborDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CGaborDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CGaborDictionary*) dicAux)->getRealAtom();
            }


            else if ( gGammaAux->dicType==5)
            {
                dicAux = new CTriangDictionary;
                ((CTriangDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CTriangDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CTriangDictionary*) dicAux)->getRealAtom();
            }

            else if ( gGammaAux->dicType==6)
            {
                dicAux = new CBatemanDictionary;
                ((CBatemanDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CBatemanDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CBatemanDictionary*) dicAux)->getRealAtom();
            }

            for (int j=0; j<dicSize; j++)
            {
               bg[j] += b[i][0]*realAtomAux[j];
            }

            // for (int j=0; j<=dicSize; j++)
            // {
            //     u[j] = u[j] - bg;
            // }
            // u = u - b[i][1]*realAtomAux;

            if (dicAux) delete dicAux;
            // cout << b[i][0] << endl;
            // cout << v[i][0] << endl;
        }
        for (int i=0; i<dicSize; i++)
        {
           u[i] = realAtom[i] - bg[i];
        }

        // u = u + realAtom;

        for (int i=0; i<dicSize; i++)
        {
            norm2u += pow(u[i],2.0);
        }

        a[0][step] = parm->innerProduct/norm2u;

        for (int i=0; i<step; i++)
        {
            a[0][i] -= b[i][0]*a[0][step];
        }

        for (int i=0; i<=step; i++)
        {
            gGammaAux->innerProduct = prevAtoms[i][0];
            gGammaAux->s            = prevAtoms[i][1];
            gGammaAux->rho          = prevAtoms[i][2];
            gGammaAux->xi           = prevAtoms[i][3];
            gGammaAux->phase        = prevAtoms[i][4];
            gGammaAux->u            = prevAtoms[i][5];
            gGammaAux->a            = prevAtoms[i][6];
            gGammaAux->b            = prevAtoms[i][7];
            gGammaAux->nextAtom     = prevAtoms[i][8];
            gGammaAux->prevAtom     = prevAtoms[i][9];
            gGammaAux->origAtomIndex= prevAtoms[i][10];
            gGammaAux->beta         = prevAtoms[i][11];
            gGammaAux->dicType      = prevAtoms[i][12];

            if ( gGammaAux->dicType==1 || gGammaAux->dicType==2 || gGammaAux->dicType==3 )
            {
                dicAux = new CExpDictionary;
                ((CExpDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CExpDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CExpDictionary*) dicAux)->getRealAtom();
            }

            else if ( gGammaAux->dicType==4)
            {
                dicAux = new CGaborDictionary;
                ((CGaborDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CGaborDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CGaborDictionary*) dicAux)->getRealAtom();
            }


            else if ( gGammaAux->dicType==5)
            {
                dicAux = new CTriangDictionary;
                ((CTriangDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CTriangDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CTriangDictionary*) dicAux)->getRealAtom();
            }

            else if ( gGammaAux->dicType==6)
            {
                dicAux = new CBatemanDictionary;
                ((CBatemanDictionary*)dicAux)-> setSignalSize(residue.getColumns());
                ((CBatemanDictionary*) dicAux)->setRealAtom(gGammaAux);
                realAtomAux = ((CBatemanDictionary*) dicAux)->getRealAtom();
            }

            for (int j=0; j<dicSize; j++)
            {
                // cout << approx[j] << endl;
                approx[j] += a[0][i]*realAtomAux[j];
            }
            // approx = approx + a[1][i]*realAtomAux;
            cgApprox.fillVector(approx,dicSize);
            if (dicAux) delete dicAux;
        }
    }

    // cout << vtb << endl;
    // cout << beta << endl;
    residue = signal - cgApprox;

    // cout << "b = [ ";

    // for (int i=0; i<step; i++)
    // {
    //     cout << b[i][0] << " ";
    // }
    // cout << "]" << endl;

    // cout << "v = [ ";

    // for (int i=0; i<step; i++)
    // {
    //     cout << v[i][0] << " ";
    // }
    // cout << "]" << endl;

    // cout << "beta = " << beta << endl;

    // cout << "Ai = [ ";
    // for (int i=0; i<=step; i++)
    // {
    //     if (i!=0) cout << "       ";
    //     for (int j=0; j<=step; j++)
    //     {
    //         cout << Ai[i][j] << " ";
    //     }

    //     if (i!=step) cout << endl;
    // }
    // cout << "]" << endl;


    // cout << "a = [ ";

    // for (int i=0; i<=step; i++)
    // {
    //     cout << a[0][i] << " ";
    // }
    // cout << "]" << endl;

    // for(int i=0;i<dicSize;i++)
    // {
    //     cout << residue[0][i] << endl;
    // }

    // residue = residue - cgRealAtom*(((strtParameter*)parm)->innerProduct);
    norm = residue.norm();
    // cout << "AQUI" << norm << endl;
    if (dic) delete dic;
    if (u) delete[] u;
    if (approx) delete[] approx;
    if (bg) delete[] bg;
    if (gGammaAux) delete gGammaAux;
    return;
}


// void updateResidue (cgMatrix<double>& residue,
//                     cgMatrix<double>& cgOrthMtx,
//                     double& norm,
//                     int dicSize,
//                     int step,
//                     strtParameter* parm)
// {
//     CDictionary* dic;
//     double* realAtom;
//     cgMatrix<double> cgRealAtom(1,dicSize,0.0);
//     cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);
//     cgMatrix<double> cgOrthRealAtom(1,dicSize,0.0);
//     cgMatrix<double> cgOrthRealAtomAux(1,dicSize,0.0);
//     double aux[dicSize];
//     double* ptr;

//     if ( parm->dicType==1 || parm->dicType==2 || parm->dicType==3 )
//     {
//         dic = new CExpDictionary;
//         ((CExpDictionary*)dic)-> setSignalSize(residue.getColumns());
//         ((CExpDictionary*) dic)->setRealAtom(parm);
//         realAtom = ((CExpDictionary*) dic)->getRealAtom();
//     }

//     else if ( parm->dicType==4)
//     {
//         dic = new CGaborDictionary;
//         ((CGaborDictionary*)dic)-> setSignalSize(residue.getColumns());
//         ((CGaborDictionary*) dic)->setRealAtom(parm);
//         realAtom = ((CGaborDictionary*) dic)->getRealAtom();
//     }


//     else if ( parm->dicType==5)
//     {
//         dic = new CTriangDictionary;
//         ((CTriangDictionary*)dic)-> setSignalSize(residue.getColumns());
//         ((CTriangDictionary*) dic)->setRealAtom(parm);
//         realAtom = ((CTriangDictionary*) dic)->getRealAtom();
//     }

//     cgRealAtom.fillVector(realAtom);
//     cgRealAtomAux = cgRealAtom*(((strtParameter*)parm)->innerProduct);

//     if (step == 0)
//     {
//         cgOrthRealAtom = cgRealAtom;
//     }
//     else
//     {
//         for (int p=0; p<step; p++)
//         {
//             ptr = cgOrthMtx.getDataLinePtr(p);
//             memcpy(aux,ptr,dicSize*sizeof(double));
//             cgOrthRealAtomAux.fillVector(aux);
//             // cout << cgRealAtom.fast_dprod(aux,0,511)[0][0] << endl;
//             cgOrthRealAtom += cgOrthRealAtomAux*(cgRealAtom.fast_dprod(aux,0,511) / pow(cgOrthRealAtomAux.norm(),2.0))[0][0];
//         }

//         cgOrthRealAtom = cgRealAtom - cgOrthRealAtom;
//     }

//     residue = residue - cgRealAtom*(((strtParameter*)parm)->innerProduct);
//     norm = residue.norm();

//     for (int j=0; j<dicSize; j++)
//     {
//         cgOrthMtx.setData(step,j,cgOrthRealAtom(0,j));
//     }

//     // ptr = cgOrthMtx.getDataLinePtr(step);
//     // memcpy(cgOrthMtx.getDataLinePtr(step),&cgOrthRealAtom,dicSize*sizeof(double));

//     // cout << cgOrthRealAtom(0,511) << endl;
//     // cout << cgOrthMtx(0,511) << endl;
//     // cgOrthMtx.fillVector(cgOrthRealAtom,step,step+dicSize);
//     if(dic) delete dic;
//     return;
// }

/*void ANN (  CFileDecomp* genData,
            cgMatrix<double>& residue,
            int dicSize,
            int& chosenDic,
            int& chosenNet,
            int step,
            int L,
            double* approxRatio)
{
    double yaux;
    double yaux2[4];
    double resAux[(int)ceil(dicSize/8.0)];
    // mclmcrInitialize();
    // const char *args[] = { "-nodisplay" };
    // mclInitializeApplication(args,1);
    // libannformpInitialize();
    mwArray x((int)ceil(dicSize/8.0), 1, mxDOUBLE_CLASS);
    mwArray y(4, 1, mxDOUBLE_CLASS);

    for (int iter = 0; iter < ceil(dicSize/8.0); ++iter)
    {
        resAux[iter] = residue.getData(0,iter*8)/residue.norm();
        cout << resAux[iter] << endl;
    }
    x.SetData(resAux,ceil(dicSize/8.0));

    if ( step == 0 || approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= -1e-1 )
    {
        if ( chosenNet == 1  || step == 0)
        {
            cout << "AQUI1" << endl;
            chosenNet = 1;
        }
    }

    else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= -5e-2 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -1e-1 )
    {
        if ( chosenNet == 1 || chosenNet == 2 || step == 0)
        {
            cout << "AQUI2" << endl;
            chosenNet = 2;
        }
    }

    // else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= -2e-2 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -5e-2 )
    // {
    //     if ( chosenNet == 1 || chosenNet == 2 || chosenNet == 3 || step == 0)
    //     {
    //         cout << "AQUI3" << endl;
    //         chosenNet = 3;
    //     }
    // }

    // else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= 0 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -2e-2 )
    // {
    //     cout << "AQUI4" << endl;
    //     chosenNet = 0;
    //     chosenDic = genData->getDicType();
    // }

    else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= 0 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -5e-2 )
    {
        if ( chosenNet == 1 || chosenNet == 2 || chosenNet == 3 || step == 0)
        {
            cout << "AQUI3" << endl;
            chosenNet = 3;
        }
    }

    if (chosenNet != 0)
    {
        if (chosenNet == 1)
        {
            net1mat(1,y,x);
            cout << "AQUI5" << endl;
        }

        else if (chosenNet == 2)
        {
            net2mat(1,y,x);
            cout << "AQUI6" << endl;
        }

        else if (chosenNet == 3)
        {
            net3mat(1,y,x);
            cout << "AQUI7" << endl;
        }

        yaux = 0.0;

        y.GetData(yaux2,3);

        for (int d = 0; d < 4; ++d)
        {
            cout << yaux2[d] << endl;
            if (yaux2[d]>yaux)
            {
                chosenDic = d+1;
                yaux = yaux2[d];
            }
        }

        if (chosenDic == 3)
        {
            chosenDic = 4; //Gabor dictionary
        }
    }

    // cout << chosenDic << endl;
    // cout << y << endl;

    // if (yaux2) delete[] yaux2;
    // if (resAux) delete[] resAux;
    return;
}*/

void DANNO (  CFileDecomp* genData,
            cgMatrix<double>& residue,
            int dicSize,
            int& chosenDic,
            int& chosenNet,
            int step,
            int L,
            double* approxRatio)
{
    // double yaux;
    // double yaux2[4];
    // double resAux[(int)ceil(dicSize/8.0)];

    double yaux = 0.0;
    double* y = new double[3];
    double* res = new double[(int)ceil(dicSize/8.0)];

    CNet* net = new CNet;

    net->setSignalSize(dicSize);
    // mclmcrInitialize();
    // const char *args[] = { "-nodisplay" };
    // mclInitializeApplication(args,1);
    // libannformpInitialize();
    // mwArray x((int)ceil(dicSize/8.0), 1, mxDOUBLE_CLASS);
    // mwArray y(4, 1, mxDOUBLE_CLASS);

    // for (int iter = 0; iter < ceil(dicSize/8.0); ++iter)
    // {
    //     resAux[iter] = residue.getData(0,iter*8)/residue.norm();
    // }

    for (int iter = 0; iter < ceil(dicSize/8.0); ++iter)
    {
        res[iter] = residue.getData(0,iter*8)/residue.norm();
        // cout << res[iter] << endl;
    }

    net->setX(res);

    // x.SetData(resAux,ceil(dicSize/8.0));

    if ( step == 0 || approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= -1e-1 )
    {
        if ( chosenNet == 1  || step == 0)
        {
            cout << "AQUI1" << endl;
            chosenNet = 1;
        }
    }

    else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= -5e-2 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -1e-1 )
    {
        if ( chosenNet == 1 || chosenNet == 2 || step == 0)
        {
            cout << "AQUI2" << endl;
            chosenNet = 2;
        }
    }

    // else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= -2e-2 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -5e-2 )
    // {
    //     if ( chosenNet == 1 || chosenNet == 2 || chosenNet == 3 || step == 0)
    //     {
    //         cout << "AQUI3" << endl;
    //         chosenNet = 3;
    //     }
    // }

    // else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= 0 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -2e-2 )
    // {
    //     cout << "AQUI4" << endl;
    //     chosenNet = 0;
    //     chosenDic = genData->getDicType();
    // }

    else if ( approxRatio[(step%L)-1]-approxRatio[(step%L)-2] <= 0 && approxRatio[(step%L)-1]-approxRatio[(step%L)-2] > -5e-2 )
    {
        if ( chosenNet == 1 || chosenNet == 2 || chosenNet == 3 || step == 0)
        {
            cout << "AQUI3" << endl;
            chosenNet = 3;
        }
    }

    if (chosenNet != 0)
    {

        // if (chosenNet == 1)
        // {
        //     net1mat(1,y,x);
        //     cout << "AQUI5" << endl;
        // }

        // else if (chosenNet == 2)
        // {
        //     net2mat(1,y,x);
        //     cout << "AQUI6" << endl;
        // }

        // else if (chosenNet == 3)
        // {
        //     net3mat(1,y,x);
        //     cout << "AQUI7" << endl;
        // }

        // y.GetData(yaux2,3);
        net->setNumber(chosenNet);
        net->net();
        y = net->getY();

        for (int d = 0; d < 3; ++d)
        {
            if (y[d]>yaux)
            {
                chosenDic = d+1;
                yaux = y[d];
            }
        }

        if (chosenDic == 3)
        {
            chosenDic = 4; //Gabor dictionary
        }
    }

    // cout << chosenDic << endl;
    // cout << y << endl;

    // if (yaux2) delete[] yaux2;
    // if (resAux) delete[] resAux;
    return;
}

// int checkPrevAtom(  int step,
//                     cgMatrix[double] prevAtoms,
//                     s,
//                     xi,
//                     tau,
//                     dicType)
// {
//     for (int j=0; j<step)
//     {
//         if prevAtoms[step][1]
//     }
// }

// int checkPrevAtom(  int step,
//                     cgMatrix[double] prevAtoms,
//                     s,
//                     xi,
//                     tau,
//                     dicType)
// {
//     for (int j=0; j<step)
//     {
//         if prevAtoms[step][1]
//     }
// }

void matchingPursuitEDA(cgMatrix<double>& residue,
                        strtParameter* chosenParm,
                        int dicSize,
                        CDataSignal* dataSignal,
                        CFileDictionary* dicData,
                        CFileDecomp* genData, int step, int chosenDic,
                        strtParameter** dicAtoms,
                        int& a0,
                        int& b0,
                        int flagOMP)
{
    int decincAsymmFlag;
    int i,j;
    int k;
    // int s;
    double s;
    int a;
    int chosenTau = 0;
    double chosenXi = 0.0;
    int delta_tau;
    double rise;
    double decay;
    double chosenOptPhase = 0.0;
    double Fs;
    double beta = 0.0;
    double rho = 0.0;
    // double Ffund;
    // double delta_f = 0;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    int dicType;
    // int signalSize;
    double  innerProd;

    double* realAtom = new double[dicSize];
    double* realAtomDec = new double [dicSize];
    double* realAtomInc = new double [dicSize];
    double* realAtomWin = new double[2*dicSize];

    CComplex* complexAtom = new CComplex[dicSize];
    CComplex* complexAtomDec = new CComplex[dicSize];
    CComplex* complexAtomInc = new CComplex[dicSize];
    CComplex* complexAtomXi = new CComplex[dicSize];
    CComplex* complexAtom2Xi = new CComplex[dicSize];


    double* conv_zxi0_w2 = new double[2*dicSize];
    double* conv_zxi0_w2_expinc = new double[2*dicSize];

    CDictionary *expDic = new CExpDictionary;
    CDictionary * gaborDic = new CGaborDictionary;
    CDictionary * triangDic = new CTriangDictionary;
    CDictionary * triangDic2 = new CTriangDictionary;
    CDictionary * bateDic = new CBatemanDictionary;


    ((CExpDictionary*)expDic)->setSignalSize(dicSize);
    ((CGaborDictionary*)gaborDic)->setSignalSize(dicSize);
    ((CTriangDictionary*)triangDic)->setSignalSize(dicSize);
    ((CTriangDictionary*)triangDic2)->setSignalSize(2*dicSize);
    ((CBatemanDictionary*)bateDic)->setSignalSize(dicSize);


    strtParameter* parm = new strtParameter;
    ((strtParameter*)chosenParm)->innerProduct=0.0;
    ((strtParameter*)chosenParm)->rho = 0.0;
    ((strtParameter*)chosenParm)->xi = 0.0;
    ((strtParameter*)chosenParm)->phase = 0.0;
    ((strtParameter*)chosenParm)->a = 0;
    ((strtParameter*)chosenParm)->b = 0;
    //((strtParameter*)chosenParm)->beta = 0;

    double xi = 0.0;
    int tau = 0;
    double maxInnerProd = 0.0;

    int d = 0;

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

     //if (MPType == 3)
     //{
       //  fileName = "fastMPKolasaModified.out";
     //}

    for(j=0;j<dicData->getNumDicBlock();j++)
    {
        // stream = fopen( fileName, "a" );
        // fprintf(stream,"%d --------------------------------------------------------------\n", step);
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
        rise = (dicData->getDicBlock())[j].rise;
        decay = (dicData->getDicBlock())[j].decay;
        int N = dicSize;

        //if(dicType==4)
        //    MPType=2;
        //else if (dicType==5)
          //  MPType=2;
        //else if(dicType==6)
          //  MPType=3;

        for (int n=0; n<dicSize; n++)
        {
            realAtom[n] = 0.0;
            realAtomDec[n] = 0.0;
            realAtomInc[n] = 0.0;
        }

        for (int n=0; n<2*dicSize; n++)
        {
            realAtomWin[n] = 0.0;
            conv_zxi0_w2[n] = 0.0;
            conv_zxi0_w2_expinc[n] = 0.0;
        }


        if (chosenDic == dicType || chosenDic == 0)
        {
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
            //else if (dataSignal->getType() == 4) // ECG
            //{
            //    Fs = ((CECGSignal*)dataSignal)->getSamplingRate();
            //}


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
                else if (dataSignal->getType() == 4) // ECG
                {
                    freqi = 0;
                    freqf = Fs;
                }

            }
            else if (freqf==9999999999 && freqi==0000000000)
            {
                freqf=Fs;
                freqi= freqf/2;
            }
            else if (freqf==9999999999)
            {
                freqf = Fs;
            }
            else if (freqi==9999999999)
            {
                freqi = freqf/s;
            }
            else if (freqi==0000000000)
            {
                freqi = freqf/2;
            }

            if (freqf>Fs)
            {
                cout << "Final frequency greater than sampling frequency !!!" << endl;
                printf("final freq. %f - Fs %f\n",freqf,Fs);
                exit(1);
            }

            if (fdiscrtype==1) // linear
            {
                // delta_f = (2*PI/freqf)* freqi;
                Nfreq = (int)(freqf/(2*freqi));
                // cout << Nfreq << endl;
                //Nfreq = (int)ceil(freqf/freqi);
                xi_vec = new double[Nfreq];
                for (i=0;i<Nfreq;i++)
                {
                    // xi_vec[i] = delta_f * i;
                    xi_vec[i] = (2*PI/Fs) * (freqi * i );
                    // cout << xi_vec[i] << endl;
                }
            }
            else if (fdiscrtype==2) // geometric with quarter-tone discretization
            {
                Nfreq = (int)ceil(24 * ( log10(freqf/freqi)/log10(2.0) ) )+1;

                xi_vec = new double[Nfreq];
                xi_vec[0] = 0.0;
                for (i=1;i<Nfreq;i++)
                {
                    xi_vec[i] = (2*PI/Fs) * (freqi * pow (2.0, (double)(i-1)/24) );
                }
            }

            if (dicType == 1) //exponential
            {

                if (MPType == 1)
                {
                    for(k=0;k<Nfreq;k++)
                    {
                        for(tau=0;tau<(double)N;tau+=delta_tau)
                        {
                            // cout << tau << endl;

                            //Decreasing
                            //xi = k * ((2*PI)/s);
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = tau;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            complexAtomDec = expDic->getComplexAtom();

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtomDec, fileName,s,tau,N);


                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                tau,
                                                N-1,
                                                0.0);
                            }
                            d++;
                        }
                        for (tau=N-1;tau>0; tau=tau-delta_tau)
                        {
                            // Increasing
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = -1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = tau;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            complexAtomInc = expDic->getComplexAtom();

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtomInc,fileName,s,tau,N);

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                -s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                0,
                                                tau,
                                                0.0);
                            }
                            d++;
                        }
                    }
                }

                else if (MPType==2)
                {
                    // int delta_tau = s;
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        //decreasing
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = tau;
                        ((strtParameter*)parm)->a = tau;
                        ((strtParameter*)parm)->b = N-1;
                        ((CExpDictionary*)expDic)->setRealAtom(parm);
                        realAtomDec = expDic->getRealAtom();
                        fastMPKolasa(   residue,
                                        maxInnerProd,
                                        chosenOptPhase,
                                        chosenXi,
                                        dicSize,
                                        tau,
                                        s,
                                        realAtomDec,
                                        fileName);

                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            tau,
                                            N-1,
                                            0.0);
                            // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                            // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                            // ((strtParameter*)chosenParm)->xi = chosenXi;
                            // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                            // ((strtParameter*)chosenParm)->a = tau;
                            // ((strtParameter*)chosenParm)->b = N-1;
                        }
                    }

                    for (tau=1; tau<=(double)N; tau=tau+delta_tau)
                    {
                        // increainsg
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = -1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->a = 0.0;
                        ((strtParameter*)parm)->b = tau;
                        ((strtParameter*)parm)->u = tau;
                        ((CExpDictionary*)expDic)->setRealAtom(parm);
                        realAtomInc = expDic->getRealAtom();

                        fastMPKolasa(   residue,
                                        maxInnerProd,
                                        chosenOptPhase,
                                        chosenXi,
                                        dicSize,
                                        tau,
                                        s,
                                        realAtomInc,
                                        fileName);
                        // cout << chosenXi << endl;
                        // cout << tau << endl;
                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << s << endl;
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            -s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            0,
                                            tau,
                                            0.0);
                            // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                            // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                            // ((strtParameter*)chosenParm)->xi = chosenXi;
                            // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                            // ((strtParameter*)chosenParm)->a = 0.0;
                            // ((strtParameter*)chosenParm)->b = tau;
                        }
                    }
                }

                else if (MPType==3)
                {
                    // cout << "AQUI" << endl;
                    ((strtParameter*)parm)->s = s;
                    ((strtParameter*)parm)->rho = 1.0/s;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = 0;
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)expDic)->setRealAtom(parm);
                    // realAtomWin = expDic->getRealAtom();
                    // memcpy(realAtomWin, (expDic->getRealAtom()), sizeof(double)*dicSize);
                    for (int k1=0; k1<dicSize;k1++)
                    {
                        realAtomWin[k1+N] = *(expDic->getRealAtom()+k1);
                    }

                    for(k=0;k<Nfreq;k++)
                    {
                        //xi = k * ((2*PI)/s);
                        xi = xi_vec[k];
                        if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                        {
                            // z1
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            // complexAtomXi = expDic->getComplexAtom();
                            // memcpy(complexAtomXi, (expDic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                            }

                            // z2
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = 2*xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            // complexAtom2Xi = expDic->getComplexAtom();
                            // memcpy(complexAtom2Xi, (expDic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                            }


                            //Decreasing exponential
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
                                                    conv_zxi0_w2,
                                                    fileName,
                                                    s);
                            // if (k == 0){
                            // }

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                chosenTau,
                                                N-1,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = xi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = chosenTau;
                                // ((strtParameter*)chosenParm)->b = N-1;
                            }

                            //Increasing exponential
                            decincAsymmFlag = -1;
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
                                                    conv_zxi0_w2_expinc,
                                                    fileName,
                                                    s);
                            // // if (k == 0){
                            // // }

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                -s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                0,
                                                chosenTau,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = -1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = xi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = 0;
                                // ((strtParameter*)chosenParm)->b = chosenTau;
                            }
                        }
                    }
                }
            }

            else if (dicType == 2) // pure sine
            {
                ((CExpDictionary*)expDic)->setSignalSize(dicSize);

                for (tau=0; tau<(double)N; tau=tau+delta_tau)
                {
                    ((strtParameter*)parm)->rho = 0.0;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = tau;
                    ((strtParameter*)parm)->a = tau;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)expDic)->setRealAtom(parm);
                    realAtom = expDic->getRealAtom();
                    fastMPKolasa(   residue,
                                    maxInnerProd,
                                    chosenOptPhase,
                                    chosenXi,
                                    dicSize,
                                    tau,
                                    s,
                                    realAtom,
                                    fileName);
                    if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    {
                        // cout << maxInnerProd << endl;
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        chosenXi,
                                        chosenOptPhase,
                                        tau,
                                        tau,
                                        N-1,
                                        0.0);
                    }
                    d++;
                }
            }

            else if (dicType == 3) // impulse
            {
                for(i=0;i<N;i++)
                {
                    innerProd = residue[0][i];
                    if (fabs(innerProd)>fabs(maxInnerProd))
                    {
                        // cout << maxInnerProd << endl;
                        maxInnerProd = innerProd;
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        0.0,
                                        0.0,
                                        i,
                                        i,
                                        i,
                                        0.0);
                        // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                        // ((strtParameter*)chosenParm)->xi = 0.0;
                        // ((strtParameter*)chosenParm)->phase = 0.0;
                        // ((strtParameter*)chosenParm)->a = i;
                        // ((strtParameter*)chosenParm)->b = i;
                    }
                    d++;
                }
            }

            else if (dicType == 4) // Gabor
            {


                if (MPType == 1)
                {
                    for(k=0;k<Nfreq;k++)
                    {
                        for(tau=0;tau<N;tau+=delta_tau)
                        {
                            //xi = k * ((2*PI)/s);
                            // cout << s << endl;
                            // cout << tau << endl;
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CGaborDictionary*)gaborDic)->setComplexAtom(parm);
                            complexAtom = gaborDic->getComplexAtom();


                            // MPTradicional(  residue,
                            //                 complexAtomDec,
                            //                 maxInnerProd,
                            //                 chosenOptPhase,
                            //                 tau,
                            //                 xi,
                            //                 N, //int dicSize
                            //                 s,
                            //                 fileName);

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                            // cout << maxInnerProd << endl;

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                0,
                                                N-1,
                                                0.0);
                            }
                            d++;
                        }
                    }
                }

                else if (MPType==2)
                {
                    int delta_tau = (int)s;
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = tau;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CGaborDictionary*)gaborDic)->setRealAtom(parm);
                        realAtom = gaborDic->getRealAtom();
                        fastMPKolasa(   residue,
                                        maxInnerProd,
                                        chosenOptPhase,
                                        chosenXi,
                                        dicSize,
                                        tau,
                                        s,
                                        realAtom,
                                        fileName);

                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            0,
                                            N-1,
                                            0.0);
                        }
                    }
                }

                else if (MPType==3)
                {


                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/((double)s);
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;//0
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CGaborDictionary*)gaborDic)->setRealAtom(parm);
                        // memcpy(realAtomWin, (gaborDic->getRealAtom()), sizeof(double)*dicSize);
                        // cout << "realAtomWin" << endl;

                        /*
                        realAtomWin[0] = *(gaborDic->getRealAtom());
                        for (int k1=1; k1<dicSize;k1++)
                        {
                            //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1);
                            realAtomWin[k1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                            realAtomWin[(2*N)-k1] = *(gaborDic->getRealAtom()+k1);
                            //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                            //realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);
                            cout<<realAtomWin[k1]<< endl;
                        }
                        */
                        for (int k1=0; k1<dicSize;k1++)
                        {
                            //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1);
                            realAtomWin[k1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                            realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);
                            //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                            //realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);

                        }

                         //cout << dicSize << endl;
                        // for (int k1=0; k1<2*dicSize; k1++)
                        // {
                        //     cout << realAtomWin[k1] << endl;
                        // }
                        // cout << "ATOM" << endl;
                        // for (int k1=0;k1<2*dicSize;k1++)
                        // {
                        //     cout << realAtomWin[k1] << endl;
                        // }

                        for(k=0;k<Nfreq;k++)
                        {
                            xi = xi_vec[k];
                            if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                            {
                                // cout << s << endl;
                                // cout << xi << endl;
                                // z1
                                ((strtParameter*)parm)->s = s;
                                ((strtParameter*)parm)->rho = 0.0; //generating complex exponential with rectangular window
                                ((strtParameter*)parm)->xi = xi;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = 0;
                                ((strtParameter*)parm)->a = 0;
                                ((strtParameter*)parm)->b = N-1;
                                ((CGaborDictionary*)gaborDic)->setComplexAtom(parm);
                                // cout << "complexAtomXi" << endl;
                                // memcpy(complexAtomXi, (gaborDic->getComplexAtom()), sizeof(CComplex)*(dicSize));

                                for (int k1=0; k1<dicSize;k1++)
                                {
                                    complexAtomXi[k1] = *(gaborDic->getComplexAtom()+k1);
                                    // cout << complexAtomXi[K1] << endl;
                                }

                                /*
                                complexAtomXi[0] = *(gaborDic->getComplexAtom());
                                for (int k1=1; k1<dicSize;k1++)
                                {
                                    //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1);
                                    complexAtomXi[k1] = *(gaborDic->getComplexAtom()+k1); //add 1 to centralize de max
                                    complexAtomXi[(2*N)-k1] = *(gaborDic->getComplexAtom()+k1);
                                    //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                                    //realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);

                                }
                                */

                                // cout << "AQUI" << endl;
                                // for (int k1=0; k1<dicSize; k1++)
                                // {
                                //     cout << complexAtomXi[k1].Real() << " " << complexAtomXi[k1].Imag() << endl;
                                // }


                                // z2
                                ((strtParameter*)parm)->s = s;
                                ((strtParameter*)parm)->rho = 0.0;  //generating complex exponential with rectangular window
                                ((strtParameter*)parm)->xi = 2*xi;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = 0;
                                ((strtParameter*)parm)->a = 0;
                                ((strtParameter*)parm)->b = N-1;
                                ((CGaborDictionary*)gaborDic)->setComplexAtom(parm);
                                // memcpy(complexAtom2Xi, (gaborDic->getComplexAtom()), sizeof(CComplex)*(dicSize));
                                // cout << "complexAtom2Xi" << endl;

                                for (int k1=0; k1<dicSize;k1++)
                                {
                                    complexAtom2Xi[k1] = *(gaborDic->getComplexAtom()+k1);
                                    // cout << complexAtom2Xi[K1] << endl;
                                }

                                /*
                                complexAtom2Xi[0] = *(gaborDic->getComplexAtom());
                                for (int k1=1; k1<dicSize;k1++)
                                {
                                    //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1);
                                    complexAtom2Xi[k1] = *(gaborDic->getComplexAtom()+k1); //add 1 to centralize de max
                                    complexAtom2Xi[(2*N)-k1] = *(gaborDic->getComplexAtom()+k1);
                                    //realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                                    //realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);

                                }
                                */
                                fastMPKolasaModified(   residue,
                                                        maxInnerProd,
                                                        chosenOptPhase,
                                                        chosenTau,
                                                        dicSize,
                                                        1,
                                                        delta_tau,
                                                        xi,
                                                        realAtomWin,
                                                        complexAtomXi,
                                                        complexAtom2Xi,
                                                        conv_zxi0_w2,
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
                                                    0,
                                                    N-1,
                                                    0.0);
                                }
                                // cout << 1/sg << endl;
                            }
                        }

                }
            }

            else if (dicType == 5) // triangular
            {
                if (MPType == 1 || MPType == 2)
                {
                    for(tau=0;tau<N;tau+=delta_tau)
                    {
                        // cout << "a: " << a << endl;
                        for(a=0;a<s;a+=delta_tau)
                        {
                            // cout << "tau: " << tau << endl;
                            //Decreasing
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = a;
                            ((strtParameter*)parm)->b = (int)s-1-a;
                            ((CTriangDictionary*)triangDic)->setComplexAtom(parm);
                            complexAtom = triangDic->getComplexAtom();

                            // cout << "AQUI" << endl;
                            // cout <<" s- " << s << "tau- " << tau << "a- " << a << endl;

                            // for (int k1=0;k1<N;k1++)
                            // {
                            //     // complexAtomXi[k1] = *(triangDic->getComplexAtom()+k1);
                            //     // cout << residue[0][k1] << endl;
                            //     cout << complexAtom[k1].Real() << endl;
                            // }

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                            // cout << maxInnerProd << endl;
                            // cout << "maxInnerProd " << maxInnerProd << endl;
                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                // cout << "s- " << s << "tau- " << tau << "a- " << a << endl;
                                // cout << maxInnerProd << endl;
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                a,
                                                (int)s-a-1,
                                                0.0);
                            }
                            d++;

                        }
                    }
                }

                else if (MPType==3)
                {

                    for(a=1;a<=((double)s/5.0);a+=1)
                    {
                        if (a>=((double)s/8.0))
                        {
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = a;//N;
                            ((strtParameter*)parm)->a = a;
                            ((strtParameter*)parm)->b = s-a-1;
                            ((CTriangDictionary*)triangDic2)->setRealAtom(parm);

                            // cout << "s " << s << " - a " << a << endl;
                            for (int k1=0; k1<2*dicSize;k1++)
                            {
                                realAtomWin[k1] = *(triangDic2->getRealAtom()+k1);
                                // cout << realAtomWin[k1] << endl;
                                //realAtomWin[(2*N-1)-k1] = *(triangDic2->getRealAtom()+k1);
                            }
                            // z1
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CTriangDictionary*)expDic)->setComplexAtom(parm);

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                            }

                            // z2
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CTriangDictionary*)expDic)->setComplexAtom(parm);

                            // cout << "AQUI1" << endl;
                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                                // cout << complexAtom2Xi << endl;
                            }
                            // cout << "AQUI2" << endl;

                            fastMPKolasaModified(   residue,
                                                    maxInnerProd,
                                                    chosenOptPhase,
                                                    chosenTau,
                                                    dicSize,
                                                    1,
                                                    delta_tau,
                                                    xi,
                                                    realAtomWin,
                                                    complexAtomXi,
                                                    complexAtom2Xi,
                                                    conv_zxi0_w2,
                                                    fileName,
                                                    s);

                            //computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);


                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                    setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            xi,
                                            chosenOptPhase,
                                            chosenTau,
                                            a,
                                            s-a-1,
                                            0.0);
                            }
                        }
                    }

                }
            }

            if (dicType == 6) //Bateman
            {

                if (MPType == 1)
                {
                    beta=rise;
                    //for (beta = 2.0; beta > 1.0/s; beta/=2.0)
                    //{
                        for(k=0;k<Nfreq;k++)
                        {
                            for(tau=0;tau<(double)N;tau+=delta_tau)
                            {
                                xi = xi_vec[k];
                                ((strtParameter*)parm)->rho = 1.0/s;
                                ((strtParameter*)parm)->xi = xi;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = tau;
                                ((strtParameter*)parm)->a = tau;
                                ((strtParameter*)parm)->b = N-1;
                                ((strtParameter*)parm)->beta = beta;
                                ((CBatemanDictionary*)bateDic)->setComplexAtom(parm);
                                complexAtom = bateDic->getComplexAtom();

                                //computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);

                                if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                                {
                                    setParameters ( chosenParm,
                                                    dicType,
                                                    maxInnerProd,
                                                    s,
                                                    xi,
                                                    chosenOptPhase,
                                                    tau,
                                                    tau,
                                                    N-1,
                                                    beta);
                                }
                                d++;
                            }
                        //}
                    }
                }

                else if (MPType==2)
                {
                    // int delta_tau = s;
                    //for (beta=1.0; beta > 1.0/s; beta/=(2.0))
                    //{
                    beta=rise;
                    rho=decay;
                        for (tau=0; tau<(double)N; tau=tau+delta_tau)
                        {
                            ((strtParameter*)parm)->rho = rho;
                            ((strtParameter*)parm)->beta = beta;
                            //((strtParameter*)parm)->rho = 1.0/s;
                            //((strtParameter*)parm)->beta = beta;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = tau; //tau
                            ((strtParameter*)parm)->b = N-1;
                            ((CBatemanDictionary*)bateDic)->setRealAtom(parm);
                            realAtom = bateDic->getRealAtom();
                            fastMPKolasa(   residue,
                                            maxInnerProd,
                                            chosenOptPhase,
                                            chosenXi,
                                            dicSize,
                                            tau,
                                            s,
                                            realAtom,
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
                                                tau,
                                                N-1,
                                                beta);

                            }
                        }
                    //}
                }

                else if (MPType==3)
                {
                    //int a=0;
                    beta=rise;
                    rho = decay;
                    ///int j = 0;
                    // for (rho =(0.7/4.0); rho <= (1.1/4.0) ; rho+=(0.1/4.0))
                    ///for (rho =0.025; rho <= 0.5 ; rho=(((1e4*rho)+250)/1e4))
                    ///{
                        // cout << "rho - " << xi << endl;
                         //cout << "k1 - " << delta_tau << endl;
                       /// for (beta=0.5; beta > rho; beta=(((1e4*beta)-500)/1e4))
                        ///{
                        // cout << "   beta - " << beta << endl;
                            // cout << int(1e3*beta) << " " << int(1e3*rho) << endl;
                            ///j = j+1;
                    // for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    //{
                            ((strtParameter*)parm)->rho = rho;
                            ((strtParameter*)parm)->beta = beta;
                            ((strtParameter*)parm)->xi = 0; //0
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0; //0
                            ((strtParameter*)parm)->a = 0; //0
                            ((strtParameter*)parm)->b = N-1;
                            ((CBatemanDictionary*)bateDic)->setRealAtom(parm);

                            //((CBatemanDictionary*)expDic)->setRealAtom(parm);

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                //realAtomWin[k1] = *(bateDic->getRealAtom()+k1); //add 1 to centralize de max
                                realAtomWin[k1] = *(bateDic->getRealAtom()+k1);
                                //realAtomWin[k1+N] = *(expDic->getRealAtom()+k1);
                                //realAtomWin[k1] = *(bateDic->getRealAtom()+k1); //realAtomWin[k1+N] = *(bateDic->getRealAtom()+k1);
                                //cout << "FES-  " <<  realAtomWin[k1] <<endl;//"   RAW-  "<< realAtomWin[k1+N] << endl;

                            }

                            for(k=0;k<Nfreq;k++)
                            {
                                //xi = k * ((2*PI)/s);
                                xi=xi_vec[k];
                                //cout << "x1-" << xi << endl;
                                if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                                {
                                   // cout << "OI" << xi << endl;
                                    // z1 
                                    ((strtParameter*)parm)->rho = rho;
                                    ((strtParameter*)parm)->beta = beta;
                                    ((strtParameter*)parm)->xi = xi;
                                    ((strtParameter*)parm)->phase = 0.0;
                                    ((strtParameter*)parm)->u = 0; //0
                                    ((strtParameter*)parm)->a = 0;//0
                                    ((strtParameter*)parm)->b = N-1;
                                    //((CBatemanDictionary*)expDic)->setComplexAtom(parm);
                                    ((CBatemanDictionary*)bateDic)->setComplexAtom(parm);

                                    for (int k1=0; k1<dicSize;k1++)
                                    {

                                        //complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                                        complexAtomXi[k1] = *(bateDic->getComplexAtom()+k1); //complexAtomXi[k1] = *(bateDic->getComplexAtom()+k1);
                                       // cout << "HeY" <<  complexAtomXi[k1].Real() << endl;
                                        //cout << "HoY" <<  complexAtomXi[k1].Imag() << endl;
                                    }

                                    // z2
                                    ((strtParameter*)parm)->rho = rho;
                                    ((strtParameter*)parm)->beta =beta;
                                    ((strtParameter*)parm)->xi = 2*xi;
                                    ((strtParameter*)parm)->phase = 0.0;
                                    ((strtParameter*)parm)->u = 0; //0
                                    ((strtParameter*)parm)->a = 0;//tau
                                    ((strtParameter*)parm)->b = N-1;
                                    ((CBatemanDictionary*)bateDic)->setComplexAtom(parm);
                                    //((CBatemanDictionary*)expDic)->setComplexAtom(parm);

                                    for (int k1=0; k1<dicSize;k1++)
                                    {


                                        //complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                                        complexAtom2Xi[k1] = *(bateDic->getComplexAtom()+k1);   //complexAtom2Xi[k1] = *(bateDic->getComplexAtom()+k1);
                                        //cout << "HOY" << complexAtom2Xi[k1].Real() << endl;



                                    }

                                      //    FILE* stream;
                                         //                           stream = fopen("BATCA.out","a");
                                        // for(int k1=0;k1<dicSize;k1++)
                                        // {
                                         //fprintf (   stream," CAR- %15.8f CAR2 - %15.8f \n",
                                          //       complexAtomXi[k1].Real(),
                                           //      complexAtom2Xi[k1].Real());
                                            //}
                                         //fflush(stream);

                                         //fclose(stream);



                                    //((CBatemanDictionary*)bateDic)->optimizeContinuousParms(residue, chosenParm);

                                    fastMPKolasaModified(   residue,
                                                            maxInnerProd,
                                                            chosenOptPhase,
                                                            chosenTau,
                                                            dicSize,
                                                            1,
                                                            delta_tau,
                                                            xi,
                                                            realAtomWin,
                                                            complexAtomXi,
                                                            complexAtom2Xi,
                                                            conv_zxi0_w2,
                                                            fileName,
                                                            (1.0/rho));

                                //cout << "AQUI" << fabs(((strtParameter*)chosenParm)->innerProduct) << endl;

                                    //((CBatemanDictionary*)bateDic)->optimizeContinuousParms(residue, chosenParm);

                                    //computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);


                                    //cout << "HERE" << fabs(((strtParameter*)chosenParm)->innerProduct) << endl;
                                    //cout << "AQUI" << fabs(maxInnerProd) << endl;
                                    //cout<<"OI"<<fabs(((strtParameter*)chosenParm)->innerProduct)<<endl;
                                    if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                                    {
                                        setParameters ( chosenParm,
                                                        dicType,
                                                        maxInnerProd,
                                                        (1.0/rho),
                                                        xi,
                                                        chosenOptPhase,
                                                        chosenTau,
                                                        chosenTau,
                                                        N-1,
                                                        beta);
                                         //((CBatemanDictionary*)bateDic)->optimizeContinuousParms(residue, chosenParm);

                                    //cout << "MIP" << maxInnerProd<< endl;
                                    //cout << "TAU:" << delta_tau << endl;
                                    //cout << "CTAU" <<chosenTau<< endl;
                                    //cout << "HEY" << innerProduct<< endl;



                                    }
                                    //((CBatemanDictionary*)bateDic)->optimizeContinuousParms(residue, chosenParm);


                           // }
                        }
                    }
                }
            }

            if (xi_vec) delete [] xi_vec;
        }
    }

    if (parm) delete parm;
    if (triangDic) delete triangDic;
    if (triangDic2) delete triangDic2;
    if (gaborDic) delete gaborDic;
    if (bateDic) delete bateDic;
    if (expDic) delete expDic;
    if (conv_zxi0_w2_expinc) delete[] conv_zxi0_w2_expinc;
    if (conv_zxi0_w2) delete[] conv_zxi0_w2;
    return;
}





// void CTriangDictionary::optimizeContinuousParms(   cgMatrix<double>& residue,
//                                 strtParameter* parm);
// {

// }

// // ---- Comtrade heuristics

// void CTriangDictionary::proceedHeuristics( cgMatrix<double>& residue,
//                         strtParameter* parm,
//                         CDataSignal* dataSignal);
// {

// }

// // with rho and xi optimization

// void CTriangDictionary::findFastBestTimeSupport(   cgMatrix<double>& residue,
//                                 strtParameter* parm);
// {

// }

// // with rho optimization

// void CTriangDictionary::findFastBestTimeSupport(   cgMatrix<double>& residue,
//                                 strtParameter* parm,
//                                 int dummy,
//                                 double coefHeur);
// {

// }

// // with rho optimization

// void CTriangDictionary::findBestTimeSupport(   cgMatrix<double>& residue,
//                             strtParameter* parm,
//                             int dummy);
// {

// }

// void CTriangDictionary::optimizeDecaying(  cgMatrix<double>& residue,
//                         strtParameter* parm_aux);
// {

// }

// void CTriangDictionary::optimizeDecayingErrorNorm(  cgMatrix<double>& residue,
//                                  strtParameter* parm_aux,
//                                  double& minResNorm,
//                                  double coef_xi0,
//                                  double ratomsample_xi0);
// {

// }

// void CTriangDictionary::quantizeFrequency( cgMatrix<double>& residue,
//                         strtParameter* parm,
//                         CDataSignal* dataSignal);
// {

// }

// void CTriangDictionary::discrimineSine(cgMatrix<double>& residue,
//                     strtParameter* parm,
//                     CDataSignal* dataSignal);
// {

// }

// void CTriangDictionary::searchSBPreviousBlock( cgMatrix<double>& residue,
//                             strtParameter* parm,
//                             CStructBook* sbPreviousBlock,
//                             int iPrevBlock,
//                             int flagFile,
//                             double coefHeur);
// {

// }

// void CTriangDictionary::evalAtomContinuity(cgMatrix<double>& residue,
//                         CStructBook* sbPreviousBlock,
//                         CStructBook* structBook,
//                         int iSignal,
//                         int iCurrBlock,
//                         int iPrevBlock,
//                         int& step,
//                         int flagFile,
//                         char* fileName,
//                         char* sbbFName,
//                         double initBlockNorm,
//                         fstream& file_stage,
//                         int flag_stage);
// {

// }

/*void matchingPursuit(   cgMatrix<double>& residue,
                        strtParameter* chosenParm,
                        int dicSize,
                        CDataSignal* dataSignal,
                        CFileDictionary* dicData,
                        CFileDecomp* genData, int step, int chosenDic)
{
    int decincAsymmFlag;
    int i,j;
    int k;
    // int s;
    double s;
    int a;
    int chosenTau = 0;cout<<"IP"<<((strtParameter*)parm)->innerProduct<<endl;
    double chosenXi = 0.0;
    int delta_tau;
    double chosenOptPhase = 0.0;
    double Fs;
    double beta = 0.0;
    double rho = 0.0;
    // double Ffund;
    // double delta_f = 0;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;
    int dicType;
    // int signalSize;
    double  innerProd;

    double* realAtom = new double[dicSize];
    double* realAtomDec = new double [dicSize];
    double* realAtomInc = new double [dicSize];
    double* realAtomWin = new double[2*dicSize];

    CComplex* complexAtom = new CComplex[dicSize];
    CComplex* complexAtomDec = new CComplex[dicSize];
    CComplex* complexAtomInc = new CComplex[dicSize];
    CComplex* complexAtomXi = new CComplex[dicSize];
    CComplex* complexAtom2Xi = new CComplex[dicSize];


    double* conv_zxi0_w2 = new double[2*dicSize];
    double* conv_zxi0_w2_expinc = new double[2*dicSize];

    CDictionary *expDic = new CExpDictionary;
    CDictionary * gaborDic = new CGaborDictionary;
    CDictionary * triangDic = new CTriangDictionary;
    CDictionary * triangDic2 = new CTriangDictionary;
    CDictionary * bateDic = new CBatemanDictionary;

    ((CExpDictionary*)expDic)->setSignalSize(dicSize);
    ((CGaborDictionary*)gaborDic)->setSignalSize(dicSize);
    ((CTriangDictionary*)triangDic)->setSignalSize(dicSize);
    ((CTriangDictionary*)triangDic2)->setSignalSize(2*dicSize);
    ((CBatemanDictionary*)bateDic)->setSignalSize(dicSize);

    strtParameter* parm = new strtParameter;
    ((strtParameter*)chosenParm)->innerProduct=0.0;
    ((strtParameter*)chosenParm)->rho = 0.0;
    ((strtParameter*)chosenParm)->xi = 0.0;
    ((strtParameter*)chosenParm)->phase = 0.0;
    ((strtParameter*)chosenParm)->a = 0;
    ((strtParameter*)chosenParm)->b = 0;
    double xi = 0.0;
    int tau = 0;
    double maxInnerProd = 0.0;


    int MPType;
    MPType = genData->getMPType();

    FILE* stream;
    char* fileName;

    int d = 0;


    for(j=0;j<dicData->getNumDicBlock();j++)
    {
        // stream = fopen( fileName, "a" );
        // fprintf(stream,"%d --------------------------------------------------------------\n", step);
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
        for (int n=0; n<dicSize; n++)
        {
            realAtom[n] = 0.0;
            realAtomDec[n] = 0.0;
            realAtomInc[n] = 0.0;
        }

        for (int n=0; n<2*dicSize; n++)
        {
            realAtomWin[n] = 0.0;
            conv_zxi0_w2[n] = 0.0;
            conv_zxi0_w2_expinc[n] = 0.0;
        }


        if (chosenDic == dicType || chosenDic == 0)
        {
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
            else if (dataSignal->getType() == 4) // ECG
            {
                Fs = ((CECGSignal*)dataSignal)->getSamplingRate();
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
                else if (dataSignal->getType() == 4) // ECG
                {
                    freqi = 0;
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
            else if (freqi==0000000000)
            {
                freqi = freqf/2;
            }

            if (freqf>Fs)
            {
                cout << "Final frequency greater than sampling frequency !!!" << endl;
                printf("final freq. %f - Fs %f\n",freqf,Fs);
                exit(1);
            }

            if (fdiscrtype==1) // linear
            {
                // delta_f = (2*PI/freqf)* freqi;
                Nfreq = (int)(freqf/(2*freqi));
                // cout << Nfreq << endl;
                //Nfreq = (int)ceil(freqf/freqi);
                xi_vec = new double[Nfreq];
                for (i=0;i<Nfreq;i++)
                {
                    // xi_vec[i] = delta_f * i;
                    xi_vec[i] = (2*PI/Fs) * (freqi * i );
                    // cout << xi_vec[i] << endl;
                }
            }
            else if (fdiscrtype==2) // geometric with quarter-tone discretization
            {
                Nfreq = (int)ceil(24 * ( log10(freqf/freqi)/log10(2.0) ) )+1;

                xi_vec = new double[Nfreq];
                xi_vec[0] = 0.0;
                for (i=1;i<Nfreq;i++)
                {
                    xi_vec[i] = (2*PI/Fs) * (freqi * pow (2.0, (double)(i-1)/24) );
                }
            }

            if (dicType == 1) //exponential
            {
                if (MPType == 1)
                {
                    for(k=0;k<Nfreq;k++)
                    {
                        for(tau=0;tau<(double)N;tau+=delta_tau)
                        {
                            // cout << tau << endl;

                            //Decreasing
                            //xi = k * ((2*PI)/s);
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = tau;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            complexAtomDec = expDic->getComplexAtom();
                            // cout << s  << ' ' << xi << ' ' << tau << endl;
                            // MPTradicional(  residue,
                            //                 complexAtomDec,
                            //                 maxInnerProd,
                            //                 chosenOptPhase,
                            //                 tau,
                            //                 xi,
                            //                 N, //int dicSize
                            //                 s,
                            //                 fileName);

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtomDec, fileName,s,tau,N);


                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                tau,
                                                N-1,
                                                0.0);
                                // cout << "AQUI" << endl;
                                // cout << ((strtParameter*)chosenParm)->s  << ' ' << ((strtParameter*)chosenParm)->xi << ' ' << ((strtParameter*)chosenParm)->u << endl;
                            }
                            d++;
                        }
                        for (tau=N-1;tau>0; tau=tau-delta_tau)
                        {
                            // Increasing
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = -1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = tau;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            complexAtomInc = expDic->getComplexAtom();

                            // for (int k1=0; k1<dicSize;k1++)
                            // {
                            //     complexAtomInc[k1] = *(expDic->getComplexAtom()+k1);
                            //     cout << complexAtomInc[k1].Real() << endl;
                            // }

                            // MPTradicional(  residue,
                            //                 complexAtomInc,
                            //                 maxInnerProd,
                            //                 chosenOptPhase,
                            //                 chosenTau,
                            //                 chosenXi,
                            //                 dic,
                            //                 N, //int dicSize
                            //                 // decincAsymmFlag,
                            //                 s,
                            //                 fileName);

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtomInc,fileName,s,tau,N);

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                -s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                0,
                                                tau,
                                                0.0);
                            }
                            d++;
                        }
                    }
                }

                else if (MPType==2)
                {
                    // int delta_tau = s;
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        //decreasing
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = tau;
                        ((strtParameter*)parm)->a = tau;
                        ((strtParameter*)parm)->b = N-1;
                        ((CExpDictionary*)expDic)->setRealAtom(parm);
                        realAtomDec = expDic->getRealAtom();
                        fastMPKolasa(   residue,
                                        maxInnerProd,
                                        chosenOptPhase,
                                        chosenXi,
                                        dicSize,
                                        tau,
                                        s,
                                        realAtomDec,
                                        fileName,
                                        d);

                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            tau,
                                            N-1,
                                            0.0);
                            // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                            // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                            // ((strtParameter*)chosenParm)->xi = chosenXi;
                            // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                            // ((strtParameter*)chosenParm)->a = tau;
                            // ((strtParameter*)chosenParm)->b = N-1;
                        }
                    }

                    for (tau=1; tau<=(double)N; tau=tau+delta_tau)
                    {
                        // increainsg
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = -1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->a = 0.0;
                        ((strtParameter*)parm)->b = tau;
                        ((strtParameter*)parm)->u = tau;
                        ((CExpDictionary*)expDic)->setRealAtom(parm);
                        realAtomInc = expDic->getRealAtom();

                        fastMPKolasa(   residue,
                                        maxInnerProd,
                                        chosenOptPhase,
                                        chosenXi,
                                        dicSize,
                                        tau,
                                        s,
                                        realAtomInc,
                                        fileName,
                                        d);
                        // cout << chosenXi << endl;
                        // cout << tau << endl;
                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << s << endl;
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            -s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            0,
                                            tau,
                                            0.0);
                            // cout << endl << "ATOM" << endl;
                            // cout << -1/double(s) << " " << tau << endl;

                            // for (int p = 0; p<512; p++)
                            // {
                            //     cout << realAtomInc[p] << endl;
                            // }

                            // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                            // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                            // ((strtParameter*)chosenParm)->xi = chosenXi;
                            // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                            // ((strtParameter*)chosenParm)->a = 0.0;
                            // ((strtParameter*)chosenParm)->b = tau;
                        }
                    }
                }

                else if (MPType==3)
                {
                    // cout << "AQUI" << endl;
                    ((strtParameter*)parm)->s = s;
                    ((strtParameter*)parm)->rho = 1.0/s;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = 0;
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)expDic)->setRealAtom(parm);
                    // realAtomWin = expDic->getRealAtom();
                    // memcpy(realAtomWin, (expDic->getRealAtom()), sizeof(double)*dicSize);
                    for (int k1=0; k1<dicSize;k1++)
                    {
                        realAtomWin[k1+N] = *(expDic->getRealAtom()+k1);
                    }

                    for(k=0;k<Nfreq;k++)
                    {
                        //xi = k * ((2*PI)/s);
                        xi = xi_vec[k];
                        if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                        {
                            // z1
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            // complexAtomXi = expDic->getComplexAtom();
                            // memcpy(complexAtomXi, (expDic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                            }

                            // z2
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;
                            ((strtParameter*)parm)->xi = 2*xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CExpDictionary*)expDic)->setComplexAtom(parm);
                            // complexAtom2Xi = expDic->getComplexAtom();
                            // memcpy(complexAtom2Xi, (expDic->getComplexAtom()), sizeof(CComplex) * (dicSize));

                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                            }


                            //Decreasing exponential
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
                                                    conv_zxi0_w2,
                                                    fileName,
                                                    s,
                                                    d);
                            // if (k == 0){
                            // }

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                chosenTau,
                                                N-1,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = xi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = chosenTau;
                                // ((strtParameter*)chosenParm)->b = N-1;
                            }

                            //Increasing exponential
                            decincAsymmFlag = -1;
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
                                                    conv_zxi0_w2_expinc,
                                                    fileName,
                                                    s,
                                                    d);
                            // // if (k == 0){
                            // // }

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                -s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                0,
                                                chosenTau,
                                                0.0);
                                // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                                // ((strtParameter*)chosenParm)->rho = -1.0/(double)s;
                                // ((strtParameter*)chosenParm)->xi = xi;
                                // ((strtParameter*)chosenParm)->phase = chosenOptPhase;
                                // ((strtParameter*)chosenParm)->a = 0;
                                // ((strtParameter*)chosenParm)->b = chosenTau;
                            }
                        }
                    }
                }
            }

            else if (dicType == 2) // pure sine
            {
                // ((CExpDictionary*)expDic)->setSignalSize(dicSize);

                for (tau=0; tau<(double)N; tau=tau+delta_tau)
                {
                    ((strtParameter*)parm)->rho = 0.0;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = tau;
                    ((strtParameter*)parm)->a = tau;
                    ((strtParameter*)parm)->b = N-1;
                    ((CExpDictionary*)expDic)->setRealAtom(parm);
                    realAtom = expDic->getRealAtom();
                    fastMPKolasa(   residue,
                                    maxInnerProd,
                                    chosenOptPhase,
                                    chosenXi,
                                    dicSize,
                                    tau,
                                    s,
                                    realAtom,
                                    fileName,
                                    d);
                    if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    {
                        // cout << maxInnerProd << endl;
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        chosenXi,
                                        chosenOptPhase,
                                        tau,
                                        tau,
                                        N-1,
                                        0.0);
                    }
                }
            }

            else if (dicType == 3) // impulse
            {
                for(i=0;i<N;i++)
                {
                    innerProd = residue[0][i];
                    if (fabs(innerProd)>fabs(maxInnerProd))
                    {
                        // cout << maxInnerProd << endl;
                        maxInnerProd = innerProd;
                        setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        0.0,
                                        0.0,
                                        i,
                                        i,
                                        i,
                                        0.0);
                        // ((strtParameter*)chosenParm)->innerProduct = maxInnerProd;
                        // ((strtParameter*)chosenParm)->rho = 1.0/(double)s;
                        // ((strtParameter*)chosenParm)->xi = 0.0;
                        // ((strtParameter*)chosenParm)->phase = 0.0;
                        // ((strtParameter*)chosenParm)->a = i;
                        // ((strtParameter*)chosenParm)->b = i;
                    }
                    d++;
                }
            }

            else if (dicType == 4) // Gabor
            {
                if (MPType == 1)
                {
                    for(k=0;k<Nfreq;k++)
                    {
                        for(tau=0;tau<N;tau+=delta_tau)
                        {
                            //Decreasing
                            //xi = k * ((2*PI)/s);
                            // cout << s << endl;
                            // cout << tau << endl;
                            xi = xi_vec[k];
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CGaborDictionary*)gaborDic)->setComplexAtom(parm);
                            complexAtom = gaborDic->getComplexAtom();


                            // MPTradicional(  residue,
                            //                 complexAtomDec,
                            //                 maxInnerProd,
                            //                 chosenOptPhase,
                            //                 tau,
                            //                 xi,
                            //                 N, //int dicSize
                            //                 s,
                            //                 fileName);

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                            // cout << maxInnerProd << endl;

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                0,
                                                N-1,
                                                0.0);
                            }
                            d++;
                        }
                    }

                    // for(k=0;k<Nfreq;k++)
                    // {
                    //     for(tau=0;tau<N;tau+=delta_tau)
                    //     {
                    //         // cout << xi /*<< " " << xi << " " << dicType << endl;
                    //         //xi = k * ((2*PI)/s);
                    //         xi = xi_vec[k];
                    //         ((strtParameter*)parm)->s = s;
                    //         ((strtParameter*)parm)->rho = (1.0/(double)s);
                    //         ((strtParameter*)parm)->xi = xi;
                    //         ((strtParameter*)parm)->phase = 0.0;
                    //         ((strtParameter*)parm)->a = tau;
                    //         ((strtParameter*)parm)->b = N-1;
                    //         ((CExpDictionary*)dic)->setComplexAtom(parm);
                    //         complexAtom = dic->getComplexAtom();
                    //         cout << xi << " " << tau << " " << &complexAtom[0] << endl;

                    //         computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                    //         // MPTradicional(  residue,
                    //         //                 complexAtom,
                    //         //                 maxInnerProd,
                    //         //                 chosenOptPhase,
                    //         //                 tau,
                    //         //                 xi,
                    //         //                 // dic,
                    //         //                 N, //int dicSize
                    //         //                 s,
                    //         //                 fileName);
                    //         // cout << maxInnerProd << endl;
                    //         if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                    //         {
                    //             setParameters ( chosenParm,
                    //                             dicType,
                    //                             maxInnerProd,
                    //                             s,
                    //                             xi,
                    //                             chosenOptPhase,
                    //                             tau,
                    //                             0,
                    //                             N-1);
                    //         }
                    //     }
                    // }
                }

                else if (MPType==2)
                {
                    int delta_tau = (int)s;
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = tau;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CGaborDictionary*)gaborDic)->setRealAtom(parm);
                        realAtom = gaborDic->getRealAtom();
                        // cout << "AQUI" << endl;
                        fastMPKolasa(   residue,
                                        maxInnerProd,
                                        chosenOptPhase,
                                        chosenXi,
                                        dicSize,
                                        tau,
                                        s,
                                        realAtom,
                                        fileName,
                                        d);
                        // cout << "AQUI2" << endl;

                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                            // cout << maxInnerProd << endl;
                            setParameters ( chosenParm,
                                            dicType,
                                            maxInnerProd,
                                            s,
                                            chosenXi,
                                            chosenOptPhase,
                                            tau,
                                            0,
                                            N-1,
                                            0.0);
                        }
                    }
                }


                else if (MPType==3)
                {
                    ((strtParameter*)parm)->s = s;
                    ((strtParameter*)parm)->rho = 1.0/s;
                    ((strtParameter*)parm)->xi = 0.0;
                    ((strtParameter*)parm)->phase = 0.0;
                    ((strtParameter*)parm)->u = N-1;
                    ((strtParameter*)parm)->a = 0;
                    ((strtParameter*)parm)->b = N-1;
                    ((CGaborDictionary*)gaborDic)->setRealAtom(parm);
                    // memcpy(realAtomWin, (gaborDic->getRealAtom()), sizeof(double)*dicSize);
                    // cout << "realAtomWin" << endl;
                    for (int k1=0; k1<dicSize;k1++)
                    {
                        realAtomWin[k1+1] = *(gaborDic->getRealAtom()+k1); //add 1 to centralize de max
                        realAtomWin[(2*N-1)-k1] = *(gaborDic->getRealAtom()+k1);
                    }

                    // cout << "AQUI" << endl;
                    // for (int k1=0; k1<2*dicSize; k1++)
                    // {
                    //     cout << realAtomWin[k1] << endl;
                    // }
                    // cout << "ATOM" << endl;
                    // for (int k1=0;k1<2*dicSize;k1++)
                    // {
                    //     cout << realAtomWin[k1] << endl;
                    // }

                    for(k=0;k<Nfreq;k++)
                    {
                        xi = xi_vec[k];
                        if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                        {
                            // cout << s << endl;
                            // cout << xi << endl;
                            // z1
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0; //generating complex exponential with rectangular window
                            ((strtParameter*)parm)->xi = xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CGaborDictionary*)expDic)->setComplexAtom(parm);
                            // cout << "complexAtomXi" << endl;
                            // memcpy(complexAtomXi, (gaborDic->getComplexAtom()), sizeof(CComplex)*(dicSize));
                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                                // cout << complexAtomXi[K1] << endl;
                            }
                            // cout << "AQUI" << endl;
                            // for (int k1=0; k1<dicSize; k1++)
                            // {
                            //     cout << complexAtomXi[k1].Real() << " " << complexAtomXi[k1].Imag() << endl;
                            // }


                            // z2
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 0.0;  //generating complex exponential with rectangular window
                            ((strtParameter*)parm)->xi = 2*xi;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = 0;
                            ((strtParameter*)parm)->a = 0;
                            ((strtParameter*)parm)->b = N-1;
                            ((CGaborDictionary*)expDic)->setComplexAtom(parm);
                            // memcpy(complexAtom2Xi, (gaborDic->getComplexAtom()), sizeof(CComplex)*(dicSize));
                            // cout << "complexAtom2Xi" << endl;
                            for (int k1=0; k1<dicSize;k1++)
                            {
                                complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                                // cout << complexAtom2Xi[K1] << endl;
                            }

                            fastMPKolasaModified(   residue,
                                                    maxInnerProd,
                                                    chosenOptPhase,
                                                    chosenTau,
                                                    dicSize,
                                                    1,
                                                    delta_tau,
                                                    xi,
                                                    realAtomWin,
                                                    complexAtomXi,
                                                    complexAtom2Xi,
                                                    conv_zxi0_w2,
                                                    fileName,
                                                    s,
                                                    d);

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                chosenTau,
                                                0,
                                                N-1,
                                                0.0);
                            }
                        }
                    }
                }
            }

            else if (dicType == 5) // triangular
            {
                if (MPType == 1 || MPType == 2)
                {
                    for(tau=0;tau<N;tau+=delta_tau)
                    {
                        // cout << "a: " << a << endl;
                        for(a=0;a<s;a+=delta_tau)
                        {
                            // cout << "tau: " << tau << endl;
                            //Decreasing
                            ((strtParameter*)parm)->s = s;
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = a;
                            ((strtParameter*)parm)->b = (int)s-1-a;
                            ((CTriangDictionary*)triangDic)->setComplexAtom(parm);
                            complexAtom = triangDic->getComplexAtom();

                            // cout << "AQUI" << endl;
                            // cout <<" s- " << s << "tau- " << tau << "a- " << a << endl;

                            // for (int k1=0;k1<N;k1++)
                            // {
                            //     // complexAtomXi[k1] = *(triangDic->getComplexAtom()+k1);
                            //     // cout << residue[0][k1] << endl;
                            //     cout << complexAtom[k1].Real() << endl;
                            // }

                            computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);
                            // cout << maxInnerProd << endl;
                            // cout << "maxInnerProd " << maxInnerProd << endl;
                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                // cout << "s- " << s << "tau- " << tau << "a- " << a << endl;
                                // cout << maxInnerProd << endl;
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                xi,
                                                chosenOptPhase,
                                                tau,
                                                a,
                                                (int)s-a-1,
                                                0.0);
                            }
                            d++;
                        }
                    }
                }

                else if (MPType==3)
                {
                    // cout << s << endl;
                    for(a=1;a<(int)s;a+=delta_tau)
                    {
                        // cout << a << endl;
                        ((strtParameter*)parm)->s = s;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = N;
                        ((strtParameter*)parm)->a = a;
                        ((strtParameter*)parm)->b = (int)s-a-1;
                        ((CTriangDictionary*)triangDic2)->setRealAtom(parm);

                        // cout << "s " << s << " - a " << a << endl;
                        for (int k1=0; k1<2*dicSize;k1++)
                        {
                            realAtomWin[k1] = *(triangDic2->getRealAtom()+k1);
                            // cout << realAtomWin[k1] << endl;
                            // realAtomWin[(2*N-1)-k1] = *(triangDic2->getRealAtom()+k1);
                        }
                        // z1
                        ((strtParameter*)parm)->rho = 0.0;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CTriangDictionary*)expDic)->setComplexAtom(parm);

                        for (int k1=0; k1<dicSize;k1++)
                        {
                            complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                        }

                        // z2
                        ((strtParameter*)parm)->rho = 0.0;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CTriangDictionary*)expDic)->setComplexAtom(parm);

                        // cout << "AQUI1" << endl;
                        for (int k1=0; k1<dicSize;k1++)
                        {
                            complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                            // cout << residue[1][k1] << endl;
                        }
                        // cout << "AQUI2" << endl;

                        fastMPKolasaModified(   residue,
                                                maxInnerProd,
                                                chosenOptPhase,
                                                chosenTau,
                                                dicSize,
                                                1,
                                                delta_tau,
                                                xi,
                                                realAtomWin,
                                                complexAtomXi,
                                                complexAtom2Xi,
                                                conv_zxi0_w2,
                                                fileName,
                                                s,
                                                d);

                        if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                        {
                                setParameters ( chosenParm,
                                        dicType,
                                        maxInnerProd,
                                        s,
                                        xi,
                                        chosenOptPhase,
                                        chosenTau,
                                        a,
                                        (int)s-a-1,
                                        0.0);
                        }
                    }
                }
            }

            if (dicType == 6) //Bateman
            {
                if (MPType == 1)
                {
                    for (beta = 2.0; beta > 1.0/s; beta/=2.0)
                    {
                        for(k=0;k<Nfreq;k++)
                        {
                            for(tau=0;tau<(double)N;tau+=delta_tau)
                            {
                                xi = xi_vec[k];
                                ((strtParameter*)parm)->rho = 1.0/s;
                                ((strtParameter*)parm)->xi = xi;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = tau;
                                ((strtParameter*)parm)->a = tau;
                                ((strtParameter*)parm)->b = N-1;
                                ((strtParameter*)parm)->beta = beta;
                                ((CBatemanDictionary*)bateDic)->setComplexAtom(parm);
                                complexAtom = bateDic->getComplexAtom();

                                computeOptimumPhase(residue,chosenOptPhase,maxInnerProd,N,xi,complexAtom, fileName,s,tau,N);

                                if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                                {
                                    setParameters ( chosenParm,
                                                    dicType,
                                                    maxInnerProd,
                                                    s,
                                                    xi,
                                                    chosenOptPhase,
                                                    tau,
                                                    tau,
                                                    N-1,
                                                    beta);
                                }
                                d++;
                            }
                        }
                    }
                }

                else if (MPType==2)
                {
                    // int delta_tau = s;
                    for (beta=1.0; beta > 1.0/s; beta/=(2.0))
                    {
                        for (tau=0; tau<(double)N; tau=tau+delta_tau)
                        {
                            ((strtParameter*)parm)->rho = 1.0/s;
                            ((strtParameter*)parm)->beta = beta;
                            ((strtParameter*)parm)->xi = 0.0;
                            ((strtParameter*)parm)->phase = 0.0;
                            ((strtParameter*)parm)->u = tau;
                            ((strtParameter*)parm)->a = tau;
                            ((strtParameter*)parm)->b = N-1;
                            ((CBatemanDictionary*)bateDic)->setRealAtom(parm);
                            realAtom = bateDic->getRealAtom();
                            fastMPKolasa(   residue,
                                            maxInnerProd,
                                            chosenOptPhase,
                                            chosenXi,
                                            dicSize,
                                            tau,
                                            s,
                                            realAtom,
                                            fileName,
                                            d);

                            if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                            {
                                setParameters ( chosenParm,
                                                dicType,
                                                maxInnerProd,
                                                s,
                                                chosenXi,
                                                chosenOptPhase,
                                                tau,
                                                tau,
                                                N-1,
                                                beta);
                            }
                        }
                    }
                }

                else if (MPType==3)
                {
                    for (beta=1.0; beta > 1.0/s; beta/=2.0)
                    {
                        ((strtParameter*)parm)->rho = 1.0/s;
                        ((strtParameter*)parm)->beta = beta;
                        ((strtParameter*)parm)->xi = 0.0;
                        ((strtParameter*)parm)->phase = 0.0;
                        ((strtParameter*)parm)->u = 0;
                        ((strtParameter*)parm)->a = 0;
                        ((strtParameter*)parm)->b = N-1;
                        ((CBatemanDictionary*)bateDic)->setRealAtom(parm);

                        for (int k1=0; k1<dicSize;k1++)
                        {
                            realAtomWin[k1+N] = *(bateDic->getRealAtom()+k1);
                        }

                        for(k=0;k<Nfreq;k++)
                        {
                            //xi = k * ((2*PI)/s);
                            xi = xi_vec[k];
                            if ( (xi==0.0) || (xi>=((2*PI)/s) ))
                            {
                                // z1
                                ((strtParameter*)parm)->rho = 0.0;
                                ((strtParameter*)parm)->beta = 0.0;
                                ((strtParameter*)parm)->xi = 0.0;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = 0;
                                ((strtParameter*)parm)->a = 0;
                                ((strtParameter*)parm)->b = N-1;
                                ((CExpDictionary*)expDic)->setComplexAtom(parm);

                                for (int k1=0; k1<dicSize;k1++)
                                {
                                    complexAtomXi[k1] = *(expDic->getComplexAtom()+k1);
                                }

                                // z2
                                ((strtParameter*)parm)->rho = 0.0;
                                ((strtParameter*)parm)->beta = 0.0;
                                ((strtParameter*)parm)->xi = 0.0;
                                ((strtParameter*)parm)->phase = 0.0;
                                ((strtParameter*)parm)->u = 0;
                                ((strtParameter*)parm)->a = 0;
                                ((strtParameter*)parm)->b = N-1;
                                ((CExpDictionary*)expDic)->setComplexAtom(parm);

                                for (int k1=0; k1<dicSize;k1++)
                                {
                                    complexAtom2Xi[k1] = *(expDic->getComplexAtom()+k1);
                                }


                                fastMPKolasaModified(   residue,
                                                        maxInnerProd,
                                                        chosenOptPhase,
                                                        chosenTau,
                                                        dicSize,
                                                        1,
                                                        delta_tau,
                                                        xi,
                                                        realAtomWin,
                                                        complexAtomXi,
                                                        complexAtom2Xi,
                                                        conv_zxi0_w2,
                                                        fileName,
                                                        s,
                                                        d);

                                if (fabs(maxInnerProd)>fabs(((strtParameter*)chosenParm)->innerProduct))
                                {
                                    setParameters ( chosenParm,
                                                    dicType,
                                                    maxInnerProd,
                                                    s,
                                                    xi,
                                                    chosenOptPhase,
                                                    chosenTau,
                                                    chosenTau,
                                                    N-1,
                                                    beta);
                                }
                            }
                        }
                    }
                }
            }

            if (xi_vec) delete [] xi_vec;
        }
    }

    if (parm) delete parm;

    if (expDic) delete expDic;
    if (gaborDic) delete gaborDic;
    if (triangDic) delete triangDic;
    if (triangDic) delete triangDic2;
    if (bateDic) delete bateDic;

    if (conv_zxi0_w2_expinc) delete[] conv_zxi0_w2_expinc;
    if (conv_zxi0_w2) delete[] conv_zxi0_w2;

    return;
}*/
