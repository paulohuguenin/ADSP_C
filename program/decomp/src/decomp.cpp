#include "decomp.h"

void decompElectric(CFileDecomp* genData,
                    CFileDecompBlockRange* blockRange,
                    CFileDictionary* dicData,
                    char* InputFile)
{
    //////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
    dataSignal = new CComtradeSignal;

    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // -----------------------------------------------------------
    //  Decomposing signal
    //
    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBook** structBook;
    structBook = new CStructBook* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookExp[numBlock];
    }

    // Setting Dictionary
    CDictionary* dic;
    // if (genData->getDicType()==1)
    // {
    dic = new CExpDictionary;
    // }
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int dicSize = (int) pow(2.0,(double)nbits);
    dic->setSignalSize(dicSize);

    // Decomposing Blocks of Signals

    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();
    int nMaxStep = genData->getNumMaxStep();


    cgMatrix<double> residue(1,dicSize,0.0);
    cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);
    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;


    // CParameter* chosenParm;

    // if (genData->getDicType()==1)
    // {
    //     chosenParm = new CExpParm;
    // }

    strtParameter* chosenParm;
    chosenParm = new strtParameter;

    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }
    // Allocate memory for candidate atoms with tume continuity
    CStructBook** sbContinuity;
    sbContinuity = new CStructBook* [numSignal];
    int k;
    for (k=0; k<numSignal; k++)
    {
        sbContinuity[k] = new CStructBookExp[numBlock];
    }

    int L = (int)ceil(((log10((double)(dicSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];
    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }
    double meanApproxRatio;

    double tolAppRatio =  genData->getApproxRatioTarget();
    double snrTarget = genData->getSNRTarget();

    double befSupInnerP,aftSupInnerP;

    char fileName[_MAX_PATH];
    strcpy(fileName, dataSignal->getFileName());
    char* pos;
    pos = strrchr( fileName, '.');
    char aux[_MAX_PATH];
    sprintf(aux,"_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    // Writing the Main Header
    FILE* iosba;
    iosba = fopen(fileName,"w");
    fprintf(iosba,"Sign. Type :          %5i\n", dataSignal->getType());
    // fprintf(iosba,"Dict. Type :          %5i\n", genData->getDicType());
    fprintf(iosba,"No. Signals:          %5i\n", dataSignal->getNumSignal());
    fprintf(iosba,"Signal Size:       %8i\n", dataSignal->getSignalSize());
    fprintf(iosba,"Block Hop:            %5i\n", dataSignal->getBlockHop());
    fprintf(iosba,"Block Size:           %5i\n", dataSignal->getBlockSize());
    if (genData->getSigType()==1)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    }
    if (genData->getSigType()==2)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    }
    fprintf(iosba,"Init. Block:          %5i\n", initBlock);
    fprintf(iosba,"Final Block:          %5i\n", finalBlock);
    fflush(iosba);
    fclose(iosba);

    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d_header.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    FILE* iosbb;
    iosbb = fopen(sbbFName,"wb");
    int dummyint;
    double dummydouble;
    dummyint = dataSignal->getType();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = genData->getDicType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getNumSignal();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getSignalSize();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getBlockHop();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getBlockSize();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    if (genData->getSigType()==1)
        dummydouble = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
    if (genData->getSigType()==2)
        dummydouble = ((CAudioSignal*)dataSignal)->getSamplingRate();
    fwrite(&dummydouble, sizeof(double), 1, iosbb);
    dummyint = initBlock;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = finalBlock;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    fclose(iosbb);

    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);
    iosbb = fopen(sbbFName,"wb");
    fclose(iosbb);

    fstream file_stage;
    if (genData->getPrintDecompStage()==1)
    {
        //file pointers
        file_stage.open("decomp_stages.out",ios::out);
        // file header
        file_stage  <<  setw (10) << setfill(' ') << "Signal" << " "
                    <<  setw (10) << setfill(' ') << "Block" << " "
                    <<  setw (10) << setfill(' ') << "No." << " "
                    <<  setw (10) << setfill(' ') << "Stage" << " "
                    <<  setw (20) << setfill(' ') << "Coef." << " "
                    <<  setw (20) << setfill(' ') << "Decay" << " "
                    <<  setw (20) << setfill(' ') << "Freq" << " "
                    <<  setw (20) << setfill(' ') << "Phase"<< " "
                    <<  setw (10) << setfill(' ') << "Ti"<< " "
                    <<  setw (10) << setfill(' ') << "Tf"<< " "
                    << endl;
    }

    int iSignal;
    double sigNorm;
    int iBlock;
    for (i=0; i < dataSignal->getNumSignal(); i++ )
    {
        // Writing the Signal Header
        iosba = fopen(fileName,"a");
        fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        fprintf(iosba,"Signal:               %5i\n",i+1);
        fprintf(iosba,"Norm:            %10.5f\n",dataSignal->getNorm(i));
        fflush(iosba);
        fclose(iosba);

        iosbb = fopen(sbbFName,"ab");
        iSignal = i+1;
        fwrite(&iSignal, sizeof(int), 1, iosbb);
        sigNorm = dataSignal->getNorm(i);
        fwrite(&sigNorm, sizeof(double), 1, iosbb);
        fclose(iosbb);

        if (dataSignal->getNorm(i)!=0.0)
        {
            for (int j=initBlock-1; j < finalBlock; j++ )
            {

                residue.zeros();
                // Loading signal into the vector
                if ((int)(j*(double)dataSignal->getBlockHop())+dicSize < dataSignal->getSignalSize())
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dicSize);
                }
                else
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dataSignal->getSignalSize() - (int)(j*(double)dataSignal->getBlockHop()));
                }
                // Normalizing initial residue (with signal norm)
                residue /= dataSignal->getNorm(i);


                norm = residue.norm();
                ((CStructBookExp*)structBook[i])[j].setNorm(norm);
                double initBlockNorm = norm;

                // Writing the Block Header
                iosba = fopen(fileName,"a");
                fprintf(iosba,"--------------------------------------------------------------\n");
                fprintf(iosba,"Block:                %5i\n",j+1);
                fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
                fprintf(iosba,"No.    Coef.           Decaying        Freq            Phase           Ti   Tf    PrevAtom AppRatio   meanAppRat befSup     aftSup     normRatio  SNR(dB)     \n");
                fflush(iosba);
                fclose(iosba);

                iosbb = fopen(sbbFName,"ab");
                iBlock = j+1;
                fwrite(&iBlock, sizeof(int), 1, iosbb);
                fwrite(&initBlockNorm, sizeof(double), 1, iosbb);
                fclose(iosbb);


                int decomp_stage;

                // Beginning of the decomposition
                if (norm!=0)
                {
                    if ((j!=0) &&
                        (genData->getFlagEvalAtomCont()==1))
                    {
                        ((CExpDictionary*)dic)->evalAtomContinuity( residue,
                                                                    sbContinuity[i],
                                                                    structBook[i],
                                                                    i,
                                                                    j,
                                                                    j-1,
                                                                    step,
                                                                    0,
                                                                    fileName,
                                                                    sbbFName,
                                                                    initBlockNorm,
                                                                    file_stage,
                                                                    genData->getPrintDecompStage());
                    }
                    int nAtomCont = step;
                    cout << "- nAtomCont: " << nAtomCont << endl;
                    do
                    {
                        if (step>=nMaxStep) break;
#ifdef DBG_WATCH_DCMP_STEP
                        residue.PrintToFile("residue.dat");
#endif
                        cout << "##########################################################################" << endl;
                        cout << "->Decomposing Signal: " << i+1 << "; Block: "<< j+1 << "; Atom: "<< step+1 << endl;
                        cout << "##########################################################################" << endl;
                        // Signal projection over a discrete dictionary
                        cout << "Signal projection over a discrete dictionary" << endl;
                        //expParm = fastMPKolasa(residue);
                        ((CExpDictionary*)dic)->fastMPKolasaModified(residue,dataSignal,dicData,chosenParm,step);
                        // ((CExpParm*)chosenParm)->printParm2Screen();
                        cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                <<"| a: "<<((strtParameter*)chosenParm)->a
                                <<"| b: "<<((strtParameter*)chosenParm)->b
                                <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                << endl;
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                        << endl;
                        }
                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        approxRatio[step%L] = fabs(((strtParameter*)chosenParm)->innerProduct)/norm;
                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }
                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;
                        cout << "Tol. Approx. Ratio: " << tolAppRatio << endl;

                        // Maximize approximation by finding optimum parameters
                        //optimizeContinuousParms(residue, chosenParm);

                        // Heuristics
                        befSupInnerP = 0.0;
                        aftSupInnerP = 0.0;

                        if (((strtParameter*)chosenParm)->a != ((strtParameter*)chosenParm)->b)  // damp and pure cases
                        {
                            if (genData->getFlagFindSupport()==1)
                            {
                                befSupInnerP = ((strtParameter*)chosenParm)->innerProduct;
                                cout<< "->Finding the best time support ..." << endl;
                                if (((strtParameter*)chosenParm)->rho==0.0)
                                {
                                    // Two-way search
                                    ((CExpDictionary*)dic)->findFastBestTimeSupport(residue,chosenParm,1,
                                                                                    genData->getCoefTempSup());
                                }
                                else
                                {
                                    // One-way search
                                    ((CExpDictionary*)dic)->findFastBestTimeSupport(residue,chosenParm,0,
                                                                                    genData->getCoefTempSup());
                                }
                                // ((CExpParm*)chosenParm)->printParm2Screen();
                                cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                        <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                        <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                        <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                        <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                        <<"| a: "<<((strtParameter*)chosenParm)->a
                                        <<"| b: "<<((strtParameter*)chosenParm)->b
                                        <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                        << endl;

                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=2;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                                <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                                <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                                << endl;
                                }
                                aftSupInnerP = ((strtParameter*)chosenParm)->innerProduct;
                            }
                            if (genData->getFlagOptDecay()==1)
                            {
                                ((CExpDictionary*)dic)->optimizeDecaying(residue,chosenParm);
                                // ((CExpParm*)chosenParm)->printParm2Screen();
                                cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                        <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                        <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                        <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                        <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                        <<"| a: "<<((strtParameter*)chosenParm)->a
                                        <<"| b: "<<((strtParameter*)chosenParm)->b
                                        <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                        << endl;
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=3;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                                <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                                <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                                <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                                << endl;
                                }
                            }

                        }
                        if (((strtParameter*)chosenParm)->rho != 0.0) // damp case
                        {
                            // Discrimine Sine
                            cout<< "->Discrimining sine..." << endl;
                            ((CExpDictionary*)dic)->discrimineSine(residue,chosenParm,dataSignal);
                            // ((CExpParm*)chosenParm)->printParm2Screen();
                            cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                        <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                        <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                        <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                        <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                        <<"| a: "<<((strtParameter*)chosenParm)->a
                                        <<"| b: "<<((strtParameter*)chosenParm)->b
                                        <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                        << endl;
                            if (genData->getPrintDecompStage()==1)
                            {
                                decomp_stage=4;
                                file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                            <<  setw (10) << setfill(' ') << j+1 << " "
                                            <<  setw (10) << setfill(' ') << step+1 << " "
                                            <<  setw (10) << setfill(' ') << decomp_stage << " "
                                            <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                            <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                            <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                            <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                            <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                            <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                            << endl;
                            }
                        }


                        if ((j!=0) && (genData->getFlagSeqBlock()==1))
                        {
                            cout<< "->Search in the Previous Block Structure Book ..." << endl;
                            ((CExpDictionary*)dic)->searchSBPreviousBlock(  residue,
                                                                            chosenParm,
                                                                            sbContinuity[i],
                                                                            j-1,
                                                                            0,
                                                                            genData->getCoefSeqBlock());
                            // ((CExpParm*)chosenParm)->printParm2Screen();
                            cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                        <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                        <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                        <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                        <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                        <<"| a: "<<((strtParameter*)chosenParm)->a
                                        <<"| b: "<<((strtParameter*)chosenParm)->b
                                        <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                        << endl;
                        }


                        // Adjusting parameters
                        ((CExpDictionary*)dic)->adjustParameters(residue,chosenParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        // ((CExpParm*)chosenParm)->printParm2Screen();
                        cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                        <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                        <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                        <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                        <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                        <<"| a: "<<((strtParameter*)chosenParm)->a
                                        <<"| b: "<<((strtParameter*)chosenParm)->b
                                        <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                        << endl;

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;

                        ((CExpDictionary*)dic)->setRealAtom(chosenParm);
                        cgRealAtom.fillVector(((CExpDictionary*)dic)->getRealAtom());
                        cgRealAtomAux = cgRealAtom*((strtParameter*)chosenParm)->innerProduct;
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        residue = residue - cgRealAtom*((strtParameter*)chosenParm)->innerProduct;
                        norm = residue.norm();
                        // ((CExpParm*)chosenParm)->printParm2Screen();
                        cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                        <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                        <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                        <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                        <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                        <<"| a: "<<((strtParameter*)chosenParm)->a
                                        <<"| b: "<<((strtParameter*)chosenParm)->b
                                        <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                        << endl;

                         // Add element to structure book
                        ((CStructBookExp*)structBook[i])[j].addElement(chosenParm);
                        int indorig = ((CStructBookExp*)structBook[i])[j].getNumElement();
                        ((CStructBookExp*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);

                        // Add element to structure book with candidate atom with continuity
                        // for the next block
                        if (((strtParameter*)chosenParm)->b==dicSize-1)
                        {
                            ((CStructBookExp*)sbContinuity[i])[j].addElement(chosenParm);
                            int ind = ((CStructBookExp*)sbContinuity[i])[j].getNumElement();
                            ((CStructBookExp*)sbContinuity[i])[j].setNextAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setPrevAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setOrigAtomIndex(ind-1,indorig-1);
                        }




                        iosba = fopen(fileName,"a");
                        ((CStructBookExp*)structBook[i])[j].saveElementASCII(   iosba,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                0,
                                                                                0,
                                                                                (norm/initBlockNorm));
                        fflush(iosba);
                        fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;
                        cout << "SNR Target: " << snrTarget<< endl;

                        iosbb = fopen(sbbFName,"ab");
                        ((CStructBookExp*)structBook[i])[j].saveElementBin(iosbb);
                        fclose(iosbb);

                        step++;
                    }
                    while(      (   (meanApproxRatio > tolAppRatio) ||
                                    (step<(L+nAtomCont)) )
                                //((norm/initBlockNorm)>1e-8)
                                && (step<nMaxStep)
                                && (20*log10(initBlockNorm/norm)<snrTarget)
                                //&& ( fabs(expParm.innerProd) > 1e-8 )
                         );
                }
                else
                {
                    cout << "  ### Block "<< j+1 <<" with null samples ### " << endl;
                    iosba = fopen(fileName,"a");
                    fprintf(iosba,"###### Block with null samples ######\n");
                    fflush(iosba);
                    fclose(iosba);
                }

                iosba = fopen(fileName,"a");
                fprintf(iosba,"99999\n");
                fflush(iosba);
                fclose(iosba);

                iosbb = fopen(sbbFName,"ab");
                dummyint = 99999;
                fwrite( &dummyint, sizeof( int ), 1, iosbb );
                fclose(iosbb);

                step = 0;
            }
            iosba = fopen(fileName,"a");
            fprintf(iosba,"88888\n");
            fflush(iosba);
            fclose(iosba);

            iosbb = fopen(sbbFName,"ab");
            dummyint = 88888;
            fwrite( &dummyint, sizeof( int ), 1, iosbb );
            fclose(iosbb);
        }
        else
        {
                cout << "  ### Signal "<< i+1 <<" with null samples ### " << endl;
                iosba = fopen(fileName,"a");
                fprintf(iosba,"###### Signal with null samples ######\n");
                fflush(iosba);
                fclose(iosba);
        }
        iosba = fopen(fileName,"a");
        fprintf(iosba,"77777\n");
        fflush(iosba);
        fclose(iosba);

        iosbb = fopen(sbbFName,"ab");
        dummyint = 77777;
        fwrite( &dummyint, sizeof( int ), 1, iosbb );
        fclose(iosbb);

        //delete [] sbPrevProj;
    }

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }


    // Deallocating objects
    delete [] approxRatio;

    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)structBook[i]);
    }
    delete [] structBook;
    // if (genData->getDicType()==1)
    // {
    delete (CExpDictionary*)dic;
    // }

    delete (CComtradeSignal*)dataSignal;

    return;
}

void decompAudio(   CFileDecomp* genData,
                    CFileDecompBlockRange* blockRange,
                    CFileDictionary* dicData,
                    char* InputFile)
{
    int chosenDic;
    int chosenNet = 0;
    int endingFl;
    // if ( genData->getFlagANN() == 1 )
    // {
    //     const char *args[] = { "-nodisplay" };
    //     mclInitializeApplication(args,1);
    //     libannformpInitialize();
    // }

    chosenDic = genData->getDicType();


    //////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
    dataSignal = new CAudioSignal;

    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBookParm** structBook;
    structBook = new CStructBookParm* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookParm[numBlock];
    }

    clock_t c2, c1; /* variáveis que contam ciclos do processador */
    float tempo;
    FILE* stream;
    stream = fopen("time.out", "w");
    fprintf(stream, "bloco           tempo \n");
    fflush(stream);
    fclose(stream);

    // Setting Dictionary
    // CDictionary* dic;
    // dic = new CExpDictionary;

    // if (genData->getDicType()==1)
    // {
        // dic = new CExpDictionary;
    // }
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int dicSize = (int) pow(2.0,(double)nbits);

    // ((CDictionary*)dic)->setSignalSize(dicSize);

    // Decomp Stage
    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();

    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }

    int nMaxStep = genData->getNumMaxStep();

    cgMatrix<double> residue(1,dicSize,0.0);
    cgMatrix<double> signal(1,dicSize,0.0);
    // mwArray residueaux;
    // cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    // cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);

    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;


    // CParameter* chosenParm;
    // if (genData->getDicType()==1)
    // {
    //     chosenParm = new CExpParm;
    // }

    strtParameter* chosenParm;
    chosenParm = new strtParameter;


    int L = (int)ceil(((log10((double)(dicSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];
    int k;
    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }

    double meanApproxRatio;

    double tolAppRatio =  genData->getApproxRatioTarget();
    double snrTarget = genData->getSNRTarget();

    double befSupInnerP,aftSupInnerP;

    strtParameter** dicAtoms = new strtParameter*[(dicData->getNumDicBlock()*dicSize)];
    for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++)
    {
        dicAtoms[kdic] = new strtParameter;
    }

    char aux[_MAX_PATH];

    // ((CStructBookParm*)structBook) -> newFileASCII(initBlock,finalBlock,dataSignal,dicData);
    char fileName[_MAX_PATH];
    strcpy(fileName, dataSignal->getFileName());
    char* pos;
    pos = strrchr( fileName, '.');
    sprintf(aux,"_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    char aux2[_MAX_PATH];
    char fileName2[_MAX_PATH];
    strcpy(fileName2, dataSignal->getFileName());
    char* pos2;
    pos2 = strrchr( fileName2, '.');
    sprintf(aux2,"_OMP_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos2[0], aux2);


    long posaux;

    // double* realAtom;

    // Writing the Main Header in ASCII

    // FILE* iosba;
    // iosba = fopen(fileName,"w");
    // fflush(iosba);
    // fclose(iosba);
    ((CStructBookParm*)structBook[0])[0].saveMainHeaderASCII(fileName,initBlock,finalBlock,dataSignal);

    if (genData->getFlagOMP() == 1)
    {
        ((CStructBookParm*)structBook[0])[0].saveMainHeaderASCII(fileName2,initBlock,finalBlock,dataSignal);
    }
    // fprintf(iosba,"Sign. Type :          %5i\n", dataSignal->getType());
    // // fprintf(iosba,"Dict. Type :          %5i\n", genData->getDicType());
    // fprintf(iosba,"No. Signals:          %5i\n", dataSignal->getNumSignal());
    // fprintf(iosba,"Signal Size:       %8i\n", dataSignal->getSignalSize());
    // fprintf(iosba,"Block Hop:            %5i\n", dataSignal->getBlockHop());
    // fprintf(iosba,"Block Size:           %5i\n", dataSignal->getBlockSize());
    // if (genData->getSigType()==1)
    // {
    //     fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    // }
    // if (genData->getSigType()==2)
    // {
    //     fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    // }
    // fprintf(iosba,"Init. Block:          %5i\n", initBlock);
    // fprintf(iosba,"Final Block:          %5i\n", finalBlock);
    // fflush(iosba);
    // fclose(iosba);

    // Writing the Main Header in Binary



    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d_header.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    // iosbb = fopen(sbbFName,"wb");
    // fclose(iosbb);
    ((CStructBookParm*)structBook[0])[0].saveMainHeaderBin(sbbFName,initBlock,finalBlock,dataSignal,dicData,genData);
    // int dummyint;
    // double dummydouble;
    // dummyint = dataSignal->getType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = genData->getDicType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getNumSignal();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getSignalSize();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getBlockHop();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getBlockSize();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummydouble = ((CAudioSignal*)dataSignal)->getSamplingRate();
    // fwrite(&dummydouble, sizeof(double), 1, iosbb);
    // dummyint = initBlock;
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = finalBlock;
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // fclose(iosbb);

    // FILE* iosbb;
    // ((CStructBookParm*)structBook)->newFileBin(initBlock,finalBlock,dataSignal);
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);
    // iosbb = fopen(sbbFName,"wb");
    // fclose(iosbb);




    // cout << "Dentro do DECOMP" << endl;
    // dicData->printToScreen();
    ofstream file_stage;
    if (genData->getPrintDecompStage()==1)
    {
        //file pointers
        file_stage.open("decomp_stages.out",ios::out);
        // file header
        file_stage  <<  setw (10) << setfill(' ') << "Signal" << " "
                    <<  setw (10) << setfill(' ') << "Block" << " "
                    <<  setw (10) << setfill(' ') << "No." << " "
                    <<  setw (10) << setfill(' ') << "Stage" << " "
                    <<  setw (20) << setfill(' ') << "Coef." << " "
                    <<  setw (20) << setfill(' ') << "Decay" << " "
                    <<  setw (20) << setfill(' ') << "Freq" << " "
                    <<  setw (20) << setfill(' ') << "Phase"<< " "
                    <<  setw (10) << setfill(' ') << "Ti"<< " "
                    <<  setw (10) << setfill(' ') << "Tf"<< " "
                    << endl;
        file_stage.close();
    }

    //int iSignal;
    //double sigNorm;
    //int iBlock;

    for (i=0; i < dataSignal->getNumSignal(); i++ )
    {
        // Writing the Signal Header]
        // iosba = fopen(fileName,"a");

        ((CStructBookParm*)structBook[0])[0].saveSignalHeaderASCII(fileName,i,dataSignal);
        if (genData->getFlagOMP() == 1)
        {
            ((CStructBookParm*)structBook[0])[0].saveSignalHeaderASCII(fileName2,i,dataSignal);
        }
        // fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        // fprintf(iosba,"Signal:               %5i\n",i+1);
        // fprintf(iosba,"Norm:            %10.5f\n",dataSignal->getNorm(i));
        // fflush(iosba);
        // fclose(iosba);


        // iosbb = fopen(sbbFName,"ab");
        ((CStructBookParm*)structBook[0])[0].saveSignalHeaderBin(sbbFName,i,dataSignal);
        // iSignal = i+1;
        // fwrite(&iSignal, sizeof(int), 1, iosbb);
        // sigNorm = dataSignal->getNorm(i);
        // fwrite(&sigNorm, sizeof(double), 1, iosbb);
        // fclose(iosbb);



        if (dataSignal->getNorm(i)!=0.0)
        {
            for (int j=initBlock-1; j < finalBlock; j++ )
            {
                endingFl = 0;

                for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++)
                {
                    dicAtoms[kdic]->innerProduct = 0.0;
                }

                int a0 = -1;
                int b0 = dicSize;

                cgMatrix<double> prevAtoms(nMaxStep,13,0.0);
                cgMatrix<double> b(nMaxStep,1,0.0);
                cgMatrix<double> v(nMaxStep,1,0.0);
                cgMatrix<double> Ai(nMaxStep,nMaxStep,0.0);
                cgMatrix<double> a(1,nMaxStep,0.0);

                streampos oldpos[nMaxStep];

                c1 = clock();
                residue.zeros();
                // Loading signal into the vector
                if ((int)(j*(double)dataSignal->getBlockHop())+dicSize < dataSignal->getSignalSize())
                {
                    // cout << "AQUI" << endl;
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dicSize);
                }
                else
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dataSignal->getSignalSize() - (int)(j*(double)dataSignal->getBlockHop()));
                }
                // Normalizing initial residue (with signal norm)
                residue /= dataSignal->getNorm(i);

                signal = residue;

                norm = residue.norm();
                ((CStructBookParm*)structBook[i])[j].setNorm(norm);
                double initBlockNorm = norm;

                // Writing the Block Header
                // iosba = fopen(fileName,"a");
                ((CStructBookParm*)structBook[0])[0].saveBlockHeaderASCII(fileName,j,initBlockNorm);
                if (genData->getFlagOMP() == 1)
                {
                    ((CStructBookParm*)structBook[0])[0].saveBlockHeaderASCII(fileName2,j,initBlockNorm);
                }
                // fprintf(iosba,"--------------------------------------------------------------\n");
                // fprintf(iosba,"Block:                %5i\n",j+1);
                // fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
                // fprintf(iosba,"No.    Coef.           Decaying        Freq            Phase           Ti   Tf    PrevAtom AppRatio   meanAppRat befSup     aftSup     normRatio  SNR(dB)     \n");
                // fflush(iosba);
                // fclose(iosba);
                // iosbb = fopen(sbbFName,"ab");
                ((CStructBookParm*)structBook[0])[0].saveBlockHeaderBin(sbbFName,j,initBlockNorm);
                // iBlock = j+1;
                // fwrite(&iBlock, sizeof(int), 1, iosbb);
                // fwrite(&initBlockNorm, sizeof(double), 1, iosbb);
                // fclose(iosbb);

                int decomp_stage;
                // Beginning of the decomposition
                if (norm!=0)
                {
                    int nAtomCont = step;
                    cout << "- nAtomCont: " << nAtomCont << endl;
                    do
                    {
                        // cout << nMaxStep << endl;
                        // cout << residue[1] << endl;
                        if (step>=nMaxStep) break;
#ifdef DBG_WATCH_DCMP_STEP
                        // residue.PrintToFile("residue.dat");
#endif
                        if ( genData->getFlagANN() == 1 )
                        {
                            DANNO(  genData,
                                    residue,
                                    dicSize,
                                    chosenDic,
                                    chosenNet,
                                    step,
                                    L,
                                    approxRatio);
                        }
                        // cout << step << endl;
                        // residue.PrintToFile("residue.dat");
                        cout << "##########################################################################" << endl;
                        cout << "->Decomposing Signal: " << i+1 << "; Block: "<< j+1 << "; Atom: "<< step+1 << endl;
                        cout << "##########################################################################" << endl;
                        // Signal projection over a discrete dictionary
                        cout << "Signal projection over a discrete dictionary" << endl;
                        //chosenParm=fastMPKolasa(residue,dataSignal,dicData,dic);
                        // ((CExpDictionary*)dic)->fastMPKolasa(residue,dataSignal,dicData,chosenParm);
                            // ((CExpDictionary*)dic)->fastMPKolasaModified(residue,dataSignal,dicData,chosenParm, step);
                        matchingPursuit(residue,chosenParm,dicSize,dataSignal,dicData,genData,step,chosenDic,dicAtoms,a0,b0,genData->getFlagOMP());
                        // net1_1(1,residue);
                            //((CExpParm*)chosenParm)->printParm2Screen();
                        cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                <<"| a: "<<((strtParameter*)chosenParm)->a
                                <<"| b: "<<((strtParameter*)chosenParm)->b
                                <<"| beta: "<<((strtParameter*)chosenParm)->beta
                                <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                <<"| chosenNet: "<<chosenNet
                                << endl;
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                        // <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                        // <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->u << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->beta << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->dicType << " "
                                        <<  setw (10) << setfill(' ') << chosenNet << " "
                                        << endl;
                        }

                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        //approxRatio[step%L] = fabs(((CExpParm*)chosenParm)->innerProd)/norm;
                        approxRatio[step%L] = fabs(((strtParameter*)chosenParm)->innerProduct)/norm;

                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }

                        // cout << approxRatio[step%L] - approxRatio[step%L - 1] << endl;

                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;
                        cout << "Tol. Approx. Ratio: " << tolAppRatio << endl;
                        // Adjusting parameters
                        // ((CExpDictionary*)dic)->adjustParameters(residue,chosenParm);
                        adjustParameters(residue, chosenParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        //((CExpParm*)chosenParm)->printParm2Screen();
                        cout << "Coef: "<<((strtParameter*)chosenParm)->innerProduct<<"| Rho: "<<((strtParameter*)chosenParm)->rho<<"| Xi: "<<((strtParameter*)chosenParm)->xi<<"| Phase: "<<((strtParameter*)chosenParm)->phase<<"| Tau: "<<((strtParameter*)chosenParm)->u<<"| a: "<<((strtParameter*)chosenParm)->a<<"| b: "<<((strtParameter*)chosenParm)->b<<"| beta: "<<((strtParameter*)chosenParm)->beta<<"| dicType: "<<((strtParameter*)chosenParm)->dicType<< endl;

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;

                        // ((CExpDictionary*)dic)->setRealAtom(chosenParm);
                        // cgRealAtom.fillVector(((CExpDictionary*)dic)->getRealAtom());
                        if (genData->getFlagOMP() == 0)
                        {
                            updateResidue(residue,norm,dicSize,chosenParm);
                        }
                        if (genData->getFlagOMP() == 1)
                        {
                            updateResidue(  residue,
                                            signal,
                                            b,
                                            v,
                                            Ai,
                                            a,
                                            norm,
                                            dicSize,
                                            step,
                                            prevAtoms,
                                            chosenParm);
                        }

                        // cout << residue.getData(0,512)<< endl;
                        // cgRealAtom.fillVector(realAtom);

                        //cgRealAtomAux = cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        // cgRealAtomAux = cgRealAtom*(((strtParameter*)chosenParm)->innerProduct);
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        //residue = residue - cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        // residue = residue - cgRealAtom*(((strtParameter*)chosenParm)->innerProduct);
                        // norm = residue.norm();
                         // Add element to structure book
                        ((CStructBookParm*)structBook[i])[j].addElement(chosenParm);
                        // int indorig = ((CStructBookParm*)structBook[i])[j].getNumElement();
                        // ((CStructBookParm*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        // ((CStructBookParm*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        // ((CStructBookParm*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);
                        int indorig = ((CStructBookParm*)structBook[i])[j].getNumElement();
                        ((CStructBookParm*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookParm*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookParm*)structBook[i])[j].setOrigAtomIndex(indorig-1,-1);


                        // Save elemnt ASCII
                        // iosba = fopen(fileName,"a");
                        ((CStructBookParm*)structBook[i])[j].saveElementASCII(  fileName,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                0,
                                                                                0,
                                                                                (norm/initBlockNorm),
                                                                                chosenNet);
                        // cout << fabs(((strtParameter*)chosenParm)->innerProduct) << endl;
                        if (    meanApproxRatio <= tolAppRatio
                                // || pow((norm/initBlockNorm),2)<=0.05
                                // || (step>=(L+nAtomCont))
                                || (step>=nMaxStep-1)
                                || (20*log10(initBlockNorm/norm)>=snrTarget)
                                || fabs(((strtParameter*)chosenParm)->innerProduct) <= 1e-12)
                            {
                                endingFl = 1;
                            }

                        if (genData->getFlagOMP() == 1)
                        {
                            if (endingFl == 1)
                            {
                                ((CStructBookParm*)structBook[i])[j].saveInnerProdASCII( fileName, fileName2, a, step, j, initBlock, posaux);
                            }
                        }

                        // if (genData->getFlagOMP() == 1)
                        // {
                        //     if (    meanApproxRatio <= tolAppRatio
                        //         // || (step>=(L+nAtomCont))
                        //         || (step>=nMaxStep-1)
                        //         || (20*log10(initBlockNorm/norm)>=snrTarget))
                        //     {
                        //         ((CStructBookParm*)structBook[i])[j].saveInnerProdASCII( fileName, fileName2, a, step, j, initBlock, posaux);
                        //     }
                        // }

                        // fflush(iosba);
                        // fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;
                        cout << "SNR Target: " << snrTarget<< endl;
                        // cout << norm << endl;
                        // iosbb = fopen(sbbFName,"ab");
                        ((CStructBookParm*)structBook[i])[j].saveElementBin(sbbFName);
                        //fclose(iosbb);
                        step++;
                    }
                    // while(      (   (meanApproxRatio > tolAppRatio) ||
                    //                 (step<(L+nAtomCont)) )
                    //             //((norm/initBlockNorm)>1e-8)
                    //             && (step<nMaxStep)
                    //             && (20*log10(initBlockNorm/norm)<snrTarget)
                    //             && fabs(((strtParameter*)chosenParm)->innerProduct) >= 1e-10
                    //             //&& ( fabs(expParm.innerProd) > 1e-8 )
                    //      );
                    while(endingFl == 0);
                }
                else
                {
                    cout << "  ### Block "<< j+1 <<" with null samples ### " << endl;
                    ((CStructBookParm*)structBook[0])[0].saveBlockNullSamplesASCII(fileName);
                    if (genData->getFlagOMP() == 1)
                    {
                        ((CStructBookParm*)structBook[0])[0].saveBlockNullSamplesASCII(fileName2);
                    }
                    // iosba = fopen(fileName,"a");
                    // fprintf(iosba,"###### Block with null samples ######\n");
                    // fflush(iosba);
                    // fclose(iosba);
                }

                ((CStructBookParm*)structBook[0])[0].saveBlockEndingASCII(fileName);
                if (genData->getFlagOMP() == 1)
                {
                    ((CStructBookParm*)structBook[0])[0].saveBlockEndingASCII(fileName2);
                }
                // iosba = fopen(fileName,"a");
                // fprintf(iosba,"99999\n");
                // fflush(iosba);
                // fclose(iosba);

                ((CStructBookParm*)structBook[0])[0].saveBlockEndingBin(sbbFName);
                // iosbb = fopen(sbbFName,"ab");
                // dummyint = 99999;
                // fwrite( &dummyint, sizeof( int ), 1, iosbb );
                // fclose(iosbb);
                step = 0;
                c2 =clock();
                tempo = ((float)(c2 - c1))/CLOCKS_PER_SEC;

                stream = fopen("time.out", "a");
                fprintf(stream, "%5d %15.8f \n", j+1, tempo);
                fflush(stream);
                fclose(stream);
            }

            ((CStructBookParm*)structBook[0])[0].saveSignalEndingASCII(fileName);
            if (genData->getFlagOMP() == 1)
            {
                ((CStructBookParm*)structBook[0])[0].saveSignalEndingASCII(fileName2);
            }
            // iosba = fopen(fileName,"a");
            // fprintf(iosba,"88888\n");
            // fflush(iosba);
            // fclose(iosba);

            ((CStructBookParm*)structBook[0])[0].saveSignalEndingBin(sbbFName);
            // iosbb = fopen(sbbFName,"ab");
            // dummyint = 88888;
            // fwrite( &dummyint, sizeof( int ), 1, iosbb );
            // fclose(iosbb);
        }
        else
        {
                cout << "  ### Signal "<< i+1 <<" with null samples ### " << endl;
                ((CStructBookParm*)structBook[0])[0].saveSignalNullSamplesASCII(fileName);
                if (genData->getFlagOMP() == 1)
                {
                    ((CStructBookParm*)structBook[0])[0].saveSignalNullSamplesASCII(fileName2);
                }
                // iosba = fopen(fileName,"a");
                // fprintf(iosba,"###### Signal with null samples ######\n");
                // fflush(iosba);
                // fclose(iosba);
        }

        ((CStructBookParm*)structBook[0])[0].saveDecompEndingASCII(fileName);
        if (genData->getFlagOMP() == 1)
        {
            ((CStructBookParm*)structBook[0])[0].saveDecompEndingASCII(fileName2);
        }
        // iosba = fopen(fileName,"a");
        // fprintf(iosba,"77777\n");
        // fflush(iosba);
        // fclose(iosba);

        ((CStructBookParm*)structBook[0])[0].saveDecompEndingBin(sbbFName);
        /*iosbb = fopen(sbbFName,"ab");
        dummyint = 77777;
        fwrite( &dummyint, sizeof( int ), 1, iosbb );
        fclose(iosbb);*/

        //delete [] sbPrevProj;

    }

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }


    // Deallocating objects
    // delete dic;
    for (i=0; i<dataSignal->getNumSignal(); i++) delete [] ((CStructBookParm*)structBook[i]);
    delete[] structBook;
    if (chosenParm!=NULL) delete chosenParm;
    if (approxRatio!=NULL) delete[] approxRatio;
    if ((CAudioSignal*)dataSignal!=NULL) delete (CAudioSignal*)dataSignal;

    for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++) delete [] dicAtoms[kdic];
    delete[] dicAtoms;

    // if ( genData->getFlagANN() == 1 )
    // {
    //     libannformpTerminate();
    //     mclTerminateApplication();
    // }
    return;
}

void decompNoise(   CFileDecomp* genData,
                    CFileDecompBlockRange* blockRange,
                    CFileDictionary* dicData)
{
    //////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
    dataSignal = new CNoiseSignal;

    dataSignal->setNumSignal(1000);
    dataSignal->setSignalSize(blockRange->getBlockSize());
    dataSignal->setSamplingRate(44100);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBook** structBook;
    structBook = new CStructBook* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookExp[numBlock];
    }

    // Setting Dictionary
    CExpDictionary* expDic = new CExpDictionary;
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int sigSize = (int) pow(2.0,(double)nbits);
    expDic->setSignalSize(sigSize);

    // Decomp Stage



    // Deallocating objects
    delete expDic;
    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)structBook[i]);
    }
    delete [] structBook;

    delete (CNoiseSignal*)dataSignal;

    return;
}

/*
void decompECG( CFileDecomp* genData,
                CFileDecompBlockRange* blockRange,
                CFileDictionary* dicData,
                char* InputFile)
{
    //////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
    dataSignal = new CECGSignal;
    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    int chosenDic;
    int chosenNet;

    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBookParm** structBook;
    structBook = new CStructBookParm* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookParm[numBlock];
    }

    clock_t c2, c1; // variáveis que contam ciclos do processador
    float tempo;

    FILE* stream;
    stream = fopen("time.out", "w");
    fprintf(stream, "bloco           tempo \n");
    fflush(stream);
    fclose(stream);

    // Setting Dictionary
    // CDictionary* dic;
    // dic = new CExpDictionary;

    // if (genData->getDicType()==1)
    // {
        // dic = new CExpDictionary;
    // }
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int dicSize = (int) pow(2.0,(double)nbits);

    // ((CDictionary*)dic)->setSignalSize(dicSize);


    // Decomp Stage
    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();

    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }

    int nMaxStep = genData->getNumMaxStep();


    cgMatrix<double> residue(1,dicSize,0.0);
    // cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    // cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);
    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;


    // CParameter* chosenParm;
    // if (genData->getDicType()==1)
    // {
    //     chosenParm = new CExpParm;
    // }

    strtParameter* chosenParm;
    chosenParm = new strtParameter;


    int L = (int)ceil(((log10((double)(dicSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];
    int k;
    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }

    double meanApproxRatio;

    double tolAppRatio =  genData->getApproxRatioTarget();
    double snrTarget = genData->getSNRTarget();

    double befSupInnerP,aftSupInnerP;

    strtParameter** dicAtoms = new strtParameter*[(dicData->getNumDicBlock()*dicSize)];
    for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++)
    {
        dicAtoms[kdic] = new strtParameter;

    }

    char aux[_MAX_PATH];

    // ((CStructBookParm*)structBook) -> newFileASCII(initBlock,finalBlock,dataSignal,dicData);
    char fileName[_MAX_PATH];
    strcpy(fileName, dataSignal->getFileName());
    char* pos;
    pos = strrchr( fileName, '.');
    sprintf(aux,"_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos[0], aux);
    // double* realAtom;

    // Writing the Main Header in ASCII

    // FILE* iosba;
    // iosba = fopen(fileName,"w");
    // fflush(iosba);
    // fclose(iosba);
    ((CStructBookParm*)structBook[0])[0].saveMainHeaderASCII(fileName,initBlock,finalBlock,dataSignal);
    // fprintf(iosba,"Sign. Type :          %5i\n", dataSignal->getType());
    // // fprintf(iosba,"Dict. Type :          %5i\n", genData->getDicType());
    // fprintf(iosba,"No. Signals:          %5i\n", dataSignal->getNumSignal());
    // fprintf(iosba,"Signal Size:       %8i\n", dataSignal->getSignalSize());
    // fprintf(iosba,"Block Hop:            %5i\n", dataSignal->getBlockHop());
    // fprintf(iosba,"Block Size:           %5i\n", dataSignal->getBlockSize());
    // if (genData->getSigType()==1)
    // {
    //     fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    // }
    // if (genData->getSigType()==2)
    // {
    //     fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    // }
    // fprintf(iosba,"Init. Block:          %5i\n", initBlock);
    // fprintf(iosba,"Final Block:          %5i\n", finalBlock);
    // fflush(iosba);
    // fclose(iosba);

    // Writing the Main Header in Binary



    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d_header.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    // iosbb = fopen(sbbFName,"wb");
    // fclose(iosbb);
    ((CStructBookParm*)structBook[0])[0].saveMainHeaderBin(sbbFName,initBlock,finalBlock,dataSignal,dicData,genData);
    // int dummyint;
    // double dummydouble;
    // dummyint = dataSignal->getType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = genData->getDicType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getNumSignal();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getSignalSize();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getBlockHop();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getBlockSize();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummydouble = ((CAudioSignal*)dataSignal)->getSamplingRate();
    // fwrite(&dummydouble, sizeof(double), 1, iosbb);
    // dummyint = initBlock;
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = finalBlock;
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // fclose(iosbb);

    // FILE* iosbb;
    // ((CStructBookParm*)structBook)->newFileBin(initBlock,finalBlock,dataSignal);
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);
    // iosbb = fopen(sbbFName,"wb");
    // fclose(iosbb);




    // cout << "Dentro do DECOMP" << endl;
    // dicData->printToScreen();
    ofstream file_stage;
    if (genData->getPrintDecompStage()==1)
    {
        //file pointers
        file_stage.open("decomp_stages.out",ios::out);
        // file header
        file_stage  <<  setw (10) << setfill(' ') << "Signal" << " "
                    <<  setw (10) << setfill(' ') << "Block" << " "
                    <<  setw (10) << setfill(' ') << "No." << " "
                    <<  setw (10) << setfill(' ') << "Stage" << " "
                    <<  setw (20) << setfill(' ') << "Coef." << " "
                    <<  setw (20) << setfill(' ') << "Decay" << " "
                    <<  setw (20) << setfill(' ') << "Freq" << " "
                    <<  setw (20) << setfill(' ') << "Phase"<< " "
                    <<  setw (10) << setfill(' ') << "Ti"<< " "
                    <<  setw (10) << setfill(' ') << "Tf"<< " "
                    << endl;
        file_stage.close();
    }

    //int iSignal;
    //double sigNorm;
    //int iBlock;

    for (i=0; i < dataSignal->getNumSignal(); i++ )
    {
        // Writing the Signal Header]
        // iosba = fopen(fileName,"a");

        ((CStructBookParm*)structBook[0])[0].saveSignalHeaderASCII(fileName,i,dataSignal);
        // fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        // fprintf(iosba,"Signal:               %5i\n",i+1);
        // fprintf(iosba,"Norm:            %10.5f\n",dataSignal->getNorm(i));
        // fflush(iosba);
        // fclose(iosba);


        // iosbb = fopen(sbbFName,"ab");
        ((CStructBookParm*)structBook[0])[0].saveSignalHeaderBin(sbbFName,i,dataSignal);
        // iSignal = i+1;
        // fwrite(&iSignal, sizeof(int), 1, iosbb);
        // sigNorm = dataSignal->getNorm(i);
        // fwrite(&sigNorm, sizeof(double), 1, iosbb);
        // fclose(iosbb);

        if (dataSignal->getNorm(i)!=0.0)
        {
            for (int j=initBlock-1; j < finalBlock; j++ )
            {
                int a0 = -1;
                int b0 = dicSize;

                for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++)
                {
                    dicAtoms[kdic]->innerProduct = 0.0;
                }
                c1 = clock();
                residue.zeros();
                // Loading signal into the vector
                if ((int)(j*(double)dataSignal->getBlockHop())+dicSize < dataSignal->getSignalSize())
                {
                    // cout << "AQUI" << endl;
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dicSize);
                }
                else
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dataSignal->getSignalSize() - (int)(j*(double)dataSignal->getBlockHop()));
                }
                // Normalizing initial residue (with signal norm)
                residue /= dataSignal->getNorm(i);
                norm = residue.norm();
                ((CStructBookParm*)structBook[i])[j].setNorm(norm);
                double initBlockNorm = norm;

                // Writing the Block Header
                // iosba = fopen(fileName,"a");
                ((CStructBookParm*)structBook[0])[0].saveBlockHeaderASCII(fileName,j,initBlockNorm);
                // fprintf(iosba,"--------------------------------------------------------------\n");
                // fprintf(iosba,"Block:                %5i\n",j+1);
                // fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
                // fprintf(iosba,"No.    Coef.           Decaying        Freq            Phase           Ti   Tf    PrevAtom AppRatio   meanAppRat befSup     aftSup     normRatio  SNR(dB)     \n");
                // fflush(iosba);
                // fclose(iosba);
                // iosbb = fopen(sbbFName,"ab");
                ((CStructBookParm*)structBook[0])[0].saveBlockHeaderBin(sbbFName,j,initBlockNorm);
                // iBlock = j+1;
                // fwrite(&iBlock, sizeof(int), 1, iosbb);
                // fwrite(&initBlockNorm, sizeof(double), 1, iosbb);
                // fclose(iosbb);
                int decomp_stage;
                // Beginning of the decomposition
                if (norm!=0)
                {

                    int nAtomCont = step;
                    cout << "- nAtomCont: " << nAtomCont << endl;
                    do
                    {
                        // cout << residue[1] << endl;
                        if (step>=nMaxStep) break;
#ifdef DBG_WATCH_DCMP_STEP
                        residue.PrintToFile("residue.dat");
#endif
                        // residue.PrintToFile("residue.dat");
                        cout << "##########################################################################" << endl;
                        cout << "->Decomposing Signal: " << i+1 << "; Block: "<< j+1 << "; Atom: "<< step+1 << endl;
                        cout << "##########################################################################" << endl;
                        // Signal projection over a discrete dictionary
                        cout << "Signal projection over a discrete dictionary" << endl;
                        //chosenParm=fastMPKolasa(residue,dataSignal,dicData,dic);
                        // ((CExpDictionary*)dic)->fastMPKolasa(residue,dataSignal,dicData,chosenParm);
                        // ((CExpDictionary*)dic)->fastMPKolasaModified(residue,dataSignal,dicData,chosenParm, step);
                        matchingPursuit(residue,chosenParm,dicSize,dataSignal,dicData,genData,step,chosenDic,dicAtoms,a0,b0,genData->getFlagOMP());
                        //((CExpParm*)chosenParm)->printParm2Screen();
                        cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                <<"| Rho: "<<((strtParameter*)chosenParm)->rho
                                <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                <<"| a: "<<((strtParameter*)chosenParm)->a
                                <<"| b: "<<((strtParameter*)chosenParm)->b
                                <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                << endl;
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                        // <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                        // <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->u << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                        << endl;
                        }

                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        //approxRatio[step%L] = fabs(((CExpParm*)chosenParm)->innerProd)/norm;
                        approxRatio[step%L] = fabs(((strtParameter*)chosenParm)->innerProduct)/norm;
                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }
                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;
                        cout << "Tol. Approx. Ratio: " << tolAppRatio << endl;
                        // Adjusting parameters
                        // ((CExpDictionary*)dic)->adjustParameters(residue,chosenParm);
                        adjustParameters(residue, chosenParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        //((CExpParm*)chosenParm)->printParm2Screen();
                        cout << "Coef: "<<((strtParameter*)chosenParm)->innerProduct<<"| Rho: "<<((strtParameter*)chosenParm)->rho<<"| Xi: "<<((strtParameter*)chosenParm)->xi<<"| Phase: "<<((strtParameter*)chosenParm)->phase<<"| Tau: "<<((strtParameter*)chosenParm)->u<<"| a: "<<((strtParameter*)chosenParm)->a<<"| b: "<<((strtParameter*)chosenParm)->b<<"| dicType: "<<((strtParameter*)chosenParm)->dicType<< endl;

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;

                        // ((CExpDictionary*)dic)->setRealAtom(chosenParm);
                        // cgRealAtom.fillVector(((CExpDictionary*)dic)->getRealAtom());

                        updateResidue(residue,norm,dicSize,chosenParm);
                        // cgRealAtom.fillVector(realAtom);

                        //cgRealAtomAux = cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        // cgRealAtomAux = cgRealAtom*(((strtParameter*)chosenParm)->innerProduct);
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        //residue = residue - cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        // residue = residue - cgRealAtom*(((strtParameter*)chosenParm)->innerProduct);
                        // norm = residue.norm();
                         // Add element to structure book
                        ((CStructBookParm*)structBook[i])[j].addElement(chosenParm);
                        // int indorig = ((CStructBookParm*)structBook[i])[j].getNumElement();
                        // ((CStructBookParm*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        // ((CStructBookParm*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        // ((CStructBookParm*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);
                        int indorig = ((CStructBookParm*)structBook[i])[j].getNumElement();
                        ((CStructBookParm*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookParm*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookParm*)structBook[i])[j].setOrigAtomIndex(indorig-1,-1);


                        // Save elemnt ASCII
                        // iosba = fopen(fileName,"a");
                        ((CStructBookParm*)structBook[i])[j].saveElementASCII(  fileName,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                0,
                                                                                0,
                                                                                (norm/initBlockNorm),
                                                                                chosenNet);
                        // fflush(iosba);
                        // fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;
                        cout << "SNR Target: " << snrTarget<< endl;
                        // iosbb = fopen(sbbFName,"ab");
                        ((CStructBookParm*)structBook[i])[j].saveElementBin(sbbFName);
                        //fclose(iosbb);
                        step++;
                    }
                    while(      (   (meanApproxRatio > tolAppRatio) ||
                                    (step<(L+nAtomCont)) )
                                //((norm/initBlockNorm)>1e-8)
                                && (step<nMaxStep)
                                && (20*log10(initBlockNorm/norm)<snrTarget)
                                //&& ( fabs(expParm.innerProd) > 1e-8 )
                         );

                }
                else
                {
                    cout << "  ### Block "<< j+1 <<" with null samples ### " << endl;
                    ((CStructBookParm*)structBook[0])[0].saveBlockNullSamplesASCII(fileName);
                    // iosba = fopen(fileName,"a");
                    // fprintf(iosba,"###### Block with null samples ######\n");
                    // fflush(iosba);
                    // fclose(iosba);
                }

                ((CStructBookParm*)structBook[0])[0].saveBlockEndingASCII(fileName);
                // iosba = fopen(fileName,"a");
                // fprintf(iosba,"99999\n");
                // fflush(iosba);
                // fclose(iosba);

                ((CStructBookParm*)structBook[0])[0].saveBlockEndingBin(sbbFName);
                // iosbb = fopen(sbbFName,"ab");
                // dummyint = 99999;
                // fwrite( &dummyint, sizeof( int ), 1, iosbb );
                // fclose(iosbb);
                step = 0;
                c2 =clock();
                tempo = ((float)(c2 - c1))/CLOCKS_PER_SEC;

                stream = fopen("time.out", "a");
                fprintf(stream, "%5d %15.8f \n", j+1, tempo);
                fflush(stream);
                fclose(stream);
            }

            ((CStructBookParm*)structBook[0])[0].saveSignalEndingASCII(fileName);
            // iosba = fopen(fileName,"a");
            // fprintf(iosba,"88888\n");
            // fflush(iosba);
            // fclose(iosba);

            ((CStructBookParm*)structBook[0])[0].saveSignalEndingBin(sbbFName);
            // iosbb = fopen(sbbFName,"ab");
            // dummyint = 88888;
            // fwrite( &dummyint, sizeof( int ), 1, iosbb );
            // fclose(iosbb);
        }
        else
        {
                cout << "  ### Signal "<< i+1 <<" with null samples ### " << endl;
                ((CStructBookParm*)structBook[0])[0].saveSignalNullSamplesASCII(fileName);
                // iosba = fopen(fileName,"a");
                // fprintf(iosba,"###### Signal with null samples ######\n");
                // fflush(iosba);
                // fclose(iosba);
        }

        ((CStructBookParm*)structBook[0])[0].saveDecompEndingASCII(fileName);
        // iosba = fopen(fileName,"a");
        // fprintf(iosba,"77777\n");
        // fflush(iosba);
        // fclose(iosba);

        ((CStructBookParm*)structBook[0])[0].saveDecompEndingBin(sbbFName);
        //iosbb = fopen(sbbFName,"ab");
        //dummyint = 77777;
        //fwrite( &dummyint, sizeof( int ), 1, iosbb );
        //fclose(iosbb);

        //delete [] sbPrevProj;

    }

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }





    // Deallocating objects
    // delete dic;
    for (i=0; i<dataSignal->getNumSignal(); i++) delete [] ((CStructBookParm*)structBook[i]);
    delete [] structBook;
    if (chosenParm!=NULL) delete chosenParm;
    if (approxRatio!=NULL) delete[] approxRatio;
    if ((CECGSignal*)dataSignal!=NULL) delete (CECGSignal*)dataSignal;

    for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++) delete [] dicAtoms[kdic];
    delete[] dicAtoms;

    return;
}
*/
void decompEDA(   CFileDecomp* genData,
                    CFileDecompBlockRange* blockRange,
                    CFileDictionary* dicData,
                    char* InputFile)
{
    int chosenDic;
    int chosenNet = 0;
    int endingFl;
    // if ( genData->getFlagANN() == 1 )
    // {
    //     const char *args[] = { "-nodisplay" };
    //     mclInitializeApplication(args,1);
    //     libannformpInitialize();
    // }

    chosenDic = genData->getDicType();


    //////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
    dataSignal = new CAudioSignal;

    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBookParm** structBook;
    structBook = new CStructBookParm* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookParm[numBlock];
    }

    clock_t c2, c1; /* variáveis que contam ciclos do processador */
    float tempo;
    FILE* stream;
    stream = fopen("time.out", "w");
    fprintf(stream, "bloco           tempo \n");
    fflush(stream);
    fclose(stream);

    // Setting Dictionary
    // CDictionary* dic;
    // dic = new CExpDictionary;

    // if (genData->getDicType()==1)
    // {
        // dic = new CExpDictionary;
    // }
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int dicSize = (int) pow(2.0,(double)nbits);

    // ((CDictionary*)dic)->setSignalSize(dicSize);

    // Decomp Stage
    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();

    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }

    int nMaxStep = genData->getNumMaxStep();

    cgMatrix<double> residue(1,dicSize,0.0);
    cgMatrix<double> signal(1,dicSize,0.0);
    // mwArray residueaux;
    // cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    // cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);

    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;


    // CParameter* chosenParm;
    // if (genData->getDicType()==1)
    // {
    //     chosenParm = new CExpParm;
    // }

    strtParameter* chosenParm;
    chosenParm = new strtParameter;


    int L = (int)ceil(((log10((double)(dicSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];
    int k;
    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }

    double meanApproxRatio;

    double tolAppRatio =  genData->getApproxRatioTarget();
    double snrTarget = genData->getSNRTarget();

    double befSupInnerP,aftSupInnerP;

    strtParameter** dicAtoms = new strtParameter*[(dicData->getNumDicBlock()*dicSize)];
    for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++)
    {
        dicAtoms[kdic] = new strtParameter;
    }

    char aux[_MAX_PATH];

    // ((CStructBookParm*)structBook) -> newFileASCII(initBlock,finalBlock,dataSignal,dicData);
    char fileName[_MAX_PATH];
    strcpy(fileName, dataSignal->getFileName());
    char* pos;
    pos = strrchr( fileName, '.');
    sprintf(aux,"_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    char aux2[_MAX_PATH];
    char fileName2[_MAX_PATH];
    strcpy(fileName2, dataSignal->getFileName());
    char* pos2;
    pos2 = strrchr( fileName2, '.');
    sprintf(aux2,"_OMP_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos2[0], aux2);


    long posaux;

    // double* realAtom;

    // Writing the Main Header in ASCII

    // FILE* iosba;
    // iosba = fopen(fileName,"w");
    // fflush(iosba);
    // fclose(iosba);
    ((CStructBookParm*)structBook[0])[0].saveMainHeaderASCII(fileName,initBlock,finalBlock,dataSignal);

    if (genData->getFlagOMP() == 1)
    {
        ((CStructBookParm*)structBook[0])[0].saveMainHeaderASCII(fileName2,initBlock,finalBlock,dataSignal);
    }
    // fprintf(iosba,"Sign. Type :          %5i\n", dataSignal->getType());
    // // fprintf(iosba,"Dict. Type :          %5i\n", genData->getDicType());
    // fprintf(iosba,"No. Signals:          %5i\n", dataSignal->getNumSignal());
    // fprintf(iosba,"Signal Size:       %8i\n", dataSignal->getSignalSize());
    // fprintf(iosba,"Block Hop:            %5i\n", dataSignal->getBlockHop());
    // fprintf(iosba,"Block Size:           %5i\n", dataSignal->getBlockSize());
    // if (genData->getSigType()==1)
    // {
    //     fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    // }
    // if (genData->getSigType()==2)
    // {
    //     fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    // }
    // fprintf(iosba,"Init. Block:          %5i\n", initBlock);
    // fprintf(iosba,"Final Block:          %5i\n", finalBlock);
    // fflush(iosba);
    // fclose(iosba);

    // Writing the Main Header in Binary



    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d_header.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    // iosbb = fopen(sbbFName,"wb");
    // fclose(iosbb);
    ((CStructBookParm*)structBook[0])[0].saveMainHeaderBin(sbbFName,initBlock,finalBlock,dataSignal,dicData,genData);
    // int dummyint;
    // double dummydouble;
    // dummyint = dataSignal->getType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = genData->getDicType();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getNumSignal();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getSignalSize();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getBlockHop();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = dataSignal->getBlockSize();
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummydouble = ((CAudioSignal*)dataSignal)->getSamplingRate();
    // fwrite(&dummydouble, sizeof(double), 1, iosbb);
    // dummyint = initBlock;
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // dummyint = finalBlock;
    // fwrite(&dummyint, sizeof(int), 1, iosbb);
    // fclose(iosbb);

    // FILE* iosbb;
    // ((CStructBookParm*)structBook)->newFileBin(initBlock,finalBlock,dataSignal);
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);
    // iosbb = fopen(sbbFName,"wb");
    // fclose(iosbb);




    // cout << "Dentro do DECOMP" << endl;
    // dicData->printToScreen();
    ofstream file_stage;
    if (genData->getPrintDecompStage()==1)
    {
        //file pointers
        file_stage.open("decomp_stages.out",ios::out);
        // file header
        file_stage  <<  setw (10) << setfill(' ') << "Signal" << " "
                    <<  setw (10) << setfill(' ') << "Block" << " "
                    <<  setw (10) << setfill(' ') << "No." << " "
                    <<  setw (10) << setfill(' ') << "Stage" << " "
                    <<  setw (20) << setfill(' ') << "Coef." << " "
                    <<  setw (20) << setfill(' ') << "Decay" << " "
                    <<  setw (20) << setfill(' ') << "Freq" << " "
                    <<  setw (20) << setfill(' ') << "Phase"<< " "
                    <<  setw (10) << setfill(' ') << "Ti"<< " "
                    <<  setw (10) << setfill(' ') << "Tf"<< " "
                    << endl;
        file_stage.close();
    }

    //int iSignal;
    //double sigNorm;
    //int iBlock;

    for (i=0; i < dataSignal->getNumSignal(); i++ )
    {
        // Writing the Signal Header]
        // iosba = fopen(fileName,"a");

        ((CStructBookParm*)structBook[0])[0].saveSignalHeaderASCII(fileName,i,dataSignal);
        if (genData->getFlagOMP() == 1)
        {
            ((CStructBookParm*)structBook[0])[0].saveSignalHeaderASCII(fileName2,i,dataSignal);
        }
        // fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        // fprintf(iosba,"Signal:               %5i\n",i+1);
        // fprintf(iosba,"Norm:            %10.5f\n",dataSignal->getNorm(i));
        // fflush(iosba);
        // fclose(iosba);


        // iosbb = fopen(sbbFName,"ab");
        ((CStructBookParm*)structBook[0])[0].saveSignalHeaderBin(sbbFName,i,dataSignal);
        // iSignal = i+1;
        // fwrite(&iSignal, sizeof(int), 1, iosbb);
        // sigNorm = dataSignal->getNorm(i);
        // fwrite(&sigNorm, sizeof(double), 1, iosbb);
        // fclose(iosbb);



        if (dataSignal->getNorm(i)!=0.0)
        {
            for (int j=initBlock-1; j < finalBlock; j++ )
            {
                int a0 = -1;
                int b0 = dicSize;

                for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++)
                {
                    dicAtoms[kdic]->innerProduct = 0.0;
                }

                endingFl = 0;
                cgMatrix<double> prevAtoms(nMaxStep,13,0.0);
                cgMatrix<double> b(nMaxStep,1,0.0);
                cgMatrix<double> v(nMaxStep,1,0.0);
                cgMatrix<double> Ai(nMaxStep,nMaxStep,0.0);
                cgMatrix<double> a(1,nMaxStep,0.0);

                streampos oldpos[nMaxStep];

                c1 = clock();
                residue.zeros();
                // Loading signal into the vector
                if ((int)(j*(double)dataSignal->getBlockHop())+dicSize < dataSignal->getSignalSize())
                {
                    // cout << "AQUI" << endl;
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dicSize);
                }
                else
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dataSignal->getSignalSize() - (int)(j*(double)dataSignal->getBlockHop()));
                }
                // Normalizing initial residue (with signal norm)
                residue /= dataSignal->getNorm(i);

                signal = residue;

                norm = residue.norm();
                ((CStructBookParm*)structBook[i])[j].setNorm(norm);
                double initBlockNorm = norm;

                // Writing the Block Header
                // iosba = fopen(fileName,"a");
                ((CStructBookParm*)structBook[0])[0].saveBlockHeaderASCII(fileName,j,initBlockNorm);
                if (genData->getFlagOMP() == 1)
                {
                    ((CStructBookParm*)structBook[0])[0].saveBlockHeaderASCII(fileName2,j,initBlockNorm);
                }
                // fprintf(iosba,"--------------------------------------------------------------\n");
                // fprintf(iosba,"Block:                %5i\n",j+1);
                // fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
                // fprintf(iosba,"No.    Coef.           Decaying        Freq            Phase           Ti   Tf    PrevAtom AppRatio   meanAppRat befSup     aftSup     normRatio  SNR(dB)     \n");
                // fflush(iosba);
                // fclose(iosba);
                // iosbb = fopen(sbbFName,"ab");
                ((CStructBookParm*)structBook[0])[0].saveBlockHeaderBin(sbbFName,j,initBlockNorm);
                // iBlock = j+1;
                // fwrite(&iBlock, sizeof(int), 1, iosbb);
                // fwrite(&initBlockNorm, sizeof(double), 1, iosbb);
                // fclose(iosbb);

                int decomp_stage;
                // Beginning of the decomposition
                if (norm!=0)
                {
                    int nAtomCont = step;
                    cout << "- nAtomCont: " << nAtomCont << endl;
                    do
                    {
                        // cout << nMaxStep << endl;
                        // cout << residue[1] << endl;
                        if (step>=nMaxStep) break;
#ifdef DBG_WATCH_DCMP_STEP
                        // residue.PrintToFile("residue.dat");
#endif

                        if ( genData->getFlagANN() == 1 )
                        {
                            DANNO(  genData,
                                    residue,
                                    dicSize,
                                    chosenDic,
                                    chosenNet,
                                    step,
                                    L,
                                    approxRatio);
                        }
                        /*if ( genData->getFlagANN() == 1 )
                        {
                            ANN(genData,
                                residue,
                                dicSize,
                                chosenDic,
                                chosenNet,
                                step,
                                L,
                                approxRatio);
                        }*/
                        // cout << step << endl;
                        // residue.PrintToFile("residue.dat");
                        //CDictionary * bateDict;
                       // bateDict = new CBatemanDictionary;

                        cout << "##########################################################################" << endl;
                        cout << "->Decomposing Signal: " << i+1 << "; Block: "<< j+1 << "; Atom: "<< step+1 << endl;
                        cout << "##########################################################################" << endl;
                        // Signal projection over a discrete dictionary
                        cout << "Signal projection over a discrete dictionary" << endl;
                        //chosenParm=fastMPKolasa(residue,dataSignal,dicData,dic);
                        // ((CExpDictionary*)dic)->fastMPKolasa(residue,dataSignal,dicData,chosenParm);
                        // ((CExpDictionary*)dic)->fastMPKolasaModified(residue,dataSignal,dicData,chosenParm, step);
                        matchingPursuitEDA(residue,chosenParm,dicSize,dataSignal,dicData,genData,step,chosenDic,dicAtoms,a0,b0,genData->getFlagOMP());
                        // net1_1(1,residue);
                            //((CBatemanDictionary*)chosenParm)->printParm2Screen();
                        cout    << "Coef: "<<((strtParameter*)chosenParm)->innerProduct
                                <<"| Decay "<<((strtParameter*)chosenParm)->rho
                                <<"| Xi: "<<((strtParameter*)chosenParm)->xi
                                <<"| Phase: "<<((strtParameter*)chosenParm)->phase
                                <<"| Tau: " <<((strtParameter*)chosenParm)->u
                                <<"| a: "<<((strtParameter*)chosenParm)->a
                                <<"| b: "<<((strtParameter*)chosenParm)->b
                                <<"| Rise: "<<((strtParameter*)chosenParm)->beta
                                <<"| dicType: "<<((strtParameter*)chosenParm)->dicType
                                <<"| chosenNet: "<<chosenNet
                                << endl;
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                        // <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                        // <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                        // <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->innerProduct << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->rho << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->xi << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->phase << " "
                                        <<  setw (20) << setfill(' ') << ((strtParameter*)chosenParm)->u << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->a << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->b << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->beta << " "
                                        <<  setw (10) << setfill(' ') << ((strtParameter*)chosenParm)->dicType << " "
                                        <<  setw (10) << setfill(' ') << chosenNet << " "
                                        << endl;
                        }

                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        //approxRatio[step%L] = fabs(((CExpParm*)chosenParm)->innerProd)/norm;
                        approxRatio[step%L] = fabs(((strtParameter*)chosenParm)->innerProduct)/norm;

                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }

                        // cout << approxRatio[step%L] - approxRatio[step%L - 1] << endl;

                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;
                        cout << "Tol. Approx. Ratio: " << tolAppRatio << endl;


                        // Maximize approximation by finding optimum parameters
                        //int dicsd;
                        //dicsd = dataSignal->getSignalSize();
                        //cout<<"vert"<<dicsd<<endl;
                        //dicsd=((CBatemanDictionary*)bateDic)->getSignalSize();
                        //cout<<"vert"<<dicsd<<endl;

                        //((CBatemanDictionary*)bateDict)->setSignalSize(dicSize);

                        //((CBatemanDictionary*)bateDic)->setSignalSize(((CBatemanDictionary*)bateDic)->getSignalSize());
                        //optimizeContinuousParms(residue, chosenParm);
                        //cout<<"*ip "<<((strtParameter*)chosenParm)->innerProduct<<endl;
                        //((CBatemanDictionary*)bateDict)->optimizeContinuousParms(residue, chosenParm);

                        //optimizeContinuousParms(residue, chosenParm);

                        //cout<<"ip* "<<((strtParameter*)chosenParm)->innerProduct<<endl;
                        //cout<<"rho "<<((strtParameter*)chosenParm)->rho<<endl;
                        // Adjusting parameters
                        // ((CExpDictionary*)dic)->adjustParameters(residue,chosenParm);
                        adjustParameters(residue, chosenParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        //((CExpParm*)chosenParm)->printParm2Screen();
                        cout << "Coef: "<<((strtParameter*)chosenParm)->innerProduct<<"| Rho: "<<((strtParameter*)chosenParm)->rho<<"| Xi: "<<((strtParameter*)chosenParm)->xi<<"| Phase: "<<((strtParameter*)chosenParm)->phase<<"| Tau: "<<((strtParameter*)chosenParm)->u<<"| a: "<<((strtParameter*)chosenParm)->a<<"| b: "<<((strtParameter*)chosenParm)->b<<"| beta: "<<((strtParameter*)chosenParm)->beta<<"| dicType: "<<((strtParameter*)chosenParm)->dicType<< endl;

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;

                        // ((CExpDictionary*)dic)->setRealAtom(chosenParm);
                        // cgRealAtom.fillVector(((CExpDictionary*)dic)->getRealAtom());
                        if (genData->getFlagOMP() == 0)
                        {
                            updateResidue(residue,norm,dicSize,chosenParm);
                        }
                        if (genData->getFlagOMP() == 1)
                        {
                            updateResidue(  residue,
                                            signal,
                                            b,
                                            v,
                                            Ai,
                                            a,
                                            norm,
                                            dicSize,
                                            step,
                                            prevAtoms,
                                            chosenParm);


                        }

                        // cout << residue.getData(0,512)<< endl;
                        // cgRealAtom.fillVector(realAtom);
//
                        //cgRealAtomAux = cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        // cgRealAtomAux = cgRealAtom*(((strtParameter*)chosenParm)->innerProduct);
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        //residue = residue - cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        // residue = residue - cgRealAtom*(((strtParameter*)chosenParm)->innerProduct);
                        // norm = residue.norm();
                         // Add element to structure book
                        ((CStructBookParm*)structBook[i])[j].addElement(chosenParm);
                        // int indorig = ((CStructBookParm*)structBook[i])[j].getNumElement();
                        // ((CStructBookParm*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        // ((CStructBookParm*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        // ((CStructBookParm*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);
                        int indorig = ((CStructBookParm*)structBook[i])[j].getNumElement();
                        ((CStructBookParm*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookParm*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookParm*)structBook[i])[j].setOrigAtomIndex(indorig-1,-1);


                        // Save elemnt ASCII
                        // iosba = fopen(fileName,"a");
                        ((CStructBookParm*)structBook[i])[j].saveElementASCII(  fileName,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                0,
                                                                                0,
                                                                                (norm/initBlockNorm),
                                                                                chosenNet);
                        // cout << fabs(((strtParameter*)chosenParm)->innerProduct) << endl;
                        cout << pow((norm/initBlockNorm),2) << endl;
                        if (    meanApproxRatio <= tolAppRatio
                                // || pow((norm/initBlockNorm),2)<=0.05
                                // || (step>=(L+nAtomCont))
                                || (step>=nMaxStep-1)
                                || (20*log10(initBlockNorm/norm)>=snrTarget)
                                || fabs(((strtParameter*)chosenParm)->innerProduct) <= 1e-12)
                            {
                                endingFl = 1;
                            }

                        if (genData->getFlagOMP() == 1)
                        {
                            if (endingFl == 1)
                            {
                                ((CStructBookParm*)structBook[i])[j].saveInnerProdASCII( fileName, fileName2, a, step, j, initBlock, posaux);
                            }
                        }

                        // if (genData->getFlagOMP() == 1)
                        // {
                        //     if (    meanApproxRatio <= tolAppRatio
                        //         // || (step>=(L+nAtomCont))
                        //         || (step>=nMaxStep-1)
                        //         || (20*log10(initBlockNorm/norm)>=snrTarget))
                        //     {
                        //         ((CStructBookParm*)structBook[i])[j].saveInnerProdASCII( fileName, fileName2, a, step, j, initBlock, posaux);
                        //     }
                        // }

                        // fflush(iosba);
                        // fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;
                        cout << "SNR Target: " << snrTarget<< endl;
                        // cout << norm << endl;
                        // iosbb = fopen(sbbFName,"ab");
                        ((CStructBookParm*)structBook[i])[j].saveElementBin(sbbFName);
                        //fclose(iosbb);
                        step++;
                       // delete bateDict;
                    }
                    // while(      (   (meanApproxRatio > tolAppRatio) ||
                    //                 (step<(L+nAtomCont)) )
                    //             //((norm/initBlockNorm)>1e-8)
                    //             && (step<nMaxStep)
                    //             && (20*log10(initBlockNorm/norm)<snrTarget)
                    //             && fabs(((strtParameter*)chosenParm)->innerProduct) >= 1e-10
                    //             //&& ( fabs(expParm.innerProd) > 1e-8 )
                    //      );
                    while(endingFl == 0);
                }
                else
                {
                    cout << "  ### Block "<< j+1 <<" with null samples ### " << endl;
                    ((CStructBookParm*)structBook[0])[0].saveBlockNullSamplesASCII(fileName);
                    if (genData->getFlagOMP() == 1)
                    {
                        ((CStructBookParm*)structBook[0])[0].saveBlockNullSamplesASCII(fileName2);
                    }
                    // iosba = fopen(fileName,"a");
                    // fprintf(iosba,"###### Block with null samples ######\n");
                    // fflush(iosba);
                    // fclose(iosba);
                }

                ((CStructBookParm*)structBook[0])[0].saveBlockEndingASCII(fileName);
                if (genData->getFlagOMP() == 1)
                {
                    ((CStructBookParm*)structBook[0])[0].saveBlockEndingASCII(fileName2);
                }
                // iosba = fopen(fileName,"a");
                // fprintf(iosba,"99999\n");
                // fflush(iosba);
                // fclose(iosba);

                ((CStructBookParm*)structBook[0])[0].saveBlockEndingBin(sbbFName);
                // iosbb = fopen(sbbFName,"ab");
                // dummyint = 99999;
                // fwrite( &dummyint, sizeof( int ), 1, iosbb );
                // fclose(iosbb);
                step = 0;
                c2 =clock();
                tempo = ((float)(c2 - c1))/CLOCKS_PER_SEC;

                stream = fopen("time.out", "a");
                fprintf(stream, "%5d %15.8f \n", j+1, tempo);
                fflush(stream);
                fclose(stream);
            }

            ((CStructBookParm*)structBook[0])[0].saveSignalEndingASCII(fileName);
            if (genData->getFlagOMP() == 1)
            {
                ((CStructBookParm*)structBook[0])[0].saveSignalEndingASCII(fileName2);
            }
            // iosba = fopen(fileName,"a");
            // fprintf(iosba,"88888\n");
            // fflush(iosba);
            // fclose(iosba);

            ((CStructBookParm*)structBook[0])[0].saveSignalEndingBin(sbbFName);
            // iosbb = fopen(sbbFName,"ab");
            // dummyint = 88888;
            // fwrite( &dummyint, sizeof( int ), 1, iosbb );
            // fclose(iosbb);

        }
        else
        {
                cout << "  ### Signal "<< i+1 <<" with null samples ### " << endl;
                ((CStructBookParm*)structBook[0])[0].saveSignalNullSamplesASCII(fileName);
                if (genData->getFlagOMP() == 1)
                {
                    ((CStructBookParm*)structBook[0])[0].saveSignalNullSamplesASCII(fileName2);
                }
                // iosba = fopen(fileName,"a");
                // fprintf(iosba,"###### Signal with null samples ######\n");
                // fflush(iosba);
                // fclose(iosba);
        }

        ((CStructBookParm*)structBook[0])[0].saveDecompEndingASCII(fileName);
        if (genData->getFlagOMP() == 1)
        {
            ((CStructBookParm*)structBook[0])[0].saveDecompEndingASCII(fileName2);
        }
        // iosba = fopen(fileName,"a");
        // fprintf(iosba,"77777\n");
        // fflush(iosba);
        // fclose(iosba);

        ((CStructBookParm*)structBook[0])[0].saveDecompEndingBin(sbbFName);
        /*iosbb = fopen(sbbFName,"ab");
        dummyint = 77777;
        fwrite( &dummyint, sizeof( int ), 1, iosbb );
        fclose(iosbb);*/

        //delete [] sbPrevProj;

    }

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }


    // Deallocating objects
    // delete dic;
    for (i=0; i<dataSignal->getNumSignal(); i++) delete [] ((CStructBookParm*)structBook[i]);
    delete[] structBook;
    if (chosenParm!=NULL) delete chosenParm;
    if (approxRatio!=NULL) delete[] approxRatio;
    if ((CAudioSignal*)dataSignal!=NULL) delete (CAudioSignal*)dataSignal;
    // if ( genData->getFlagANN() == 1 )
    // {
    //     libannformpTerminate();
    //     mclTerminateApplication();
    // }
    for (int kdic = 0; kdic < (dicData->getNumDicBlock()*dicSize); kdic++) delete [] dicAtoms[kdic];
    delete[] dicAtoms;
    return;
}
