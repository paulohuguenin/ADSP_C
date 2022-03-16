#include "encode.h"

void optimizeRDExp(	CFileRDBitRange* rdBitRange,
					CDataSignal* dataSignal,
					strtSBBHeader sbbHeader,
					CStructBook** structBook,
					int Nfreq,
					char* InputFile)
{
	int i, j, k, t;
    //////////////////////////////////
    // Defining the quantizers range

    int init_nbit_amp, end_nbit_amp, delta_nbit_amp;
    int init_nbit_rho, end_nbit_rho, delta_nbit_rho;
    int init_nbit_phase, end_nbit_phase, delta_nbit_phase;

    init_nbit_amp   =  rdBitRange->getInitNbitAmp();
    end_nbit_amp    =  rdBitRange->getEndNbitAmp();
    delta_nbit_amp  =  rdBitRange->getDeltaNbitAmp();
    init_nbit_rho   =  rdBitRange->getInitNbitRho();
    end_nbit_rho    =  rdBitRange->getEndNbitRho();
    delta_nbit_rho  =  rdBitRange->getDeltaNbitRho();
    init_nbit_phase =  rdBitRange->getInitNbitPhase();
    end_nbit_phase  =  rdBitRange->getEndNbitPhase();
    delta_nbit_phase  = rdBitRange->getDeltaNbitPhase();

    int numQuant=  (((end_nbit_amp - init_nbit_amp)/delta_nbit_amp) +1) *
                   (((end_nbit_rho - init_nbit_rho)/delta_nbit_rho) +1) *
                   (((end_nbit_phase - init_nbit_phase)/delta_nbit_phase) +1)
                   + 1;

    // Setting fixed quantizers
    int nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
    int nb_sample = (int)ceil( log(sbbHeader.subBlockSize) / log(2.0) );
    cout << "nb_xi" << " " << "nb_sample" << endl;
    cout << nb_xi << " " << nb_sample << endl;

    //int nb_block = (int)ceil( log(sbbHeader.numBlock) / log(2.0) );
    // Setting variable quantizers
    strtQuantExp* quantExp;
    quantExp = new strtQuantExp[numQuant];
    quantExp[0].nb_amp = 0;
    quantExp[0].nb_rho = 0;
    quantExp[0].nb_phase = 0;
    quantExp[0].nb_xi = 0;
    quantExp[0].nb_sample = 0;
    t=1;
    for (i=init_nbit_amp;i<=end_nbit_amp;i=i+delta_nbit_amp)
    {
        for (j=init_nbit_rho;j<=end_nbit_rho;j=j+delta_nbit_rho)
        {
            for (k=init_nbit_phase;k<=end_nbit_phase;k=k+delta_nbit_phase)
            {
                quantExp[t].nb_amp = i;
                quantExp[t].nb_rho = j;
                quantExp[t].nb_phase = k;
                quantExp[t].nb_xi = nb_xi;
                quantExp[t].nb_sample = nb_sample;
                t++;
            }
        }
    }

	// variables
	int deltaSupMax;
	strtRD* rd = new strtRD[numQuant];
    strtRDIndex rdIndex;
    strtContinuousExp* pSB;
    int numElement;
	strtContinuousExp* pSBQ;
	int numElementQ;
    double 	min_amp,max_amp,
    		min_rho,max_rho,
    		min_phase,max_phase;
    double* recSignal;
    double* blockSignal;
    recSignal = new double[sbbHeader.blockSize];
    blockSignal = new double[sbbHeader.blockSize];
    double** origSignal;
    origSignal = dataSignal->getSignal();
    //file pointers
    char FName[_MAX_PATH];
    char aux[_MAX_PATH];
    char* pos;
    strcpy(FName, InputFile);
    pos = strrchr( FName, '.');
    sprintf(aux,"_rdpoint.out");
    strcpy( &pos[0], aux);
    fstream file_rdpoint(FName,ios::out);
    strcpy(FName, InputFile);
    pos = strrchr( FName, '.');
    sprintf(aux,"_rdopcurve.out");
    strcpy( &pos[0], aux);
    fstream file_rdopcurve(FName,ios::out);
    strcpy(FName, InputFile);
    pos = strrchr( FName, '.');
    sprintf(aux,"_rdreport.out");
    strcpy( &pos[0], aux);
    fstream file_report(FName,ios::out);
    // file header
    file_rdpoint << setw (10) << setfill(' ') << "Signal" << " "
		            <<  setw (10) << setfill(' ') << "Block" << " "
		            <<  setw (10) << setfill(' ') << "nb_amp" << " "
					<<  setw (10) << setfill(' ') << "nb_rho" << " "
			        <<  setw (10) << setfill(' ') << "nb_phase"<< " "
			        <<  setw (10) << setfill(' ') << "nb_xi"<< " "
			        <<  setw (10) << setfill(' ') << "nb_sample"<< " "
			        <<  setw (10) << setfill(' ') << "nb_dsup"<< " "
			        <<  setw (20) << setfill(' ') << "Rate"<< " "
			        <<  setw (20) << setfill(' ') << "Distortion"<< " "
			        <<  setw (20) << setfill(' ') << "NumElementQ"<< " "
			        << endl;
	file_rdopcurve << setw (10) << setfill(' ') << "Signal"<< " "
		            <<  setw (10) << setfill(' ') << "Block"<< " "
		            <<  setw (10) << setfill(' ') << "nb_amp"<< " "
					<<  setw (10) << setfill(' ') << "nb_rho"<< " "
			        <<  setw (10) << setfill(' ') << "nb_phase"<< " "
			        <<  setw (10) << setfill(' ') << "nb_xi"<< " "
			        <<  setw (10) << setfill(' ') << "nb_sample"<< " "
			        <<  setw (10) << setfill(' ') << "nb_dsup"<< " "
			        <<  setw (20) << setfill(' ') << "Rate"<< " "
			        <<  setw (20) << setfill(' ') << "Distortion"<< " "
                    <<  setw (20) << setfill(' ') << "Lambda"<< " "
                    <<  setw (20) << setfill(' ') << "NumElementQ"<< " "
			        << endl;
	file_report << "=====================" << endl;
	file_report << "Parameters range" << endl;
	file_report << "=====================" << endl;
	file_report << setw (10) << setfill(' ') << "Signal"<< " "
		            <<  setw (10) << setfill(' ') << "Block"<< " "
		            <<  setw (20) << setfill(' ') << "min_amp"<< " "
					<<  setw (20) << setfill(' ') << "max_amp"<< " "
			        <<  setw (20) << setfill(' ') << "min_rho"<< " "
			        <<  setw (20) << setfill(' ') << "max_rho"<< " "
			        <<  setw (20) << setfill(' ') << "min_phase"<< " "
			        <<  setw (20) << setfill(' ') << "max_phase"<< " "
			        <<  setw (20) << setfill(' ') << "deltaSupMax"<< " "
			        << endl;
    cout << "numQuant: " << numQuant << endl;
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		double norm = dataSignal->getNorm(i);
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			cout << "Signal: " << i+1 << " -- " << "Block: " << j+1 << endl;
			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			numElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			pSBQ = new strtContinuousExp[numElement];

			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp,
															  max_amp,
															  min_rho,
															  max_rho,
															  min_phase,
															  max_phase);

			deltaSupMax = ((CStructBookExp*)structBook[i])[j].findDeltaSupMax();

			int nb_deltasup = (int)ceil( log2(deltaSupMax+1) );


			file_report << setw (10) << setfill(' ') << i+1<< " "
		            <<  setw (10) << setfill(' ') << j+1<< " "
		            <<  setw (20) << setfill(' ') << min_amp<< " "
					<<  setw (20) << setfill(' ') << max_amp<< " "
			        <<  setw (20) << setfill(' ') << min_rho<< " "
			        <<  setw (20) << setfill(' ') << max_rho<< " "
			        <<  setw (20) << setfill(' ') << min_phase<< " "
			        <<  setw (20) << setfill(' ') << max_phase<< " "
			        <<  setw (20) << setfill(' ') << deltaSupMax << " "
			        << endl;

			// init vector
			for (k=0;k<sbbHeader.blockSize;k++)
				blockSignal[k] = 0.0;
			//----------------------------
			int initBlock = j*sbbHeader.blockSize;
			int copylen = min(dataSignal->getSignalSize()-initBlock, sbbHeader.blockSize);
			memcpy(blockSignal,&origSignal[i][initBlock],copylen*sizeof(double));

			printf(" *** Obtain RD points\n");
			k=0;
			rd[0].dist = computeSqrNorm(blockSignal,sbbHeader.blockSize);
			rd[0].rate = 0;
			quantExp[k].nb_deltasup = nb_deltasup;
			file_rdpoint << setw (10) << setfill(' ') << i+1<< " "
						<<  setw (10) << setfill(' ') << j+1<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_amp<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_rho<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_phase<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_xi<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_sample<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_deltasup<< " "
						<<  setw (20) << setfill(' ') << rd[k].rate<< " "
						<<  setw (20) << setfill(' ') << rd[k].dist<< " "
						<<  setw (20) << setfill(' ') << 0 << " "
						<< endl;
			for (k=1;k<numQuant;k++)
			{
				quantExp[k].nb_deltasup = nb_deltasup;
				cout << "Quantizer: " << k+1 << endl;
				printf(" ****** Quantize Structure Book \n");
				quantizeStructBookExp(min_amp,max_amp,
									  min_rho,max_rho,
									  min_phase,max_phase,
									  quantExp[k],
									  pSB,numElement,
									  pSBQ,numElementQ); // output elements

				// fstream file_quantstrbook("quantstrbook.out",ios::out);
// 				file_quantstrbook << "===============================================" << endl;
// 				file_quantstrbook << setw (10) << setfill(' ') << "Signal"<< " "
// 		            <<  setw (10) << setfill(' ') << "Block"<< " "
// 		            <<  setw (10) << setfill(' ') << "nb_amp"<< " "
// 					<<  setw (10) << setfill(' ') << "nb_rho"<< " "
// 			        <<  setw (10) << setfill(' ') << "nb_phase"<< " "
// 			        <<  setw (10) << setfill(' ') << "nb_xi"<< " "
// 			        <<  setw (10) << setfill(' ') << "nb_sample"<< " "
// 			        << endl;
// 				file_quantstrbook << setw (10) << setfill(' ') << i+1<< " "
// 		            <<  setw (10) << setfill(' ') << j+1<< " "
// 		            <<  setw (10) << setfill(' ') << quantExp[k].nb_amp<< " "
// 					<<  setw (10) << setfill(' ') << quantExp[k].nb_rho<< " "
// 			        <<  setw (10) << setfill(' ') << quantExp[k].nb_phase<< " "
// 			        <<  setw (10) << setfill(' ') << quantExp[k].nb_xi<< " "
// 			        <<  setw (10) << setfill(' ') << quantExp[k].nb_sample<< " "
// 			        << endl;
// 			    file_quantstrbook << "NumElementQ: " << numElementQ << endl;
// 			    file_quantstrbook << "-------------------------------------------------"<< endl;
// 			    file_quantstrbook <<  setw (10) << setfill(' ') << "amp"<< " "
// 					<<  setw (10) << setfill(' ') << "rho"<< " "
// 			        <<  setw (10) << setfill(' ') << "phase"<< " "
// 			        <<  setw (10) << setfill(' ') << "xi"<< " "
// 			        <<  setw (10) << setfill(' ') << "a"<< " "
// 			        <<  setw (10) << setfill(' ') << "b"<< " "
// 			        << endl;
// 			    int j2;
// 			    for( j2 = 0; j2< numElementQ ; j2++)
// 				{
// 			    	file_quantstrbook <<  setw (10) << setfill(' ') << pSBQ[j2].innerProduct << " "
// 					<<  setw (10) << setfill(' ') << pSBQ[j2].rho<< " "
// 			        <<  setw (10) << setfill(' ') << pSBQ[j2].phase<< " "
// 			        <<  setw (10) << setfill(' ') << pSBQ[j2].xi<< " "
// 			        <<  setw (10) << setfill(' ') << pSBQ[j2].a<< " "
// 			        <<  setw (10) << setfill(' ') << pSBQ[j2].b<< " "
// 			        << endl;
// 			    }
// 			    file_quantstrbook.close();

				printf(" ****** Synthesize signal \n");
			    synthSignalSBExp(	recSignal,
			    					norm,
			    					sbbHeader.blockSize,
			    					pSBQ,
			    					numElementQ);

// 				fstream file_sigcomp("sigcomp.out",ios::out);
// 				for(j2 = 0; j2< sbbHeader.blockSize ; j2++)
// 				{
// 					file_sigcomp << blockSignal[j2] << " " << recSignal[j2] << endl;
// 				}
// 				file_sigcomp.close();

				//printf(" ****** Compute MSE \n");
				//rd[k].dist = computeMSE(blockSignal,recSignal,sbbHeader.blockSize);
				printf(" ****** Compute Square Error \n");
				rd[k].dist = computeSqrError(blockSignal,recSignal,0,sbbHeader.blockSize-1);
				printf(" ****** Compute rate \n");
				rd[k].rate = computeRateSBExp(quantExp[k],numElementQ);

				// cout << "Dist/rate: "<< rd[k].dist << " " << rd[k].rate << endl;
// 				cout << "Press enter" << endl;
// 				int lala;
// 				cin >> lala;

				file_rdpoint << setw (10) << setfill(' ') << i+1<< " "
		            <<  setw (10) << setfill(' ') << j+1<< " "
		            <<  setw (10) << setfill(' ') << quantExp[k].nb_amp<< " "
					<<  setw (10) << setfill(' ') << quantExp[k].nb_rho<< " "
			        <<  setw (10) << setfill(' ') << quantExp[k].nb_phase<< " "
			        <<  setw (10) << setfill(' ') << quantExp[k].nb_xi<< " "
			        <<  setw (10) << setfill(' ') << quantExp[k].nb_sample<< " "
			        <<  setw (10) << setfill(' ') << quantExp[k].nb_deltasup<< " "
			        <<  setw (20) << setfill(' ') << rd[k].rate<< " "
			        <<  setw (20) << setfill(' ') << rd[k].dist<< " "
			        <<  setw (20) << setfill(' ') << numElementQ<< " "
			        << endl;

			}

			printf(" *** Obtain RD operation curve\n");
			rdIndex.theta_vec = NULL;
			rdIndex.index_vec = NULL;
			computeRDOpCurve(numQuant,rd,rdIndex);

			int ind;
			strtOpCurveExp* opCurve;
			opCurve = new strtOpCurveExp[rdIndex.numElement];
			for (k=0;k<rdIndex.numElement;k++)
			{
				ind = rdIndex.index_vec[k];
				opCurve[k].nb_amp = quantExp[ind].nb_amp;
	    		opCurve[k].nb_rho = quantExp[ind].nb_rho;
	    		opCurve[k].nb_phase = quantExp[ind].nb_phase;
	    		opCurve[k].nb_xi = quantExp[ind].nb_xi;
	    		opCurve[k].nb_sample = quantExp[ind].nb_sample;
	    		opCurve[k].nb_deltasup = quantExp[ind].nb_deltasup;
	    		opCurve[k].dist = rd[ind].dist;
	    		opCurve[k].rate = rd[ind].rate;
        		opCurve[k].lambda = rdIndex.theta_vec[k];

	    		file_rdopcurve << setw (10) << setfill(' ') << i+1<< " "
		            <<  setw (10) << setfill(' ') << j+1<< " "
		            <<  setw (10) << setfill(' ') << quantExp[ind].nb_amp<< " "
					<<  setw (10) << setfill(' ') << quantExp[ind].nb_rho<< " "
			        <<  setw (10) << setfill(' ') << quantExp[ind].nb_phase<< " "
			        <<  setw (10) << setfill(' ') << quantExp[ind].nb_xi<< " "
			        <<  setw (10) << setfill(' ') << quantExp[ind].nb_sample<< " "
			        <<  setw (10) << setfill(' ') << quantExp[ind].nb_deltasup<< " "
			        <<  setw (20) << setfill(' ') << rd[ind].rate<< " "
			        <<  setw (20) << setfill(' ') << rd[ind].dist<< " "
                    <<  setw (20) << setfill(' ') << rdIndex.theta_vec[k]<< " "
                    <<  setw (20) << setfill(' ') << numElementQ<< " "
			        << endl;
			}

			((CStructBookExp*)structBook[i])[j].setOpCurve(opCurve,rdIndex.numElement);

			delete [] opCurve;
			if (rdIndex.theta_vec != NULL) delete [] rdIndex.theta_vec;
			if (rdIndex.index_vec != NULL) delete [] rdIndex.index_vec;

			delete [] pSBQ;
		}
	}

	// close files
	file_rdpoint.close();
	file_rdopcurve.close();
	file_report.close();

	delete [] rd;
	delete [] recSignal;
	delete [] blockSignal;
}

void optimizeRDExpAmpRangeOpCurve(	CFileRDBitRange* rdBitRange,
									CDataSignal* dataSignal,
									strtSBBHeader sbbHeader,
									CStructBook** structBook,
									int Nfreq,
									char* InputFile)
{
	int i, j, k, t, iAmpRange;
    //////////////////////////////////
    // Defining the quantizers range

    int init_nbit_amp, end_nbit_amp, delta_nbit_amp;
    int init_nbit_rho, end_nbit_rho, delta_nbit_rho;
    int init_nbit_phase, end_nbit_phase, delta_nbit_phase;

    init_nbit_amp   =  rdBitRange->getInitNbitAmp();
    end_nbit_amp    =  rdBitRange->getEndNbitAmp();
    delta_nbit_amp  =  rdBitRange->getDeltaNbitAmp();
    init_nbit_rho   =  rdBitRange->getInitNbitRho();
    end_nbit_rho    =  rdBitRange->getEndNbitRho();
    delta_nbit_rho  =  rdBitRange->getDeltaNbitRho();
    init_nbit_phase =  rdBitRange->getInitNbitPhase();
    end_nbit_phase  =  rdBitRange->getEndNbitPhase();
    delta_nbit_phase  = rdBitRange->getDeltaNbitPhase();

    int numQuant=  (((end_nbit_rho - init_nbit_rho)/delta_nbit_rho) +1) *
                   (((end_nbit_phase - init_nbit_phase)/delta_nbit_phase) +1)
                   + 1;

    // Setting fixed quantizers
    int nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
    int nb_sample = (int)ceil( log(sbbHeader.blockSize) / log(2.0) );
    cout << "nb_xi" << " " << "nb_sample" << endl;
    cout << nb_xi << " " << nb_sample << endl;

    //int nb_block = (int)ceil( log(sbbHeader.numBlock) / log(2.0) );
    // Setting variable quantizers
    strtQuantExp* quantExp;
    quantExp = new strtQuantExp[numQuant];
    quantExp[0].nb_amp = 0;
    quantExp[0].nb_rho = 0;
    quantExp[0].nb_phase = 0;
    quantExp[0].nb_xi = 0;
    quantExp[0].nb_sample = 0;
    t=1;
	for (j=init_nbit_rho;j<=end_nbit_rho;j=j+delta_nbit_rho)
	{
		for (k=init_nbit_phase;k<=end_nbit_phase;k=k+delta_nbit_phase)
		{
			quantExp[t].nb_amp = 0;
			quantExp[t].nb_rho = j;
			quantExp[t].nb_phase = k;
			quantExp[t].nb_xi = nb_xi;
			quantExp[t].nb_sample = nb_sample;
			t++;
		}
	}
	// variables
	int deltaSupMax;
	strtRD* rd = new strtRD[numQuant];
    strtRDIndex rdIndex;
    strtContinuousExp* pSB;
    int numElement;
    double 	min_amp,max_amp,
    		min_rho,max_rho,
    		min_phase,max_phase;
    ////
    double* recSignal;
    double* residueSignal;
    double* blockSignal;
    recSignal = new double[sbbHeader.blockSize];
    residueSignal = new double[sbbHeader.blockSize];
    blockSignal = new double[sbbHeader.blockSize];
    double** origSignal;
    origSignal = dataSignal->getSignal();
    ////
    int nbamp;
    ////
    strtOpCurveQuantExp**** opCurve;

	opCurve = new strtOpCurveQuantExp***[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		opCurve[i] =  new strtOpCurveQuantExp**[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			opCurve[i][j] = new strtOpCurveQuantExp*[end_nbit_amp];
			for (nbamp=init_nbit_amp;nbamp<=end_nbit_amp;nbamp=nbamp+delta_nbit_amp)
			{
				opCurve[i][j][nbamp-1] = new strtOpCurveQuantExp[nbamp];
			}
		}
	}

   	//file pointers
   	char FName[_MAX_PATH];
    char aux[_MAX_PATH];
    char* pos;
    strcpy(FName, InputFile);
    pos = strrchr( FName, '.');
    sprintf(aux,"_rdbyamp_rdpoint.out");
    strcpy( &pos[0], aux);
    fstream file_rdpoint(FName,ios::out);
    strcpy(FName, InputFile);
    pos = strrchr( FName, '.');
    sprintf(aux,"_rdbyamp_rdopcurve.out");
    strcpy( &pos[0], aux);
    fstream file_rdopcurve(FName,ios::out);
    // file header
    file_rdpoint << setw (10) << setfill(' ') << "Signal" << " "
		            <<  setw (10) << setfill(' ') << "Block" << " "
		            <<  setw (10) << setfill(' ') << "nb_amp" << " "
			        <<  setw (10) << setfill(' ') << "iAmpRange" << " "
			        <<  setw (10) << setfill(' ') << "nb_rho" << " "
			        <<  setw (10) << setfill(' ') << "nb_phase"<< " "
			        <<  setw (10) << setfill(' ') << "nb_xi"<< " "
			        <<  setw (10) << setfill(' ') << "nb_sample"<< " "
			        <<  setw (10) << setfill(' ') << "nb_dsup"<< " "
			        <<  setw (20) << setfill(' ') << "Rate"<< " "
			        <<  setw (20) << setfill(' ') << "Distortion"<< " "
			        << endl;
	file_rdopcurve << setw (10) << setfill(' ') << "Signal"<< " "
		            <<  setw (10) << setfill(' ') << "Block"<< " "
		            <<  setw (10) << setfill(' ') << "nb_amp"<< " "
		            <<  setw (10) << setfill(' ') << "iAmpRange" << " "
					<<  setw (10) << setfill(' ') << "nb_rho"<< " "
			        <<  setw (10) << setfill(' ') << "nb_phase"<< " "
			        <<  setw (10) << setfill(' ') << "nb_xi"<< " "
			        <<  setw (10) << setfill(' ') << "nb_sample"<< " "
                    <<  setw (10) << setfill(' ') << "nb_dsup"<< " "
			        <<  setw (20) << setfill(' ') << "Rate"<< " "
			        <<  setw (20) << setfill(' ') << "Distortion"<< " "
                    <<  setw (20) << setfill(' ') << "Lambda"<< " "
			        << endl;

    cout << "numQuant: " << numQuant << endl;
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		double norm = dataSignal->getNorm(i);
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			cout << "Signal: " << i+1 << " -- " << "Block: " << j+1 << endl;
			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			numElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp,
															  max_amp,
															  min_rho,
															  max_rho,
															  min_phase,
															  max_phase);
            deltaSupMax = ((CStructBookExp*)structBook[i])[j].findDeltaSupMax();

			int nb_deltasup = (int)ceil( log2(deltaSupMax+1) );

			double L_amp = max_amp - min_amp;

			// init vector
			for (k=0;k<sbbHeader.blockSize;k++)
				blockSignal[k] = 0.0;
			//----------------------------
			// int initBlock = j*sbbHeader.blockSize;
// 			int copylen = min(sbbHeader.signalSize-initBlock, sbbHeader.blockSize);
// 			memcpy(blockSignal,&origSignal[i][initBlock],copylen*sizeof(double));

			for (nbamp=init_nbit_amp;nbamp<=end_nbit_amp;nbamp=nbamp+delta_nbit_amp)
			{
				int NAmpRange = nbamp;

				CStructBook* structBookByAmp;
				structBookByAmp = new CStructBookExp[NAmpRange];

				double* ampRangeLimit;
				ampRangeLimit = new double[NAmpRange+1];

				double* ampBar;
				ampBar = new double[NAmpRange];

				double* ampRangeNumElement;
				ampRangeNumElement = new double[NAmpRange];

				ampRangeLimit[0]=-0.001;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				    ampBar[NAmpRange-1-iAmpRange] =
				    	(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				    ((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSB,
																			numElement,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

                    strtContinuousExp* pSBAmp =  ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
                    int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();
					strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElementAmp];
					int numElementAmpQ;
//
					synthSignalSBExp( 	blockSignal,
										norm,
										sbbHeader.blockSize,
										pSBAmp,
										numElementAmp);

// 				 	synthSignalSBExp( 	recSignal,
// 										norm,
// 										sbbHeader.blockSize,
// 										pSBAmp,
// 										numElementAmp);
//
// 					for (k=0;k<sbbHeader.blockSize;k++)
// 					{
// 						residueSignal[k] = blockSignal[k] - recSignal[k];
// 					}

					k=0;
					rd[0].dist = computeSqrNorm(blockSignal,sbbHeader.blockSize);
					rd[0].rate = 0;
					quantExp[0].nb_amp = nbamp;
					quantExp[k].nb_deltasup = nb_deltasup;
					file_rdpoint << setw (10) << setfill(' ') << i+1 << " "
						<<  setw (10) << setfill(' ') << j+1 << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_amp << " "
						<<  setw (10) << setfill(' ') << iAmpRange+1 << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_rho << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_phase<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_xi<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_sample << " "
						<<  setw (20) << setfill(' ') << rd[k].rate << " "
						<<  setw (20) << setfill(' ') << rd[k].dist << " "
						<< endl;
					for (k=1;k<numQuant;k++)
					{
						quantExp[k].nb_amp = nbamp;
                        quantExp[k].nb_deltasup = nb_deltasup;
						cout << "Quantizer: " << k+1 << endl;
						printf(" ****** Quantize Structure Book \n");
						quantizeStructBookExp(min_amp,max_amp,
											  min_rho,max_rho,
											  min_phase,max_phase,
											  quantExp[k],
											  pSBAmp,numElementAmp,
											  pSBAmpQ,numElementAmpQ); // output elements

						printf(" ****** Synthesize signal \n");
						synthSignalSBExp(	recSignal,
											norm,
											sbbHeader.blockSize,
											pSBAmpQ,
											numElementAmpQ);

					   	// for (t=0;t<sbbHeader.blockSize;t++)
// 					   	{
//  							recSignal[t] = residueSignal[t] + recSignal[t];
//  						}

						printf(" ****** Compute Square Error \n");
						rd[k].dist = computeSqrError(blockSignal,recSignal,0,sbbHeader.blockSize-1);
						printf(" ****** Compute rate \n");
						rd[k].rate = computeRateSBExp(quantExp[k],numElementAmpQ);


						file_rdpoint << setw (10) << setfill(' ') << i+1 << " "
						<<  setw (10) << setfill(' ') << j+1 << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_amp << " "
						<<  setw (10) << setfill(' ') << iAmpRange+1 << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_rho << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_phase<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_xi<< " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_sample << " "
						<<  setw (10) << setfill(' ') << quantExp[k].nb_deltasup << " "
						<<  setw (20) << setfill(' ') << rd[k].rate << " "
						<<  setw (20) << setfill(' ') << rd[k].dist << " "
						<< endl;

					}
					printf(" *** Obtain RD operation curve\n");
					rdIndex.theta_vec = NULL;
					rdIndex.index_vec = NULL;
					computeRDOpCurve(numQuant,rd,rdIndex);

					int ind;

					opCurve[i][j][nbamp-1][iAmpRange].numElement = rdIndex.numElement;
					opCurve[i][j][nbamp-1][iAmpRange].quantExp = new strtQuantExp[rdIndex.numElement];
					opCurve[i][j][nbamp-1][iAmpRange].rate = new double[rdIndex.numElement];
					opCurve[i][j][nbamp-1][iAmpRange].dist = new double[rdIndex.numElement];
					opCurve[i][j][nbamp-1][iAmpRange].lambda = new double[rdIndex.numElement];
					for (k=0;k<rdIndex.numElement;k++)
					{
						ind = rdIndex.index_vec[k];
						opCurve[i][j][nbamp-1][iAmpRange].quantExp[k] = quantExp[ind];
						opCurve[i][j][nbamp-1][iAmpRange].dist[k] = rd[ind].dist;
						opCurve[i][j][nbamp-1][iAmpRange].rate[k] = rd[ind].rate;
						opCurve[i][j][nbamp-1][iAmpRange].lambda[k] = rdIndex.theta_vec[k];

						file_rdopcurve << setw (10) << setfill(' ') << i+1<< " "
						<<  setw (10) << setfill(' ') << j+1<< " "
						<<  setw (10) << setfill(' ') << quantExp[ind].nb_amp << " "
						<<  setw (10) << setfill(' ') << iAmpRange+1 << " "
						<<  setw (10) << setfill(' ') << quantExp[ind].nb_rho<< " "
						<<  setw (10) << setfill(' ') << quantExp[ind].nb_phase<< " "
						<<  setw (10) << setfill(' ') << quantExp[ind].nb_xi  << " "
						<<  setw (10) << setfill(' ') << quantExp[ind].nb_sample << " "
						<<  setw (10) << setfill(' ') << quantExp[ind].nb_deltasup << " "
						<<  setw (20) << setfill(' ') << rd[ind].rate<< " "
						<<  setw (20) << setfill(' ') << rd[ind].dist << " "
						<<  setw (20) << setfill(' ') << opCurve[i][j][nbamp-1][iAmpRange].lambda[k] << " "
						<< endl;
					}


					if (rdIndex.theta_vec != NULL) delete [] rdIndex.theta_vec;
					if (rdIndex.index_vec != NULL) delete [] rdIndex.index_vec;

					delete [] pSBAmpQ;
				}
				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;

			}
		}
	}
	file_rdpoint.close();
    file_rdopcurve.close();

	/////////////////////////
	/// SAVE FILE
	cout << "Saving file opcurve_amprange.bin" << endl;
	double dummydouble;
	int dummyint;
	strcpy(FName, InputFile);
    pos = strrchr( FName, '.');
    sprintf(aux,"_rdbyamp_opcurve.bin");
    strcpy( &pos[0], aux);
	FILE* io_opcurve;
	io_opcurve = fopen(FName,"wb");
	dummyint =  sbbHeader.numSignal;
	fwrite ( &dummyint, sizeof ( int ), 1, io_opcurve );
	dummyint =  sbbHeader.numBlock;
	fwrite ( &dummyint, sizeof ( int ), 1, io_opcurve );
	dummyint =  init_nbit_amp;
	fwrite ( &dummyint, sizeof ( int ), 1, io_opcurve );
	dummyint =  end_nbit_amp;
	fwrite ( &dummyint, sizeof ( int ), 1, io_opcurve );
	dummyint =  delta_nbit_amp;
	fwrite ( &dummyint, sizeof ( int ), 1, io_opcurve );
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		for (j=0; j<sbbHeader.numBlock; j++)
		{
			dummydouble =  ((CStructBookExp*)structBook[i])[j].getMinAmp();
			fwrite ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			dummydouble =  ((CStructBookExp*)structBook[i])[j].getMaxAmp();
			fwrite ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			dummydouble =  ((CStructBookExp*)structBook[i])[j].getMinRho();
			fwrite ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			dummydouble =  ((CStructBookExp*)structBook[i])[j].getMaxRho();
			fwrite ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			dummydouble =  ((CStructBookExp*)structBook[i])[j].getMinPhase();
			fwrite ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			dummydouble =  ((CStructBookExp*)structBook[i])[j].getMaxPhase();
			fwrite ( &dummydouble, sizeof ( double ), 1, io_opcurve );

			for (nbamp=init_nbit_amp;nbamp<=end_nbit_amp;nbamp=nbamp+delta_nbit_amp)
			{
				for (iAmpRange=0;iAmpRange<nbamp;iAmpRange++)
				{
					int nels = opCurve[i][j][nbamp-1][iAmpRange].numElement;
					fwrite ( &nels, sizeof ( int ), 1, io_opcurve );
					fwrite ( opCurve[i][j][nbamp-1][iAmpRange].quantExp, sizeof ( strtQuantExp ), nels, io_opcurve);
					fwrite ( opCurve[i][j][nbamp-1][iAmpRange].dist, sizeof ( double ), nels, io_opcurve);
					fwrite ( opCurve[i][j][nbamp-1][iAmpRange].rate, sizeof ( double ), nels, io_opcurve);
					fwrite ( opCurve[i][j][nbamp-1][iAmpRange].lambda, sizeof (double), nels, io_opcurve);
				}
			}
		}
	}
	fclose(io_opcurve);

	//////////////
	// DEALLOCATE
	cout << "deallocate" << endl;

	delete [] rd;
	delete [] recSignal;
	delete [] residueSignal;
	delete [] blockSignal;

	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (nbamp=init_nbit_amp;nbamp<=end_nbit_amp;nbamp=nbamp+delta_nbit_amp)
			{
    			delete [] opCurve[i][j][nbamp-1];

    		}
			delete [] opCurve[i][j];
		}
		delete [] opCurve[i];
	}
}

void quantizeStructBookExp(   double min_amp, double max_amp,
							  double min_rho,double max_rho,
							  double min_phase, double max_phase,
							  strtQuantExp quantExp,
							  strtContinuousExp* pSB, int numElement,
							  strtContinuousExp* pSBQ,int &numElementQ)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(quantExp.nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp -1.0) );

    double step_rho, nlevel_rho;
    nlevel_rho = round (pow(2.0, static_cast<double>(quantExp.nb_rho)) );
    if (nlevel_rho==0.0)
        step_rho=1e+30;
    else
        step_rho = fabs( (max_rho - min_rho) / (nlevel_rho - 1.0) );

    double step_phase, nlevel_phase;
    nlevel_phase = round(pow(2.0, static_cast<double>(quantExp.nb_phase)) );
    if (nlevel_phase==0.0)
        step_phase =1e+30;
    else
        step_phase = fabs( (max_phase - min_phase) / (nlevel_phase -1.0) );

    double coef;
    numElementQ=0;
    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct - min_amp,step_amp);
        if (coef!=0.0)
        {
        	pSBQ[numElementQ].innerProduct = coef;
        	pSBQ[numElementQ].phase = linearQuant(pSB[i].phase - min_phase,step_phase);
        	if (pSB[i].rho<0)
        		pSBQ[numElementQ].rho = - linearQuant(fabs(pSB[i].rho)- min_rho,step_rho);
        	else
        		pSBQ[numElementQ].rho = linearQuant(fabs(pSB[i].rho)- min_rho,step_rho);

        	pSBQ[numElementQ].xi = pSB[i].xi;
        	pSBQ[numElementQ].a = pSB[i].a;
        	pSBQ[numElementQ].b = pSB[i].b;

        	numElementQ++;
        }
    }
}

double quantizeStructBookExpDistAtom( double min_amp, double max_amp,
									  double min_rho,double max_rho,
									  double min_phase, double max_phase,
									  strtQuantExp quantExp,
									  strtContinuousExp* pSB, int numElement,
									  strtContinuousExp* pSBQ,int &numElementQ,
									  int sigSize, double signorm)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(quantExp.nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp -1.0) );

    double step_rho, nlevel_rho;
    nlevel_rho = round (pow(2.0, static_cast<double>(quantExp.nb_rho)) );
    if (nlevel_rho==0.0)
        step_rho=1e+30;
    else
        step_rho = fabs( (max_rho - min_rho) / (nlevel_rho - 1.0) );

    double step_phase, nlevel_phase;
    nlevel_phase = round(pow(2.0, static_cast<double>(quantExp.nb_phase)) );
    if (nlevel_phase==0.0)
        step_phase =1e+30;
    else
        step_phase = fabs( (max_phase - min_phase) / (nlevel_phase -1.0) );


    double distAtom=0.0;
    double* origAtom = new double[sigSize];
    double* recAtom = new double[sigSize];

    for (j=0; j < sigSize; j++)
	{
		origAtom[j] = 0.0;
		recAtom[j]  = 0.0;
	}

    double coef;
    numElementQ=0;
    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct - min_amp,step_amp);
        if (coef!=0.0)
        {
        	pSBQ[numElementQ].innerProduct = coef;
        	pSBQ[numElementQ].phase = linearQuant(pSB[i].phase - min_phase,step_phase);
        	if (pSB[i].rho<0)
        		pSBQ[numElementQ].rho = - linearQuant(fabs(pSB[i].rho)- min_rho,step_rho);
        	else
        		pSBQ[numElementQ].rho = linearQuant(fabs(pSB[i].rho)- min_rho,step_rho);

        	pSBQ[numElementQ].xi = pSB[i].xi;
        	pSBQ[numElementQ].a = pSB[i].a;
        	pSBQ[numElementQ].b = pSB[i].b;




        	synthSignalAtomExp(	origAtom,
								sigSize,
								pSB[i]);

			synthSignalAtomExp(	recAtom,
								sigSize,
								pSBQ[numElementQ]);

			for (j=pSB[i].a; j<=pSB[i].b; j++)
			{
				origAtom[j] = signorm*pSB[i].innerProduct *origAtom[j];
				recAtom[j] = signorm*coef*recAtom[j];
			}

			distAtom += computeSqrError(origAtom,
										recAtom,
										pSB[i].a,
										pSB[i].b);

        	numElementQ++;

        }
        else
        {
        	synthSignalAtomExp(	origAtom,
								sigSize,
								pSB[i]);
			for (j=pSB[i].a; j<=pSB[i].b; j++)
			{
				origAtom[j] = signorm*pSB[i].innerProduct *origAtom[j];
				distAtom+= origAtom[j]*origAtom[j];
			}
        }
    }

    delete [] origAtom;
    delete [] recAtom;

    return distAtom;
}



void quantizeStructBookExpRecUnquantNoReset( double min_amp, double max_amp,
									  double min_rho,double max_rho,
									  double min_phase, double max_phase,
									  strtQuantExp quantExp,
									  strtContinuousExp* pSB, int numElement,
									  strtContinuousExp* pSBQ,int &numElementQ,
									  double* recSignal, double signorm, int signalSize)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(quantExp.nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp -1.0) );

    double step_rho, nlevel_rho;
    nlevel_rho = round (pow(2.0, static_cast<double>(quantExp.nb_rho)) );
    if (nlevel_rho==0.0)
        step_rho=1e+30;
    else
        step_rho = fabs( (max_rho - min_rho) / (nlevel_rho - 1.0) );

    double step_phase, nlevel_phase;
    nlevel_phase = round(pow(2.0, static_cast<double>(quantExp.nb_phase)) );
    if (nlevel_phase==0.0)
        step_phase =1e+30;
    else
        step_phase = fabs( (max_phase - min_phase) / (nlevel_phase - 1.0) );

    double coef;
    numElementQ=0;

    double* recAtom;
    recAtom = new double[signalSize];
    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct - min_amp,step_amp);
        if (coef!=0.0)
        {
        	pSBQ[numElementQ].innerProduct = coef;
        	pSBQ[numElementQ].phase = linearQuant(pSB[i].phase - min_phase,step_phase);
        	if (pSB[i].rho<0)
        		pSBQ[numElementQ].rho = - linearQuant(fabs(pSB[i].rho)- min_rho,step_rho);
        	else
        		pSBQ[numElementQ].rho = linearQuant(fabs(pSB[i].rho)- min_rho,step_rho);

        	pSBQ[numElementQ].xi = pSB[i].xi;
        	pSBQ[numElementQ].a = pSB[i].a;
        	pSBQ[numElementQ].b = pSB[i].b;

        	numElementQ++;

        	synthSignalAtomExp(	recAtom,
								signalSize,
								pSB[i]);

			for (j=pSB[i].a; j<=pSB[i].b; j++)
				recSignal[j] = recSignal[j] + signorm*pSB[i].innerProduct*recAtom[j];
        }
    }
    delete [] recAtom;
}

void quantizeStructBookExpRecUnquantNoResetArithCod(  double min_amp, double max_amp,
													  strtQuantExp quantExp,
													  strtContinuousExp* pSB, int numElement,
													  strtContinuousExp* pSBQ,int &numElementQ,
													  double* recSignal, double signorm, int signalSize,
													  int* ampIndex, int* rhoIndex, int* xiIndex,
													  int* phiIndex, int* aIndex, int* bIndex,
													  int fdiscrtype, double step_xi)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(quantExp.nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp - 1.0) );

    double coef;
    numElementQ=0;

    double step_rho = quantExp.qrho;
    double step_phase = quantExp.qphi;

    double* recAtom;
    recAtom = new double[signalSize];
    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct,step_amp);
        if (coef!=0.0)
        {
        	pSBQ[numElementQ].innerProduct = coef;
        	pSBQ[numElementQ].phase = linearQuant(pSB[i].phase,step_phase);
        	pSBQ[numElementQ].rho = linearQuant(pSB[i].rho,step_rho);
        	pSBQ[numElementQ].xi = pSB[i].xi;
        	pSBQ[numElementQ].a = pSB[i].a;
         	pSBQ[numElementQ].b = pSB[i].b;

//         	cout << ">> Coef: " << pSB[i].innerProduct << " -- " <<  pSBQ[numElementQ].innerProduct << " // " <<  step_amp << endl;
//         	cout << ">> Phi : " << pSB[i].phase << " -- " <<  pSBQ[numElementQ].phase << " // " <<  step_phase << endl;
//         	cout << ">> Rho : " << pSB[i].rho << " -- " <<  pSBQ[numElementQ].rho << " // " <<  step_rho << endl;
        	/////
        	ampIndex[numElementQ] = linearQuantIndex(pSB[i].innerProduct,step_amp);
        	rhoIndex[numElementQ] = linearQuantIndex(pSB[i].rho,step_rho);
        	phiIndex[numElementQ] = linearQuantIndex(pSB[i].phase,step_phase);
        	if (fdiscrtype==1)
        	{
        		xiIndex[numElementQ] = linearQuantIndex(pSB[i].xi, step_xi);
			}
       		else if (fdiscrtype==2)
       		{
       			xiIndex[numElementQ] = geometricQuantIndex(pSB[i].xi, step_xi, 24);
       		}
        	aIndex[numElementQ] = pSB[i].a;
        	bIndex[numElementQ] = pSB[i].b;

        	numElementQ++;

        	synthSignalAtomExp(	recAtom,
								signalSize,
								pSB[i]);

			for (j=pSB[i].a; j<=pSB[i].b; j++)
				recSignal[j] = recSignal[j] + signorm*pSB[i].innerProduct*recAtom[j];
        }
    }
    delete [] recAtom;
}

int computeNumElementAmpQ(double min_amp, double max_amp,
						  int nb_amp,
						  strtContinuousExp* pSB, int numElement)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp - 1.0) );

    double coef;
    int numElementQ=0;

    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct,step_amp);
        if (coef!=0.0)
        {
			numElementQ++;
		}
	}

	return numElementQ;
}


void quantizeStructBookExpAmp(  double min_amp, double max_amp,
							  	int nb_amp,
							  	strtContinuousExp* pSB, int numElement,
							  	strtContinuousExp* pSBQ,int &numElementQ)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp - 1.0) );

    double coef;
    numElementQ=0;
    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct - min_amp,step_amp);
        if (coef!=0.0)
        {
        	pSBQ[numElementQ].innerProduct = coef;
        	pSBQ[numElementQ].phase = pSB[i].phase;
        	pSBQ[numElementQ].rho = pSB[i].rho;
        	pSBQ[numElementQ].xi = pSB[i].xi;
        	pSBQ[numElementQ].a = pSB[i].a;
        	pSBQ[numElementQ].b = pSB[i].b;

        	numElementQ++;
        }
    }
}

void dequantizeStructBookExpArithCodAmpVec(   int nb_amp, double min_amp, double max_amp,
											  unsigned* ampIndex, int numElement,
											  double* ampvec)
{
    int i;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp- 1.0) );

    for (i=0; i< numElement; i++)
    {
		ampvec[i] = linearDeQuantIndex(static_cast<int>(ampIndex[i]), step_amp);
    }
}

void dequantizeStructBookExpArithCod( int nb_amp, double min_amp, double max_amp,
									  int fdiscrtype, double step_xi,
									  double step_rho, double step_phase,
									  unsigned* ampIndex, int* rhoIndex, unsigned* xiIndex,
									  unsigned* phiIndex, unsigned* aIndex, unsigned* bIndex,
									  strtContinuousExp* pSB, int numElement)
{
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp- 1.0) );

    for (i=0; i< numElement; i++)
    {
		pSB[i].innerProduct = linearDeQuantIndex(static_cast<int>(ampIndex[i]), step_amp);

		pSB[i].rho = linearDeQuantIndex(rhoIndex[i],step_rho);
		pSB[i].phase = linearDeQuantIndex(static_cast<int>(phiIndex[i]),step_phase);
		if (fdiscrtype==1)
		{
			pSB[i].xi = linearDeQuantIndex(static_cast<int>(xiIndex[i]), step_xi);
		}
		else if (fdiscrtype==2)
		{
			pSB[i].xi = geometricDeQuantIndex(static_cast<int>(xiIndex[i]), step_xi, 24);
		}
		pSB[i].a = static_cast<int>(aIndex[i]);
		pSB[i].b = static_cast<int>(bIndex[i]);


// 		cout 	<<  pSB[i].innerProduct << " "
// 				<<  pSB[i].rho << " "
// 				<<  pSB[i].phase << " "
// 				<<  pSB[i].xi << " "
// 				<<  pSB[i].a << " "
// 				<<  pSB[i].b << " " << endl;
    }
}

double sumSqrErrorStructBookExpAmp( double min_amp, double max_amp,
									int nb_amp,
									strtContinuousExp* pSB, int numElement)
{
	double sumSqrError=0.0;
    int i,j;
    double step_amp, nlevel_amp;
    nlevel_amp = round(pow(2.0, static_cast<double>(nb_amp) ) );
    if (nlevel_amp==0.0)
        step_amp = 1e+30;
    else
        step_amp = fabs( (max_amp - min_amp) / (nlevel_amp - 1.0) );

    double coef;
    for (i=0; i< numElement; i++)
    {
        coef = linearQuant(pSB[i].innerProduct - min_amp,step_amp);
		sumSqrError +=  (pSB[i].innerProduct - min_amp - coef)*
						(pSB[i].innerProduct - min_amp - coef);
    }
    return sumSqrError;
}



void synthSignalSBExp(   double* recSignal,
						 double norm,
						 int sigSize,
						 strtContinuousExp* sb,
						 int numElement)
{
    CExpDictionary* expDic = new CExpDictionary;
    expDic->setSignalSize(sigSize);
    CParameter* parm;
    parm = new CExpParm;
    int i,j;
    double* sigAux;
//      cout << "norm: "<< norm << endl;
//      cout << "sigSize: "<< sigSize << endl;
//  	cout << "numElement: "<< numElement << endl;
    for (i=0;i<sigSize;i++)
    {
        recSignal[i] = 0.0;
    }
    for (i=0;i<numElement;i++)
    {
        ((CExpParm*)parm)->setParm(sb[i]);
//        ((CExpParm*)parm)->printParm2Screen();
        expDic->setRealAtomOnSupport(parm);
        sigAux = expDic->getRealAtom();
        for (j=sb[i].a;j<=sb[i].b;j++)
        {
            recSignal[j] += norm*sb[i].innerProduct*sigAux[j];
        }
    }
    delete expDic;
    delete parm;
}

void synthSignalSBExpNoReset( double* recSignal,
							  double norm,
							  int sigSize,
							  strtContinuousExp* sb,
							  int numElement)
{
    CExpDictionary* expDic = new CExpDictionary;
    expDic->setSignalSize(sigSize);
    CParameter* parm;
    parm = new CExpParm;
    int i,j;
    double* sigAux;
//     cout << "norm: "<< norm << endl;
//     cout << "sigSize: "<< sigSize << endl;
//     cout << "numElement: "<< numElement << endl;
    for (i=0;i<numElement;i++)
    {
        ((CExpParm*)parm)->setParm(sb[i]);
        //((CExpParm*)parm)->printParm2Screen();
        expDic->setRealAtomOnSupport(parm);
        sigAux = expDic->getRealAtom();
//        cout << i << "; a: " << sb[i].a << "; b: " << sb[i].b << endl;
        for (j=sb[i].a;j<=sb[i].b;j++)
        {
            recSignal[j] += norm*sb[i].innerProduct*sigAux[j];
        }
    }
    delete expDic;
    delete parm;
}

double computeRateSBExp(strtQuantExp quantExp,int numElementQ)
{
//	double rate = numElementQ*( quantExp.nb_amp + quantExp.nb_rho + quantExp.nb_phase +
//	                            quantExp.nb_xi + 2*quantExp.nb_sample + 1 /* rho sign bit*/);
    double rate = numElementQ*( quantExp.nb_amp + quantExp.nb_rho + quantExp.nb_phase +
	                            quantExp.nb_xi + quantExp.nb_sample + 1 /* rho sign bit*/ +
	                            quantExp.nb_deltasup);

	return rate;
}

double computeRateParameter(double nb,int numElementQ)
{
	double rate = nb*numElementQ;
	return rate;
}

double computeRateSBExpAmp(double nb,int numElementQ)
{
	double rate = nb*numElementQ;
	return rate;
}

double computeRateSBExpFreq(double nb,int numElementQ)
{
	double rate = nb*numElementQ;
	return rate;
}

double computeRateSBExpSample(double nb,int numElementQ)
{
	double rate = nb*numElementQ;
	return rate;
}

double computeRateUnifPdf(int nQStep,
						  int numElementQ)
{
	double rate = numElementQ*( -log2( 1.0/ static_cast<double>( nQStep) ) );
	return rate;
}

double computeRateSBExpPhaseUnifPdf(int nQPhiStep,
									int numElementQ)
{
	double rate = numElementQ*( -log2( 1.0/ static_cast<double>( nQPhiStep) ) );
	return rate;
}

double computeRateSBExpDecayUnifPdf(int nQRhoStep,
									int numElementQ)
{
	double rate = numElementQ*( -log2( 1.0/ static_cast<double>( nQRhoStep) ) );
	return rate;
}

double computeRatePdf(	double* vec, int numElement,
						double* qstepEdge, double* qstepProb, int numQStep)
{
	double rate=0.0;
	int i,j;

	int* qstepCount = new int[numQStep];
	for(j=0;j<numQStep;j++)
		qstepCount[j] = 0;

	int count=0;
	j=0;
	i=0;
	while(i<numElement)
	{
		qstepCount[j] = 0;
		while( vec[i] < qstepEdge[j+1] )
		{
			qstepCount[j]++;
			i++;
			if (i==numElement) break;
		}
		j++;
		if (j==numQStep) break;
	}

	for(j=0;j<numQStep;j++)
	{
		rate +=  static_cast<double>(qstepCount[j]) * (-log2(qstepProb[j]));
	}
	delete [] qstepCount;

	return rate;
}

double computeRateSBExpDecayGGDPdf(	double* rhovec, int numElement,
									double* qstepEdge, double* qstepProb, int numQStep)
{
	double rate=0.0;
	int i,j;

	int* qstepCount = new int[numQStep];
	for(j=0;j<numQStep;j++)
		qstepCount[j] = 0;

	int count=0;
	j=0;
	i=0;
	while(i<numElement)
	{
		qstepCount[j] = 0;
		while( rhovec[i] < qstepEdge[j+1] )
		{
			qstepCount[j]++;
			i++;
			if (i==numElement) break;
		}
		j++;
		if (j==numQStep) break;
	}

	for(j=0;j<numQStep;j++)
	{
		rate +=  static_cast<double>(qstepCount[j]) * (-log2(qstepProb[j]));
	}
	delete [] qstepCount;

	return rate;
}

double computeRateSBExpDecayGGDPdfPrint(	double* rhovec, int numElement,
									double* qstepEdge, double* qstepProb, int numQStep,
									int flag)
{
	double rate=0.0;
	int i,j;

	int* qstepCount = new int[numQStep];
	for(j=0;j<numQStep;j++)
		qstepCount[j] = 0;

	int count=0;
	j=0;
	i=0;
	while(i<numElement)
	{
		qstepCount[j] = 0;
		while( rhovec[i] < qstepEdge[j+1] )
		{
			qstepCount[j]++;
			i++;
			if (i==numElement) break;
		}
		j++;
		if (j==numQStep) break;
	}

	fstream file("rdbyamp_rhoQStepCountProb.out",ios::out);

	for(j=0;j<numQStep;j++)
	{
		rate +=  static_cast<double>(qstepCount[j]) * (-log2(qstepProb[j]));

		if (flag ==1)
		{
			file << setw (20) << setfill(' ') << qstepCount[j] <<" "
				 << setw (20) << setfill(' ') << qstepProb[j]  <<" "<< endl;
		}

	}
	file.close();
	delete [] qstepCount;

	return rate;
}

double computeEntropyUnifPdf(int nQStep)
{
	double entropy=0.0;
	int i;
	for(i=0;i<nQStep;i++)
	{
		entropy +=  -  ( 1.0/ static_cast<double>( nQStep) ) * log2( 1.0/ static_cast<double>( nQStep) ) ;
	}
	return entropy;
}

double computeEntropySBExpPhaseUnifPdf(int nQPhiStep)
{
	double entropy=0.0;
	int i;
	for(i=0;i<nQPhiStep;i++)
	{
		entropy +=  -  ( 1.0/ static_cast<double>( nQPhiStep) ) * log2( 1.0/ static_cast<double>( nQPhiStep) ) ;
	}
	return entropy;
}

double computeEntropySBExpDecayUnifPdf(int nQRhoStep)
{
	double entropy=0.0;
	int i;
	for(i=0;i<nQRhoStep;i++)
	{
		entropy +=  -  ( 1.0/ static_cast<double>( nQRhoStep) ) * log2( 1.0/ static_cast<double>( nQRhoStep) ) ;
	}
	return entropy;
}

double computeEntropyPdf( double* qstepProb, int numQStep)
{
	double entropy=0.0;
	int j;

	for(j=0;j<numQStep;j++)
	{
		entropy +=  - (qstepProb[j]) * log2(qstepProb[j]);
	}

	return entropy;
}

double computeEntropySBExpDecayGGDPdf( double* qstepProb, int numQStep)
{
	double entropy=0.0;
	int j;


	for(j=0;j<numQStep;j++)
	{
		entropy +=  - (qstepProb[j]) * log2(qstepProb[j]);
	}

	return entropy;
}

void computeRDOpCurve(	int numQuant,
						strtRD* rd,
						strtRDIndex& rdIndex)
{
	double* theta_vec;
	double* theta_vec_aux;
	int* index_vec;
	int* index_vec_aux;
	int numElement = 0;
	theta_vec = NULL;
	theta_vec_aux = NULL;
	index_vec = NULL;
	index_vec_aux = NULL;

	int firstPoint=0;
	int secondPoint=0;

	int i;

	double max_rate=0;
	// Find first point and max rate
	//printf("- Find first point and max rate.\n");
	for(i=0; i<numQuant; i++)
	{
		if (rd[firstPoint].rate > rd[i].rate)
		{
			firstPoint = i;
		}
		if (max_rate < rd[i].rate)
		{
			max_rate = rd[i].rate;
		}
	}

	//printf("- Building Convex hull.\n");
	double min_theta;
	double theta=0;
	double delta_rate=0;
	double delta_dist=0;
	double min_delta_rate=0;

	while(1)
	{
		min_theta = pi / 2.0;
		min_delta_rate = max_rate;
		for(i=0; i<numQuant; i++)
		{
			if ((i!=firstPoint) &&
				(rd[i].rate > rd[firstPoint].rate) )
			{
				delta_rate = rd[i].rate - rd[firstPoint].rate;
				delta_dist = rd[i].dist - rd[firstPoint].dist;

				if (delta_rate==0)
				{
					if(delta_dist>=0)
						theta = pi/2.0;
					else
						theta = -pi/2.0;
				}
				else
				{
					theta = atan(delta_dist/delta_rate);
				}

				if ( (theta < min_theta) && (theta <0))
				{
					min_theta = theta;
					secondPoint = i;
				}
			}
		}

		if (numElement!=0)
		{
			if (index_vec_aux!=NULL)
			{
				delete [] index_vec_aux;
			}
			index_vec_aux = new int[numElement];
			memcpy(index_vec_aux,index_vec,sizeof(int)*numElement);

			if (theta_vec_aux!=NULL)
			{
				delete [] theta_vec_aux;
			}
			theta_vec_aux = new double[numElement];
			memcpy(theta_vec_aux,theta_vec,sizeof(double)*numElement);
		}
		if (index_vec!=NULL)
		{
			delete [] index_vec;
		}
		if (theta_vec!=NULL)
		{
			delete [] theta_vec;
		}
		numElement++;
		index_vec = new int[numElement];
		theta_vec = new double[numElement];
		if ((numElement-1)!=0)
		{
			memcpy(index_vec,index_vec_aux,sizeof(int)*(numElement-1));
			memcpy(theta_vec,theta_vec_aux,sizeof(double)*(numElement-1));
		}
		index_vec[numElement-1] = firstPoint;
		theta_vec[numElement-1] = min_theta;

		if (min_theta>0)
		{
			theta_vec[numElement-1] = 0;
			break;
		}

		firstPoint = secondPoint;

		if (rd[firstPoint].rate == max_rate)
		{
			if (index_vec_aux!=NULL)
			{
				delete [] index_vec_aux;
			}
			index_vec_aux = new int[numElement];
			memcpy(index_vec_aux,index_vec,sizeof(int)*numElement);

			if (theta_vec_aux!=NULL)
			{
				delete [] theta_vec_aux;
			}
			theta_vec_aux = new double[numElement];
			memcpy(theta_vec_aux,theta_vec,sizeof(double)*numElement);

			if (index_vec!=NULL)
			{
				delete [] index_vec;
			}
			if (theta_vec!=NULL)
			{
				delete [] theta_vec;
			}
			numElement++;

			index_vec = new int[numElement];
			memcpy(index_vec,index_vec_aux,sizeof(int)*(numElement-1));
			index_vec[numElement-1] = firstPoint;

			theta_vec = new double[numElement];
			memcpy(theta_vec,theta_vec_aux,sizeof(double)*(numElement-1));
			theta_vec[numElement-1] = 0;

			break;
		}
	}

	//printf("- Load convex hull.\n");
	//printf("- numElement %d nChannel %d\n",numElement,nChannel);
	//ofstream outFile("bestq_ind.dat",ios::out);

	if (rdIndex.theta_vec == NULL)
	{
		rdIndex.theta_vec = new double[numElement];
	}
	if (rdIndex.index_vec == NULL)
	{
		rdIndex.index_vec = new int[numElement];
	}
	rdIndex.numElement = numElement;
	printf("      Index | Rate| Theta\n");
	for(int g=0;g<numElement;g++)
	{
		//outFile << index_vec[g] << ' ';
		rdIndex.theta_vec[g] = theta_vec[g];
		rdIndex.index_vec[g] = index_vec[g];

		printf("      %d %f %f\n",rdIndex.index_vec[g],
							rd[index_vec[g]].rate,
							rdIndex.theta_vec[g]);
	}
	cout << endl;

	if (index_vec!=NULL) delete [] index_vec;
	if (index_vec_aux!=NULL) delete [] index_vec_aux;
	if (theta_vec!=NULL) delete [] theta_vec;
	if (theta_vec_aux!=NULL) delete [] theta_vec_aux;
}

void optimizeRDExpByAmp(CFileRDBitRange* rdBitRange,
						CDataSignal* dataSignal,
						strtSBBHeader sbbHeader,
						CStructBook** structBook,
						int Nfreq)
{
    //////////////////////////////////
    // Defining the quantizers range

    int init_nbit_amp, end_nbit_amp, delta_nbit_amp;
    int init_nbit_rho, end_nbit_rho, delta_nbit_rho;
    int init_nbit_phase, end_nbit_phase, delta_nbit_phase;

    init_nbit_amp   =  rdBitRange->getInitNbitAmp();
    end_nbit_amp    =  rdBitRange->getEndNbitAmp();
    delta_nbit_amp  =  rdBitRange->getDeltaNbitAmp();
    init_nbit_rho   =  rdBitRange->getInitNbitRho();
    end_nbit_rho    =  rdBitRange->getEndNbitRho();
    delta_nbit_rho  =  rdBitRange->getDeltaNbitRho();
    init_nbit_phase =  rdBitRange->getInitNbitPhase();
    end_nbit_phase  =  rdBitRange->getEndNbitPhase();
    delta_nbit_phase  = rdBitRange->getDeltaNbitPhase();

    int numQuant=  (((end_nbit_amp - init_nbit_amp)/delta_nbit_amp) +1) *
                   (((end_nbit_rho - init_nbit_rho)/delta_nbit_rho) +1) *
                   (((end_nbit_phase - init_nbit_phase)/delta_nbit_phase) +1);

    int numQuantByAmp =  (((end_nbit_rho - init_nbit_rho)/delta_nbit_rho) +1) *
                   (((end_nbit_phase - init_nbit_phase)/delta_nbit_phase) +1);


    // Setting fixed quantizers
    int nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
    int nb_sample = (int)ceil( log(sbbHeader.blockSize) / log(2.0) );

	// variables
    strtContinuousExp* pSB;
    int numElement;
    double 	min_amp,max_amp,
    		min_rho,max_rho,
    		min_phase,max_phase;
    double 	step_amp, step_rho, step_phase;
    double* recSignal;
    double* blockSignal;
    recSignal = new double[sbbHeader.blockSize];
    blockSignal = new double[sbbHeader.blockSize];
    //file pointers
    fstream file_rdbyamp("rdbyamp.out",ios::out);
    fstream file_report("rdbyamp_report.out",ios::out);
    // file header
    file_rdbyamp << setw (10) << setfill(' ') << "Signal" << " "
		            <<  setw (10) << setfill(' ') << "Block" << " "
		            <<  setw (10) << setfill(' ') << "indAmp" << " "
		            <<  setw (10) << setfill(' ') << "nb_amp" << " "
					<<  setw (10) << setfill(' ') << "nb_rho" << " "
			        <<  setw (10) << setfill(' ') << "nb_phase"<< " "
			        <<  setw (10) << setfill(' ') << "nb_xi"<< " "
			        <<  setw (10) << setfill(' ') << "nb_sample"<< " "
			        <<  setw (20) << setfill(' ') << "NumAmpRange"<< " "
			        <<  setw (20) << setfill(' ') << "LowAmpLimit"<< " "
			        <<  setw (20) << setfill(' ') << "UpperAmpLimit"<< " "
			        <<  setw (20) << setfill(' ') << "AmpBar"<< " "
			        <<  setw (20) << setfill(' ') << "MeanCoef"<< " "
			        <<  setw (20) << setfill(' ') << "MeanAtomDist"<< " "
			        <<  setw (20) << setfill(' ') << "MeanCoefDist"<< " "
			        << endl;

	int i, j;

	FILE* iosba;
    iosba = fopen("rdbyamp_header.out","w");
    fprintf(iosba,"Sign. Type :          %5i\n", sbbHeader.type);
    fprintf(iosba,"Dict. Type :          %5i\n", sbbHeader.dicType);
    fprintf(iosba,"No. Signals:          %5i\n", sbbHeader.numSignal);
    fprintf(iosba,"Signal Size:       %8i\n", sbbHeader.signalSize);
    fprintf(iosba,"Block Hop:            %5i\n", sbbHeader.blockHop);
    fprintf(iosba,"Block Size :          %5i\n", sbbHeader.blockSize);
    fprintf(iosba,"No. Blocks :          %5i\n", sbbHeader.numBlock);
    fprintf(iosba,"Samp. Freq :     %10.2f\n", sbbHeader.Fs);
    for (i=0;i<sbbHeader.numSignal;i++)
	{
    	fprintf(iosba,"NormSig. %d :     %10.2f\n",i+1, sbbHeader.norm[i]);
    }
    fflush(iosba);
    fclose(iosba);

    FILE* iosbb;
    iosbb = fopen("headergrp.sbb","wb");
    int dummyint;
    double dummydouble;
    dummyint = sbbHeader.type;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = sbbHeader.dicType;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = sbbHeader.numSignal;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = sbbHeader.signalSize;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = sbbHeader.blockHop;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = sbbHeader.blockSize;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummydouble = sbbHeader.Fs;
    fwrite(&dummydouble, sizeof(double), 1, iosbb);
    fclose(iosbb);


    int nb_amp, nb_rho, nb_phase, iAmpRange, iAtom;

    end_nbit_amp = 4;
    end_nbit_rho = 16;
    end_nbit_phase = 16;

    double*** DfTable[end_nbit_amp][end_nbit_rho][end_nbit_phase];

    for (nb_phase=0;nb_phase<end_nbit_phase;nb_phase++)
	{
		for (nb_rho=0;nb_rho<end_nbit_rho;nb_rho++)
		{
			for (nb_amp=0;nb_amp<end_nbit_amp;nb_amp++)
			{
				DfTable[nb_amp][nb_rho][nb_phase] = new double**[sbbHeader.numSignal];
				for (i=0;i<sbbHeader.numSignal;i++)
				{
					DfTable[nb_amp][nb_rho][nb_phase][i] = new double*[sbbHeader.numBlock];
					for (j=0;j<sbbHeader.numBlock;j++)
					{
						DfTable[nb_amp][nb_rho][nb_phase][i][j] = new double[nb_amp+1];
					}
				}
			}
		}
	}

    int** statFreqTable[end_nbit_amp][end_nbit_amp];
    int** statFreqTableQ[end_nbit_amp][end_nbit_amp];
    for (nb_amp=0;nb_amp<end_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<end_nbit_amp;iAmpRange++)
		{
			statFreqTable[nb_amp][iAmpRange] = new int*[sbbHeader.numSignal];
			statFreqTableQ[nb_amp][iAmpRange] = new int*[sbbHeader.numSignal];
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				statFreqTable[nb_amp][iAmpRange][i] = new int[sbbHeader.numBlock];
				statFreqTableQ[nb_amp][iAmpRange][i] = new int[sbbHeader.numBlock];
			}
		}
	}


    for (i=0;i<sbbHeader.numSignal;i++)
	{
		double norm = dataSignal->getNorm(i);
		for (j=0;j<sbbHeader.numBlock;j++)
		{

			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			numElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElement];

			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp,
															  max_amp,
															  min_rho,
															  max_rho,
															  min_phase,
															  max_phase);
			// To keep the same range of rho and phase for all blocks
			// min_rho = 0.0;
// 			max_rho = 4.0;
// 			min_phase = 0.0,
// 			max_phase = 2*pi;
			double L_amp = max_amp - min_amp;
			double L_rho = max_rho - min_rho;
			double L_phase = max_phase - min_phase;

			file_report << endl;
			file_report << "=====================" << endl;
			file_report << "Parameters range" << endl;
			file_report << "=====================" << endl;
			file_report << setw (10) << setfill(' ') << "Signal"<< " "
							<<  setw (10) << setfill(' ') << "Block"<< " "
							<<  setw (20) << setfill(' ') << "min_amp"<< " "
							<<  setw (20) << setfill(' ') << "max_amp"<< " "
							<<  setw (20) << setfill(' ') << "min_rho"<< " "
							<<  setw (20) << setfill(' ') << "max_rho"<< " "
							<<  setw (20) << setfill(' ') << "min_phase"<< " "
							<<  setw (20) << setfill(' ') << "max_phase"<< " "
							<< endl;

			file_report << setw (10) << setfill(' ') << i+1<< " "
		            <<  setw (10) << setfill(' ') << j+1<< " "
		            <<  setw (20) << setfill(' ') << min_amp<< " "
					<<  setw (20) << setfill(' ') << max_amp<< " "
			        <<  setw (20) << setfill(' ') << min_rho<< " "
			        <<  setw (20) << setfill(' ') << max_rho<< " "
			        <<  setw (20) << setfill(' ') << min_phase<< " "
			        <<  setw (20) << setfill(' ') << max_phase<< " "
			        << endl;

			file_report << endl;
			file_report << "================================" << endl;
			file_report << "Parameters range by amplitude" << endl;
			file_report << "================================" << endl;
			file_report << setw (10) << setfill(' ') << "Signal"<< " "
							<<  setw (10) << setfill(' ') << "Block"<< " "
							<<  setw (20) << setfill(' ') << "nb_amp" << " "
							<<  setw (20) << setfill(' ') << "numAmpRange"<< " "
							<<  setw (20) << setfill(' ') << "min_rho"<< " "
							<<  setw (20) << setfill(' ') << "max_rho"<< " "
							<<  setw (20) << setfill(' ') << "min_phase"<< " "
							<<  setw (20) << setfill(' ') << "max_phase"<< " "
							<< endl;

			int t=1 ;
			int indAmp=0;
			for (nb_amp=1;nb_amp<=end_nbit_amp;nb_amp++)
    		{
				// cout << "Press enter" << endl;
// 				int lala;
// 				cin >> lala;

				indAmp++;

    			int NAmpRange;
				double* ampRangeLimit;
				double* ampBar;
				double* meanCoefByAmp;
				double* min_rho_vec;
				double* max_rho_vec;
				double* min_phase_vec;
				double* max_phase_vec;
				CStructBook* structBookByAmp;
    			// Set amplitude range
				//NAmpRange = pow(2.0, static_cast<double>(nb_amp));
				NAmpRange = nb_amp;
				structBookByAmp = new CStructBookExp[NAmpRange];
				ampBar = new double[NAmpRange];
				meanCoefByAmp = new double[NAmpRange];
				ampRangeLimit = new double[NAmpRange+1];
				ampRangeLimit[0]=-0.001;

				min_rho_vec = new double[NAmpRange];
				max_rho_vec = new double[NAmpRange];
				min_phase_vec = new double[NAmpRange];
				max_phase_vec = new double[NAmpRange];


				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				    ampBar[NAmpRange-1-iAmpRange] =
				    	(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				    ((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSB,
																			numElement,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit
					double dummy;
					((CStructBookExp*)structBookByAmp)[iAmpRange].findParmRange(  dummy,
																				  dummy,
																				  min_rho_vec[iAmpRange],
																				  max_rho_vec[iAmpRange],
																				  min_phase_vec[iAmpRange],
																				  max_phase_vec[iAmpRange]);
					file_report << setw (10) << setfill(' ') << i+1<< " "
								<<  setw (10) << setfill(' ') << j+1<< " "
								<<  setw (20) << setfill(' ') << nb_amp<< " "
								<<  setw (20) << setfill(' ') << iAmpRange+1 << " "
								<<  setw (20) << setfill(' ') << min_rho_vec[iAmpRange]<< " "
								<<  setw (20) << setfill(' ') << max_rho_vec[iAmpRange]<< " "
								<<  setw (20) << setfill(' ') << min_phase_vec[iAmpRange]<< " "
								<<  setw (20) << setfill(' ') << max_phase_vec[iAmpRange]<< " "
								<< endl;
					int ntotal =  ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();
					strtContinuousExp* pSBAux = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
					meanCoefByAmp[iAmpRange] = 0.0;
					if (ntotal!=0)
					{
						for (int nel=0;nel<ntotal;nel++)
						{
							meanCoefByAmp[iAmpRange] += pSBAux[nel].innerProduct;
						}

						meanCoefByAmp[iAmpRange] = meanCoefByAmp[iAmpRange]/ntotal;
					}
				}

				step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp), min_amp, max_amp);

				for (nb_rho=1;nb_rho<=end_nbit_rho;nb_rho++)
				{
					for (nb_phase=1;nb_phase<=end_nbit_phase;nb_phase++)
					{

						// quantExp[t].nb_amp = nb_amp;
// 						quantExp[t].nb_rho = nb_rho;
// 						quantExp[t].nb_phase = nb_phase;
// 						quantExp[t].nb_xi = nb_xi;
// 						quantExp[t].nb_sample = nb_sample;
						cout << "Signal: " << i+1 << " -- " << "Block: " << j+1 << endl;
						cout << "Quantizer: " << t << " " << nb_amp << " " << nb_rho << " " << nb_phase << endl;
						cout << "step_rho: " << step_rho << "; step_phase: " << step_phase << endl;
 						t++;
						for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
						{
							cout << "- AmpRange: " << iAmpRange+1 << endl;
							step_rho = computeMidRiseQuantStep(  	static_cast<double>(nb_rho),
																	min_rho_vec[iAmpRange],
																	max_rho_vec[iAmpRange]);
							step_phase = computeMidRiseQuantStep( static_cast<double>(nb_phase),
																	min_phase_vec[iAmpRange],
																	max_phase_vec[iAmpRange]);

							strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
							int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

							statFreqTable[nb_amp-1][iAmpRange][i][j] = numElementAmp;

							int numElementAmpQ=0;
							double distAtom=0.0;
							double meanDistAtom=0.0;
							double distCoef=0.0;
							double meanDistCoef=0.0;
							if (numElementAmp!=0)
							{
								for (iAtom=0;iAtom<numElementAmp;iAtom++)
								{
									pSBAmpQ[iAtom] =  quantizeAtomExp(	step_amp,step_rho, step_phase,
																		min_amp, min_rho_vec[iAmpRange], min_phase_vec[iAmpRange],
																		pSBAmp[iAtom]);

									if (pSBAmpQ[iAtom].innerProduct!=0.0)
									{
										cout << "DRho: " << fabs(pSBAmp[iAtom].rho - pSBAmpQ[iAtom].rho) << endl;
										synthSignalAtomExp(	blockSignal,
															sbbHeader.blockSize,
															pSBAmp[iAtom]);

										synthSignalAtomExp(	recSignal,
															sbbHeader.blockSize,
															pSBAmpQ[iAtom]);

										distAtom += sbbHeader.blockSize*computeMSE(blockSignal,recSignal,sbbHeader.blockSize);

										distCoef += (pSBAmp[iAtom].innerProduct-pSBAmpQ[iAtom].innerProduct)*
													(pSBAmp[iAtom].innerProduct-pSBAmpQ[iAtom].innerProduct);
										numElementAmpQ++;
									}
								}
								meanDistAtom =  distAtom/numElementAmp;
								meanDistCoef =  distCoef/numElementAmp;
							}

							statFreqTableQ[nb_amp-1][iAmpRange][i][j] = numElementAmpQ;

							// cout << "numElementAmp" << ";" <<  "distAtom" << ";" << "meanDistAtom" << endl;
                            // cout << numElementAmp << ";" <<  distAtom << ";" << meanDistAtom << endl;

							DfTable[nb_amp-1][nb_rho-1][nb_phase-1][i][j][iAmpRange] = meanDistAtom;

							file_rdbyamp 	<< setw (10) << setfill(' ') << i+1 << " "
											<<  setw (10) << setfill(' ') << j+1 << " "
											<<  setw (10) << setfill(' ') << indAmp << " "
											<<  setw (10) << setfill(' ') << nb_amp << " "
											<<  setw (10) << setfill(' ') << nb_rho << " "
											<<  setw (10) << setfill(' ') << nb_phase<< " "
											<<  setw (10) << setfill(' ') << nb_xi << " "
											<<  setw (10) << setfill(' ') << nb_sample << " "
											<<  setw (20) << setfill(' ') << iAmpRange+1<< " "
											<<  setw (20) << setfill(' ') << ampRangeLimit[iAmpRange] << " "
											<<  setw (20) << setfill(' ') << ampRangeLimit[iAmpRange+1]<< " "
											<<  setw (20) << setfill(' ') << ampBar[iAmpRange]<< " "
											<<  setw (20) << setfill(' ') << meanCoefByAmp[iAmpRange] << " "
											<<  setw (20) << setfill(' ') << meanDistAtom << " "
											<<  setw (20) << setfill(' ') << meanDistCoef << " "
											<< endl;
						}
					} // nb_phase
				} // nb_rho
				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;
				delete [] meanCoefByAmp;

				delete [] min_rho_vec;
				delete [] max_rho_vec;
				delete [] min_phase_vec;
				delete [] max_phase_vec;
			} // nb_amp
			delete [] pSBAmpQ;
		} // block
	} // signal


	file_report << endl << endl;
	file_report << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	file_report << "Frequency Tables - Non quantized atoms" << endl;
	file_report << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			file_report << "--------------------------------------------------------------" << endl;
			file_report << setw (20) << setfill(' ') << "Signal"<< " :"
						<< setw (10) << setfill(' ') << i+1 << endl
						<< setw (20) << setfill(' ') << "Block"<< " :"
						<< setw (10) << setfill(' ') << j+1 << endl;
			file_report << "--------------------------------------------------------------" << endl;
			file_report << setw (20) << setfill(' ') << "nbamp | nAmpRange" << " " ;
			for (nb_amp=1;nb_amp<=end_nbit_amp;nb_amp++)
    		{
    			file_report << setw (20) << setfill(' ') <<  nb_amp << " ";
    		}
    		file_report << endl;
    		for (nb_amp=1;nb_amp<=end_nbit_amp;nb_amp++)
    		{
    			file_report << setw (20) << setfill(' ') <<  nb_amp << " ";
    			for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
				{
    				file_report << setw (20) << setfill(' ') <<  statFreqTable[nb_amp-1][iAmpRange][i][j] << " ";
    			}
    			file_report << endl;
    		}
    	}
    }

    file_report << endl << endl;
	file_report << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	file_report << "Frequency Tables - Quantized atoms" << endl;
	file_report << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			file_report << "--------------------------------------------------------------" << endl;
			file_report << setw (20) << setfill(' ') << "Signal"<< " :"
						<< setw (10) << setfill(' ') << i+1 << endl
						<< setw (20) << setfill(' ') << "Block"<< " :"
						<< setw (10) << setfill(' ') << j+1 << endl;
			file_report << "--------------------------------------------------------------" << endl;
			file_report << setw (20) << setfill(' ') << "nbamp | nAmpRange" << " " ;
			for (nb_amp=1;nb_amp<=end_nbit_amp;nb_amp++)
    		{
    			file_report << setw (20) << setfill(' ') <<  nb_amp << " ";
    		}
    		file_report << endl;
    		for (nb_amp=1;nb_amp<=end_nbit_amp;nb_amp++)
    		{
    			file_report << setw (20) << setfill(' ') <<  nb_amp << " ";
    			for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
				{
    				file_report << setw (20) << setfill(' ') <<  statFreqTableQ[nb_amp-1][iAmpRange][i][j] << " ";
    			}
    			file_report << endl;
    		}
    	}
    }



	file_report << endl << endl;
	file_report << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	file_report << "Atom Distortion Tables" << endl;
	file_report << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

	iosbb = fopen("rdbyamp_atomdist.sbb","wb");
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (nb_amp=1;nb_amp<=end_nbit_amp;nb_amp++)
    		{
    			for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
				{
    				file_report << "--------------------------------------------------------------" << endl;
    				file_report << setw (20) << setfill(' ') << "Signal"<< " :"
    							<< setw (10) << setfill(' ') << i+1 << endl
    							<< setw (20) << setfill(' ') << "Block"<< " :"
    							<< setw (10) << setfill(' ') << j+1 << endl
    							<< setw (20) << setfill(' ') << "nb_amp"<< " :"
    							<< setw (10) << setfill(' ') << nb_amp << endl
    							<< setw (20) << setfill(' ') << "iAmpRange"<< " :"
    							<< setw (10) << setfill(' ') << iAmpRange+1 << endl;
    				file_report << "--------------------------------------------------------------" << endl;
    				file_report << setw (20) << setfill(' ') << "R_rho | R_phi" << " " ;
    				for (nb_phase=1;nb_phase<=end_nbit_phase;nb_phase++)
					{
						file_report << setw (20) << setfill(' ') << nb_phase << " " ;
					}
					file_report << endl;
					for (nb_rho=1;nb_rho<=end_nbit_rho;nb_rho++)
					{
						file_report << setw (20) << setfill(' ') << nb_rho << " ";
						for (nb_phase=1;nb_phase<=end_nbit_phase;nb_phase++)
						{
							file_report << setw (20) << setfill(' ') << DfTable[nb_amp-1][nb_rho-1][nb_phase-1][i][j][iAmpRange] << " ";
							dummydouble = DfTable[nb_amp-1][nb_rho-1][nb_phase-1][i][j][iAmpRange];
						    fwrite(&dummydouble, sizeof(double), 1, iosbb);
						}
						file_report << endl;
					}
				}
    		}
		}
	}
	fclose(iosbb);
	// close files
	file_rdbyamp.close();
	file_report.close();


	delete [] recSignal;
	delete [] blockSignal;

	for (nb_phase=0;nb_phase<end_nbit_phase;nb_phase++)
	{
		for (nb_rho=0;nb_rho<end_nbit_rho;nb_rho++)
		{
			for (nb_amp=0;nb_amp<end_nbit_amp;nb_amp++)
			{
				for (i=0;i<sbbHeader.numSignal;i++)
				{
					for (j=0;j<sbbHeader.numBlock;j++)
					{
						delete [] DfTable[nb_amp][nb_rho][nb_phase][i][j];
					}
					delete [] DfTable[nb_amp][nb_rho][nb_phase][i];
				}
				delete [] DfTable[nb_amp][nb_rho][nb_phase];
			}
		}
	}

	for (nb_amp=0;nb_amp<end_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<end_nbit_amp;iAmpRange++)
		{
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				delete [] statFreqTable[nb_amp][iAmpRange][i];
				delete [] statFreqTableQ[nb_amp][iAmpRange][i];
			}
			delete [] statFreqTable[nb_amp][iAmpRange];
			delete [] statFreqTableQ[nb_amp][iAmpRange];
		}
	}
}

void calcTableRDExpByAmp(	CFileRDBitRange* rdBitRange,
							CDataSignal* dataSignal,
							strtSBBHeader sbbHeader,
							CStructBook** structBook,
							int Nfreq,
							char* InputFile)
{

    // Setting fixed quantizers
    int nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
    int nb_sample = (int)ceil( log(sbbHeader.blockSize) / log(2.0) );
    int nb_rhoSign = 1;

	// variables
    double* recSignal;
    double* blockSignal;
    recSignal = new double[sbbHeader.blockSize];
    blockSignal = new double[sbbHeader.blockSize];

	// counters
	int i, j, k, i_qrho, i_qphi;
	/////////////////


    strtContinuousExp* pSB;
    int numElement;

    double** min_amp= new double*[sbbHeader.numSignal];
    double** max_amp= new double*[sbbHeader.numSignal];
    double** min_rho= new double*[sbbHeader.numSignal];
    double** max_rho= new double*[sbbHeader.numSignal];
    double** min_phase= new double*[sbbHeader.numSignal];
    double** max_phase= new double*[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		min_amp[i]= new double[sbbHeader.numBlock];
    	max_amp[i]= new double[sbbHeader.numBlock];
    	min_rho[i]= new double[sbbHeader.numBlock];
    	max_rho[i]= new double[sbbHeader.numBlock];
    	min_phase[i]= new double[sbbHeader.numBlock];
    	max_phase[i]= new double[sbbHeader.numBlock];
	}

    int nb_amp, nb_rho, nb_phase, iAmpRange, iAtom;

    double 	step_amp;
    ///////////////////
    int max_nbit_amp = rdBitRange->getEndNbitAmp(); // <<<<<<<<    change here
    int min_nbit_amp = rdBitRange->getInitNbitAmp();
	///////////////////
	double delta_nb;
    double min_rho_ref = 0.0;
    double max_rho_ref = 4.0;

    double NbRhoMax = 22.0;
    delta_nb = 0.25; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    change here
    int N_qrho = static_cast<int>(round((NbRhoMax-1.0)/delta_nb)) + 1;
    double* step_rho = new double[N_qrho];

    j=0;
    double nb = 0.0;
    for (nb=NbRhoMax;nb>=1.0;nb=nb-delta_nb)
	{

		step_rho[j] = computeMidRiseQuantStep( nb,
												 min_rho_ref,
												 max_rho_ref);
		j++;
	}
	///////////======

    ////
    double min_phase_ref = 0;
    double max_phase_ref = 2*pi;
    double NbPhiMax = 18.0;
    delta_nb = 0.25; //	 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    change here
    int N_qphi = static_cast<int>(round((NbPhiMax-1.0)/delta_nb)) + 1;
    double* step_phase = new double[N_qphi];

    j=0;
    for (nb=NbPhiMax;nb>=1.0;nb=nb-delta_nb)
	{
		step_phase[j] = computeMidRiseQuantStep( nb,
												   min_phase_ref,
												   max_phase_ref);
		j++;
	}


	///////////==============================================================
    double**** DfTable;

    DfTable = new double***[max_nbit_amp-min_nbit_amp+1];
    for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		DfTable[nb_amp-min_nbit_amp] = new double**[N_qrho];
		for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
		{
			DfTable[nb_amp-min_nbit_amp][i_qrho] = new double*[N_qphi];
			for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
			{
				DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi]= new double[nb_amp];
			}
		}
	}
	for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
		{
			for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
			{
				for (i=0;i<nb_amp;i++)
				{
					DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi][i]= 0.0;
				}
			}
		}
	}

    int**** statFreqTable = new int***[max_nbit_amp-min_nbit_amp+1];
    int**** statFreqTableQ = new int***[max_nbit_amp-min_nbit_amp+1];
    for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		statFreqTable[nb_amp-min_nbit_amp] =  new int**[nb_amp];
		statFreqTableQ[nb_amp-min_nbit_amp]=  new int**[nb_amp];
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			statFreqTable[nb_amp-min_nbit_amp][iAmpRange] = new int*[sbbHeader.numSignal];
			statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange] = new int*[sbbHeader.numSignal];
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				statFreqTable[nb_amp-min_nbit_amp][iAmpRange][i] = new int[sbbHeader.numBlock];
				statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange][i] = new int[sbbHeader.numBlock];
			}
		}
	}


    for (i=0;i<sbbHeader.numSignal;i++)
	{
		double norm = dataSignal->getNorm(i);
		for (j=0;j<sbbHeader.numBlock;j++)
		{

			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			numElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElement];

			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp[i][j],
															  max_amp[i][j],
															  min_rho[i][j],
															  max_rho[i][j],
															  min_phase[i][j],
															  max_phase[i][j]);

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phase = max_phase[i][j] - min_phase[i][j];

			int t=1 ;
			int indAmp=0;
			for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
    		{
    			cout 	<< "Signal: " << i+1 << " -- "
    					<< "Block: " << j+1 << " -- "
    					<< "nbamp: " << nb_amp <<endl;
				// cout << "Press enter" << endl;
// 				int lala;
// 				cin >> lala;

				indAmp++;

    			int NAmpRange;
				double* ampRangeLimit;
				double* ampBar;
				double* meanCoefByAmp;
				double* min_rho_vec;
				double* max_rho_vec;
				double* min_phase_vec;
				double* max_phase_vec;
				CStructBook* structBookByAmp;
    			// Set amplitude range
				//NAmpRange = pow(2.0, static_cast<double>(nb_amp));
				NAmpRange = nb_amp;
				structBookByAmp = new CStructBookExp[NAmpRange];
				ampBar = new double[NAmpRange];
				meanCoefByAmp = new double[NAmpRange];
				ampRangeLimit = new double[NAmpRange+1];
				ampRangeLimit[0]=-0.001;

				min_rho_vec = new double[NAmpRange];
				max_rho_vec = new double[NAmpRange];
				min_phase_vec = new double[NAmpRange];
				max_phase_vec = new double[NAmpRange];


				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				    ampBar[NAmpRange-1-iAmpRange] =
				    	(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				    ((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSB,
																			numElement,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit
					///////////////
					double dummy;
					((CStructBookExp*)structBookByAmp)[iAmpRange].findParmRange(  dummy,
																				  dummy,
																				  min_rho_vec[iAmpRange],
																				  max_rho_vec[iAmpRange],
																				  min_phase_vec[iAmpRange],
																				  max_phase_vec[iAmpRange]);
					////////////////
					int ntotal =  ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();
					strtContinuousExp* pSBAux = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
					meanCoefByAmp[iAmpRange] = 0.0;
					if (ntotal!=0)
					{
						for (int nel=0;nel<ntotal;nel++)
						{
							meanCoefByAmp[iAmpRange] += pSBAux[nel].innerProduct;
						}

						meanCoefByAmp[iAmpRange] = meanCoefByAmp[iAmpRange]/ntotal;
					}
				}

				step_amp = computeMidRiseQuantStep(  static_cast<double>(nb_amp), min_amp[i][j], max_amp[i][j]);
				for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
				{
					for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
					{
 						t++;
						for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
						{
							strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
							int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

							statFreqTable[nb_amp-min_nbit_amp][iAmpRange][i][j] = numElementAmp;

							int numElementAmpQ=0;
							double distAtom=0.0;
							double meanDistAtom=0.0;
							double distCoef=0.0;
							double meanDistCoef=0.0;
							if (numElementAmp!=0)
							{
								for (iAtom=0;iAtom<numElementAmp;iAtom++)
								{
									pSBAmpQ[iAtom] =  quantizeAtomExp(	step_amp,step_rho[i_qrho], step_phase[i_qphi],
																		min_amp[i][j], min_rho[i][j], min_phase[i][j],
																		pSBAmp[iAtom]);

									if (pSBAmpQ[iAtom].innerProduct!=0.0)
									{

										double timeconstant;
										if  (fabs(pSBAmp[iAtom].rho) > 0)
											timeconstant = 1.0/fabs(pSBAmp[iAtom].rho);
										else
											timeconstant = 1e100;
										// =====================
										// If the timeconstant  or time support less than 500 samples uses
										// the numerical distortion computation
										// Otherwise use Closed form of the atom distortion
										// if (  ( (pSBAmp[iAtom].b -pSBAmp[iAtom].a)<=500) || (timeconstant<=500) )
// 										{
// 											synthSignalAtomExp(	blockSignal,
// 															sbbHeader.blockSize,
// 															pSBAmp[iAtom]);
//
// 											synthSignalAtomExp(	recSignal,
// 																sbbHeader.blockSize,
// 																pSBAmpQ[iAtom]);
//
// 											distAtom += computeSqrError(blockSignal,
// 																		recSignal,
// 																		pSBAmp[iAtom].a,
// 																		pSBAmp[iAtom].b);
// 										}
// 										else
// 										{
											distAtom += calcAtomDistFunction(pSBAmp[iAtom],pSBAmpQ[iAtom],sbbHeader.blockSize);
										//}

										distCoef += (pSBAmp[iAtom].innerProduct-pSBAmpQ[iAtom].innerProduct)*
													(pSBAmp[iAtom].innerProduct-pSBAmpQ[iAtom].innerProduct);

										numElementAmpQ++;


									}
								}
								meanDistAtom =  distAtom/numElementAmpQ;
								meanDistCoef =  distCoef/numElementAmpQ;
							}

							statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange][i][j] = numElementAmpQ;
							DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi][iAmpRange] += distAtom;

						}
					} // nb_phase
				} // nb_rho
				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;
				delete [] meanCoefByAmp;

				delete [] min_rho_vec;
				delete [] max_rho_vec;
				delete [] min_phase_vec;
				delete [] max_phase_vec;
			} // nb_amp
			delete [] pSBAmpQ;
		} // block
	} // signal

	// Mean DFTable	among signals and blocks
	for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			//////
			int ntotal = 0;
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				for (j=0;j<sbbHeader.numBlock;j++)
				{
					ntotal += statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange][i][j];
				}
			}
			/////
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
				{
					DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi][iAmpRange] =
											DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi][iAmpRange]/
											static_cast<double>(ntotal);
				}
			}
		}
	}


	//=======================================================================
	//
	// Save table information in binary file to be loaded through MATLAB
	//
	//

    int dummyint;
    double dummydouble;

	char fileName[_MAX_PATH];
	char aux[_MAX_PATH];
	char* pos;
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
	sprintf(aux, "_dftable_byamp_nbamp%d-%d.bin",min_nbit_amp,max_nbit_amp);
	strcpy( &pos[0], aux);
	cout << fileName << endl;
	FILE* iobin = fopen(fileName,"wb");
	// Save header information
	//  - number of signals
	// 	- number of blocks
	//	- signal size
	//	- block hop
	//	- block size
	//	- maximum number of bits for amplitude
	//	- number of rho (decaying) quantization steps
	//	- number of phase quantization steps
	//
	dummyint = sbbHeader.numSignal;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = sbbHeader.numBlock;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = sbbHeader.signalSize;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = sbbHeader.blockHop;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = sbbHeader.blockSize;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = max_nbit_amp;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = min_nbit_amp;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = N_qrho;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	dummyint = N_qphi;
	fwrite(&dummyint, sizeof(int), 1, iobin);

	// Save parameters ranges
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			dummydouble =  min_amp[i][j];
			fwrite(&dummydouble, sizeof(double), 1, iobin);
			dummydouble =  max_amp[i][j];
			fwrite(&dummydouble, sizeof(double), 1, iobin);
			dummydouble =  min_rho[i][j];
			fwrite(&dummydouble, sizeof(double), 1, iobin);
			dummydouble =  max_rho[i][j];
			fwrite(&dummydouble, sizeof(double), 1, iobin);
			dummydouble = min_phase[i][j];
			fwrite(&dummydouble, sizeof(double), 1, iobin);
			dummydouble = max_phase[i][j];
			fwrite(&dummydouble, sizeof(double), 1, iobin);
		}
	}
	// Save rho and phase quatization steps
	for (i=0;i<N_qrho;i++)
	{
		dummydouble = step_rho[i];
		fwrite(&dummydouble, sizeof(double), 1, iobin);

	}
	for (i=0;i<N_qphi;i++)
	{
		dummydouble = step_phase[i];
		fwrite(&dummydouble, sizeof(double), 1, iobin);
	}

	// Save frequency table of non quantized atoms
    for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			//////
			int ntotal = 0;
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				for (j=0;j<sbbHeader.numBlock;j++)
				{
					ntotal += statFreqTable[nb_amp-min_nbit_amp][iAmpRange][i][j];
				}
			}
			/////
			fwrite(&ntotal, sizeof(int), 1, iobin);
		}
	}

    // Save frequency table of quantized atoms (remind the dead zone)
    for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			//////
			int ntotal = 0;
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				for (j=0;j<sbbHeader.numBlock;j++)
				{
					ntotal += statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange][i][j];
				}
			}
			/////
			fwrite(&ntotal, sizeof(int), 1, iobin);
		}
	}

    // Save the mean atom distortion atom for each amplitude range
	for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
				{
					dummydouble = DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi][iAmpRange];
					fwrite(&dummydouble, sizeof(double), 1, iobin);
				}
			}
		}
	}
	fclose(iobin);

	/////////////////////////////////
	// DEALLOCATE
	delete [] recSignal;
	delete [] blockSignal;
	for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
		{
			for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
			{
				delete [] DfTable[nb_amp-min_nbit_amp][i_qrho][i_qphi];
			}
			delete [] DfTable[nb_amp-min_nbit_amp][i_qrho];
		}
		delete [] DfTable[nb_amp-min_nbit_amp];
	}
	delete [] DfTable;
	for (nb_amp=min_nbit_amp;nb_amp<=max_nbit_amp;nb_amp++)
	{
		for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
		{
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				delete [] statFreqTable[nb_amp-min_nbit_amp][iAmpRange][i];
				delete [] statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange][i];
			}
			delete [] statFreqTable[nb_amp-min_nbit_amp][iAmpRange];
			delete [] statFreqTableQ[nb_amp-min_nbit_amp][iAmpRange];
		}
		delete [] statFreqTable[nb_amp-min_nbit_amp];
		delete [] statFreqTableQ[nb_amp-min_nbit_amp];
	}
	delete [] statFreqTable;
	delete [] statFreqTableQ;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] min_amp[i];
    	delete [] max_amp[i];
    	delete [] min_rho[i];
    	delete [] max_rho[i];
    	delete [] min_phase[i];
    	delete [] max_phase[i];
	}
	delete [] min_amp;
	delete [] max_amp;
	delete [] min_rho;
	delete [] max_rho;
	delete [] min_phase;
	delete [] max_phase;
	///======================================
	delete [] step_rho;
	delete [] step_phase;
}

double calcAtomDistFunction(strtContinuousExp sb,
							strtContinuousExp sbQ,
							int sigSize)
{
	double df, df1,df2, sqrnorm;
	double innerProd = sb.innerProduct;
    double rho = sb.rho;
    double xi = sb.xi;
    double phase = sb.phase;
    int a = sb.a;
    int b = sb.b;

    double rhoq = sbQ.rho;
    double phaseq = sbQ.phase;

    double kappa2, kappaq2;

	if (a!=b)
	{

		if (rho<0.0)
		{
			rho = -rho;
			phase = - phase;
			rhoq = -rhoq;
			phaseq = - phaseq;
			a = sigSize-1 - sb.b;
			b = sigSize-1 - sb.a;
		}

		kappa2 =  indefIntegralAtomSqrNorm(rho,xi,phase,static_cast<double>(b-a)) -
              indefIntegralAtomSqrNorm(rho,xi,phase,0);
		kappaq2 =  indefIntegralAtomSqrNorm(rhoq,xi,phaseq,static_cast<double>(b-a)) -
				   indefIntegralAtomSqrNorm(rhoq,xi,phaseq,0);

		df = (1/kappa2)  * ( indefIntegralDfa(rho, xi, phase, static_cast<double>(b-a)) -
							 indefIntegralDfa(rho, xi, phase,  0   )) +
			 (1/kappaq2) *(  indefIntegralDfa(rhoq, xi, phaseq, static_cast<double>(b-a)) -
							 indefIntegralDfa(rhoq, xi, phaseq,  0 )) -
			 (1/(sqrt(kappaq2)* sqrt(kappa2)))*
							(   indefIntegralDfb(rho,rhoq, xi, phase, phaseq,static_cast<double>(b-a)) -
								indefIntegralDfb(rho,rhoq, xi, phase, phaseq, 0   )) -
			 (1/(sqrt(kappaq2)* sqrt(kappa2)))*
							(   indefIntegralDfc(rho,rhoq, phase, phaseq,static_cast<double>(b-a)) -
								indefIntegralDfc(rho,rhoq, phase, phaseq, 0   ));
	}
	else
	{
		if ( ( (cos(phase) > 0.0) && ( cos(phaseq) >0.0) ) ||
        	 ( (cos(phase) < 0.0) && ( cos(phaseq) <0.0) ) )
        {
        	df = 0.0;
        }
    	else if ( ( ( cos(phase) == 0.0) && ( cos(phaseq) != 0.0) ) ||
             	  ( ( cos(phase) != 0.0) && ( cos(phaseq) == 0.0) ) )
        {
        	df = 1.0;
        }
    	else
    	{
        	df = 4.0;
		}
	}
	return df;
}

double indefIntegralAtomSqrNorm(double rho, double xi, double phase, double t)
{
	double f;
	if ( (xi!=0) && (rho!=0) )
	{
		f =- exp(-2*rho*t)*(rho*rho*(cos(2*phase + 2*t*xi) + 1) + (xi*xi) - rho*xi*sin(2*phase + 2*t*xi)) /
		     (4*(rho*rho*rho) + 4*rho*(xi*xi));
	}
	else if ( (xi!=0) && (rho==0) )
	{
		f = t/2 + sin(2*phase + 2*t*xi)/(4*xi);
	}
	else if ( (xi==0) && (rho!=0) )
	{
		f = -exp(-2*rho*t)*(cos(phase)*cos(phase))/(2*rho);
	}
	else
	{
		f = t*cos(phase)*cos(phase);
	}
	return f;
}

double indefIntegralDfa(double rho, double xi, double phi, double t)
{
	double f;
	if ( (xi!=0) && (rho!=0) )
	{
		f = - exp(-2*rho*t)*((rho*rho)*(cos(2*phi + 2*t*xi) + 1) + (xi*xi) - rho*xi*sin(2*phi + 2*t*xi)) /
              (4*(rho*rho*rho) + 4*rho*(xi*xi));
	}
	else if ( (xi!=0) && (rho==0) )
	{
		f = t/2 + sin(2*phi + 2*t*xi)/(4*xi);
	}
	else if ( (xi==0) && (rho!=0) )
	{
		f = -(cos(phi)*cos(phi))*exp(-2*rho*t)/(2*rho);
	}
	else
	{
		f = t*(cos(phi)*cos(phi));
	}
	return f;
}

double indefIntegralDfb(double rho, double rhoq, double xi, double phi, double phiq, double t)
{
	double f;
	if ( (xi!=0) && ((rho+rhoq)!=0) )
	{
		f = - exp(-rho*t - rhoq*t)*(rho*cos(phi + phiq + 2*t*xi) + rhoq*cos(phi + phiq + 2*t*xi) - 2*xi*sin(phi + phiq + 2*t*xi)) /
              (((rho*rho) + 2*rho*rhoq + (rhoq*rhoq) + 4*(xi*xi)));
	}
	else if ( (xi!=0) && ((rho+rhoq)==0) )
	{
		f = sin(phi + phiq + 2*t*xi)/(2*xi);
	}
	else if ( (xi==0) && ((rho+rhoq)!=0) )
	{
		f = - exp(-rho*t - rhoq*t)*cos(phi + phiq)/(rho + rhoq);
	}
	else
	{
		f = t*cos(phi + phiq);
	}
	return f;
}

double indefIntegralDfc(double rho, double rhoq, double phi, double phiq, double t)
{
	double f;
	if ( (rho+rhoq)!=0)
	{
		f = - exp(-rho*t - rhoq*t)*cos(phi - phiq)/(rho + rhoq);
	}
	else
	{
		f = t*cos(phi - phiq);
	}
	return f;
}


double computeMidThreadQuantStep( double nb, double min, double max)
{
	double step;
    double nlevel = (pow(2.0, nb ) - 1.0);
    if (nlevel==0.0)
        step = 1e+30;
    else
        step = fabs( (max - min) / (nlevel-1.0) );

    return step;
}

double computeMidRiseQuantStep( double nb, double min, double max)
{
	double step;
    double nlevel = (pow(2.0, nb ) );
    if (nlevel==0.0)
        step = 1e+30;
    else
        step = fabs( (max - min) / (nlevel-1.0) );

    return step;
}


strtContinuousExp quantizeAtomExp(	double step_amp,double step_rho, double step_phase,
									double min_amp,double min_rho, double min_phase,
									strtContinuousExp pSB)
{
	strtContinuousExp pSBQ;

	pSBQ.innerProduct = linearQuant(pSB.innerProduct-min_amp,step_amp) + min_amp;
	pSBQ.phase = linearQuant(pSB.phase-min_phase,step_phase) + min_phase;
	if (pSB.rho<0)
		pSBQ.rho = - ( linearQuant(fabs(pSB.rho) - min_rho,step_rho) + min_rho );
	else
		pSBQ.rho = linearQuant(fabs(pSB.rho)-min_rho,step_rho) + min_rho;

	pSBQ.xi = pSB.xi;
	pSBQ.a = pSB.a;
	pSBQ.b = pSB.b;

	return pSBQ;
}

void synthSignalAtomExp( double* recSignal,
						 int sigSize,
						 strtContinuousExp sb)
{
    CExpDictionary* expDic = new CExpDictionary;
    expDic->setSignalSize(sigSize);
    CParameter* parm;
    parm = new CExpParm;
    int i,j;
    double* sigAux;
 //    for (i=0;i<sigSize;i++)
//     {
//         recSignal[i] = 0.0;
//     }
	((CExpParm*)parm)->setParm(sb);
	// ((CExpParm*)parm)->printParm2Screen();
	expDic->setRealAtomOnSupport(parm);
	sigAux = expDic->getRealAtom();
	for (j=sb.a;j<=sb.b;j++)
	{
		recSignal[j] = sigAux[j];
	}
    delete expDic;
    delete parm;
}

void groupBlockExp(	strtSBBHeader sbbHeader, CStructBook** structBook,
					strtSBBHeader sbbHeaderGroup, CStructBook** structBookGroup,
					int numBlockPerGroup, int iSignal)
{
	int i,j,k,ind, ishift;
	strtContinuousExp* sb,*sb_aux;
	int sbNumElement;
	int grpBlock;
	int prevAtom,nextAtom;
	CExpDictionary expDic;
	CExpDictionary expDicGrp;
	strtContinuousExp sbgrp;
	double* rAtomBlock;
	double prevRho;

	cgMatrix<double> cgSignal(1,sbbHeaderGroup.blockSize,0.0);
	cgMatrix<double> innerProd(1,1,0.0);
	expDic.setSignalSize(sbbHeader.blockSize);
	expDicGrp.setSignalSize(sbbHeaderGroup.blockSize);

#ifdef DEBUG_GRPBLOCK
	cout << "- Saving atom grouping report ... " << endl;
	int numGrpEle =0;
	fstream file_grprep("strbookgrp_rep.out",ios::out);
	file_grprep  << setw (20) << setfill(' ') << "numGrpEle"<<" "
				 << setw (20) << setfill(' ') << "Block" <<" "
				 << setw (20) << setfill(' ') << "Coef" <<" "
				 << setw (20) << setfill(' ') << "Rho" <<" "
				 << setw (20) << setfill(' ') << "Freq" <<" "
				 << setw (20) << setfill(' ') << "Phi" <<" "
				 << setw (20) << setfill(' ') << "A" <<" "
				 << setw (20) << setfill(' ') << "B" <<" " << endl;
    cout  << setw (20) << setfill(' ') << "numGrpEle"<<" "
				 << setw (20) << setfill(' ') << "Block" <<" "
				 << setw (20) << setfill(' ') << "Coef" <<" "
				 << setw (20) << setfill(' ') << "Rho" <<" "
				 << setw (20) << setfill(' ') << "Freq" <<" "
				 << setw (20) << setfill(' ') << "Phi" <<" "
				 << setw (20) << setfill(' ') << "A" <<" "
				 << setw (20) << setfill(' ') << "B" <<" " << endl;
	for (i=sbbHeader.numBlock-1; i>=0; i--)
    {
        sb = ((CStructBookExp*)structBook[iSignal])[i].getStructBook();
        sbNumElement = ((CStructBookExp*)structBook[iSignal])[i].getNumElement();

        cout << iSignal << ";"<< i << ";" << sbNumElement << endl;

        grpBlock = (int) (floor(static_cast<double>(i)/static_cast<double>(numBlockPerGroup)) );

        for (j=0;j<sbNumElement;j++)
        {
			//cout << "Element: " << j+1 << endl;
            prevAtom = sb[j].prevAtom;
            nextAtom = sb[j].nextAtom;

            if (nextAtom==-1)
            {
            	if (prevAtom==-1)
            	{

					file_grprep  << setw (20) << setfill(' ') << numGrpEle <<" "
								 << setw (20) << setfill(' ') << i <<" "
								 << setw (20) << setfill(' ') << sb[j].innerProduct <<" "
								 << setw (20) << setfill(' ') << sb[j].rho <<" "
								 << setw (20) << setfill(' ') << sb[j].xi <<" "
								 << setw (20) << setfill(' ') << sb[j].phase <<" "
								 << setw (20) << setfill(' ') << sb[j].a <<" "
								 << setw (20) << setfill(' ') << sb[j].b <<" " << endl;
					cout  << setw (20) << setfill(' ') << numGrpEle <<" "
								 << setw (20) << setfill(' ') << i <<" "
								 << setw (20) << setfill(' ') << sb[j].innerProduct <<" "
								 << setw (20) << setfill(' ') << sb[j].rho <<" "
								 << setw (20) << setfill(' ') << sb[j].xi <<" "
								 << setw (20) << setfill(' ') << sb[j].phase <<" "
								 << setw (20) << setfill(' ') << sb[j].a <<" "
								 << setw (20) << setfill(' ') << sb[j].b <<" " << endl;

					numGrpEle++;


            	}
            	else
            	{


					file_grprep  << setw (20) << setfill(' ') << numGrpEle<<" "
								 << setw (20) << setfill(' ') << i <<" "
								 << setw (20) << setfill(' ') << sb[j].innerProduct <<" "
								 << setw (20) << setfill(' ') << sb[j].rho <<" "
								 << setw (20) << setfill(' ') << sb[j].xi <<" "
								 << setw (20) << setfill(' ') << sb[j].phase <<" "
								 << setw (20) << setfill(' ') << sb[j].a <<" "
								 << setw (20) << setfill(' ') << sb[j].b <<" " << endl;
					cout  << setw (20) << setfill(' ') << numGrpEle<<" "
								 << setw (20) << setfill(' ') << i <<" "
								 << setw (20) << setfill(' ') << sb[j].innerProduct <<" "
								 << setw (20) << setfill(' ') << sb[j].rho <<" "
								 << setw (20) << setfill(' ') << sb[j].xi <<" "
								 << setw (20) << setfill(' ') << sb[j].phase <<" "
								 << setw (20) << setfill(' ') << sb[j].a <<" "
								 << setw (20) << setfill(' ') << sb[j].b <<" " << endl;


					ishift = 1;
					while ( (prevAtom!=-1) && ((i-ishift)>= grpBlock*numBlockPerGroup ) )
					{
						sb_aux = ((CStructBookExp*)structBook[iSignal])[i-ishift].getStructBook();
						sb_aux[prevAtom].nextAtom = i-ishift+1;


						file_grprep  << setw (20) << setfill(' ') << numGrpEle <<" "
									 << setw (20) << setfill(' ') << i-ishift <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].innerProduct <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].rho <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].xi <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].phase <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].a <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].b <<" " << endl;
			             cout  << setw (20) << setfill(' ') << numGrpEle <<" "
									 << setw (20) << setfill(' ') << i-ishift <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].innerProduct <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].rho <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].xi <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].phase <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].a <<" "
									 << setw (20) << setfill(' ') << sb_aux[prevAtom].b <<" " << endl;

						ishift++;
						prevAtom = sb_aux[prevAtom].prevAtom;
					}
					numGrpEle++;
            	}
         	}
        }
    }
    file_grprep.close();
    cout << "- Atom grouping report saved ... " << endl;
#endif

	for (i=sbbHeader.numBlock-1; i>=0; i--)
    {
		cout << "Block: " << i+1 << endl;
        sb = ((CStructBookExp*)structBook[iSignal])[i].getStructBook();
        sbNumElement = ((CStructBookExp*)structBook[iSignal])[i].getNumElement();

        grpBlock = (int) (floor(static_cast<double>(i)/static_cast<double>(numBlockPerGroup)) );

        for (j=0;j<sbNumElement;j++)
        {
//			cout << "Element: " << j+1 << endl;
            prevAtom = sb[j].prevAtom;
            nextAtom = sb[j].nextAtom;

            if (nextAtom==-1)
            {
            	if (prevAtom==-1)
            	{

            		sbgrp.innerProduct = sb[j].innerProduct;
            		sbgrp.rho = sb[j].rho;
            		sbgrp.xi = sb[j].xi;
            		/////
            		// Tuning the phase reference
            		sbgrp.phase =  	sb[j].phase -
            						(sb[j].xi *
            						(i % numBlockPerGroup) *
            						sbbHeader.blockSize);
    				sbgrp.phase = sbgrp.phase - ceil(sbgrp.phase/(2*pi))*(2*pi);
    				if (sbgrp.phase < 0) sbgrp.phase = sbgrp.phase + 2*pi;
    				///////
					sbgrp.a = (i % numBlockPerGroup) * sbbHeader.blockSize +  sb[j].a;
					sbgrp.b = (i % numBlockPerGroup) * sbbHeader.blockSize +  sb[j].b;

					sbgrp.prevAtom = 0;
					sbgrp.nextAtom = 0;
					sbgrp.origAtomIndex = 0;

					//cout << "Block: "<< i << "; Element: " << j << endl;
// 					cout << "grpBlock: " << grpBlock <<  endl;
// 					cout << sbgrp.innerProduct << " " << sbgrp.rho << " " << sbgrp.xi << " "
// 						<< sbgrp.phase << " " << sbgrp.a << " " << sbgrp.b << " " << endl;
					((CStructBookExp*)structBookGroup[iSignal])[grpBlock].addElement(sbgrp);
            	}
            	else
            	{
            		// cout << "Block: "<< i << "; Element: " << j << endl;
// 					cout << "grpBlock: " << grpBlock <<  endl;

					sbgrp.innerProduct = sb[j].innerProduct;
            		sbgrp.rho = sb[j].rho;
            		sbgrp.xi = sb[j].xi;
            		/////
            		sbgrp.phase =  	sb[j].phase -
            						(sb[j].xi *
            						(i % numBlockPerGroup) *
            						sbbHeader.blockSize);
    				sbgrp.phase = sbgrp.phase - ceil(sbgrp.phase/(2*pi))*(2*pi);
    				if (sbgrp.phase < 0) sbgrp.phase = sbgrp.phase + 2*pi;
    				///////
					sbgrp.a = (i % numBlockPerGroup) * sbbHeader.blockSize +  sb[j].a;
					sbgrp.b = (i % numBlockPerGroup) * sbbHeader.blockSize +  sb[j].b;

// 					cout << "(i % numBlockPerGroup): " << (i % numBlockPerGroup) <<
//  						"; sbgrp.a: " << sbgrp.a << "; sbgrp.b: " << sbgrp.b << endl;
					// ----
					expDic.setRealAtomOnSupport(sb[j]);
					rAtomBlock = expDic.getRealAtom();
					for (k=sb[j].a; k<=sb[j].b;k++)
					{
						ind = (i % numBlockPerGroup) * sbbHeader.blockSize + k;
						cgSignal[0][ind] = sb[j].innerProduct*rAtomBlock[k];
					}

					prevRho = sbgrp.rho;

					ishift = 1;
					while ( (prevAtom!=-1) && ((i-ishift)>= grpBlock*numBlockPerGroup ) )
					{
						sb_aux = ((CStructBookExp*)structBook[iSignal])[i-ishift].getStructBook();


						if ( fabs(sb_aux[prevAtom].rho-prevRho) > DECAYVARIATIONFACTOR *fabs(prevRho)) break;

						sb_aux[prevAtom].nextAtom = i-ishift+1;

//						cout << "sb_aux[prevAtom].nextAtom: " << sb_aux[prevAtom].nextAtom << endl;


						prevRho = sb_aux[prevAtom].rho;
						////////
						sbgrp.phase =  	sb_aux[prevAtom].phase -
            						(sb_aux[prevAtom].xi *
            						((i-ishift) % numBlockPerGroup) *
            						sbbHeader.blockSize);
						sbgrp.phase = sbgrp.phase - ceil(sbgrp.phase/(2*pi))*(2*pi);
    					if (sbgrp.phase < 0) sbgrp.phase = sbgrp.phase + 2*pi;
						///////
						sbgrp.a = ((i-ishift) % numBlockPerGroup) * sbbHeader.blockSize +  sb_aux[prevAtom].a;
// 						cout << "((i-ishift) % numBlockPerGroup): " << ((i-ishift) % numBlockPerGroup) <<
//  						"; sbgrp.a: " << sbgrp.a << "; sbgrp.b: " << sbgrp.b << endl;
						///////
						// ----
						expDic.setRealAtomOnSupport(sb_aux[prevAtom]);
						rAtomBlock = expDic.getRealAtom();
						for (k=sb_aux[prevAtom].a; k<=sb_aux[prevAtom].b;k++)
						{
							ind = ((i-ishift) % numBlockPerGroup) * sbbHeader.blockSize + k;
							cgSignal[0][ind] = sb_aux[prevAtom].innerProduct*rAtomBlock[k];
						}

						ishift++;
						prevAtom = sb_aux[prevAtom].prevAtom;
					}

					if (ishift > 1)
					{
		    			expDicGrp.setRealAtomOnSupport(sbgrp);
						rAtomBlock = expDicGrp.getRealAtom();
						//innerProd = cgSignal*rAtomBlock;
						innerProd = cgSignal.fast_dprod(rAtomBlock,sbgrp.a,sbgrp.b);
						sbgrp.innerProduct = innerProd[0][0];
					}


// 					cout << sbgrp.innerProduct << " " << sbgrp.rho << " " << sbgrp.xi << " "
//  						<< sbgrp.phase << " " << sbgrp.a << " " << sbgrp.b << " " << endl;
					sbgrp.prevAtom = 0;
					sbgrp.nextAtom = 0;
					sbgrp.origAtomIndex = 0;
					((CStructBookExp*)structBookGroup[iSignal])[grpBlock].addElement(sbgrp);

            	}
         	}
        }
    }
}


void encodeSBExp(	char* InputFile,
					double rateTarget,
					strtSBBHeader sbbHeader,
                	CStructBook** structBook,
                	int flFixedLambdaVSUnifBlockRate)
{
	/////////////////////////////////
	// Load opcurve.bin file
	strtOpCurveExp* opCurve;
	int numElement, totalNumElement;
	double dummydouble;
	FILE* io_opcurve;
	io_opcurve = fopen("opcurve.bin","rb");
	if (io_opcurve ==NULL) cout << "Could not open opcurve.bin file" << endl;

	int i,j,k;
	totalNumElement=0;
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		for (j=0; j<sbbHeader.numBlock; j++)
		{
			///
			fread ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			((CStructBookExp*)structBook[i])[j].setMinAmp(dummydouble);
			///
			fread ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			((CStructBookExp*)structBook[i])[j].setMaxAmp(dummydouble);
			///
			fread ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			((CStructBookExp*)structBook[i])[j].setMinRho(dummydouble);
			///
			fread ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			((CStructBookExp*)structBook[i])[j].setMaxRho(dummydouble);
			///
			fread ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			((CStructBookExp*)structBook[i])[j].setMinPhase(dummydouble);
			///
			fread ( &dummydouble, sizeof ( double ), 1, io_opcurve );
			((CStructBookExp*)structBook[i])[j].setMaxPhase(dummydouble);
			///
			fread ( &numElement, sizeof ( int ), 1, io_opcurve );
			opCurve = new strtOpCurveExp[numElement];
			fread (	opCurve,sizeof ( strtOpCurveExp ), numElement, io_opcurve);
			((CStructBookExp*)structBook[i])[j].setOpCurve(opCurve,numElement);
			delete [] opCurve;
			totalNumElement += numElement;
		}
	}
	fclose(io_opcurve);

	double min_amp, min_rho, min_phase;
	double max_amp, max_rho, max_phase;

	double *lambda_vec;
	lambda_vec = new double[totalNumElement];
	int t = 0;
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		for (j=0; j<sbbHeader.numBlock; j++)
		{
			cout << "--------------------------------------------------------------" << endl;
			cout << setw (20) << setfill(' ') << "Signal"<< " :"
						<< setw (10) << setfill(' ') << i+1 << endl
						<< setw (20) << setfill(' ') << "Block"<< " :"
						<< setw (10) << setfill(' ') << j+1 << endl;
			cout << "--------------------------------------------------------------" << endl;
			cout << "- Ranges" << endl;

			min_amp =  ((CStructBookExp*)structBook[i])[j].getMinAmp();
			max_amp =  ((CStructBookExp*)structBook[i])[j].getMaxAmp();
			min_rho =  ((CStructBookExp*)structBook[i])[j].getMinRho();
			max_rho =  ((CStructBookExp*)structBook[i])[j].getMaxRho();
			min_phase =  ((CStructBookExp*)structBook[i])[j].getMinPhase();
			max_phase =  ((CStructBookExp*)structBook[i])[j].getMaxPhase();

			cout << setw (20) << setfill(' ') << "min_amp:"<< " "
				 << setw (20) << setfill(' ') << min_amp << endl
				 << setw (20) << setfill(' ') << "max_amp:"<< " "
				 << setw (20) << setfill(' ') << max_amp << endl
				 << setw (20) << setfill(' ') << "min_rho:"<< " "
				 << setw (20) << setfill(' ') << min_rho << endl
				 << setw (20) << setfill(' ') << "max_rho:"<< " "
				 << setw (20) << setfill(' ') << max_rho << endl
				 << setw (20) << setfill(' ') << "min_phase:"<< " "
				 << setw (20) << setfill(' ') << min_phase << endl
				 << setw (20) << setfill(' ') << "max_phase:"<< " "
				 << setw (20) << setfill(' ') << max_phase << endl;
		    cout << "--------------------------------------------------------------" << endl;
		    cout << "- Operational curve" << endl;
		    cout << setw (20) << setfill(' ') << "Rate"<< " "
				 << setw (20) << setfill(' ') << "Dist"<< " "
				 << setw (20) << setfill(' ') << "Lambda"<< endl;
			opCurve = ((CStructBookExp*)structBook[i])[j].getOpCurve();
			numElement = ((CStructBookExp*)structBook[i])[j].getNumElemOpCurve();
			for (k=0; k<numElement; k++)
			{
				cout << setw (20) << setfill(' ') << opCurve[k].rate   << " "
				     << setw (20) << setfill(' ') << opCurve[k].dist   << " "
				     << setw (20) << setfill(' ') << opCurve[k].lambda << endl;
				lambda_vec[t] = opCurve[k].lambda;
				t++;
			}
		}
	}

	bubble_srtdouble( lambda_vec, totalNumElement );

	double *rate_vec;
	rate_vec = new double[totalNumElement];

	int chosenOpCurvePointMat[totalNumElement][sbbHeader.numSignal][sbbHeader.numBlock];


	fstream file_rdresult("rdresult.out",ios::out);
	file_rdresult << setw (25) << setfill(' ') << "Lambda"<<" "
	              << setw (25) << setfill(' ') << "Bitrate (b/s)"<<" "
	              << setw (25) << setfill(' ') << "Distortion (Sqr. Error)"<< " " << endl;



	for (t=0; t<totalNumElement; t++)
	{
		double sqrerror =0.0;
		rate_vec[t] =0.0;
		for (i=0; i<sbbHeader.numSignal; i++)
		{
			for (j=0; j<sbbHeader.numBlock; j++)
			{
				opCurve = ((CStructBookExp*)structBook[i])[j].getOpCurve();
				numElement = ((CStructBookExp*)structBook[i])[j].getNumElemOpCurve();
				for (k=0; k<numElement; k++)
				{
					if (lambda_vec[t] <= opCurve[k].lambda)
					{
						sqrerror += opCurve[k].dist;
						rate_vec[t] += opCurve[k].rate;
						chosenOpCurvePointMat[t][i][j] = k;
						break;
					}
				}
			}
		}
		file_rdresult << setw (25) << setfill(' ') << setprecision(10) << lambda_vec[t] <<" "
	              << setw (30) << setfill(' ') << setprecision(10) << rate_vec[t]/
	              													  (sbbHeader.numSignal*sbbHeader.signalSize) <<" "
	              << setw (25) << setfill(' ') << setprecision(10) << sqrerror << " " << endl;
	}


	file_rdresult.close();

	int chosenLambdaInd;
	double minRateDiff = 1e8;
	for (t=0; t<totalNumElement; t++)
	{
		if (fabs(rate_vec[t] - (rateTarget*sbbHeader.signalSize)) < minRateDiff)
		{
			chosenLambdaInd = t;
			minRateDiff = fabs(rate_vec[t] - (rateTarget*sbbHeader.signalSize));
		}
	}
	cout << "--- Ordered lambda and rate---- " << endl;
	for (i=0; i<totalNumElement; i++)
	{
		cout << setw (20) << setfill(' ') << lambda_vec[i] << " "
			 << setw (20) << setfill(' ') << rate_vec[i] << endl;
	}
	cout << "--------------------------------" << endl;

	cout << "flFixedLambdaVSUnifBlockRate: " << flFixedLambdaVSUnifBlockRate << endl;
	cout << "--------------------------------" << endl;
	///////////////

	strtOpCurveExp chosenOpCurvePoint;
	strtContinuousExp* pSB;
	strtContinuousExp* pSBQ;
	double** recSignal;
	recSignal = new double*[sbbHeader.numSignal];

	double* recBlockSignal;
    recBlockSignal = new double[sbbHeader.blockSize];


    double nbits_signal_total=0;
    double sqrerror_total=0.0;
    double distAtom =0.0;

	for (i=0; i<sbbHeader.numSignal; i++)
	{
		recSignal[i] = new double[sbbHeader.signalSize];
		double rateTargetBlock = rateTarget/static_cast<double>(sbbHeader.numBlock);
		for (j=0; j<sbbHeader.numBlock; j++)
		{
			opCurve = ((CStructBookExp*)structBook[i])[j].getOpCurve();
			numElement = ((CStructBookExp*)structBook[i])[j].getNumElemOpCurve();

			cout << "Signal: " << i+1 << "; Block: " << j+1 << endl;
			///////////////////////////////////////////////
			// Choose RD optimum quantizer for each signal
			printf(" *** Choose RD optimum quantizer\n");
			minRateDiff=1e8;
			int chosenInd=0;

			if (flFixedLambdaVSUnifBlockRate == 0)
			{
				////////////////////////////////////////////////////////////////
				// Chose regarding bitrate distributed uniformily among blocks
				for (k=0;k<numElement;k++)
				{
					//cout << "- rateTarget*sbbHeader.blockSize: " << rateTarget*sbbHeader.blockSize
					//     << "; opCurve rate: "<< opCurve[k].rate
					//	 << "; minRateDiff: " << minRateDiff << endl;
					if (fabs(opCurve[k].rate - (rateTarget*sbbHeader.blockSize)) < minRateDiff)
					{
						chosenOpCurvePoint = opCurve[k];
						minRateDiff = fabs(opCurve[k].rate - (rateTarget*sbbHeader.blockSize));
					}
				}
			}
			else if (flFixedLambdaVSUnifBlockRate == 1)
			{
				///////////////////////////////////////////////////////////////
				// Chose regarding fixed lambda
				chosenOpCurvePoint = opCurve[ chosenOpCurvePointMat[chosenLambdaInd][i][j] ];
			}

			cout << "- nb_amp: " << chosenOpCurvePoint.nb_amp << endl
			     << "- nb_rho: " << chosenOpCurvePoint.nb_rho << endl
			     << "- nb_phase: " << chosenOpCurvePoint.nb_phase << endl
			     << "- nb_xi: "  << chosenOpCurvePoint.nb_xi << endl
			     << "- nb_sample: " << chosenOpCurvePoint.nb_sample<< endl;
			     //<< "- nb_block: " << chosenOpCurvePoint.nb_block<< endl;

			//////
			min_amp =  ((CStructBookExp*)structBook[i])[j].getMinAmp();
			max_amp =  ((CStructBookExp*)structBook[i])[j].getMaxAmp();
			min_rho =  ((CStructBookExp*)structBook[i])[j].getMinRho();
			max_rho =  ((CStructBookExp*)structBook[i])[j].getMaxRho();
			min_phase =  ((CStructBookExp*)structBook[i])[j].getMinPhase();
			max_phase =  ((CStructBookExp*)structBook[i])[j].getMaxPhase();
			///////
			printf(" ****** Quantize Structure Book \n");
			strtQuantExp chosenQuantExp;
			chosenQuantExp.nb_amp = chosenOpCurvePoint.nb_amp;
			chosenQuantExp.nb_rho = chosenOpCurvePoint.nb_rho;
			chosenQuantExp.nb_phase = chosenOpCurvePoint.nb_phase;
			chosenQuantExp.nb_xi = chosenOpCurvePoint.nb_xi;
			chosenQuantExp.nb_sample = chosenOpCurvePoint.nb_sample;
			chosenQuantExp.nb_block = chosenOpCurvePoint.nb_block;

			nbits_signal_total+= chosenOpCurvePoint.rate;
    		sqrerror_total+= chosenOpCurvePoint.dist;
			/////
			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
			int sbNumElementQ;

			pSBQ = new strtContinuousExp[sbNumElement];

			distAtom = quantizeStructBookExpDistAtom( min_amp,max_amp,
													  min_rho,max_rho,
													  min_phase,max_phase,
													  chosenQuantExp,
													  pSB,sbNumElement,
													  pSBQ,sbNumElementQ, // output elements
									  				  sbbHeader.blockSize, sbbHeader.norm[i]);
			quantizeStructBookExp(min_amp,max_amp,
								  min_rho,max_rho,
								  min_phase,max_phase,
								  chosenQuantExp,
								  pSB,sbNumElement,
								  pSBQ,sbNumElementQ); // output elements
			cout << "sbNumElementQ: " << sbNumElementQ << endl;



			printf(" ****** Synthesize block signal \n");
			synthSignalSBExp(	recBlockSignal,
								sbbHeader.norm[i],
								sbbHeader.blockSize,
								pSBQ,
								sbNumElementQ);
			int initBlockSample = j*sbbHeader.blockSize;
			if ((initBlockSample + sbbHeader.blockSize) >= sbbHeader.signalSize)
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
			}
			else
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*sbbHeader.blockSize);
			}

			delete [] pSBQ;
		}
	}

	double bitrate = static_cast<double>(nbits_signal_total)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Final rate: %f \n",nbits_signal_total);
    printf("Signal size: %d \n",sbbHeader.signalSize);
    printf("Bitrate: %f bit/sample\n",bitrate);
    //printf("Total MSE: %f \n",mse_total);
    //printf("Total MSE per Atom: %f \n",mse_total_per_atom);
    printf("Total Square Error: %f \n",sqrerror_total);
    //printf("Total Square Error per Atom: %f \n",mse_total_per_atom*
    //                                            static_cast<double>(sbbHeader.signalSize));
    printf("Total Atom Square Error: %f \n",distAtom);

    char fileName[_MAX_PATH];
    strcpy(fileName, InputFile);
    char* pos;
    pos = strrchr( fileName, '.');
    char aux[_MAX_PATH];
    sprintf(aux,"_rec%09.6fbps.wav",bitrate);
    strcpy( &pos[0], aux);

    CDataSignal* outSignal;
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSamplingRate(sbbHeader.Fs);
    outSignal->setSignal(recSignal);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;


	/////////
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		delete [] recSignal[i];
	}
	delete [] recSignal;
	delete [] recBlockSignal;
	///////
	delete [] lambda_vec;
	delete [] rate_vec;

}



void encodeSBExpAmpRange(	char*InputFile,
							double rateTarget,
							strtSBBHeader sbbHeader,
							CStructBook** structBook,
							CFileDfTable* DfTableFile,
							int Nfreq,
							CDataSignal* dataSignal)
{

    // Setting fixed quantizers
    int nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
    int nb_sample = (int)ceil( log(sbbHeader.blockSize) / log(2.0) );
    int nb_rhoSign = 1;

	////////////////
	int i, j, k, i_qrho, i_qphi;
	int iAmpRange;
	int t;
	////////////////

    double** min_amp= new double*[sbbHeader.numSignal];
    double** max_amp= new double*[sbbHeader.numSignal];
    double** min_rho= new double*[sbbHeader.numSignal];
    double** max_rho= new double*[sbbHeader.numSignal];
    double** min_phase= new double*[sbbHeader.numSignal];
    double** max_phase= new double*[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		min_amp[i]= new double[sbbHeader.numBlock];
    	max_amp[i]= new double[sbbHeader.numBlock];
    	min_rho[i]= new double[sbbHeader.numBlock];
    	max_rho[i]= new double[sbbHeader.numBlock];
    	min_phase[i]= new double[sbbHeader.numBlock];
    	max_phase[i]= new double[sbbHeader.numBlock];
	}

	///////////======
	strtContinuousExp* pSB;
    int sbNumElement;

    double**** DfTable = DfTableFile->getDfTable();
    int N_qrho = DfTableFile->getNumQRho();
    int N_qphi = DfTableFile->getNumQPhi();
    double* q_rho = DfTableFile->getQRho();
    double* q_phi = DfTableFile->getQPhi();

    int nb_amp;
    double nb_rho, nb_phi;

    cout << "================================================" << endl;
    cout << "- Computing rate for different lambda/R_amp" << endl;

    int lambda_vec_length = sbbHeader.numSignal*sbbHeader.numBlock*16;
    double* lambda_vec = new double[lambda_vec_length];
    t=0;
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			cout << i <<" " << j << endl;
			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp[i][j],
																  max_amp[i][j],
																  min_rho[i][j],
																  max_rho[i][j],
																  min_phase[i][j],
																  max_phase[i][j]);

			double L_amp = max_amp[i][j] - min_amp[i][j];
    		for (nb_amp=1;nb_amp<=16;nb_amp++)
    		{
    			// Save the available lambdas
				lambda_vec[t] = (L_amp*L_amp/12) * 2 * log(2) * pow(2.0,-2.0*static_cast<double>(nb_amp));
    			t++;
    		}
		}
	}

	int numNoDuplicateLambda;
	bubble_srtdouble_noduplicate(lambda_vec, lambda_vec_length, numNoDuplicateLambda);

	strtRDByAmpQuant* rdByAmpQuant = new strtRDByAmpQuant[numNoDuplicateLambda];
	for (t=0;t<numNoDuplicateLambda;t++)
    {
    	rdByAmpQuant[t].rate = new double**[sbbHeader.numSignal];
    	rdByAmpQuant[t].quantExp = new strtQuantExp**[sbbHeader.numSignal];
    	rdByAmpQuant[t].nAtom = new int**[sbbHeader.numSignal];
    	rdByAmpQuant[t].numAmpRange = new int*[sbbHeader.numSignal];
		for (i=0;i<sbbHeader.numSignal;i++)
		{
			rdByAmpQuant[t].rate[i] = new double*[sbbHeader.numBlock];
			rdByAmpQuant[t].quantExp[i] = new strtQuantExp*[sbbHeader.numBlock];
			rdByAmpQuant[t].numAmpRange[i] = new int[sbbHeader.numBlock];
			rdByAmpQuant[t].nAtom[i] = new int*[sbbHeader.numBlock];
		}
	}

    for (t=0;t<numNoDuplicateLambda;t++)
    {
    	rdByAmpQuant[t].totalRate = 0;
    	double totalRate = 0.0;
    	for (i=0;i<sbbHeader.numSignal;i++)
		{
			for (j=0;j<sbbHeader.numBlock;j++)
			{
				cout << "Lambda: " << lambda_vec[t] <<  " -- Signal: " << i+1 << " -- " << "Block: " << j+1 << endl;

				pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
				sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();


				strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];

				int sbNumElementQ;

				double L_amp = max_amp[i][j] - min_amp[i][j];
				double L_rho = max_rho[i][j] - min_rho[i][j];
				double L_phi = max_phase[i][j] - min_phase[i][j];

				nb_amp = static_cast<int>( round( 0.5* log2(L_amp*L_amp*log(2)/(6*lambda_vec[t])) ) );
				if (nb_amp<1) nb_amp=1;
				if (nb_amp>16) nb_amp=16;

				quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
											pSB, sbNumElement,
											pSBQ,sbNumElementQ);

    			int NAmpRange = nb_amp;

    			rdByAmpQuant[t].numAmpRange[i][j] = NAmpRange;
    			rdByAmpQuant[t].rate[i][j] = new double[NAmpRange];
				rdByAmpQuant[t].quantExp[i][j] = new strtQuantExp[NAmpRange];
				rdByAmpQuant[t].nAtom[i][j] = new int[NAmpRange];

				CStructBook* structBookByAmp;
				structBookByAmp = new CStructBookExp[NAmpRange];

				double* ampRangeLimit;
				ampRangeLimit = new double[NAmpRange+1];

				double* ampBar;
				ampBar = new double[NAmpRange];

				double* ampRangeNumElement;
				ampRangeNumElement = new double[NAmpRange];

				ampRangeLimit[0]=-0.001;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				    ampBar[NAmpRange-1-iAmpRange] =
				    	(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				    ((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
																			sbNumElementQ,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

					ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();
				}

				int rate=0;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					// Find nb_rho and nb_phi leading minimum Df Lagragian
					double minLagrangianDf = 1e30;
					double LagrangianDf;

					for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
					{
						if (log2( 1 + (L_rho/q_rho[i_qrho]) ) < 1) break;
						for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
						{
							if (log2( 1 + (L_phi/q_phi[i_qphi]) ) < 1) break;
							LagrangianDf =  (ampBar[iAmpRange]*ampBar[iAmpRange]) *
												DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] +
											lambda_vec[t] * ( log2( 1 + (L_rho/q_rho[i_qrho]) ) + log2( 1 + (L_phi/q_phi[i_qphi]) ) );
							if (LagrangianDf < minLagrangianDf )
							{
								minLagrangianDf = LagrangianDf;
								nb_rho = log2( 1 + (L_rho/q_rho[i_qrho]) );
								nb_phi = log2( 1 + (L_phi/q_phi[i_qphi]) );
							}
						} // nb_phase
					}  // nb_rho

					strtQuantExp quantExp;
					quantExp.nb_amp = nb_amp;
					quantExp.nb_rho = static_cast<int>(round(nb_rho));
					quantExp.nb_phase = static_cast<int>(round(nb_phi));
					quantExp.nb_xi = nb_xi;
					quantExp.nb_sample = nb_sample;

					rate += computeRateSBExp(quantExp,ampRangeNumElement[iAmpRange]);


					rdByAmpQuant[t].rate[i][j][iAmpRange] = rate;
					rdByAmpQuant[t].quantExp[i][j][iAmpRange] = quantExp;
					rdByAmpQuant[t].nAtom[i][j][iAmpRange] = ampRangeNumElement[iAmpRange];

				} // iAmpRange
				totalRate += static_cast<double>(rate);

				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;
				delete [] ampRangeNumElement;
				delete [] pSBQ;
			} // block
		} // signal
		rdByAmpQuant[t].lambda = lambda_vec[t];
		rdByAmpQuant[t].totalRate = totalRate;
    } // lambda

    fstream file_rdbyamp("rdbyamp_quant.out",ios::out);

	file_rdbyamp << setw (20) << setfill(' ') << "IndLambda"<<" "
	             << setw (20) << setfill(' ') << "Lambda"<<" "
	             << setw (20) << setfill(' ') << "Total Rate"<<" "
	             << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "NAmpRange/nbamp"<<" "
	             << setw (20) << setfill(' ') << "NAtomsInRange"<<" "
	             << setw (20) << setfill(' ') << "iAmpRange"<<" "
	             << setw (20) << setfill(' ') << "Rate"<<" "
	             << setw (20) << setfill(' ') << "nb_rho"<<" "
	             << setw (20) << setfill(' ') << "nb_phi"<<" "
	             << setw (20) << setfill(' ') << "nb_xi"<<" "
	             << setw (20) << setfill(' ') << "nb_sample"<<" " << endl;



	for (t=0;t<numNoDuplicateLambda;t++)
    {
		for (i=0;i<sbbHeader.numSignal;i++)
		{
			for (j=0;j<sbbHeader.numBlock;j++)
			{
    			for (iAmpRange=0;iAmpRange < rdByAmpQuant[t].numAmpRange[i][j];iAmpRange++)
				{
					file_rdbyamp << setw (20) << setfill(' ') << t+1 << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].lambda << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].totalRate << " "
								 << setw (20) << setfill(' ') << i+1 << " "
								 << setw (20) << setfill(' ') << j+1 << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].numAmpRange[i][j] << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].nAtom[i][j][iAmpRange] << " "
								 << setw (20) << setfill(' ') << iAmpRange+1 << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].rate[i][j][iAmpRange] << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].quantExp[i][j][iAmpRange].nb_rho << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].quantExp[i][j][iAmpRange].nb_phase << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].quantExp[i][j][iAmpRange].nb_xi << " "
								 << setw (20) << setfill(' ') << rdByAmpQuant[t].quantExp[i][j][iAmpRange].nb_sample << " " << endl;
				}
    		}
		}
	}
	file_rdbyamp.close();

	// Obtain the chosen lamda given a desired rate
	double minRateDiff = 1e30;
	int chosenLambdaInd;
	for (t=0;t<numNoDuplicateLambda;t++)
    {
    	if (fabs(rdByAmpQuant[t].totalRate - (rateTarget*sbbHeader.signalSize)) < minRateDiff)
		{
			chosenLambdaInd = t;
			minRateDiff = fabs( rdByAmpQuant[t].totalRate - (rateTarget*sbbHeader.signalSize));
		}
    }


    t = chosenLambdaInd;

    double** origSignal;
    origSignal = dataSignal->getSignal();

    double** recSignal;
	recSignal = new double*[sbbHeader.numSignal];

	double* recBlockSignal;
    recBlockSignal = new double[sbbHeader.blockSize];

    double finalRate=0;
    double sqrerror_total=0.0;

    for (i=0;i<sbbHeader.numSignal;i++)
	{
		recSignal[i] = new double[sbbHeader.signalSize];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

			int NAmpRange = rdByAmpQuant[t].numAmpRange[i][j];

			CStructBook* structBookByAmp;
			structBookByAmp = new CStructBookExp[NAmpRange];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];

			double* ampRangeNumElement;
			ampRangeNumElement = new double[NAmpRange];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSB,
																		sbNumElement,
																		ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																		ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit
			}


			// Initialize block reconstruction vector (RESET)
			for (k=0;k <sbbHeader.blockSize ;k++)
			{
				recBlockSignal[k] = 0.0;
			}

			for (iAmpRange=0;iAmpRange < rdByAmpQuant[t].numAmpRange[i][j];iAmpRange++)
			{
				strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
				int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElementAmp];
				int numElementAmpQ;

				strtQuantExp chosenQuantExp = rdByAmpQuant[t].quantExp[i][j][iAmpRange];

				quantizeStructBookExp(min_amp[i][j],max_amp[i][j],
									  min_rho[i][j],max_rho[i][j],
									  min_phase[i][j],max_phase[i][j],
									  chosenQuantExp,
									  pSBAmp,numElementAmp,
									  pSBAmpQ,numElementAmpQ); // output elements

				printf(" ****** Synthesize block signal \n");
				synthSignalSBExpNoReset(recBlockSignal,
										sbbHeader.norm[i],
										sbbHeader.blockSize,
										pSBAmpQ,
										numElementAmpQ);

				finalRate += computeRateSBExp(chosenQuantExp,numElementAmpQ);



				delete [] pSBAmpQ;
			}

			int initBlockSample = j*sbbHeader.blockSize;
			if ((initBlockSample + sbbHeader.blockSize) >= sbbHeader.signalSize)
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
			}
			else
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*sbbHeader.blockSize);
			}

			delete [] (CStructBookExp*)structBookByAmp;
			delete [] ampBar;
			delete [] ampRangeLimit;
			delete [] ampRangeNumElement;
		} // block
		sqrerror_total += computeSqrError(origSignal[i],recSignal[i],0,sbbHeader.signalSize-1);
	} //signal

	double bitrate = static_cast<double>(finalRate)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Bitrate: %f bit/sample\n",bitrate);
    printf("Total Square Error: %f \n",sqrerror_total);

    char fileName[_MAX_PATH];
    strcpy(fileName, InputFile);
    char* pos;
    pos = strrchr( fileName, '.');
    char aux[_MAX_PATH];
    sprintf(aux,"_rec%09.6fbps.wav",bitrate);
    strcpy( &pos[0], aux);

    CDataSignal* outSignal;
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSignal(recSignal);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;


	for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] min_amp[i];
    	delete [] max_amp[i];
    	delete [] min_rho[i];
    	delete [] max_rho[i];
    	delete [] min_phase[i];
    	delete [] max_phase[i];
	}
	delete [] min_amp;
	delete [] max_amp;
	delete [] min_rho;
	delete [] max_rho;
	delete [] min_phase;
	delete [] max_phase;



    for (t=0;t<numNoDuplicateLambda;t++)
    {
		for (i=0;i<sbbHeader.numSignal;i++)
		{
			for (j=0;j<sbbHeader.numBlock;j++)
			{
				delete [] rdByAmpQuant[t].rate[i][j];
				delete [] rdByAmpQuant[t].quantExp[i][j];
				delete [] rdByAmpQuant[t].nAtom[i][j];
			}
			delete [] rdByAmpQuant[t].rate[i];
			delete [] rdByAmpQuant[t].quantExp[i];
			delete [] rdByAmpQuant[t].numAmpRange[i];
			delete [] rdByAmpQuant[t].nAtom[i];
		}
		delete [] rdByAmpQuant[t].rate;
    	delete [] rdByAmpQuant[t].quantExp;
    	delete [] rdByAmpQuant[t].numAmpRange;
    	delete [] rdByAmpQuant[t].nAtom;
	}
	delete [] rdByAmpQuant;
	/////////
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		delete [] recSignal[i];
	}
	delete [] recSignal;
	delete [] recBlockSignal;
}


// 	strtRDByAmpQuant*** rdByAmpQuant = new strtRDByAmpQuant**[sbbHeader.numSignal];
//     for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		rdByAmpQuant[i] = new strtRDByAmpQuant*[sbbHeader.numBlock];
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			rdByAmpQuant[i][j] = new strtRDByAmpQuant[16];
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
//     		{
//     			rdByAmpQuant[i][j][nb_amp-1].quantExp = new strtQuantExp[nb_amp];
//     			rdByAmpQuant[i][j][nb_amp-1].rate = new double[nb_amp];
//     		}
// 		}
// 	}
//     for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			cout << "Signal: " << i+1 << " -- " << "Block: " << j+1 << endl;
//
// 			pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
// 			sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
//
// 			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];
//
// 			int sbNumElementQ;
//
// 			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp[i][j],
// 															  max_amp[i][j],
// 															  min_rho[i][j],
// 															  max_rho[i][j],
// 															  min_phase[i][j],
// 															  max_phase[i][j]);
//
// 			double L_amp = max_amp[i][j] - min_amp[i][j];
// 			double L_rho = max_rho[i][j] - min_rho[i][j];
// 			double L_phi = max_phase[i][j] - min_phase[i][j];
//
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
//     		{
//     			double lambda = (L_amp*L_amp/12) * 2 * log(2) * pow(2.0,-2.0*static_cast<double>(nb_amp));
//
//
// 				quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
// 											pSB, sbNumElement,
// 											pSBQ,sbNumElementQ);
//
//     			int NAmpRange = nb_amp;
//
// 				CStructBook* structBookByAmp;
// 				structBookByAmp = new CStructBookExp[NAmpRange];
//
// 				double* ampRangeLimit;
// 				ampRangeLimit = new double[NAmpRange+1];
//
// 				double* ampBar;
// 				ampBar = new double[NAmpRange];
//
// 				double* ampRangeNumElement;
// 				ampRangeNumElement = new double[NAmpRange];
//
// 				ampRangeLimit[0]=-0.001;
// 				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
// 				{
// 					ampRangeLimit[iAmpRange+1] =
// 						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
// 				    ampBar[NAmpRange-1-iAmpRange] =
// 				    	(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));
//
// 				    ((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
// 																			sbNumElementQ,
// 																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
// 																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit
//
// 					ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();
// 				}
//
// 				double totalRate = 0;
// 				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
// 				{
// 					// Find nb_rho and nb_phi leading minimum Df Lagragian
// 					double minLagrangianDf = 1e30;
// 					double LagrangianDf;
// 					for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
// 					{
// 						if (log2( 1 + (L_rho/q_rho[i_qrho]) ) < 1) break;
// 						for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
// 						{
// 							if (log2( 1 + (L_phi/q_phi[i_qphi]) ) < 1) break;
// 							LagrangianDf =  (ampBar[iAmpRange]*ampBar[iAmpRange]) *
// 												DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] +
// 											lambda * ( log2( 1 + (L_rho/q_rho[i_qrho]) ) + log2( 1 + (L_phi/q_phi[i_qphi]) ) );
// 							if (LagrangianDf < minLagrangianDf )
// 							{
// 								minLagrangianDf = LagrangianDf;
// 								nb_rho = log2( 1 + (L_rho/q_rho[i_qrho]) );
// 								nb_phi = log2( 1 + (L_phi/q_phi[i_qphi]) );
// 							}
// 						} // nb_phase
// 					}  // nb_rho
//
// 					strtQuantExp quantExp;
// 					quantExp.nb_amp = nb_amp;
// 					quantExp.nb_rho = static_cast<int>(round(nb_rho));
// 					quantExp.nb_phase = static_cast<int>(round(nb_phi));
// 					quantExp.nb_xi = nb_xi;
// 					quantExp.nb_sample = nb_sample;
// 					int rate = computeRateSBExp(quantExp,ampRangeNumElement[iAmpRange]);
// 					totalRate += rate;
//
// 					rdByAmpQuant[i][j][nb_amp-1].rate[iAmpRange] = rate;
// 					rdByAmpQuant[i][j][nb_amp-1].quantExp[iAmpRange] = quantExp;
//
// 				} // iAmpRange
//
// 				rdByAmpQuant[i][j][nb_amp-1].totalRate = totalRate;
// 				rdByAmpQuant[i][j][nb_amp-1].lambda = lambda;
//
// 				delete [] (CStructBookExp*)structBookByAmp;
// 				delete [] ampBar;
// 				delete [] ampRangeLimit;
// 				delete [] ampRangeNumElement;
// 			} // nb_amp
// 		} // block
// 	} // signal


//     fstream file_rdbyamp("rdbyampquant.csv",ios::out);
//
// 	file_rdbyamp << setw (20) << setfill(' ') << "Signal"<<";"
// 	             << setw (20) << setfill(' ') << "Block"<<";"
// 	             << setw (20) << setfill(' ') << "Total Rate"<<";"
// 	             << setw (20) << setfill(' ') << "Lambda"<<";"
// 	             << setw (20) << setfill(' ') << "Rate"<<";"
// 	             << setw (20) << setfill(' ') << "nb_amp"<<";"
// 	             << setw (20) << setfill(' ') << "iAmpRange"<<";"
// 	             << setw (20) << setfill(' ') << "nb_rho"<<";"
// 	             << setw (20) << setfill(' ') << "nb_phi"<<";"
// 	             << setw (20) << setfill(' ') << "nb_xi"<<";"
// 	             << setw (20) << setfill(' ') << "nb_sample"<<";" << endl;
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
//     		{
//     			for (iAmpRange=0;iAmpRange<nb_amp;iAmpRange++)
// 				{
// 					file_rdbyamp << setw (20) << setfill(' ') << i+1 << ";"
// 								 << setw (20) << setfill(' ') << j+1 << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].totalRate << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].lambda << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].rate[iAmpRange] << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].quantExp[iAmpRange].nb_amp << ";"
// 								 << setw (20) << setfill(' ') << iAmpRange+1 << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].quantExp[iAmpRange].nb_rho << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].quantExp[iAmpRange].nb_phase << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].quantExp[iAmpRange].nb_xi << ";"
// 								 << setw (20) << setfill(' ') << rdByAmpQuant[i][j][nb_amp-1].quantExp[iAmpRange].nb_sample << ";" << endl;
// 				}
//     		}
// 		}
// 	}
// 	file_rdbyamp.close();
//
// 	cout << "================================================" << endl;
//     cout << "- Optimal Bit Allocation for a given rate target" << endl;
//
//
//     for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
// 			{
// 				lambda_vec[t] =
// 				t++;
// 			}
// 		}
// 	}
//
//
    //
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
// 			{
// 				if (fabs(rate_vec[t] - (rateTarget*sbbHeader.signalSize)) < minRateDiff)
// 				{
//
// 					minRateDiff = fabs(rate_vec[t] - (rateTarget*sbbHeader.signalSize));
// 				}
// 				rdByAmpQuant[i][j][nb_amp-1].totalRate
// 			}
// 		}
// 	}
//     for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
//     		{
//     			delete [] rdByAmpQuant[i][j][nb_amp-1].quantExp;
//     			delete [] rdByAmpQuant[i][j][nb_amp-1].rate;
//     		}
// 			delete [] rdByAmpQuant[i][j];
// 		}
// 		delete [] rdByAmpQuant[i];
// 	}
// 	delete [] rdByAmpQuant;

// void encodeSBExpAmpRangeOpcurve( char* InputFile,
// 								double rateTarget,
// 								strtSBBHeader sbbHeader,
// 								CStructBook** structBook)
// {
//
// 	strtOpCurveExp**** opCurve;
//
// 	opCurve = new strtOpCurveExp***[sbbHeader.numSignal];
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		opCurve[i] =  new strtOpCurveExp**[sbbHeader.numBlock];
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			opCurve[i][j] = new strtOpCurveExp*[16];
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
//     		{
// 				opCurve[i][j][nbamp-1] = new strtOpCurveExp[nbamp];
// 			}
// 		}
// 	}
//
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			for (nb_amp=1;nb_amp<=16;nb_amp++)
//     		{
//     			delete [] opCurve[i][j][nbamp-1];
//
//     		}
// 			delete [] opCurve[i][j];
// 		}
// 		delete [] opCurve[i];
// 	}
//
//
// }

void encodeSBExpAmpRangeBisecSearch(char*InputFile,
									double rateTarget,
									strtSBBHeader sbbHeader,
									CStructBook** structBook,
									CFileDfTable* DfTableFile,
									int Nfreq,
									CDataSignal* dataSignal,
									int flLambdaBisecSearch,
									double initLambda)
{

    // Setting fixed quantizers
    int nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
    int nb_sample = (int)ceil( log(sbbHeader.blockSize) / log(2.0) );
    int nb_rhoSign = 1;

	////////////////
	int i, j, k, i_qrho, i_qphi;
	int iAmpRange;
	int t;
	////////////////

    double** min_amp= new double*[sbbHeader.numSignal];
    double** max_amp= new double*[sbbHeader.numSignal];
    double** min_rho= new double*[sbbHeader.numSignal];
    double** max_rho= new double*[sbbHeader.numSignal];
    double** min_phase= new double*[sbbHeader.numSignal];
    double** max_phase= new double*[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		min_amp[i]= new double[sbbHeader.numBlock];
    	max_amp[i]= new double[sbbHeader.numBlock];
    	min_rho[i]= new double[sbbHeader.numBlock];
    	max_rho[i]= new double[sbbHeader.numBlock];
    	min_phase[i]= new double[sbbHeader.numBlock];
    	max_phase[i]= new double[sbbHeader.numBlock];
	}

	///////////======
	strtContinuousExp* pSB;
    int sbNumElement;

    double**** DfTable = DfTableFile->getDfTable();
    int N_qrho = DfTableFile->getNumQRho();
    int N_qphi = DfTableFile->getNumQPhi();
    double* q_rho = DfTableFile->getQRho();
    double* q_phi = DfTableFile->getQPhi();

    int nb_amp;
    double nb_rho, nb_phi;

	// Find the range of the parameters for each Signal/Block
	double max_max_amp=0.0;
	double min_max_amp=1e30;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp[i][j],
																  max_amp[i][j],
																  min_rho[i][j],
																  max_rho[i][j],
																  min_phase[i][j],
																  max_phase[i][j]);
			if (max_max_amp<max_amp[i][j]) max_max_amp = max_amp[i][j];
			if (min_max_amp>max_amp[i][j]) min_max_amp = max_amp[i][j];
		}
	}

	// Compute the Amplitude Square Error for each Signal/Block/nb_amp
	double*** distAmp;
	distAmp = new double**[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		distAmp[i] = new double*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			distAmp[i][j] = new double[16];
			////
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
			////
			for (nb_amp=1;nb_amp<16;nb_amp++)
			{
				distAmp[i][j][nb_amp-1] = sumSqrErrorStructBookExpAmp(  min_amp[i][j], max_amp[i][j], nb_amp,
																		pSB, sbNumElement);
			}
		}
	}

	////////////////////
	cout << "- Computing rate for lower bound lambda" << endl;
	double leftLambda = ((min_max_amp*min_max_amp)/12.0)* (2*log(2)*pow(2.0,-2*30));
	double leftRate;
	double leftNbAmp;
	getRateNbAmpGivenLambda(    leftLambda,
								sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
							    DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								distAmp,
								leftRate, // output
								leftNbAmp); // output

	///////////////////////////
	cout << "- Computing rate for upper bound lambda" << endl;
	double rightLambda = ((max_max_amp*max_max_amp)/12.0)* (2*log(2)*pow(2.0,-2*1));
	double rightRate;
	double rightNbAmp;
	getRateNbAmpGivenLambda(    rightLambda,
								sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
							    DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								distAmp,
								rightRate, // output
								rightNbAmp); // output

	double epsilon = leftLambda*0.5;
	double chosenLambda;
	double chosenNbAmp, chosen_nb_amp_aux;
	double totalRate, chosenRate;

	cout << "- leftLambda" << " " << "/ leftRate" << endl;
	cout << leftLambda << " / " << leftRate << endl;
	cout << "- rightLambda" << " " << "/ rightRate" << endl;
	cout << rightLambda << " / " << rightRate << endl;

	double desiredLambda;
	if (initLambda==0.0)
	{
		desiredLambda = (rightLambda + leftLambda)/2.0;
	}
	else
	{
		cout << "- The initial lambda is given by USER: " << initLambda << endl;
		desiredLambda = initLambda;
	}

	chosenLambda = initLambda;

	int countTol = 0;
	while(1)
	{
		if (flLambdaBisecSearch!=1)
		{
			cout << "- The lambda bisectionsearch is not used." << endl;
			getRateNbAmpGivenLambda(chosenLambda,
								sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
							    DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								distAmp,
								totalRate, // output
								chosenNbAmp);// output);
			break;
		}
		//////////////////////////////////////
		getRateNbAmpGivenLambda(desiredLambda,
								sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
							    DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								distAmp,
								totalRate, // output
								chosen_nb_amp_aux); // output

		cout << "- leftLambda" << " " << "/ leftRate" << endl;
		cout << leftLambda << " / " << leftRate << endl;

		cout << "- rightLambda" << " " << "/ rightRate" << endl;
		cout << rightLambda << " / " << rightRate << endl;

		cout << "* desiredLambda: " << desiredLambda << endl;
		cout << "* totalRate: " << totalRate << endl;
		cout << "* TARGETRate: " << rateTarget*sbbHeader.signalSize << endl;

		if ( totalRate <= rateTarget*sbbHeader.signalSize)
		{
			rightLambda = desiredLambda;
			rightRate = totalRate;
			rightNbAmp = chosen_nb_amp_aux;
		}
		else
		{
			leftLambda = desiredLambda;
			leftRate = totalRate;
			leftNbAmp = chosen_nb_amp_aux;
		}
		//////////////////////////////////////////
		if (fabs(rightLambda - leftLambda) < epsilon)
		{
			if ( fabs(leftRate-rateTarget*sbbHeader.signalSize) <
				 fabs(rightRate-rateTarget*sbbHeader.signalSize) )
			{
				chosenLambda = leftLambda;
				chosenNbAmp = leftNbAmp;
				chosenRate = leftRate;
			}
			else
			{
				chosenLambda = rightLambda;
				chosenNbAmp = rightNbAmp;
				chosenRate = rightRate;
			}
			/////
			break;
		}

		desiredLambda = (rightLambda + leftLambda)/2.0;
	}

	strtRDByAmpQuant* rdByAmpQuant = new strtRDByAmpQuant;
	int NAmpRange = static_cast<int>(chosenNbAmp);
	//////
	rdByAmpQuant->rate = new double**[sbbHeader.numSignal];
	rdByAmpQuant->quantExp = new strtQuantExp**[sbbHeader.numSignal];
	rdByAmpQuant->nAtom = new int**[sbbHeader.numSignal];
	rdByAmpQuant->numAmpRange = new int*[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		rdByAmpQuant->rate[i] = new double*[sbbHeader.numBlock];
		rdByAmpQuant->quantExp[i] = new strtQuantExp*[sbbHeader.numBlock];
		rdByAmpQuant->numAmpRange[i] = new int[sbbHeader.numBlock];
		rdByAmpQuant->nAtom[i] = new int*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			rdByAmpQuant->numAmpRange[i][j] = NAmpRange;
			rdByAmpQuant->rate[i][j] = new double[NAmpRange];
			rdByAmpQuant->quantExp[i][j] = new strtQuantExp[NAmpRange];
			rdByAmpQuant->nAtom[i][j] = new int[NAmpRange];
		}
	}
	rdByAmpQuant->lambda = chosenLambda;
	computeRDByAmpQuant(	sbbHeader,structBook,
							min_amp, max_amp,
							min_rho, max_rho,
							min_phase, max_phase,
							DfTable,
							N_qrho, N_qphi,
							q_rho, q_phi,
							nb_xi, nb_sample,
							chosenNbAmp,
							rdByAmpQuant);
	totalRate = rdByAmpQuant->totalRate;

	cout << "#######################" << endl;
	cout << "# Chosen lambda: " << chosenLambda << endl;
	cout << "# total rate: " << totalRate << endl;
	cout << "#######################" << endl;

    ////////////////////////////////////////////////////////
    fstream file_rdbyamp("rdbyamp_quant.out",ios::out);
	file_rdbyamp << setw (20) << setfill(' ') << "Lambda"<<" "
	             << setw (20) << setfill(' ') << "Total Rate"<<" "
	             << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "NAmpRange/nbamp"<<" "
	             << setw (20) << setfill(' ') << "NAtomsInRange"<<" "
	             << setw (20) << setfill(' ') << "iAmpRange"<<" "
	             << setw (20) << setfill(' ') << "Rate"<<" "
	             << setw (20) << setfill(' ') << "nb_rho"<<" "
	             << setw (20) << setfill(' ') << "nb_phi"<<" "
	             << setw (20) << setfill(' ') << "nb_xi"<<" "
	             << setw (20) << setfill(' ') << "nb_sample"<<" " << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (iAmpRange=0;iAmpRange < rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			{
				file_rdbyamp << setw (20) << setfill(' ') << rdByAmpQuant->lambda << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->totalRate << " "
							 << setw (20) << setfill(' ') << i+1 << " "
							 << setw (20) << setfill(' ') << j+1 << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->numAmpRange[i][j] << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
							 << setw (20) << setfill(' ') << iAmpRange+1 << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->rate[i][j][iAmpRange] << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_rho << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_phase << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_xi << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_sample << " " << endl;
			}
		}
	}
	file_rdbyamp.close();

	/////////////////////////////////////////////////////
	// ENCODE WITH THE CHOSEN LAMBDA
	double** origSignal;
    origSignal = dataSignal->getSignal();

    double** recSignal;
	recSignal = new double*[sbbHeader.numSignal];

	double* recBlockSignal;
    recBlockSignal = new double[sbbHeader.blockSize];

    double** recSignalNoQuant;
	recSignalNoQuant = new double*[sbbHeader.numSignal];

    double* recBlockSignalNoQuant;
    recBlockSignalNoQuant = new double[sbbHeader.blockSize];

    double finalRate=0;
    double sqrerror_total=0.0;

    for (i=0;i<sbbHeader.numSignal;i++)
	{
		recSignal[i] = new double[sbbHeader.signalSize];
		recSignalNoQuant[i] = new double[sbbHeader.signalSize];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

			int NAmpRange = rdByAmpQuant->numAmpRange[i][j];

			CStructBook* structBookByAmp;
			structBookByAmp = new CStructBookExp[NAmpRange];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];

			double* ampRangeNumElement;
			ampRangeNumElement = new double[NAmpRange];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSB,
																		sbNumElement,
																		ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																		ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit
			}


			// Initialize block reconstruction vector (RESET)
			for (k=0;k <sbbHeader.blockSize ;k++)
			{
				recBlockSignal[k] = 0.0;
				recBlockSignalNoQuant[k] = 0.0;
			}

			for (iAmpRange=0;iAmpRange < rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			{
				strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
				int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElementAmp];
				int numElementAmpQ;

				strtQuantExp chosenQuantExp = rdByAmpQuant->quantExp[i][j][iAmpRange];

				quantizeStructBookExpRecUnquantNoReset(  min_amp[i][j],max_amp[i][j],
													  min_rho[i][j],max_rho[i][j],
													  min_phase[i][j],max_phase[i][j],
													  chosenQuantExp,
													  pSBAmp,numElementAmp,
													  pSBAmpQ,numElementAmpQ,
													  recBlockSignalNoQuant,
													  sbbHeader.norm[i],
													  sbbHeader.blockSize );

				printf(" ****** Synthesize signal %d block %d amp_range %d **** \n",i+1,j+1,iAmpRange+1);
				synthSignalSBExpNoReset(recBlockSignal,
										sbbHeader.norm[i],
										sbbHeader.blockSize,
										pSBAmpQ,
										numElementAmpQ);

				finalRate += computeRateSBExp(chosenQuantExp,numElementAmpQ);



				delete [] pSBAmpQ;
			}

			int initBlockSample = j*sbbHeader.blockSize;
			if ((initBlockSample + sbbHeader.blockSize) >= sbbHeader.signalSize)
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
				memcpy(	&recSignalNoQuant[i][initBlockSample],
						recBlockSignalNoQuant,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
			}
			else
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*sbbHeader.blockSize);
				memcpy(	&recSignalNoQuant[i][initBlockSample],
						recBlockSignalNoQuant,
						sizeof(double)*sbbHeader.blockSize);
			}

			delete [] (CStructBookExp*)structBookByAmp;
			delete [] ampBar;
			delete [] ampRangeLimit;
			delete [] ampRangeNumElement;
		} // block
		sqrerror_total += computeSqrError(origSignal[i],recSignal[i],0,sbbHeader.signalSize-1);
	} //signal

	double bitrate = static_cast<double>(finalRate)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Final rate: %20.10f \n",finalRate);
    printf("Signal size: %d \n",sbbHeader.signalSize);
    printf("Bitrate: %14.10f bit/sample\n",bitrate);
    printf("Total Square Error: %14.10f \n",sqrerror_total);

    //////// Save to file the reconstructed signals
    char fileName[_MAX_PATH];
    char* pos;
    char aux[_MAX_PATH];
    CDataSignal* outSignal;
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_rec%09.6fbps.wav",bitrate);
    strcpy( &pos[0], aux);
    //
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSignal(recSignal);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;

    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_rec%09.6fbps_noquant.wav",bitrate);
    strcpy( &pos[0], aux);
    //
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSignal(recSignalNoQuant);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;

	for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] min_amp[i];
    	delete [] max_amp[i];
    	delete [] min_rho[i];
    	delete [] max_rho[i];
    	delete [] min_phase[i];
    	delete [] max_phase[i];
	}
	delete [] min_amp;
	delete [] max_amp;
	delete [] min_rho;
	delete [] max_rho;
	delete [] min_phase;
	delete [] max_phase;

	//////////////////////////////////////////
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			delete [] rdByAmpQuant->rate[i][j];
			delete [] rdByAmpQuant->quantExp[i][j];
			delete [] rdByAmpQuant->nAtom[i][j];
			delete [] distAmp[i][j];
		}
		delete [] distAmp[i];
		delete [] rdByAmpQuant->rate[i];
		delete [] rdByAmpQuant->quantExp[i];
		delete [] rdByAmpQuant->numAmpRange[i];
		delete [] rdByAmpQuant->nAtom[i];
	}
	delete [] distAmp;
	delete [] rdByAmpQuant->rate;
	delete [] rdByAmpQuant->quantExp;
	delete [] rdByAmpQuant->numAmpRange;
	delete [] rdByAmpQuant->nAtom;
	delete rdByAmpQuant;
	/////////
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		delete [] recSignal[i];
		delete [] recSignalNoQuant[i];
	}
	delete [] recSignal;
	delete [] recSignalNoQuant;

	delete [] recBlockSignal;
	delete [] recBlockSignalNoQuant;
	/////

}


void  getRateNbAmpGivenLambda(  double lambda,
								strtSBBHeader sbbHeader,CStructBook** structBook,
								double** min_amp, double** max_amp,
								double** min_rho,double** max_rho,
								double** min_phase,double** max_phase,
								double**** DfTable,
								int N_qrho,int N_qphi,
								double* q_rho,double* q_phi,
								int nb_xi,int nb_sample,
								double*** distAmp,
								double& totalRate,
								double& chosen_nb_amp)
{
	int i,j, iAmpRange, i_qrho,i_qphi;

    int nb_amp;
    double nb_rho, nb_phi;

	totalRate = 0.0;
	double chosenRate;

	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{

			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();


			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];

			int sbNumElementQ;

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

// 			cout << "L_amp: " << L_amp << endl;
// 			cout << "L_rho: " << L_rho << endl;
// 			cout << "L_phi: " << L_phi << endl;

			double minLagrangian = 1e100;
			double LagrangianDAmp[16];
			double Lagrangian[16];
			for (nb_amp=1;nb_amp<=16;nb_amp++)
			{

				quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
											pSB, sbNumElement,
											pSBQ,sbNumElementQ);

				LagrangianDAmp[nb_amp-1] = distAmp[i][j][nb_amp-1] +
										  ((lambda) *
										   static_cast<double>(nb_amp + nb_xi + 2*nb_sample + 1 /* rho sign bit*/) *
										   sbNumElementQ);
			    Lagrangian[nb_amp-1] = LagrangianDAmp[nb_amp-1];


				int NAmpRange = nb_amp;

				CStructBook* structBookByAmp;
				structBookByAmp = new CStructBookExp[NAmpRange];

				double* ampRangeLimit;
				ampRangeLimit = new double[NAmpRange+1];

				double* ampBar;
				ampBar = new double[NAmpRange];

				double* ampRangeNumElement;
				ampRangeNumElement = new double[NAmpRange];

				// Separate atoms among amplitude ranges and return the number of atoms in each range
				ampRangeLimit[0]=-0.001;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
					ampBar[NAmpRange-1-iAmpRange] =
						(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

					((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
																			sbNumElementQ,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

					ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				}

				// Find the minimum Df Lagrangian and return the nb_rho and nb_phi
				int rate=0;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					// Find nb_rho and nb_phi leading minimum Df Lagragian

					double LagrangianDf;
					double minLagrangianDf=1e100;
					for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
					{
						if (log2( 1 + (L_rho/q_rho[i_qrho]) ) < 1) break;
						for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
						{
							if (log2( 1 + (L_phi/q_phi[i_qphi]) ) < 1) break;
							LagrangianDf =  (ampBar[iAmpRange]*ampBar[iAmpRange]) * DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] +
											lambda * ( log2( 1 + (L_rho/q_rho[i_qrho]) ) + log2( 1 + (L_phi/q_phi[i_qphi]) ) );

							if (LagrangianDf <= minLagrangianDf )
							{
								minLagrangianDf = LagrangianDf;
								nb_rho = log2( 1 + (L_rho/q_rho[i_qrho]) );
								nb_phi = log2( 1 + (L_phi/q_phi[i_qphi]) );
							}
						} // nb_phase
					}  // nb_rho

					Lagrangian[nb_amp-1] += ampRangeNumElement[iAmpRange] * minLagrangianDf;

// 					cout << " nb_amp: " << nb_amp << " "
// 						 << " iAmpRange: " << iAmpRange << " "
// 						 << " ampRangeNumElement[iAmpRange]: " << ampRangeNumElement[iAmpRange] << " "
// 						 << " nb_rho: " <<  nb_rho << " "
// 						 << " nb_phi: " << nb_phi << endl;

					strtQuantExp quantExp;
					quantExp.nb_amp = nb_amp;
					quantExp.nb_rho = static_cast<int>(round(nb_rho));
					quantExp.nb_phase = static_cast<int>(round(nb_phi));
					quantExp.nb_xi = nb_xi;
					quantExp.nb_sample = nb_sample;
					rate+= computeRateSBExp(quantExp,ampRangeNumElement[iAmpRange]);

				} // iAmpRange

				if (Lagrangian[nb_amp-1] < minLagrangian)
				{
					minLagrangian = Lagrangian[nb_amp-1];
					chosen_nb_amp = static_cast<double>(nb_amp);
					chosenRate = rate;
				}

				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;
				delete [] ampRangeNumElement;
			} // nb_amp
			delete [] pSBQ;
			totalRate += static_cast<double>(chosenRate);
		} // block
	} // signal
}

void computeRDByAmpQuant(   strtSBBHeader sbbHeader, CStructBook** structBook,
							double** min_amp, double** max_amp,
							double** min_rho, double** max_rho,
							double** min_phase, double** max_phase,
							double**** DfTable,
							int N_qrho, int N_qphi,
							double* q_rho, double* q_phi,
							int nb_xi, int nb_sample,
							double chosenNbAmp,
							strtRDByAmpQuant* rdByAmpQuant)
{
	int i,j, iAmpRange, i_qrho,i_qphi;

    int nb_amp;
    double nb_rho, nb_phi;

	double totalRate = 0.0;
	cout << setw (10) << setfill(' ') << " signal" << " "
		 << setw (10) << setfill(' ') << " block" << " "
		 << setw (10) << setfill(' ') << " nb_amp" << " "
		 << setw (10) << setfill(' ') << " iAmpRange"  << " "
		 << setw (10) << setfill(' ') << " NumEl" << " "
		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (10) << setfill(' ') << " R_xi"  << " "
		 << setw (10) << setfill(' ') << " R_sample"  << " "
		 << setw (10) << setfill(' ') << " R_rhosign"  << " "
		 << setw (10) << setfill(' ') << " iQRho" << " "
		 << setw (10) << setfill(' ') << " R_rho"  << " "
		 << setw (10) << setfill(' ') << " iQPhi" << " "
		 << setw (10) << setfill(' ') << " R_phi" << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{

			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();


			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];

			int sbNumElementQ;

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

			cout << "L_amp: " << L_amp << endl;
			cout << "L_rho: " << L_rho << endl;
			cout << "L_phi: " << L_phi << endl;

			nb_amp = static_cast<int>( chosenNbAmp );


			quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
										pSB, sbNumElement,
										pSBQ,sbNumElementQ);

			int NAmpRange = nb_amp;

			CStructBook* structBookByAmp;
			structBookByAmp = new CStructBookExp[NAmpRange];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];

			double* ampRangeNumElement;
			ampRangeNumElement = new double[NAmpRange];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
																		sbNumElementQ,
																		ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																		ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

				ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

			}

			int rate=0;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				// Find nb_rho and nb_phi leading minimum Df Lagragian
				double minLagrangianDf = 1e100;
				double LagrangianDf;
				int chosenQRhoIndex;
				int chosenQPhiIndex;
				for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
				{
					if (log2( 1 + (L_rho/q_rho[i_qrho]) ) < 1) break;
					for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
					{
						if (log2( 1 + (L_phi/q_phi[i_qphi]) ) < 1) break;
						LagrangianDf =  (ampBar[iAmpRange]*ampBar[iAmpRange]) *
											DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] +
										rdByAmpQuant->lambda * ( log2( 1 + (L_rho/q_rho[i_qrho]) ) + log2( 1 + (L_phi/q_phi[i_qphi]) ) );
						if (LagrangianDf <= minLagrangianDf )
						{
							minLagrangianDf = LagrangianDf;
							chosenQRhoIndex = i_qrho;
							chosenQPhiIndex = i_qphi;
							nb_rho = log2( 1 + (L_rho/q_rho[i_qrho]) );
							nb_phi = log2( 1 + (L_phi/q_phi[i_qphi]) );
						}
					} // nb_phase
				}  // nb_rho

				cout << setw (10) << setfill(' ') << i+1 << " "
					 << setw (10) << setfill(' ') << j+1 << " "
					 << setw (10) << setfill(' ') << nb_amp << " "
					 << setw (10) << setfill(' ') << iAmpRange+1 << " "
					 << setw (10) << setfill(' ') << ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') <<  nb_amp*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') <<  nb_xi*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') <<  2*nb_sample*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') <<  1*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') <<  chosenQRhoIndex<< " "
					 << setw (10) << setfill(' ') <<  static_cast<int>(round(nb_rho))*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') <<  chosenQPhiIndex << " "
					 << setw (10) << setfill(' ') << static_cast<int>(round(nb_phi))*ampRangeNumElement[iAmpRange] << endl;

				strtQuantExp quantExp;
				quantExp.nb_amp = nb_amp;
				quantExp.nb_rho = static_cast<int>(round(nb_rho));
				quantExp.nb_phase = static_cast<int>(round(nb_phi));
				quantExp.nb_xi = nb_xi;
				quantExp.nb_sample = nb_sample;

				rdByAmpQuant->rate[i][j][iAmpRange] = computeRateSBExp(quantExp,ampRangeNumElement[iAmpRange]);
				rdByAmpQuant->quantExp[i][j][iAmpRange] = quantExp;
				rdByAmpQuant->nAtom[i][j][iAmpRange] = ampRangeNumElement[iAmpRange];

				rate+=rdByAmpQuant->rate[i][j][iAmpRange];

			} // iAmpRange
			totalRate += static_cast<double>(rate);

			delete [] (CStructBookExp*)structBookByAmp;
			delete [] ampBar;
			delete [] ampRangeLimit;
			delete [] ampRangeNumElement;
			delete [] pSBQ;
		} // block
	} // signal
	rdByAmpQuant->totalRate = totalRate;
}

void encodeSBExpAmpRangeArithCod(char*InputFile,
								double rateTarget,
								strtSBBHeader sbbHeader,
								CStructBook** structBook,
								CFileDfTable* DfTableFile,
								CFileDictionary* dicData,
								CDataSignal* dataSignal,
								double initLambda)
{

	int Nfreq = dicData->getNumFreq();
	int fdiscrtype = dicData->getFDiscrType(0);
	double freqi = dicData->getFreqi(0);
	double step_xi = (2*pi/sbbHeader.Fs)*freqi;

    // Setting fixed quantizers
    double nb_xi = log2((double)Nfreq);
    double nb_sample = log2((double)sbbHeader.subBlockSize);
    int nb_rhoSign = 1;

	////////////////
	int i, j, k, i_qrho, i_qphi;
	int iAmpRange;
	int t;
	////////////////

    double** min_amp= new double*[sbbHeader.numSignal];
    double** max_amp= new double*[sbbHeader.numSignal];
    double** min_rho= new double*[sbbHeader.numSignal];
    double** max_rho= new double*[sbbHeader.numSignal];
    double** min_phase= new double*[sbbHeader.numSignal];
    double** max_phase= new double*[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		min_amp[i]= new double[sbbHeader.numBlock];
    	max_amp[i]= new double[sbbHeader.numBlock];
    	min_rho[i]= new double[sbbHeader.numBlock];
    	max_rho[i]= new double[sbbHeader.numBlock];
    	min_phase[i]= new double[sbbHeader.numBlock];
    	max_phase[i]= new double[sbbHeader.numBlock];
	}

	///////////======
	strtContinuousExp* pSB;
    int sbNumElement;

    double**** DfTable = DfTableFile->getDfTable();
    int N_qrho = DfTableFile->getNumQRho();
    int N_qphi = DfTableFile->getNumQPhi();
    double* q_rho = DfTableFile->getQRho();
    double* q_phi = DfTableFile->getQPhi();

    int nb_amp;
    double nb_rho, nb_phi;
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Find the range of the parameters for each Signal/Block
	double max_max_amp=0.0;
	double min_max_amp=1e30;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp[i][j],
																  max_amp[i][j],
																  min_rho[i][j],
																  max_rho[i][j],
																  min_phase[i][j],
																  max_phase[i][j]);
			if (max_max_amp<max_amp[i][j]) max_max_amp = max_amp[i][j];
			if (min_max_amp>max_amp[i][j]) min_max_amp = max_amp[i][j];
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Compute the Amplitude Square Error for each Signal/Block/nb_amp
	double*** distAmp;
	distAmp = new double**[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		distAmp[i] = new double*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			distAmp[i][j] = new double[16];
			////
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
			////
			for (nb_amp=1;nb_amp<16;nb_amp++)
			{
				distAmp[i][j][nb_amp-1] = sumSqrErrorStructBookExpAmp(  min_amp[i][j], max_amp[i][j], nb_amp,
																		pSB, sbNumElement);
			}
		}
	}

#ifdef DEBUG
	fstream file_deltasupqstep("rdbyamp_deltasupQStep.out",ios::out);
	file_deltasupqstep << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" " << endl;
#endif
	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
	cout << "- Compute probabilities with respect to TIME SUPPORT step quantization using" << endl;
	cout << "  Generalized Gaussian Distribution" << endl;
	int deltaSupMax=0;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			int deltaMaxAux = ((CStructBookExp*)structBook[i])[j].findDeltaSupMax();
			if (deltaMaxAux > deltaSupMax) deltaSupMax = deltaMaxAux;
		}
	}
	//int nQDeltaSupStep = sbbHeader.blockSize;
	int nQDeltaSupStep = deltaSupMax+1;
	double* deltaSupQStepCenter = new double[nQDeltaSupStep];
	double* deltaSupQStepEdge = new double[nQDeltaSupStep];
	double* deltaSupQStepProb = new double[nQDeltaSupStep];
	double* deltaSupQStepCumProb = new double[nQDeltaSupStep];
	double step_deltasup = 1.0;
// 	setQStepFeatureGGD( deltaSupQStepCenter,
// 						deltaSupQStepEdge,
// 						deltaSupQStepProb,
// 						deltaSupQStepCumProb,
// 						step_deltasup,
// 						nQDeltaSupStep,
// 						0,
// 						DELTASUPGGDMEAN, DELTASUPGGDSCALE, DELTASUPGGDSHAPE);
	setQStepFeatureGGDDeltaSup( deltaSupQStepCenter,
								deltaSupQStepEdge,
								deltaSupQStepProb,
								deltaSupQStepCumProb,
								step_deltasup,
								nQDeltaSupStep,
								0);
#ifdef DEBUG
	file_deltasupqstep<< setw (20) << setfill(' ') << i+1<<" "
				 << setw (20) << setfill(' ') << j+1<<" "
				 << setw (20) << setfill(' ') << deltaSupQStepCenter[0]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepCenter[static_cast<int>(nQDeltaSupStep/2)]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepCenter[nQDeltaSupStep-1]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepEdge[0]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepEdge[static_cast<int>(nQDeltaSupStep/2)]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepEdge[nQDeltaSupStep]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepProb[0]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepProb[static_cast<int>(nQDeltaSupStep/2)]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepProb[nQDeltaSupStep-1]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepCumProb[0]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepCumProb[static_cast<int>(nQDeltaSupStep/2)]<< " "
				 << setw (20) << setfill(' ') << deltaSupQStepCumProb[nQDeltaSupStep-1]<< " " << endl;
	for(i=0;i<nQDeltaSupStep;i++)
	{
		file_deltasupqstep << deltaSupQStepProb[i] <<endl;
	}

	file_deltasupqstep.close();
#endif
	double H_deltasup;
#ifdef USEGGDDSUP
	H_deltasup = computeEntropyPdf(deltaSupQStepProb,nQDeltaSupStep);
#else
	H_deltasup = computeEntropyUnifPdf(nQDeltaSupStep);
#endif

	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
	cout << "- Compute probabilities with respect to AMPLITUDE step quantization using" << endl;
	cout << "  Generalized Gaussian Distribution" << endl;
#ifdef DEBUG
	fstream file_ampqstep("rdbyamp_ampQStep.out",ios::out);
	file_ampqstep << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "nb_amp"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" " << endl;
#endif
	double**** ampQStepCenter = new double***[sbbHeader.numSignal];
	double**** ampQStepEdge = new double***[sbbHeader.numSignal];
	double**** ampQStepProb = new double***[sbbHeader.numSignal];
	double**** ampQStepCumProb = new double***[sbbHeader.numSignal];
	int nQAmpStep;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		ampQStepCenter[i] = new double**[sbbHeader.numBlock];
		ampQStepEdge[i] = new double**[sbbHeader.numBlock];
		ampQStepProb[i] = new double**[sbbHeader.numBlock];
	    ampQStepCumProb[i] = new double**[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			ampQStepCenter[i][j] = new double*[16];
			ampQStepEdge[i][j] = new double*[16];
			ampQStepProb[i][j] = new double*[16];
			ampQStepCumProb[i][j] = new double*[16];
			for (nb_amp=16; nb_amp>0; nb_amp--)
			{
				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp), min_amp[i][j], max_amp[i][j]);

				nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

				//cout << "-nb_amp: " << nb_amp << "; step_amp: " << step_amp << "; nQAmpStep: " << nQAmpStep << endl;

				// normalize step_amp
				step_amp = step_amp/max_amp[i][j];
				//cout << "step_amp normalized: " << step_amp << endl;
				ampQStepCenter[i][j][nb_amp-1] = new double[nQAmpStep];
				ampQStepEdge[i][j][nb_amp-1] = new double[nQAmpStep+1];
				ampQStepProb[i][j][nb_amp-1] = new double[nQAmpStep];
				ampQStepCumProb[i][j][nb_amp-1] = new double[nQAmpStep];

				setQStepFeatureGGD( ampQStepCenter[i][j][nb_amp-1],
									ampQStepEdge[i][j][nb_amp-1],
									ampQStepProb[i][j][nb_amp-1],
									ampQStepCumProb[i][j][nb_amp-1],
									step_amp,
									nQAmpStep,
									0,
									AMPGGDMEAN, AMPGGDSCALE, AMPGGDSHAPE);

#ifdef DEBUG
				file_ampqstep<< setw (20) << setfill(' ') << i+1<<" "
							 << setw (20) << setfill(' ') << j+1<<" "
							 << setw (20) << setfill(' ') << nb_amp<<" "
							 << setw (20) << setfill(' ') << ampQStepCenter[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepCenter[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepCenter[i][j][nb_amp-1][nQAmpStep-1]<< " "
							 << setw (20) << setfill(' ') << ampQStepEdge[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepEdge[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepEdge[i][j][nb_amp-1][nQAmpStep]<< " "
							 << setw (20) << setfill(' ') << ampQStepProb[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepProb[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepProb[i][j][nb_amp-1][nQAmpStep-1]<< " "
							 << setw (20) << setfill(' ') << ampQStepCumProb[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepCumProb[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepCumProb[i][j][nb_amp-1][nQAmpStep-1]<< " " << endl;
#endif

			}
		}
	}
#ifdef DEBUG
	file_ampqstep.close();
#endif

    int nQPhiStep;

	// Compute Entropy
	cout << "- Compute entropy table for Amp and Phase" << endl;
	double*** entropyAmpTable = new double**[sbbHeader.numSignal];
	double*** entropyPhiTable = new double**[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		entropyAmpTable[i] = new double*[sbbHeader.numBlock];
		entropyPhiTable[i] = new double*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			entropyAmpTable[i][j] = new double[16];
			entropyPhiTable[i][j] = new double[N_qphi];
			for (nb_amp=0;nb_amp<16;nb_amp++)
			{
				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp+1), min_amp[i][j], max_amp[i][j]);
				nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

#ifdef USEGGD
				entropyAmpTable[i][j][nb_amp] = computeEntropyPdf(	ampQStepProb[i][j][nb_amp],
																	nQAmpStep);
#else
				entropyAmpTable[i][j][nb_amp] = computeEntropySBExpDecayUnifPdf(nQAmpStep);
#endif
				//cout << "nb_amp: " << nb_amp << "; H_amp: " <<  entropyAmpTable[i][j][nb_amp] << endl;
			}
			for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
			{
				if ( (max_phase[i][j]/q_phi[i_qphi]) < 1.0) break;
				nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;

				entropyPhiTable[i][j][i_qphi] = computeEntropySBExpPhaseUnifPdf(nQPhiStep);
			}
		}
	}


	////////////////////



	double chosenLambda = initLambda;


	chosenLambda = initLambda;

	strtRDByAmpQuant* rdByAmpQuant = new strtRDByAmpQuant;

	rdByAmpQuant->rate = new double**[sbbHeader.numSignal];
	rdByAmpQuant->quantExp = new strtQuantExp**[sbbHeader.numSignal];
	rdByAmpQuant->nAtom = new int**[sbbHeader.numSignal];
	rdByAmpQuant->numAmpRange = new int*[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		rdByAmpQuant->rate[i] = new double*[sbbHeader.numBlock];
		rdByAmpQuant->quantExp[i] = new strtQuantExp*[sbbHeader.numBlock];
		rdByAmpQuant->numAmpRange[i] = new int[sbbHeader.numBlock];
		rdByAmpQuant->nAtom[i] = new int*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			int NAmpRange = 16;
			rdByAmpQuant->numAmpRange[i][j] = 0;
			rdByAmpQuant->rate[i][j] = new double[NAmpRange];
			rdByAmpQuant->quantExp[i][j] = new strtQuantExp[NAmpRange];
			rdByAmpQuant->nAtom[i][j] = new int[NAmpRange];
		}
	}
	cout << "- Proceed the bit allocation for each amplitude range." << endl;
	rdByAmpQuant->lambda = chosenLambda;
	cout << "chosenLambda: " << chosenLambda << endl;
	computeFullRDByAmpQuantArithCod(sbbHeader,structBook,
									min_amp, max_amp,
									min_rho, max_rho,
									min_phase, max_phase,
									DfTable,
									N_qrho, N_qphi,
									q_rho, q_phi,
									nb_xi, nb_sample,
									distAmp,
									rdByAmpQuant,
									ampQStepEdge,
									ampQStepProb,
									entropyAmpTable,
									entropyPhiTable,
									deltaSupQStepEdge,
									deltaSupQStepProb,
									H_deltasup,
									deltaSupMax);
	double totalRate = rdByAmpQuant->totalRate;

	cout << "#######################" << endl;
	cout << "# Chosen lambda: " << chosenLambda << endl;
	cout << "# total rate: " << totalRate << endl;
	cout << "#######################" << endl;

    ////////////////////////////////////////////////////////
// #ifdef DEBUG
//     fstream file_rdbyamp("rdbyamp_quant_arithcod.out",ios::out);
// 	file_rdbyamp << setw (20) << setfill(' ') << "Lambda"<<" "
// 	             << setw (20) << setfill(' ') << "Total Rate"<<" "
// 	             << setw (20) << setfill(' ') << "Signal"<<" "
// 	             << setw (20) << setfill(' ') << "Block"<<" "
// 	             << setw (20) << setfill(' ') << "NAmpRange/nbamp"<<" "
// 	             << setw (20) << setfill(' ') << "NAtomsInRange"<<" "
// 	             << setw (20) << setfill(' ') << "iAmpRange"<<" "
// 	             << setw (20) << setfill(' ') << "Rate"<<" "
// 	             << setw (20) << setfill(' ') << "qrho"<<" "
// 	             << setw (20) << setfill(' ') << "qphi"<<" "
// 	             << setw (20) << setfill(' ') << "nb_xi"<<" "
// 	             << setw (20) << setfill(' ') << "nb_sample"<<" " << endl;
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		for (j=0;j<sbbHeader.numBlock;j++)
// 		{
// 			for (iAmpRange=0;iAmpRange < rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
// 			{
// 				file_rdbyamp << setw (20) << setfill(' ') << rdByAmpQuant->lambda << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->totalRate << " "
// 							 << setw (20) << setfill(' ') << i+1 << " "
// 							 << setw (20) << setfill(' ') << j+1 << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->numAmpRange[i][j] << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
// 							 << setw (20) << setfill(' ') << iAmpRange+1 << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->rate[i][j][iAmpRange] << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].qrho << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].qphi << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_xi << " "
// 							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_sample << " " << endl;
// 			}
// 		}
// 	}
// 	file_rdbyamp.close();
// #endif
	/////////////////////////////////////////////////////
	// ENCODE WITH THE CHOSEN LAMBDA
	/////

	cout << " >>> ENCODE WITH THE CHOSEN LAMBDA <<<<<" << endl;

    int iSubBlock;
    int subBlockSize = sbbHeader.subBlockSize;
    int numSubBlock = static_cast<int>(sbbHeader.blockSize / subBlockSize);

	double** origSignal;
    origSignal = dataSignal->getSignal();

    double** recSignal;
	recSignal = new double*[sbbHeader.numSignal];

	double* recBlockSignal;
    recBlockSignal = new double[sbbHeader.blockSize];

    double** recSignalNoQuant;
	recSignalNoQuant = new double*[sbbHeader.numSignal];

    double* recBlockSignalNoQuant;
    recBlockSignalNoQuant = new double[sbbHeader.blockSize];

    double finalRate=0;
    double sqrerror_total=0.0;


    /////
#ifdef DEBUG
    fstream file_encIndex("rdbyamp_encodeIndex.out",ios::out);
	file_encIndex << setw (14) << setfill(' ') << "Signal"<<" "
	              << setw (14) << setfill(' ') << "Block"<<" "
	              << setw (14) << setfill(' ') << "subBlock"<<" "
	              << setw (14) << setfill(' ') << "iAmpRange"<<" "
	              << setw (14) << setfill(' ') << "Element"<<" "
	              << setw (14) << setfill(' ') << "ampIndex"<<" "
	              << setw (14) << setfill(' ') << "rhoIndex"<<" "
	              << setw (14) << setfill(' ') << "xiIndex"<<" "
	              << setw (14) << setfill(' ') << "phiIndex"<<" "
	              << setw (14) << setfill(' ') << "initSamp"<<" "
	              << setw (14) << setfill(' ') << "deltaSamp"<<" "
	              << setw (14) << setfill(' ') << "aIndex"<<" "
	              << setw (14) << setfill(' ') << "bIndex"<<" " << endl;
#endif

    FILE* code_file;
    //////// Save encoded file
    char fileName[_MAX_PATH];
    char* pos;
    char aux[_MAX_PATH];
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,".mpz");
    strcpy( &pos[0], aux);
    //
    code_file = fopen (fileName,"wb");


	int numBytesHeader=0;
	int numBytes=0;

	unsigned short int dummyushint;
	unsigned int dummyuint;
	int dummyint;
	float dummyfloat;
	double dummydouble;
    /// FILL HEADER //////////////////
    // numSignal
    dummyint =sbbHeader.numSignal;
    fwrite (&dummyint, 1 , sizeof(unsigned int) , code_file );
    numBytesHeader += sizeof(int);
    // numBlock
    //dummyint =sbbHeader.numBlock;
    //fwrite (&dummyint, 1 , sizeof(int) , code_file );
    //numBytesHeader += sizeof(int);
    // signalSize
    dummyint = sbbHeader.signalSize;
    fwrite (&dummyint, 1 , sizeof(unsigned int) , code_file );
    numBytesHeader += sizeof(int);
    // blockSize
    dummyint =sbbHeader.blockSize;
    fwrite (&dummyint, 1 , sizeof(unsigned int) , code_file );
    numBytesHeader += sizeof(int);
    // subBlockSize
    dummyint = subBlockSize;
    fwrite (&dummyint, 1 , sizeof(unsigned int) , code_file );
    numBytesHeader += sizeof(unsigned int);
    // Fs - sample rate
    dummyuint = static_cast<unsigned int>(sbbHeader.Fs);
    fwrite (&dummyuint, 1 , sizeof(unsigned int) , code_file );
    numBytesHeader += sizeof(unsigned int);
    // lambda
    dummyfloat = static_cast<float>(rdByAmpQuant->lambda);
    fwrite (&dummyfloat, 1 , sizeof(float) , code_file );
    numBytesHeader += sizeof(float);
    // signal norm
    float* signalNorm = new float[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
    {
    	signalNorm[i] = static_cast<float>(sbbHeader.norm[i]);
    	cout << "Signal norm - " << i+1 << " : " << signalNorm[i] << endl;
    }
    fwrite (signalNorm, sbbHeader.numSignal  , sizeof(float) , code_file );
    numBytesHeader += sizeof(float) * sbbHeader.numSignal;
    delete [] signalNorm;

    for (i=0;i<sbbHeader.numSignal;i++)
	{
    	for (j=0;j<sbbHeader.numBlock;j++)
		{
    		/// FILL HEADER //////////////////
			// nbamp/NAmpRAnge
			dummyuint = static_cast<unsigned int>(rdByAmpQuant->numAmpRange[i][j]);
			cout << "nb_amp: " << dummyuint << endl;
    		fwrite (&dummyuint, 1 , sizeof(unsigned int) , code_file );
    		numBytesHeader += sizeof(unsigned int);
			// max_amp
			dummyfloat = static_cast<float>(max_amp[i][j]);
			cout << "max_amp: " << dummyfloat << endl;
    		fwrite (&dummyfloat, 1 , sizeof(float) , code_file );
    		numBytesHeader += sizeof(float);
			// max_rho
			dummyfloat = static_cast<float>(max_rho[i][j]);
			cout << "max_rho: " << dummyfloat << endl;
    		fwrite (&dummyfloat, 1 , sizeof(float) , code_file );
    		numBytesHeader += sizeof(float);
    		// deltaSupMax
			dummyuint = static_cast<unsigned int>(deltaSupMax);
			cout << "deltaSupMax: " << dummyuint << endl;
    		fwrite (&dummyuint, 1 , sizeof(unsigned int) , code_file );
    		numBytesHeader += sizeof(unsigned int);
			// OBS:
			//   min_amp = 0
			//   min_rho =0
			//   min_phi = 0
			//   max_phi = 2*pi
		}
	}

#ifdef USEGGD
	cout << N_qrho << endl;
	double** rhoQStepCenter = new double*[N_qrho];
	double** rhoQStepEdge = new double*[N_qrho];
	double** rhoQStepProb = new double*[N_qrho];
	double** rhoQStepCumProb = new double*[N_qrho];
#endif


    for (i=0;i<sbbHeader.numSignal;i++)
	{
		recSignal[i] = new double[sbbHeader.signalSize];
		recSignalNoQuant[i] = new double[sbbHeader.signalSize];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			cout << "-->>> Encoding Signal " << i+1<< "; Block " << j+1 << endl;
			////////////////////////////////////////////////////////////////////////
			// Compute probabilities with respect to decay step quantization
			cout << "- Compute probabilities and entropies with respect to DECAY step quantization using" << endl;
			cout << "  Generalized Gaussian Distribution " << endl;

			int nQRhoStep;
#ifdef USEGGD
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{

				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;

				nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;
// 				cout << " i_qrho: " << i_qrho << endl ;
// 				cout << " nQRhoStep: " << nQRhoStep << " -- "<< (max_rho[i][j]/q_rho[i_qrho]) << endl ;

				rhoQStepCenter[i_qrho] = new double[nQRhoStep];
				rhoQStepEdge[i_qrho] = new double[nQRhoStep+1];
				rhoQStepProb[i_qrho] = new double[nQRhoStep];
				rhoQStepCumProb[i_qrho] = new double[nQRhoStep];

				setQStepFeatureGGD( rhoQStepCenter[i_qrho],
									rhoQStepEdge[i_qrho],
									rhoQStepProb[i_qrho],
									rhoQStepCumProb[i_qrho],
									q_rho[i_qrho],
									nQRhoStep,
									static_cast<int>(nQRhoStep/2),
									DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);
			}
#endif

			/////

			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			///////
			double L_amp = max_amp[i][j] - min_amp[i][j];

			int NAmpRange = rdByAmpQuant->numAmpRange[i][j];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];

			double* ampRangeNumElement;
			ampRangeNumElement = new double[NAmpRange];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));
			}
			/////////

			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];
			int sbNumElementQ;

			quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], NAmpRange,
										pSB, sbNumElement,
										pSBQ,sbNumElementQ);

			cout << "sbNumElementQ: " << sbNumElementQ << endl;

			/////
			int* ampIndex = new int[sbNumElementQ];
			int* rhoIndex = new int[sbNumElementQ];
			int* xiIndex  = new int[sbNumElementQ];
			int* phiIndex = new int[sbNumElementQ];
			int* aIndex   = new int[sbNumElementQ];
			int* bIndex   = new int[sbNumElementQ];
			//////
			int numIndex=0;
			int** tabNumIndex = new int*[numSubBlock];
			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				tabNumIndex[iSubBlock] = new int[NAmpRange+1];
			}
			/////////
			for (k=0;k<sbbHeader.blockSize;k++)
			{
				recBlockSignal[k] = 0.0;
			}
			////////////////////////////////////
			// Arithmetic Encoder Varaibles
			unsigned char* compressed_data=0;
			unsigned max_code_bytes = 0xFFFFFFE0U;

			Arithmetic_Codec ace;
			ace.set_buffer(max_code_bytes, compressed_data);

			// START ENCODER
			ace.start_encoder();
			////////////////////////////////////
			// MODEL AMP
			Static_Data_Model ampModel;
			double step_amp = computeMidRiseQuantStep( static_cast<double>(NAmpRange), min_amp[i][j], max_amp[i][j]);
			unsigned nlevel_amp = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;
			//unsigned nlevel_amp = static_cast<unsigned>(round(pow(2.0, static_cast<double>(chosenQuantExp.nb_amp) ) ));
			cout << "- nlevel_amp: " << nlevel_amp << endl;
#ifdef USEGGD
			ampModel.set_distribution(static_cast<unsigned>(nlevel_amp),
									  ampQStepProb[i][j][NAmpRange - 1]);
#else
			ampModel.set_distribution(static_cast<unsigned>(nlevel_amp));
#endif

			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				CStructBook* structBookByBlock = new CStructBookExp;

				((CStructBookExp*)structBookByBlock)->sepBySubBlock(pSBQ,
																	sbNumElementQ,
																	iSubBlock*subBlockSize,      //lowerSubBlockLimit
																	(iSubBlock+1)*subBlockSize); //upperSubBlockLimit

				strtContinuousExp* pSBSubBlock = ((CStructBookExp*)structBookByBlock)->getStructBook();
				int numElementSubBlock = ((CStructBookExp*)structBookByBlock)->getNumElement();

				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					CStructBook* structBookByAmp = new CStructBookExp;

					((CStructBookExp*)structBookByAmp)->sepByAmp(	pSBSubBlock,
																	numElementSubBlock,
																	ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																	ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

					strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)->getStructBook();
					int numElementAmp = ((CStructBookExp*)structBookByAmp)->getNumElement();

					strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElementAmp];
					int numElementAmpQ;

					strtQuantExp chosenQuantExp = rdByAmpQuant->quantExp[i][j][iAmpRange];

					quantizeStructBookExpRecUnquantNoResetArithCod(	  min_amp[i][j],max_amp[i][j],
																	  chosenQuantExp,
																	  pSBAmp,numElementAmp,
																	  pSBAmpQ,numElementAmpQ,
																	  recBlockSignalNoQuant,
																	  sbbHeader.norm[i],
																	  sbbHeader.blockSize,
																	  &ampIndex[numIndex], &rhoIndex[numIndex], &xiIndex[numIndex],
																	  &phiIndex[numIndex], &aIndex[numIndex], &bIndex[numIndex],
																	  fdiscrtype, step_xi);


					cout << "- chosenQuantExp.i_qrho: " << chosenQuantExp.i_qrho << endl;
					cout << "- chosenQuantExp.i_qphi: " << chosenQuantExp.i_qphi << endl;

					printf(" ****** Synthesize signal %d block %d subBlock %d amp_range %d **** \n",i+1,j+1,iSubBlock+1,iAmpRange+1);
					synthSignalSBExpNoReset(recBlockSignal,
											sbbHeader.norm[i],
											sbbHeader.blockSize,
											pSBAmpQ,
											numElementAmpQ);


					tabNumIndex[iSubBlock][iAmpRange] = numIndex;
					//cout << "tabNumIndex[iSubBlock][iAmpRange]: "<< tabNumIndex[iSubBlock][iAmpRange] <<endl;
					//cout << "numElementAmpQ: " << numElementAmpQ << endl;
					numIndex += numElementAmpQ;
					//cout << "numIndex: " << numIndex << endl;


					/// DEALLOCATE /////
					delete structBookByAmp;
					delete [] pSBAmpQ;
				} // iAmpRange

				tabNumIndex[iSubBlock][NAmpRange] = numIndex;
				//cout << "tabNumIndex[iSubBlock][NAmpRange]: "<< tabNumIndex[iSubBlock][NAmpRange] <<endl;

				// ENCODE AMP
				cout << tabNumIndex[iSubBlock][0] << " " <<  tabNumIndex[iSubBlock][NAmpRange] << endl;
				for (k=tabNumIndex[iSubBlock][0];k < tabNumIndex[iSubBlock][NAmpRange];k++)
				{
					ace.encode(static_cast<unsigned>(ampIndex[k]),ampModel);
				}
				// ESCAPE SYMBOL: 0
				ace.encode(static_cast<unsigned>(0),ampModel);

				/// DEALLOCATE /////
				delete structBookByBlock;
				////

			}// iSubBlock

			unsigned nbyte_amp = ace.write_to_file(code_file); // Alternalively => ace.stop_encoder();
			numBytes += static_cast<int>(nbyte_amp);

			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				// START ENCODER
				ace.start_encoder();

				/// MODEL FREQUENCY
				Static_Data_Model xiModel;
				cout << "- Nfreq: " << Nfreq << endl;
				xiModel.set_distribution(static_cast<unsigned>(Nfreq));

				/// MODEL TIME SUPPORT
				/// A
				Static_Data_Model initTimeSupModel;
				cout << "- sbbHeader.subBlockSize: " << subBlockSize << endl;
				initTimeSupModel.set_distribution(static_cast<unsigned>(subBlockSize));

				/// B/ DELTA
				Static_Data_Model timeSupModel;
				//cout << "- sbbHeader.blockSize: " << sbbHeader.blockSize << endl;
				cout << "- deltaSupMax: " << deltaSupMax << endl;
#ifdef USEGGDDSUP
				timeSupModel.set_distribution(static_cast<unsigned>(deltaSupMax+1), deltaSupQStepProb);
#else
				timeSupModel.set_distribution(static_cast<unsigned>(deltaSupMax+1));
#endif

				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					strtQuantExp chosenQuantExp = rdByAmpQuant->quantExp[i][j][iAmpRange];

					/// MODEL DECAY
					nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/chosenQuantExp.qrho))*2 + 1;
					Static_Data_Model decayModel;
					cout << "- nQRhoStep: " << nQRhoStep << endl;
					cout << "   max_rho: " << max_rho[i][j] << endl;
					cout << "   chosenQuantExp.qrho: " << chosenQuantExp.qrho << endl;
#ifdef USEGGD
					decayModel.set_distribution(static_cast<unsigned>(nQRhoStep),
												rhoQStepProb[chosenQuantExp.i_qrho]);
#else
					decayModel.set_distribution(static_cast<unsigned>(nQRhoStep));
#endif
					/// MODEL PHASE
					int nQPhiStep = static_cast<int>(ceil(max_phase[i][j]/chosenQuantExp.qphi)) + 1;
					Static_Data_Model phiModel;
					cout << "- nQPhiStep: " << nQPhiStep << endl;
					phiModel.set_distribution(static_cast<unsigned>(nQPhiStep));

					for (k=tabNumIndex[iSubBlock][iAmpRange];k<tabNumIndex[iSubBlock][iAmpRange+1];k++)
					{
						// ENCODE DECAY
						ace.encode(static_cast<unsigned>(rhoIndex[k] + static_cast<int>(nQRhoStep/2)),decayModel);
						// ENCODE PHASE
						ace.encode(static_cast<unsigned>(phiIndex[k]),phiModel);
						// ENCODE FREQUENCY
						ace.encode(static_cast<unsigned>(xiIndex[k]),xiModel);
						// ENCODE A
						ace.encode(static_cast<unsigned>(aIndex[k]%subBlockSize),initTimeSupModel);
						// ENCODE B/DELTA
						ace.encode(static_cast<unsigned>(bIndex[k]-aIndex[k]),timeSupModel);
					}
#ifdef DEBUG
					for (k=tabNumIndex[iSubBlock][iAmpRange];k<tabNumIndex[iSubBlock][iAmpRange+1];k++)
					{
						file_encIndex << setw (14) << setfill(' ') << i+1<<" "
									  << setw (14) << setfill(' ') << j+1<<" "
									  << setw (14) << setfill(' ') << iSubBlock+1<<" "
									  << setw (14) << setfill(' ') << iAmpRange+1<<" "
									  << setw (14) << setfill(' ') << k+1<<" "
									  << setw (14) << setfill(' ') << ampIndex[k] <<" "
									  << setw (14) << setfill(' ') << rhoIndex[k] <<" "
									  << setw (14) << setfill(' ') << xiIndex[k]  <<" "
									  << setw (14) << setfill(' ') << phiIndex[k] <<" "
									  << setw (14) << setfill(' ') << aIndex[k]%subBlockSize   <<" "
									  << setw (14) << setfill(' ') << bIndex[k]-aIndex[k]   <<" "
									  << setw (14) << setfill(' ') << aIndex[k]   <<" "
									  << setw (14) << setfill(' ') << bIndex[k]   <<" " << endl;
					}
#endif
				}// iAmpRange
			 	nbyte_amp = ace.write_to_file(code_file); // Alternalively => ace.stop_encoder();
				numBytes += static_cast<int>(nbyte_amp);

			} // iSubBlock


			int initBlockSample = j*sbbHeader.blockSize;
			if ((initBlockSample + sbbHeader.blockSize) >= sbbHeader.signalSize)
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
				memcpy(	&recSignalNoQuant[i][initBlockSample],
						recBlockSignalNoQuant,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
			}
			else
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*sbbHeader.blockSize);
				memcpy(	&recSignalNoQuant[i][initBlockSample],
						recBlockSignalNoQuant,
						sizeof(double)*sbbHeader.blockSize);
			}

			delete [] ampBar;
			delete [] ampRangeLimit;
			delete [] ampRangeNumElement;
			////
			delete [] ampIndex;
			delete [] rhoIndex;
			delete [] xiIndex;
			delete [] phiIndex;
			delete [] aIndex;
			delete [] bIndex;

#ifdef USEGGD
			for (i_qrho=0; i_qrho<N_qrho; i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
				delete [] rhoQStepCenter[i_qrho];
				delete [] rhoQStepEdge[i_qrho];
				delete [] rhoQStepProb[i_qrho];
				delete [] rhoQStepCumProb[i_qrho];
			}
#endif
			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				delete [] tabNumIndex[iSubBlock];
			}
			delete [] tabNumIndex;

			/////
			delete [] pSBQ;

		} // block
		sqrerror_total += computeSqrError(origSignal[i],recSignal[i],0,sbbHeader.signalSize-1);
	} //signal

#ifdef USEGGD
	delete [] rhoQStepCenter;
	delete [] rhoQStepEdge;
	delete [] rhoQStepProb;
	delete [] rhoQStepCumProb;
#endif

	fclose(code_file);
#ifdef DEBUG
	file_encIndex.close();
#endif
	finalRate = totalRate;
	double bitrate = static_cast<double>(finalRate)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Final rate calculated: %20.10f \n",finalRate);
    printf("Signal size: %d \n",sbbHeader.signalSize);
    printf("Bitrate calculated: %14.10f bit/sample\n",bitrate);
    printf("Total Square Error: %14.10f \n",sqrerror_total);
    /////////
    printf("Final rate from arithmetic codec: %20.10f bits \n",static_cast<double>(numBytes*8));
    printf("Header rate from arithmetic codec: %20.10f bits\n",static_cast<double>(numBytesHeader*8));
    double bitrateAC = static_cast<double>(numBytes*8)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
	double bitrateACHeader = static_cast<double>(numBytesHeader*8)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Bitrate from arithmetic codec: %14.10f bit/sample\n",bitrateAC);
    printf("Bitrate header from arithmetic codec: %14.10f bit/sample\n",bitrateACHeader);

    fstream file_enc_report("rdbyamp_encode_report.out",ios::out);
	file_enc_report << setw (25) << setfill(' ') << "Lambda"<<" "
	              << setw (25) << setfill(' ') << "Bitrate Calc. (b/s)"<<" "
	              << setw (30) << setfill(' ') << "Bitrate Arith. Cod. (b/s)"<<" "
	              << setw (25) << setfill(' ') << "Header bitrate (b/s)"<<" "
	              << setw (25) << setfill(' ') << "Distortion (Sqr. Error)"<< " " << endl;

	file_enc_report << setw (25) << setfill(' ') << setprecision(10) << chosenLambda <<" "
	              << setw (30) << setfill(' ') << setprecision(10) << bitrate<<" "
	              << setw (25) << setfill(' ') << setprecision(10) << bitrateAC<<" "
	              << setw (25) << setfill(' ') << setprecision(10) << bitrateACHeader <<" "
	              << setw (25) << setfill(' ') << setprecision(10) << sqrerror_total << " " << endl;
	file_enc_report.close();

    //////// Save to file the reconstructed signals
    CDataSignal* outSignal;
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_rec%09.6fbps.wav",bitrateAC);
    strcpy( &pos[0], aux);
    //
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSamplingRate(sbbHeader.Fs);
    outSignal->setSignal(recSignal);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;

    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_rec%09.6fbps_noquant.wav",bitrateAC);
    strcpy( &pos[0], aux);
    //
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSamplingRate(sbbHeader.Fs);
    outSignal->setSignal(recSignalNoQuant);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;



    //////////////////////
    /// DEALLOCATE ///////
    delete [] deltaSupQStepCenter;
	delete [] deltaSupQStepEdge;
	delete [] deltaSupQStepProb;
	delete [] deltaSupQStepCumProb;
	///////////
    for (i=0;i<sbbHeader.numSignal;i++)
	{
    	for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (nb_amp=0;nb_amp<16;nb_amp++)
			{
				delete [] ampQStepCenter[i][j][nb_amp];
				delete [] ampQStepEdge[i][j][nb_amp];
				delete [] ampQStepProb[i][j][nb_amp];
				delete [] ampQStepCumProb[i][j][nb_amp];
			}
			delete [] ampQStepCenter[i][j];
			delete [] ampQStepEdge[i][j];
			delete [] ampQStepProb[i][j];
			delete [] ampQStepCumProb[i][j];
		}
		delete [] ampQStepCenter[i];
		delete [] ampQStepEdge[i];
		delete [] ampQStepProb[i];
		delete [] ampQStepCumProb[i];
	}
	delete [] ampQStepCenter;
	delete [] ampQStepEdge;
	delete [] ampQStepProb;
	delete [] ampQStepCumProb;
    ///////////////////////
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] min_amp[i];
    	delete [] max_amp[i];
    	delete [] min_rho[i];
    	delete [] max_rho[i];
    	delete [] min_phase[i];
    	delete [] max_phase[i];
	}
	delete [] min_amp;
	delete [] max_amp;
	delete [] min_rho;
	delete [] max_rho;
	delete [] min_phase;
	delete [] max_phase;

	//////////////////////////////////////////
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			delete [] rdByAmpQuant->rate[i][j];
			delete [] rdByAmpQuant->quantExp[i][j];
			delete [] rdByAmpQuant->nAtom[i][j];
			delete [] distAmp[i][j];
		}
		delete [] distAmp[i];
		delete [] rdByAmpQuant->rate[i];
		delete [] rdByAmpQuant->quantExp[i];
		delete [] rdByAmpQuant->numAmpRange[i];
		delete [] rdByAmpQuant->nAtom[i];
	}
	delete [] distAmp;
	delete [] rdByAmpQuant->rate;
	delete [] rdByAmpQuant->quantExp;
	delete [] rdByAmpQuant->numAmpRange;
	delete [] rdByAmpQuant->nAtom;
	delete rdByAmpQuant;
	/////////
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		delete [] recSignal[i];
		delete [] recSignalNoQuant[i];
	}
	delete [] recSignal;
	delete [] recSignalNoQuant;

	delete [] recBlockSignal;
	delete [] recBlockSignalNoQuant;
	////////////
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			delete [] entropyAmpTable[i][j];
			delete [] entropyPhiTable[i][j];
		}
		delete [] entropyAmpTable[i];
		delete [] entropyPhiTable[i];
	}
	delete [] entropyAmpTable;
	delete [] entropyPhiTable;
}

void  computeFullRDByAmpQuantArithCod(  strtSBBHeader sbbHeader,CStructBook** structBook,
										double** min_amp, double** max_amp,
										double** min_rho,double** max_rho,
										double** min_phase,double** max_phase,
										double**** DfTable,
										int N_qrho,int N_qphi,
										double* q_rho,double* q_phi,
										double nb_xi,double nb_sample,
										double*** distAmp,
										strtRDByAmpQuant* rdByAmpQuant,
										double**** ampQStepEdge,
										double**** ampQStepProb,
										double*** entropyAmpTable,
										double*** entropyPhiTable,
										double* deltaSupQStepEdge,
										double* deltaSupQStepProb,
										double H_deltasup,
										int deltaSupMax)
{
	int i,j, iAmpRange, i_qrho,i_qphi, k;

    int nb_amp;
    double nb_rho, nb_phi;

	rdByAmpQuant->totalRate = 0.0;
	double chosenRate;

#ifdef USEGGD
	double** rhoQStepCenter = new double*[N_qrho];
	double** rhoQStepEdge = new double*[N_qrho];
	double** rhoQStepProb = new double*[N_qrho];
	double** rhoQStepCumProb = new double*[N_qrho];
#endif
	double* entropyRhoTable = new double[N_qrho];

#ifdef DEBUG
	fstream file_rhoqstep("rdbyamp_rhoQStep.out",ios::out);
	file_rhoqstep << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "i_qrho"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" " << endl;
#endif


	cout << setw (5) << setfill(' ') << " sig" << " "
		 << setw (4) << setfill(' ') << " bl" << " "
		 << setw (7) << setfill(' ') << " nbamp" << " "
		 << setw (7) << setfill(' ') << " iAmpR"  << " "
		 << setw (10) << setfill(' ') << " NumEl" << " "
		 << setw (10) << setfill(' ') << " R_xi"  << " "
		 << setw (10) << setfill(' ') << " R_sample"  << " "
		 << setw (10) << setfill(' ') << " R_deltasup"  << " "
		 << setw (10) << setfill(' ') << " H_deltasup"  << " "
		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (14) << setfill(' ') << " H_amp" << " "
		 << setw (10) << setfill(' ') << " iQRho" << " "
		 << setw (10) << setfill(' ') << " R_rho"  << " "
		 << setw (10) << setfill(' ') << " H_rho" << " "
		 << setw (10) << setfill(' ') << " iQPhi" << " "
		 << setw (10) << setfill(' ') << " R_phi" << " "
		 << setw (10) << setfill(' ') << " H_phi" << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			cout << "-->>> Compute RDByAmpQuant for Signal " << i+1<< "; Block " << j+1 << endl;
			////////////////////////////////////////////////////////////////////////
			// Compute probabilities with respect to decay step quantization
			cout << "- Compute probabilities and entropies with respect to DECAY step quantization using" << endl;
			cout << "  Generalized Gaussian Distribution " << endl;
			//cout << N_qrho << endl;

			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;

 				int nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

#ifdef USEGGD
				rhoQStepCenter[i_qrho] = new double[nQRhoStep];
				rhoQStepEdge[i_qrho] = new double[nQRhoStep+1];
				rhoQStepProb[i_qrho] = new double[nQRhoStep];
				rhoQStepCumProb[i_qrho] = new double[nQRhoStep];

				setQStepFeatureGGD( rhoQStepCenter[i_qrho],
									rhoQStepEdge[i_qrho],
									rhoQStepProb[i_qrho],
									rhoQStepCumProb[i_qrho],
									q_rho[i_qrho],
									nQRhoStep,
									static_cast<int>(nQRhoStep/2),
									DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);

				entropyRhoTable[i_qrho] = computeEntropySBExpDecayGGDPdf(	rhoQStepProb[i_qrho],
																			nQRhoStep);
#else
				entropyRhoTable[i_qrho] = computeEntropySBExpDecayUnifPdf(nQRhoStep);
#endif

#ifdef USEGGD
#ifdef DEBUG
				file_rhoqstep<< setw (20) << setfill(' ') << i+1<<" "
							 << setw (20) << setfill(' ') << j+1<<" "
							 << setw (20) << setfill(' ') << i_qrho+1<<" "
							 << setw (20) << setfill(' ') << rhoQStepCenter[i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCenter[i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCenter[i_qrho][nQRhoStep-1]<< " "
							 << setw (20) << setfill(' ') << rhoQStepEdge[i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepEdge[i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepEdge[i_qrho][nQRhoStep]<< " "
							 << setw (20) << setfill(' ') << rhoQStepProb[i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepProb[i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepProb[i_qrho][nQRhoStep-1]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCumProb[i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCumProb[i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCumProb[i_qrho][nQRhoStep-1]<< " " << endl;
#endif
#endif
			}

			////////////////////////////////////////////////////////////////////////////////
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];

			int sbNumElementQ;

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

			double minLagrangian = 1e100;
			double Lagrangian;
			double R_amp, chosenR_amp;
			double H_amp, chosenH_amp;
			double R_deltasup, chosenR_deltasup;
			double chosenH_deltasup;
			double chosen_nb_amp=0;
			for (nb_amp=1;nb_amp<=16;nb_amp++)
			{

				quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
											pSB, sbNumElement,
											pSBQ,sbNumElementQ);

				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp), min_amp[i][j], max_amp[i][j]);
				int nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

#ifdef USEGGD
				double* ampvec = new double[sbNumElementQ];

				// normalize amp by max_amp (max_amp=1)
				for (k=0;k<sbNumElementQ;k++)
					ampvec[k] = ( pSBQ[k].innerProduct / max_amp[i][j]);

				bubble_srtdouble(ampvec, sbNumElementQ);

				R_amp = computeRatePdf( ampvec,
										sbNumElementQ,
										ampQStepEdge[i][j][nb_amp-1],
										ampQStepProb[i][j][nb_amp-1],
										nQAmpStep);
				delete [] ampvec;
#else
				R_amp = computeRateUnifPdf(nQAmpStep, sbNumElementQ);
#endif
				H_amp = entropyAmpTable[i][j][nb_amp-1];

				////////////////////
				int nQDeltaSupStep = deltaSupMax+1;
#ifdef USEGGDDSUP
				double* deltasupvec = new double[sbNumElementQ];

				// normalize amp by max_amp (max_amp=1)
				for (k=0;k<sbNumElementQ;k++)
					deltasupvec[k] = ( pSBQ[k].b - pSBQ[k].a);

				bubble_srtdouble(deltasupvec, sbNumElementQ);

				R_deltasup = computeRatePdf(deltasupvec,
											sbNumElementQ,
											deltaSupQStepEdge,
											deltaSupQStepProb,
											nQDeltaSupStep);
				delete [] deltasupvec;
#else
				R_deltasup = computeRateUnifPdf(nQDeltaSupStep, sbNumElementQ);
#endif

				////////////////////

				Lagrangian = distAmp[i][j][nb_amp-1] +
										   ((rdByAmpQuant->lambda) *
										    static_cast<double>(H_amp + nb_xi + nb_sample + H_deltasup) *
										    sbNumElementQ );

				int NAmpRange = nb_amp;

				CStructBook* structBookByAmp;
				structBookByAmp = new CStructBookExp[NAmpRange];

				double* ampRangeLimit;
				ampRangeLimit = new double[NAmpRange+1];

				double* ampBar;
				ampBar = new double[NAmpRange];

				double* ampRangeNumElement;
				ampRangeNumElement = new double[NAmpRange];

				// Separate atoms among amplitude ranges and return the number of atoms in each range
				ampRangeLimit[0]=-0.001;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
					ampBar[NAmpRange-1-iAmpRange] =
						(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

					((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
																			sbNumElementQ,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

					ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				}

				// Find the minimum Df Lagrangian and return the nb_rho and nb_phi
				double* rateAmpRange = new double[NAmpRange];
				strtQuantExp* quantExp = new strtQuantExp[NAmpRange];
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					///////
					strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
					int sbNumElementAmp = ampRangeNumElement[iAmpRange];
					double* rhovec = new double[sbNumElementAmp];

					for (k=0;k<sbNumElementAmp;k++)
						rhovec[k] = pSBAmp[k].rho;

					bubble_srtdouble(rhovec,sbNumElementAmp );

					// Find nb_rho and nb_phi leading minimum Df Lagragian
					double chosenR_rho, chosenR_phi;
					double chosenH_rho, chosenH_phi;
					int chosenQRhoIndex =0;
					int chosenQPhiIndex =0;
					double R_rho, R_phi;
					double H_rho, H_phi;
					int nQPhiStep=0;
					int nQRhoStep=0;
					double LagrangianDf=0.0;
					double minLagrangianDf=1e100;
					for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
					{
						if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
						nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

						H_rho = entropyRhoTable[i_qrho];

						for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
						{
							if ( (max_phase[i][j]/q_phi[i_qphi]) < 1.0) break;
							nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;

							H_phi = entropyPhiTable[i][j][i_qphi];

							LagrangianDf =  ( ampRangeNumElement[iAmpRange] *
											  (ampBar[iAmpRange]*ampBar[iAmpRange]) *
								 			  DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] ) +
											rdByAmpQuant->lambda *
											  ampRangeNumElement[iAmpRange] * (H_rho + H_phi);

 							if (LagrangianDf <= minLagrangianDf)
 							{
 								//cout << "minLagrangianDf: " << minLagrangianDf << endl;
 								//cout << "LagrangianDf: " << LagrangianDf << endl;
 								minLagrangianDf = LagrangianDf;
  								chosenQRhoIndex = i_qrho;
  								chosenQPhiIndex = i_qphi;
  								chosenH_rho = H_rho;
  								chosenH_phi = H_phi;
 							}

						} // nb_phase
					}  // nb_rho

					Lagrangian += minLagrangianDf;

// 					cout << "chosenQRhoIndex: " << chosenQRhoIndex << endl;
// 					cout << "chosenQPhiIndex: " << chosenQPhiIndex << endl;
// 					cout << "minLagrangianDf: " << minLagrangianDf << endl;

					nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[chosenQRhoIndex]))*2 + 1;
#ifdef USEGGD

					R_rho = computeRateSBExpDecayGGDPdf( 	rhovec,
															sbNumElementAmp,
															rhoQStepEdge[chosenQRhoIndex],
															rhoQStepProb[chosenQRhoIndex],
															nQRhoStep);
#else
					R_rho = computeRateSBExpDecayUnifPdf(nQRhoStep,
														 ampRangeNumElement[iAmpRange]);
#endif
					chosenR_rho = R_rho;
					///////
					nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[chosenQPhiIndex])) + 1;

					R_phi = computeRateSBExpPhaseUnifPdf(	nQPhiStep,
															ampRangeNumElement[iAmpRange]);
					chosenR_phi = R_phi;

					cout << setw (5) << setfill(' ') << i+1 << " "
						 << setw (4) << setfill(' ') << j+1 << " "
						 << setw (7) << setfill(' ') << nb_amp << " "
						 << setw (7) << setfill(' ') << iAmpRange+1 << " "
						 << setw (10) << setfill(' ') << ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << nb_xi*ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << nb_sample*ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << (R_deltasup/sbNumElementQ) * ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << H_deltasup << " "
						 << setw (10) << setfill(' ') << (R_amp/sbNumElementQ) * ampRangeNumElement[iAmpRange] << " "
						 << setw (14) << setfill(' ') << H_amp << " "
						 << setw (10) << setfill(' ') << chosenQRhoIndex << " "
						 << setw (10) << setfill(' ') << chosenR_rho << " "
						 << setw (10) << setfill(' ') << chosenH_rho << " "
						 << setw (10) << setfill(' ') << chosenQPhiIndex << " "
						 << setw (10) << setfill(' ') << chosenR_phi << " "
						 << setw (10) << setfill(' ') << chosenH_phi << " " << endl;

					rateAmpRange[iAmpRange] =  (R_amp/sbNumElementQ) * ampRangeNumElement[iAmpRange];
					rateAmpRange[iAmpRange] += computeRateSBExpFreq(static_cast<double>(nb_xi),ampRangeNumElement[iAmpRange]);
					rateAmpRange[iAmpRange] += computeRateSBExpSample(static_cast<double>(nb_sample),ampRangeNumElement[iAmpRange]);
					rateAmpRange[iAmpRange] += (R_deltasup/sbNumElementQ) * ampRangeNumElement[iAmpRange];
					rateAmpRange[iAmpRange] += chosenR_rho;
					rateAmpRange[iAmpRange] += chosenR_phi;

					quantExp[iAmpRange].nb_amp = nb_amp;
					quantExp[iAmpRange].nb_rho = 0;
					quantExp[iAmpRange].nb_phase = 0;
					quantExp[iAmpRange].nb_xi = 0;
					quantExp[iAmpRange].nb_sample = 0;
					quantExp[iAmpRange].i_qrho = chosenQRhoIndex;
					quantExp[iAmpRange].qrho = q_rho[chosenQRhoIndex];
					quantExp[iAmpRange].R_rho = chosenR_rho;
					quantExp[iAmpRange].H_rho = chosenH_rho;
					quantExp[iAmpRange].i_qphi = chosenQPhiIndex;
					quantExp[iAmpRange].qphi = q_phi[chosenQPhiIndex];
					quantExp[iAmpRange].R_phi = chosenR_phi;
					quantExp[iAmpRange].H_phi = chosenH_phi;
					quantExp[iAmpRange].R_xi = nb_xi;
					quantExp[iAmpRange].R_sample = nb_sample;
					quantExp[iAmpRange].R_amp = R_amp;
					quantExp[iAmpRange].H_amp = H_amp;
					quantExp[iAmpRange].R_deltasup = R_deltasup;
					quantExp[iAmpRange].H_deltasup = H_deltasup;

					delete [] rhovec;

				} // iAmpRange

				if (Lagrangian < minLagrangian)
				{
					minLagrangian = Lagrangian;
					chosen_nb_amp = static_cast<double>(nb_amp);
					chosenR_amp = R_amp;
					chosenH_amp = H_amp;
					chosenR_deltasup = R_deltasup;
					chosenH_deltasup = H_deltasup;
					chosenRate = 0.0;
					rdByAmpQuant->numAmpRange[i][j] = static_cast<int>(chosen_nb_amp);
					for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
					{
						chosenRate+= rateAmpRange[iAmpRange];
						rdByAmpQuant->rate[i][j][iAmpRange] = rateAmpRange[iAmpRange];
						rdByAmpQuant->quantExp[i][j][iAmpRange] = quantExp[iAmpRange];
						rdByAmpQuant->nAtom[i][j][iAmpRange] = ampRangeNumElement[iAmpRange];
					}
				}

				delete [] rateAmpRange;
				delete [] quantExp;
				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;
				delete [] ampRangeNumElement;
			} // nb_amp

			rdByAmpQuant->totalRate += static_cast<double>(chosenRate);

			delete [] pSBQ;
			//////////////////////
#ifdef USEGGD
			for (i_qrho=0; i_qrho<N_qrho; i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;

				if (rhoQStepCenter[i_qrho]!=NULL)
				{
					delete [] rhoQStepCenter[i_qrho];
					rhoQStepCenter[i_qrho] =NULL;
				}
				if (rhoQStepEdge[i_qrho]!=NULL)
				{
					delete [] rhoQStepEdge[i_qrho];
					rhoQStepEdge[i_qrho]=NULL;
				}
				if ( rhoQStepProb[i_qrho]!=NULL)
				{
					delete [] rhoQStepProb[i_qrho];
					rhoQStepProb[i_qrho]=NULL;
				}
				if (rhoQStepCumProb[i_qrho]!=NULL)
				{
					delete [] rhoQStepCumProb[i_qrho];
					rhoQStepCumProb[i_qrho]=NULL;
				}
			}
#endif

		} // block
	} // signal

	delete [] entropyRhoTable;
#ifdef USEGGD
	if (rhoQStepCenter!=NULL) delete [] rhoQStepCenter;
	if (rhoQStepEdge!=NULL) delete [] rhoQStepEdge;
	if (rhoQStepProb!=NULL) delete [] rhoQStepProb;
	if (rhoQStepCumProb!=NULL) delete [] rhoQStepCumProb;
#endif

	cout << "######################################################" << endl;
	cout << "# The chosen parameters                               " << endl;
	cout << "######################################################" << endl;
	cout << setw (5) << setfill(' ') << " sig" << " "
		 << setw (4) << setfill(' ') << " bl" << " "
		 << setw (7) << setfill(' ') << " nbamp" << " "
		 << setw (7) << setfill(' ') << " iAmpR"  << " "
		 << setw (10) << setfill(' ') << " NumEl" << " "
		 << setw (10) << setfill(' ') << " R_xi"  << " "
		 << setw (10) << setfill(' ') << " R_sample"  << " "
		 << setw (10) << setfill(' ') << " R_deltasup"  << " "
		 << setw (10) << setfill(' ') << " H_deltasup"  << " "
		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (14) << setfill(' ') << " H_amp" << " "
		 << setw (10) << setfill(' ') << " iQRho" << " "
		 << setw (10) << setfill(' ') << " R_rho"  << " "
		 << setw (10) << setfill(' ') << " H_rho" << " "
		 << setw (10) << setfill(' ') << " iQPhi" << " "
		 << setw (10) << setfill(' ') << " R_phi" << " "
		 << setw (10) << setfill(' ') << " H_phi" << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			double totalNAtom =0;
			for (iAmpRange=0;iAmpRange<rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			{
				totalNAtom += static_cast<double>(rdByAmpQuant->nAtom[i][j][iAmpRange]);
			}
			for (iAmpRange=0;iAmpRange<rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			{
				cout << setw (5) << setfill(' ') << i+1 << " "
					 << setw (4) << setfill(' ') << j+1 << " "
					 << setw (7) << setfill(' ') << rdByAmpQuant->numAmpRange[i][j] << " "
					 << setw (7) << setfill(' ') << iAmpRange+1 << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].R_xi*rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].R_sample*rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
					 << setw (10) << setfill(' ') << (rdByAmpQuant->quantExp[i][j][iAmpRange].R_deltasup/totalNAtom) *rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].H_deltasup<< " "
					 << setw (10) << setfill(' ') << (rdByAmpQuant->quantExp[i][j][iAmpRange].R_amp/totalNAtom) *rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
					 << setw (14) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].H_amp << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].i_qrho << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].R_rho << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].H_rho << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].i_qphi << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].R_phi << " "
					 << setw (10) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].H_phi << " " << endl;
			 }
		}
	}



#ifdef DEBUG
	file_rhoqstep.close();
#endif
}


void encodeSBExpAmpRangeBisecSearchArithCod(char*InputFile,
											double rateTarget,
											strtSBBHeader sbbHeader,
											CStructBook** structBook,
											CFileDfTable* DfTableFile,
											CFileDictionary* dicData,
											CDataSignal* dataSignal,
											int flLambdaBisecSearch,
											double initLambda)
{

	int Nfreq = dicData->getNumFreq();
	int fdiscrtype = dicData->getFDiscrType(0);
	double freqi = dicData->getFreqi(0);
	double step_xi = (2*pi/sbbHeader.Fs)*freqi;

    // Setting fixed quantizers
    double nb_xi = log2((double)Nfreq);
    double nb_sample = log2((double)sbbHeader.blockSize);
    int nb_rhoSign = 1;

	////////////////
	int i, j, k, i_qrho, i_qphi;
	int iAmpRange;
	int t;
	////////////////

    double** min_amp= new double*[sbbHeader.numSignal];
    double** max_amp= new double*[sbbHeader.numSignal];
    double** min_rho= new double*[sbbHeader.numSignal];
    double** max_rho= new double*[sbbHeader.numSignal];
    double** min_phase= new double*[sbbHeader.numSignal];
    double** max_phase= new double*[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		min_amp[i]= new double[sbbHeader.numBlock];
    	max_amp[i]= new double[sbbHeader.numBlock];
    	min_rho[i]= new double[sbbHeader.numBlock];
    	max_rho[i]= new double[sbbHeader.numBlock];
    	min_phase[i]= new double[sbbHeader.numBlock];
    	max_phase[i]= new double[sbbHeader.numBlock];
	}

	///////////======
	strtContinuousExp* pSB;
    int sbNumElement;

    double**** DfTable = DfTableFile->getDfTable();
    int N_qrho = DfTableFile->getNumQRho();
    int N_qphi = DfTableFile->getNumQPhi();
    double* q_rho = DfTableFile->getQRho();
    double* q_phi = DfTableFile->getQPhi();

    int nb_amp;
    double nb_rho, nb_phi;
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Find the range of the parameters for each Signal/Block
	double max_max_amp=0.0;
	double min_max_amp=1e30;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			((CStructBookExp*)structBook[i])[j].findParmRange(min_amp[i][j],
																  max_amp[i][j],
																  min_rho[i][j],
																  max_rho[i][j],
																  min_phase[i][j],
																  max_phase[i][j]);
			if (max_max_amp<max_amp[i][j]) max_max_amp = max_amp[i][j];
			if (min_max_amp>max_amp[i][j]) min_max_amp = max_amp[i][j];
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Compute the Amplitude Square Error for each Signal/Block/nb_amp
	double*** distAmp;
	distAmp = new double**[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		distAmp[i] = new double*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			distAmp[i][j] = new double[16];
			////
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
			////
			for (nb_amp=1;nb_amp<16;nb_amp++)
			{
				distAmp[i][j][nb_amp-1] = sumSqrErrorStructBookExpAmp(  min_amp[i][j], max_amp[i][j], nb_amp,
																		pSB, sbNumElement);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
	cout << "- Compute probabilities with respect to AMPLITUDE step quantization using" << endl;
	cout << "  Generalized Gaussian Distribution" << endl;
#ifdef DEBUG
	fstream file_ampqstep("rdbyamp_ampQStep.out",ios::out);
	file_ampqstep << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "nb_amp"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" " << endl;
#endif
	double**** ampQStepCenter = new double***[sbbHeader.numSignal];
	double**** ampQStepEdge = new double***[sbbHeader.numSignal];
	double**** ampQStepProb = new double***[sbbHeader.numSignal];
	double**** ampQStepCumProb = new double***[sbbHeader.numSignal];
	int nQAmpStep;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		ampQStepCenter[i] = new double**[sbbHeader.numBlock];
		ampQStepEdge[i] = new double**[sbbHeader.numBlock];
		ampQStepProb[i] = new double**[sbbHeader.numBlock];
	    ampQStepCumProb[i] = new double**[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			ampQStepCenter[i][j] = new double*[16];
			ampQStepEdge[i][j] = new double*[16];
			ampQStepProb[i][j] = new double*[16];
			ampQStepCumProb[i][j] = new double*[16];
			for (nb_amp=16; nb_amp>0; nb_amp--)
			{
				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp), min_amp[i][j], max_amp[i][j]);

				nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

				cout << "-nb_amp: " << nb_amp << "; step_amp: " << step_amp << "; nQAmpStep: " << nQAmpStep << endl;

				// normalize step_amp
				step_amp = step_amp/max_amp[i][j];
				cout << "step_amp normalized: " << step_amp << endl;
				ampQStepCenter[i][j][nb_amp-1] = new double[nQAmpStep];
				ampQStepEdge[i][j][nb_amp-1] = new double[nQAmpStep+1];
				ampQStepProb[i][j][nb_amp-1] = new double[nQAmpStep];
				ampQStepCumProb[i][j][nb_amp-1] = new double[nQAmpStep];

				setQStepFeatureGGD( ampQStepCenter[i][j][nb_amp-1],
									ampQStepEdge[i][j][nb_amp-1],
									ampQStepProb[i][j][nb_amp-1],
									ampQStepCumProb[i][j][nb_amp-1],
									step_amp,
									nQAmpStep,
									0,
									AMPGGDMEAN, AMPGGDSCALE, AMPGGDSHAPE);
									//DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);

#ifdef DEBUG
				file_ampqstep<< setw (20) << setfill(' ') << i+1<<" "
							 << setw (20) << setfill(' ') << j+1<<" "
							 << setw (20) << setfill(' ') << nb_amp<<" "
							 << setw (20) << setfill(' ') << ampQStepCenter[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepCenter[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepCenter[i][j][nb_amp-1][nQAmpStep-1]<< " "
							 << setw (20) << setfill(' ') << ampQStepEdge[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepEdge[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepEdge[i][j][nb_amp-1][nQAmpStep]<< " "
							 << setw (20) << setfill(' ') << ampQStepProb[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepProb[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepProb[i][j][nb_amp-1][nQAmpStep-1]<< " "
							 << setw (20) << setfill(' ') << ampQStepCumProb[i][j][nb_amp-1][0]<< " "
							 << setw (20) << setfill(' ') << ampQStepCumProb[i][j][nb_amp-1][static_cast<int>(nQAmpStep/2)]<< " "
							 << setw (20) << setfill(' ') << ampQStepCumProb[i][j][nb_amp-1][nQAmpStep-1]<< " " << endl;
#endif

			}
		}
	}
#ifdef DEBUG
	file_ampqstep.close();
#endif

	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
	cout << "- Compute probabilities with respect to DECAY step quantization using" << endl;
	cout << "  Generalized Gaussian Distribution" << endl;
#ifdef DEBUG
	fstream file_rhoqstep("rdbyamp_rhoQStep.out",ios::out);
	file_rhoqstep << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "i_qrho"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "center"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "edge"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "prob"<<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" "
	             << setw (20) << setfill(' ') << "cumProb" <<" " << endl;
#endif
    int nQRhoStep;
    int nQPhiStep;
	double**** rhoQStepCenter = new double***[sbbHeader.numSignal];
	double**** rhoQStepEdge = new double***[sbbHeader.numSignal];
	double**** rhoQStepProb = new double***[sbbHeader.numSignal];
	double**** rhoQStepCumProb = new double***[sbbHeader.numSignal];

	//cout << N_qrho << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		rhoQStepCenter[i] = new double**[sbbHeader.numBlock];
		rhoQStepEdge[i] = new double**[sbbHeader.numBlock];
		rhoQStepProb[i] = new double**[sbbHeader.numBlock];
		rhoQStepCumProb[i] = new double**[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			rhoQStepCenter[i][j] = new double*[N_qrho];
			rhoQStepEdge[i][j] = new double*[N_qrho];
			rhoQStepProb[i][j] = new double*[N_qrho];
			rhoQStepCumProb[i][j] = new double*[N_qrho];
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;

				nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

				//cout << "- nQRhoStep: " << nQRhoStep << endl;

				rhoQStepCenter[i][j][i_qrho] = new double[nQRhoStep];
				rhoQStepEdge[i][j][i_qrho] = new double[nQRhoStep+1];
				rhoQStepProb[i][j][i_qrho] = new double[nQRhoStep];
				rhoQStepCumProb[i][j][i_qrho] = new double[nQRhoStep];

				setQStepFeatureGGD( rhoQStepCenter[i][j][i_qrho],
									rhoQStepEdge[i][j][i_qrho],
									rhoQStepProb[i][j][i_qrho],
									rhoQStepCumProb[i][j][i_qrho],
									q_rho[i_qrho],
									nQRhoStep,
									static_cast<int>(nQRhoStep/2),
									DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);

#ifdef DEBUG
				file_rhoqstep<< setw (20) << setfill(' ') << i+1<<" "
							 << setw (20) << setfill(' ') << j+1<<" "
							 << setw (20) << setfill(' ') << i_qrho+1<<" "
							 << setw (20) << setfill(' ') << rhoQStepCenter[i][j][i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCenter[i][j][i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCenter[i][j][i_qrho][nQRhoStep-1]<< " "
							 << setw (20) << setfill(' ') << rhoQStepEdge[i][j][i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepEdge[i][j][i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepEdge[i][j][i_qrho][nQRhoStep]<< " "
							 << setw (20) << setfill(' ') << rhoQStepProb[i][j][i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepProb[i][j][i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepProb[i][j][i_qrho][nQRhoStep-1]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCumProb[i][j][i_qrho][0]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCumProb[i][j][i_qrho][static_cast<int>(nQRhoStep/2)]<< " "
							 << setw (20) << setfill(' ') << rhoQStepCumProb[i][j][i_qrho][nQRhoStep-1]<< " " << endl;
#endif
			}
		}
	}
#ifdef DEBUG
	file_rhoqstep.close();
#endif
	// Compute Entropy
	cout << "- Compute entropy table for Amp, Decay and Phase" << endl;
	double*** entropyAmpTable = new double**[sbbHeader.numSignal];
	double*** entropyRhoTable = new double**[sbbHeader.numSignal];
	double*** entropyPhiTable = new double**[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		entropyAmpTable[i] = new double*[sbbHeader.numBlock];
		entropyRhoTable[i] = new double*[sbbHeader.numBlock];
		entropyPhiTable[i] = new double*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			entropyAmpTable[i][j] = new double[16];
			entropyRhoTable[i][j] = new double[N_qrho];
			entropyPhiTable[i][j] = new double[N_qphi];
			for (nb_amp=0;nb_amp<16;nb_amp++)
			{
				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp+1), min_amp[i][j], max_amp[i][j]);
				nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

#ifdef USEGGD
				entropyAmpTable[i][j][nb_amp] = computeEntropyPdf(	ampQStepProb[i][j][nb_amp],
																	nQAmpStep);
#else
				entropyAmpTable[i][j][nb_amp] = computeEntropySBExpDecayUnifPdf(nQAmpStep);
#endif
				cout << "nb_amp: " << nb_amp << "; H_amp: " <<  entropyAmpTable[i][j][nb_amp] << endl;
			}
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
				nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

#ifdef USEGGD
				entropyRhoTable[i][j][i_qrho] = computeEntropySBExpDecayGGDPdf(	rhoQStepProb[i][j][i_qrho],
																				nQRhoStep);
#else
				entropyRhoTable[i][j][i_qrho] = computeEntropySBExpDecayUnifPdf(nQRhoStep);
#endif
			}
			for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
			{
				if ( (max_phase[i][j]/q_rho[i_qphi]) < 1.0) break;
				nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;

				entropyPhiTable[i][j][i_qphi] = computeEntropySBExpPhaseUnifPdf(nQPhiStep);
			}
		}
	}


	////////////////////
	cout << "- Computing rate for lower bound lambda" << endl;
	double leftLambda = ((min_max_amp*min_max_amp)/12.0)* (2*log(2)*pow(2.0,-2*30));
	double leftRate;
	double** leftNbAmp = new double*[sbbHeader.numSignal];
	for (i=0;i< sbbHeader.numSignal;i++)
	{
		leftNbAmp[i] = new double[sbbHeader.numBlock];
	}

	if (flLambdaBisecSearch==1)
	{
		getRateNbAmpGivenLambdaArithCod(    leftLambda,
									sbbHeader,structBook,
									min_amp, max_amp,
									min_rho, max_rho,
									min_phase, max_phase,
									DfTable,
									N_qrho, N_qphi,
									q_rho, q_phi,
									nb_xi, nb_sample,
									distAmp,
									leftRate, // output
									leftNbAmp,// output
									ampQStepEdge,
									ampQStepProb,
									rhoQStepEdge,
									rhoQStepProb,
									entropyAmpTable,
									entropyRhoTable,
									entropyPhiTable);
	}

	///////////////////////////
	cout << "- Computing rate for upper bound lambda" << endl;
	double rightLambda = ((max_max_amp*max_max_amp)/12.0)* (2*log(2)*pow(2.0,-2*1));
	double rightRate;
	double** rightNbAmp = new double*[sbbHeader.numSignal];
	for (i=0;i< sbbHeader.numSignal;i++)
	{
		rightNbAmp[i] = new double[sbbHeader.numBlock];
	}

	if (flLambdaBisecSearch==1)
	{
		getRateNbAmpGivenLambdaArithCod(    rightLambda,
									sbbHeader,structBook,
									min_amp, max_amp,
									min_rho, max_rho,
									min_phase, max_phase,
									DfTable,
									N_qrho, N_qphi,
									q_rho, q_phi,
									nb_xi, nb_sample,
									distAmp,
									rightRate, // output
									rightNbAmp,// output
									ampQStepEdge,
									ampQStepProb,
									rhoQStepEdge,
									rhoQStepProb,
									entropyAmpTable,
									entropyRhoTable,
									entropyPhiTable);
	}

	double epsilon = leftLambda*0.5;
	double chosenLambda;
	double** chosenNbAmp = new double*[sbbHeader.numSignal];
	double** chosen_nb_amp_aux = new double*[sbbHeader.numSignal];
	for (i=0;i< sbbHeader.numSignal;i++)
	{
		chosenNbAmp[i] = new double[sbbHeader.numBlock];
		chosen_nb_amp_aux[i] = new double[sbbHeader.numBlock];
	}

	double totalRate, chosenRate;
	cout << "- leftLambda" << " " << "/ leftRate" << endl;
	cout << "  " << leftLambda << " / " << leftRate << endl;
	cout << "- rightLambda" << " " << "/ rightRate" << endl;
	cout << "  " << rightLambda << " / " << rightRate << endl;

	double desiredLambda;
	if (initLambda==0.0)
	{
		desiredLambda = (rightLambda + leftLambda)/2.0;
	}
	else
	{
		cout << "- The initial lambda is given by USER: " << initLambda << endl;
		desiredLambda = initLambda;
	}

	chosenLambda = initLambda;

	int countTol = 0;
	while(1)
	{
		if (flLambdaBisecSearch!=1)
		{
			cout << "- The lambda bisectionsearch is not used." << endl;
			cout << "- Get the number of bits allocated to amplitude." << endl;
			getRateNbAmpGivenLambdaArithCod(chosenLambda,
								sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
							    DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								distAmp,
								totalRate, // output
								chosenNbAmp,// output
								ampQStepEdge,
								ampQStepProb,
								rhoQStepEdge,
								rhoQStepProb,
								entropyAmpTable,
								entropyRhoTable,
								entropyPhiTable);
			break;
		}
		//////////////////////////////////////

		getRateNbAmpGivenLambdaArithCod(desiredLambda,
								sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
							    DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								distAmp,
								totalRate, // output
								chosen_nb_amp_aux,// output
								ampQStepEdge,
								ampQStepProb,
								rhoQStepEdge,
								rhoQStepProb,
								entropyAmpTable,
								entropyRhoTable,
								entropyPhiTable);

		cout << "- leftLambda" << " " << "/ leftRate" << endl;
		cout << leftLambda << " / " << leftRate << endl;

		cout << "- rightLambda" << " " << "/ rightRate" << endl;
		cout << rightLambda << " / " << rightRate << endl;

		cout << "* desiredLambda: " << desiredLambda << endl;
		cout << "* totalRate: " << totalRate << endl;
		cout << "* TARGETRate: " << rateTarget*sbbHeader.signalSize << endl;


		if ( (totalRate==rightRate) || (totalRate==leftRate))
			countTol++;
		else
			countTol=0;

		if ( totalRate <= rateTarget*sbbHeader.signalSize)
		{
			rightLambda = desiredLambda;
			rightRate = totalRate;

			for (i=0;i<sbbHeader.numSignal;i++)
			{
				memcpy(rightNbAmp[i],chosen_nb_amp_aux[i], sbbHeader.numBlock*sizeof(double));
			}
		}
		else
		{
			leftLambda = desiredLambda;
			leftRate = totalRate;
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				memcpy(leftNbAmp[i],chosen_nb_amp_aux[i], sbbHeader.numBlock*sizeof(double));
			}
		}
		//////////////////////////////////////////
		if ( (fabs(rightLambda - leftLambda) < epsilon) ||
			 (countTol==6) )
		{
			if ( fabs(leftRate-rateTarget*sbbHeader.signalSize) <
				 fabs(rightRate-rateTarget*sbbHeader.signalSize) )
			{
				chosenLambda = leftLambda;
				chosenRate = leftRate;
				for (i=0;i<sbbHeader.numSignal;i++)
				{
					memcpy(chosenNbAmp[i],leftNbAmp[i], sbbHeader.numBlock*sizeof(double));
				}
			}
			else
			{
				chosenLambda = rightLambda;
				chosenRate = rightRate;
				for (i=0;i<sbbHeader.numSignal;i++)
				{
					memcpy(chosenNbAmp[i],rightNbAmp[i], sbbHeader.numBlock*sizeof(double));
				}
			}
			/////
			break;
		}

		desiredLambda = (rightLambda + leftLambda)/2.0;
	}

	strtRDByAmpQuant* rdByAmpQuant = new strtRDByAmpQuant;

	rdByAmpQuant->rate = new double**[sbbHeader.numSignal];
	rdByAmpQuant->quantExp = new strtQuantExp**[sbbHeader.numSignal];
	rdByAmpQuant->nAtom = new int**[sbbHeader.numSignal];
	rdByAmpQuant->numAmpRange = new int*[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		rdByAmpQuant->rate[i] = new double*[sbbHeader.numBlock];
		rdByAmpQuant->quantExp[i] = new strtQuantExp*[sbbHeader.numBlock];
		rdByAmpQuant->numAmpRange[i] = new int[sbbHeader.numBlock];
		rdByAmpQuant->nAtom[i] = new int*[sbbHeader.numBlock];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			int NAmpRange = static_cast<int>(chosenNbAmp[i][j]);
			rdByAmpQuant->numAmpRange[i][j] = NAmpRange;
			rdByAmpQuant->rate[i][j] = new double[NAmpRange];
			rdByAmpQuant->quantExp[i][j] = new strtQuantExp[NAmpRange];
			rdByAmpQuant->nAtom[i][j] = new int[NAmpRange];
		}
	}
	cout << "- Proceed the remaining bit allocation for each amplitude range." << endl;
	rdByAmpQuant->lambda = chosenLambda;
	computeRDByAmpQuantArithCod(sbbHeader,structBook,
								min_amp, max_amp,
								min_rho, max_rho,
								min_phase, max_phase,
								DfTable,
								N_qrho, N_qphi,
								q_rho, q_phi,
								nb_xi, nb_sample,
								chosenNbAmp,
								rdByAmpQuant,
								ampQStepEdge,
								ampQStepProb,
								rhoQStepEdge,
								rhoQStepProb,
								entropyAmpTable,
								entropyRhoTable,
								entropyPhiTable);
	totalRate = rdByAmpQuant->totalRate;

	cout << "#######################" << endl;
	cout << "# Chosen lambda: " << chosenLambda << endl;
	cout << "# total rate: " << totalRate << endl;
	cout << "#######################" << endl;

    ////////////////////////////////////////////////////////
#ifdef DEBUG
    fstream file_rdbyamp("rdbyamp_quant_arithcod.out",ios::out);
	file_rdbyamp << setw (20) << setfill(' ') << "Lambda"<<" "
	             << setw (20) << setfill(' ') << "Total Rate"<<" "
	             << setw (20) << setfill(' ') << "Signal"<<" "
	             << setw (20) << setfill(' ') << "Block"<<" "
	             << setw (20) << setfill(' ') << "NAmpRange/nbamp"<<" "
	             << setw (20) << setfill(' ') << "NAtomsInRange"<<" "
	             << setw (20) << setfill(' ') << "iAmpRange"<<" "
	             << setw (20) << setfill(' ') << "Rate"<<" "
	             << setw (20) << setfill(' ') << "qrho"<<" "
	             << setw (20) << setfill(' ') << "qphi"<<" "
	             << setw (20) << setfill(' ') << "nb_xi"<<" "
	             << setw (20) << setfill(' ') << "nb_sample"<<" " << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (iAmpRange=0;iAmpRange < rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			{
				file_rdbyamp << setw (20) << setfill(' ') << rdByAmpQuant->lambda << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->totalRate << " "
							 << setw (20) << setfill(' ') << i+1 << " "
							 << setw (20) << setfill(' ') << j+1 << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->numAmpRange[i][j] << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->nAtom[i][j][iAmpRange] << " "
							 << setw (20) << setfill(' ') << iAmpRange+1 << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->rate[i][j][iAmpRange] << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].qrho << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].qphi << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_xi << " "
							 << setw (20) << setfill(' ') << rdByAmpQuant->quantExp[i][j][iAmpRange].nb_sample << " " << endl;
			}
		}
	}
	file_rdbyamp.close();
#endif
	/////////////////////////////////////////////////////
	// ENCODE WITH THE CHOSEN LAMBDA
	double** origSignal;
    origSignal = dataSignal->getSignal();

    double** recSignal;
	recSignal = new double*[sbbHeader.numSignal];

	double* recBlockSignal;
    recBlockSignal = new double[sbbHeader.blockSize];

    double** recSignalNoQuant;
	recSignalNoQuant = new double*[sbbHeader.numSignal];

    double* recBlockSignalNoQuant;
    recBlockSignalNoQuant = new double[sbbHeader.blockSize];

    double finalRate=0;
    double sqrerror_total=0.0;


    /////
#ifdef DEBUG
    fstream file_encIndex("rdbyamp_encodeIndex.out",ios::out);
	file_encIndex << setw (14) << setfill(' ') << "Signal"<<" "
	              << setw (14) << setfill(' ') << "Block"<<" "
	              << setw (14) << setfill(' ') << "iAmpRange"<<" "
	              << setw (14) << setfill(' ') << "Element"<<" "
	              << setw (14) << setfill(' ') << "ampIndex"<<" "
	              << setw (14) << setfill(' ') << "rhoIndex"<<" "
	              << setw (14) << setfill(' ') << "xiIndex"<<" "
	              << setw (14) << setfill(' ') << "phiIndex"<<" "
	              << setw (14) << setfill(' ') << "aIndex"<<" "
	              << setw (14) << setfill(' ') << "bIndex"<<" " << endl;
#endif

    FILE* code_file;
    //////// Save encoded file
    char fileName[_MAX_PATH];
    char* pos;
    char aux[_MAX_PATH];
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,".mpz");
    strcpy( &pos[0], aux);
    //
    code_file = fopen (fileName,"wb");


	int numBytesHeader=0;
	int numBytes=0;

	unsigned short int dummyushint;
	unsigned int dummyuint;
	int dummyint;
	float dummyfloat;
	double dummydouble;
    /// FILL HEADER //////////////////
    // numSignal
    dummyint =sbbHeader.numSignal;
    fwrite (&dummyint, 1 , sizeof(int) , code_file );
    numBytesHeader += sizeof(int);
    // numBlock
    dummyint =sbbHeader.numBlock;
    fwrite (&dummyint, 1 , sizeof(int) , code_file );
    numBytesHeader += sizeof(int);
    // signalSize
    dummyint = sbbHeader.signalSize;
    fwrite (&dummyint, 1 , sizeof(int) , code_file );
    numBytesHeader += sizeof(int);
    // blockSize
    dummyint =sbbHeader.blockSize;
    fwrite (&dummyint, 1 , sizeof(int) , code_file );
    numBytesHeader += sizeof(int);
    // Fs - sample rate
    dummyuint = static_cast<unsigned int>(sbbHeader.Fs);
    fwrite (&dummyuint, 1 , sizeof(unsigned int) , code_file );
    numBytesHeader += sizeof(unsigned int);
    // lambda
    dummyfloat = static_cast<float>(rdByAmpQuant->lambda);
    fwrite (&dummyfloat, 1 , sizeof(float) , code_file );
    numBytesHeader += sizeof(float);
    // signal norm
    float* signalNorm = new float[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
    {
    	signalNorm[i] = static_cast<float>(sbbHeader.norm[i]);
    	cout << "Signal norm - " << i+1 << " : " << signalNorm[i] << endl;
    }
    fwrite (signalNorm, sbbHeader.numSignal  , sizeof(float) , code_file );
    numBytesHeader += sizeof(float) * sbbHeader.numSignal;

    for (i=0;i<sbbHeader.numSignal;i++)
	{
    	for (j=0;j<sbbHeader.numBlock;j++)
		{
    		/// FILL HEADER //////////////////
			// nbamp/NAmpRAnge
			dummyuint = static_cast<unsigned int>(rdByAmpQuant->numAmpRange[i][j]);
			cout << "nb_amp: " << dummyuint << endl;
    		fwrite (&dummyuint, 1 , sizeof(unsigned int) , code_file );
    		numBytesHeader += sizeof(unsigned int);
			// max_amp
			dummyfloat = static_cast<float>(max_amp[i][j]);
			cout << "max_amp: " << dummyfloat << endl;
    		fwrite (&dummyfloat, 1 , sizeof(float) , code_file );
    		numBytesHeader += sizeof(float);
			// max_rho
			dummyfloat = static_cast<float>(max_rho[i][j]);
			cout << "max_rho: " << dummyfloat << endl;
    		fwrite (&dummyfloat, 1 , sizeof(float) , code_file );
    		numBytesHeader += sizeof(float);
			// OBS:
			//   min_amp = 0
			//   min_rho =0
			//   min_phi = 0
			//   max_phi = 2*pi
		}
	}

    for (i=0;i<sbbHeader.numSignal;i++)
	{
		recSignal[i] = new double[sbbHeader.signalSize];
		recSignalNoQuant[i] = new double[sbbHeader.signalSize];
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();

			double L_amp = max_amp[i][j] - min_amp[i][j];

			int NAmpRange = rdByAmpQuant->numAmpRange[i][j];

			CStructBook* structBookByAmp;
			structBookByAmp = new CStructBookExp[NAmpRange];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];

			double* ampRangeNumElement;
			ampRangeNumElement = new double[NAmpRange];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSB,
																		sbNumElement,
																		ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																		ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit
			}
			////
			for (iAmpRange=0;iAmpRange < rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			//for (iAmpRange=rdByAmpQuant->numAmpRange[i][j]-1;iAmpRange >= 0;iAmpRange--)
			{
				strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
				int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElementAmp];
				int numElementAmpQ;

				strtQuantExp chosenQuantExp = rdByAmpQuant->quantExp[i][j][iAmpRange];

				numElementAmpQ = computeNumElementAmpQ(	  min_amp[i][j],max_amp[i][j],
														  chosenQuantExp.nb_amp,
														  pSBAmp,numElementAmp);

				/// FILL HEADER //////////
				dummyuint = static_cast<unsigned int>(numElementAmpQ);
    			fwrite (&dummyuint, 1 , sizeof(unsigned int) , code_file );
    			numBytesHeader += sizeof(unsigned int);

				delete [] pSBAmpQ;
			}


			unsigned char* compressed_data=0;
			unsigned max_code_bytes = 0xFFFFFFE0U;
			////////////////////////////////////
			// Arithmetic Encoder Varaibles
			Arithmetic_Codec ace;
			ace.set_buffer(max_code_bytes, compressed_data);

			// START ENCODER
			ace.start_encoder();

			for (k=0;k<sbbHeader.blockSize;k++)
			{
				recBlockSignal[k] = 0.0;
			}

			for (iAmpRange=0;iAmpRange < rdByAmpQuant->numAmpRange[i][j];iAmpRange++)
			//for (iAmpRange=rdByAmpQuant->numAmpRange[i][j]-1;iAmpRange >= 0;iAmpRange--)
			{
				strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
				int numElementAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				strtContinuousExp* pSBAmpQ = new strtContinuousExp[numElementAmp];
				int numElementAmpQ;

				strtQuantExp chosenQuantExp = rdByAmpQuant->quantExp[i][j][iAmpRange];

				int* ampIndex = new int[numElementAmp];
				int* rhoIndex = new int[numElementAmp];
				int* xiIndex  = new int[numElementAmp];
				int* phiIndex = new int[numElementAmp];
				int* aIndex   = new int[numElementAmp];
				int* bIndex   = new int[numElementAmp];

				quantizeStructBookExpRecUnquantNoResetArithCod(	  min_amp[i][j],max_amp[i][j],
																  chosenQuantExp,
																  pSBAmp,numElementAmp,
																  pSBAmpQ,numElementAmpQ,
																  recBlockSignalNoQuant,
																  sbbHeader.norm[i],
																  sbbHeader.blockSize,
																  ampIndex, rhoIndex, xiIndex,
													 			  phiIndex, aIndex, bIndex,
													 			  fdiscrtype, step_xi);

				cout << "- chosenQuantExp.i_qrho: " << chosenQuantExp.i_qrho << endl;
				cout << "- chosenQuantExp.i_qphi: " << chosenQuantExp.i_qphi << endl;

				// AMP
				Static_Data_Model ampModel;
				double step_amp = computeMidRiseQuantStep( static_cast<double>(chosenQuantExp.nb_amp), min_amp[i][j], max_amp[i][j]);
				unsigned nlevel_amp = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;
				//unsigned nlevel_amp = static_cast<unsigned>(round(pow(2.0, static_cast<double>(chosenQuantExp.nb_amp) ) ));
				cout << "- nlevel_amp: " << nlevel_amp << endl;
#ifdef USEGGD
				ampModel.set_distribution(static_cast<unsigned>(nlevel_amp),
										  ampQStepProb[i][j][chosenQuantExp.nb_amp - 1]);
#else
				ampModel.set_distribution(static_cast<unsigned>(nlevel_amp));
#endif
				for (k=0;k<numElementAmpQ;k++)
				{
					//cout << " Passou: " << k << endl;
					//int ind = ampIndex[k];
					//cout << "ampIndex[k]: " << ampIndex[k] << "; prob: " << ampQStepProb[i][j][chosenQuantExp.nb_amp - 1][ind]<< endl;
					ace.encode(static_cast<unsigned>(ampIndex[k]),ampModel);
				}

				/// DECAY
				int nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/chosenQuantExp.qrho))*2 + 1;
				Static_Data_Model decayModel;
				cout << "- nQRhoStep: " << nQRhoStep << endl;
#ifdef USEGGD
 				decayModel.set_distribution(static_cast<unsigned>(nQRhoStep),
 											rhoQStepProb[i][j][chosenQuantExp.i_qrho]);
#else
				decayModel.set_distribution(static_cast<unsigned>(nQRhoStep));
#endif
				for (k=0;k<numElementAmpQ;k++)
				{

					//cout << "- rhoIndex: " <<  rhoIndex[k] << ";shift: " << static_cast<int>(nQRhoStep/2) <<  endl;

					ace.encode(static_cast<unsigned>(rhoIndex[k] + static_cast<int>(nQRhoStep/2)),decayModel);
				}

				/// PHASE
				int nQPhiStep = static_cast<int>(ceil(max_phase[i][j]/chosenQuantExp.qphi)) + 1;
				Static_Data_Model phiModel;
				cout << "- nQPhiStep: " << nQPhiStep << endl;
				phiModel.set_distribution(static_cast<unsigned>(nQPhiStep));
				for (k=0;k<numElementAmpQ;k++)
				{
					//cout << "- phiIndex: " <<  phiIndex[k] <<  endl;
					ace.encode(static_cast<unsigned>(phiIndex[k]),phiModel);
				}

				/// FREQUENCY
				Static_Data_Model xiModel;
				cout << "- Nfreq: " << Nfreq << endl;
				xiModel.set_distribution(static_cast<unsigned>(Nfreq));
				for (k=0;k<numElementAmpQ;k++)
				{
					ace.encode(static_cast<unsigned>(xiIndex[k]),xiModel);
				}

				/// TIME SUPPORT
				/// A
				Static_Data_Model timeSupModel;
				cout << "- sbbHeader.blockSize: " << sbbHeader.blockSize << endl;
				timeSupModel.set_distribution(static_cast<unsigned>(sbbHeader.blockSize));
				for (k=0;k<numElementAmpQ;k++)
				{
					ace.encode(static_cast<unsigned>(aIndex[k]),timeSupModel);
				}
				/// B
				for (k=0;k<numElementAmpQ;k++)
				{
					ace.encode(static_cast<unsigned>(bIndex[k]),timeSupModel);
				}

#ifdef DEBUG
				for (k=0;k<numElementAmpQ;k++)
				{
					file_encIndex << setw (14) << setfill(' ') << i+1<<" "
								  << setw (14) << setfill(' ') << j+1<<" "
								  << setw (14) << setfill(' ') << iAmpRange+1<<" "
								  << setw (14) << setfill(' ') << k+1<<" "
								  << setw (14) << setfill(' ') << ampIndex[k] <<" "
								  << setw (14) << setfill(' ') << rhoIndex[k] <<" "
								  << setw (14) << setfill(' ') << xiIndex[k]  <<" "
								  << setw (14) << setfill(' ') << phiIndex[k] <<" "
								  << setw (14) << setfill(' ') << aIndex[k]   <<" "
								  << setw (14) << setfill(' ') << bIndex[k]   <<" " << endl;
				}
#endif
				printf(" ****** Synthesize signal %d block %d amp_range %d **** \n",i+1,j+1,iAmpRange+1);
				synthSignalSBExpNoReset(recBlockSignal,
										sbbHeader.norm[i],
										sbbHeader.blockSize,
										pSBAmpQ,
										numElementAmpQ);

				delete [] ampIndex;
				delete [] rhoIndex;
				delete [] xiIndex;
				delete [] phiIndex;
				delete [] aIndex;
				delete [] bIndex;

				delete [] pSBAmpQ;
			}

			unsigned nbyte_amp = ace.write_to_file(code_file); // Alternalively => ace.stop_encoder();
			numBytes += static_cast<int>(nbyte_amp);

			int initBlockSample = j*sbbHeader.blockSize;
			if ((initBlockSample + sbbHeader.blockSize) >= sbbHeader.signalSize)
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
				memcpy(	&recSignalNoQuant[i][initBlockSample],
						recBlockSignalNoQuant,
						sizeof(double)*(sbbHeader.signalSize-initBlockSample));
			}
			else
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*sbbHeader.blockSize);
				memcpy(	&recSignalNoQuant[i][initBlockSample],
						recBlockSignalNoQuant,
						sizeof(double)*sbbHeader.blockSize);
			}

			delete [] (CStructBookExp*)structBookByAmp;
			delete [] ampBar;
			delete [] ampRangeLimit;
			delete [] ampRangeNumElement;
		} // block
		sqrerror_total += computeSqrError(origSignal[i],recSignal[i],0,sbbHeader.signalSize-1);
	} //signal

	fclose(code_file);
#ifdef DEBUG
	file_encIndex.close();
#endif
	finalRate = totalRate;
	double bitrate = static_cast<double>(finalRate)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Final rate calculated: %20.10f \n",finalRate);
    printf("Signal size: %d \n",sbbHeader.signalSize);
    printf("Bitrate calculated: %14.10f bit/sample\n",bitrate);
    printf("Total Square Error: %14.10f \n",sqrerror_total);
    /////////
    printf("Final rate from arithmetic codec: %20.10f bits \n",static_cast<double>(numBytes*8));
    printf("Header rate from arithmetic codec: %20.10f bits\n",static_cast<double>(numBytesHeader*8));
    double bitrateAC = static_cast<double>(numBytes*8)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
	double bitrateACHeader = static_cast<double>(numBytesHeader*8)/
                     static_cast<double>(static_cast<double>(sbbHeader.signalSize)*
                              static_cast<double>(sbbHeader.numSignal) );
    printf("Bitrate from arithmetic codec: %14.10f bit/sample\n",bitrateAC);
    printf("Bitrate header from arithmetic codec: %14.10f bit/sample\n",bitrateACHeader);

    fstream file_enc_report("rdbyamp_encode_report.out",ios::out);
	file_enc_report << setw (25) << setfill(' ') << "Lambda"<<" "
	              << setw (25) << setfill(' ') << "Bitrate Calc. (b/s)"<<" "
	              << setw (30) << setfill(' ') << "Bitrate Arith. Cod. (b/s)"<<" "
	              << setw (25) << setfill(' ') << "Header bitrate (b/s)"<<" "
	              << setw (25) << setfill(' ') << "Distortion (Sqr. Error)"<< " " << endl;

	file_enc_report << setw (25) << setfill(' ') << setprecision(10) << chosenLambda <<" "
	              << setw (30) << setfill(' ') << setprecision(10) << bitrate<<" "
	              << setw (25) << setfill(' ') << setprecision(10) << bitrateAC<<" "
	              << setw (25) << setfill(' ') << setprecision(10) << bitrateACHeader <<" "
	              << setw (25) << setfill(' ') << setprecision(10) << sqrerror_total << " " << endl;
	file_enc_report.close();

    //////// Save to file the reconstructed signals
    CDataSignal* outSignal;
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_rec%09.6fbps.wav",bitrateAC);
    strcpy( &pos[0], aux);
    //
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSignal(recSignal);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;

    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_rec%09.6fbps_noquant.wav",bitrateAC);
    strcpy( &pos[0], aux);
    //
    outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(sbbHeader.numSignal);
    outSignal->setSignalSize(sbbHeader.signalSize);
    outSignal->setSignal(recSignalNoQuant);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete outSignal;



    //////////////////////
    /// DEALLOCATE ///////
    for (i=0;i<sbbHeader.numSignal;i++)
	{
    	for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (nb_amp=0;nb_amp<16;nb_amp++)
			{
				delete [] ampQStepCenter[i][j][nb_amp];
				delete [] ampQStepEdge[i][j][nb_amp];
				delete [] ampQStepProb[i][j][nb_amp];
				delete [] ampQStepCumProb[i][j][nb_amp];
			}
			delete [] ampQStepCenter[i][j];
			delete [] ampQStepEdge[i][j];
			delete [] ampQStepProb[i][j];
			delete [] ampQStepCumProb[i][j];
		}
		delete [] ampQStepCenter[i];
		delete [] ampQStepEdge[i];
		delete [] ampQStepProb[i];
		delete [] ampQStepCumProb[i];
	}
	delete [] ampQStepCenter;
	delete [] ampQStepEdge;
	delete [] ampQStepProb;
	delete [] ampQStepCumProb;
    ///////////////////////
    for (i=0;i<sbbHeader.numSignal;i++)
	{
    	for (j=0;j<sbbHeader.numBlock;j++)
		{
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
				delete [] rhoQStepCenter[i][j][i_qrho];
				delete [] rhoQStepEdge[i][j][i_qrho];
				delete [] rhoQStepProb[i][j][i_qrho];
				delete [] rhoQStepCumProb[i][j][i_qrho];
			}
			delete [] rhoQStepCenter[i][j];
			delete [] rhoQStepEdge[i][j];
			delete [] rhoQStepProb[i][j];
			delete [] rhoQStepCumProb[i][j];
		}
		delete [] rhoQStepCenter[i];
		delete [] rhoQStepEdge[i];
		delete [] rhoQStepProb[i];
		delete [] rhoQStepCumProb[i];
	}
	delete [] rhoQStepCenter;
	delete [] rhoQStepEdge;
	delete [] rhoQStepProb;
	delete [] rhoQStepCumProb;

	for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] min_amp[i];
    	delete [] max_amp[i];
    	delete [] min_rho[i];
    	delete [] max_rho[i];
    	delete [] min_phase[i];
    	delete [] max_phase[i];
	}
	delete [] min_amp;
	delete [] max_amp;
	delete [] min_rho;
	delete [] max_rho;
	delete [] min_phase;
	delete [] max_phase;

	//////////////////////////////////////////
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			delete [] rdByAmpQuant->rate[i][j];
			delete [] rdByAmpQuant->quantExp[i][j];
			delete [] rdByAmpQuant->nAtom[i][j];
			delete [] distAmp[i][j];
		}
		delete [] distAmp[i];
		delete [] rdByAmpQuant->rate[i];
		delete [] rdByAmpQuant->quantExp[i];
		delete [] rdByAmpQuant->numAmpRange[i];
		delete [] rdByAmpQuant->nAtom[i];
	}
	delete [] distAmp;
	delete [] rdByAmpQuant->rate;
	delete [] rdByAmpQuant->quantExp;
	delete [] rdByAmpQuant->numAmpRange;
	delete [] rdByAmpQuant->nAtom;
	delete rdByAmpQuant;
	/////////
	for (i=0; i<sbbHeader.numSignal; i++)
	{
		delete [] recSignal[i];
		delete [] recSignalNoQuant[i];
	}
	delete [] recSignal;
	delete [] recSignalNoQuant;

	delete [] recBlockSignal;
	delete [] recBlockSignalNoQuant;
	/////
	for (i=0;i< sbbHeader.numSignal;i++)
	{
		delete [] chosenNbAmp[i];
		delete [] leftNbAmp[i];
		delete [] rightNbAmp[i];
		delete [] chosen_nb_amp_aux[i];
	}
	delete [] chosenNbAmp;
	delete [] leftNbAmp;
	delete [] rightNbAmp;
	delete [] chosen_nb_amp_aux;
	////////////
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			delete [] entropyAmpTable[i][j];
			delete [] entropyRhoTable[i][j];
			delete [] entropyPhiTable[i][j];
		}
		delete [] entropyAmpTable[i];
		delete [] entropyRhoTable[i];
		delete [] entropyPhiTable[i];
	}
	delete [] entropyAmpTable;
	delete [] entropyRhoTable;
	delete [] entropyPhiTable;
}

void  getRateNbAmpGivenLambdaArithCod(  double lambda,
										strtSBBHeader sbbHeader,CStructBook** structBook,
										double** min_amp, double** max_amp,
										double** min_rho,double** max_rho,
										double** min_phase,double** max_phase,
										double**** DfTable,
										int N_qrho,int N_qphi,
										double* q_rho,double* q_phi,
										double nb_xi,double nb_sample,
										double*** distAmp,
										double& totalRate,
										double** chosen_nb_amp,
										double**** ampQStepEdge,
										double**** ampQStepProb,
										double**** rhoQStepEdge,
										double**** rhoQStepProb,
										double*** entropyAmpTable,
										double*** entropyRhoTable,
										double*** entropyPhiTable)
{
	int i,j, iAmpRange, i_qrho,i_qphi, k;

    int nb_amp;
    double nb_rho, nb_phi;

	totalRate = 0.0;
	double chosenRate;

	cout << setw (10) << setfill(' ') << " signal" << " "
		 << setw (10) << setfill(' ') << " block" << " "
		 << setw (10) << setfill(' ') << " nb_amp" << " "
		 << setw (10) << setfill(' ') << " iAmpRange"  << " "
		 << setw (10) << setfill(' ') << " NumEl" << " "
//		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (10) << setfill(' ') << " R_xi"  << " "
		 << setw (10) << setfill(' ') << " R_sample"  << " "
		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (14) << setfill(' ') << " H_amp" << " "
		 << setw (10) << setfill(' ') << " iQRho" << " "
		 << setw (10) << setfill(' ') << " R_rho"  << " "
		 << setw (10) << setfill(' ') << " H_rho" << " "
		 << setw (10) << setfill(' ') << " iQPhi" << " "
		 << setw (10) << setfill(' ') << " R_phi" << " "
		 << setw (10) << setfill(' ') << " H_phi" << endl;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{
			/////
			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();


			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];

			int sbNumElementQ;

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

// 			cout << "L_amp: " << L_amp << endl;
// 			cout << "L_rho: " << L_rho << endl;
// 			cout << "L_phi: " << L_phi << endl;

			double minLagrangian = 1e100;
			double LagrangianDAmp[16];
			double Lagrangian[16];
			double R_amp, chosenR_amp;
			double H_amp, chosenH_amp;
			for (nb_amp=1;nb_amp<=16;nb_amp++)
			{

				quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
											pSB, sbNumElement,
											pSBQ,sbNumElementQ);

				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp), min_amp[i][j], max_amp[i][j]);
				int nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

#ifdef USEGGD
				double* ampvec = new double[sbNumElementQ];

				// normalize amp by max_amp (max_amp=1)
				for (k=0;k<sbNumElementQ;k++)
					ampvec[k] = ( pSBQ[k].innerProduct / max_amp[i][j]);

				bubble_srtdouble(ampvec, sbNumElementQ);

				R_amp = computeRatePdf( ampvec,
										sbNumElementQ,
										ampQStepEdge[i][j][nb_amp-1],
										ampQStepProb[i][j][nb_amp-1],
										nQAmpStep);
				delete [] ampvec;
#else
				R_amp = computeRateUnifPdf(nQAmpStep, sbNumElementQ);
#endif
				H_amp = entropyAmpTable[i][j][nb_amp-1];


// 				LagrangianDAmp[nb_amp-1] = distAmp[i][j][nb_amp-1] +
// 										  ((lambda) *
// 										   static_cast<double>(nb_amp + nb_xi + 2*nb_sample) *
// 										   sbNumElementQ);
				LagrangianDAmp[nb_amp-1] = distAmp[i][j][nb_amp-1] +
										   ((lambda) *
										    static_cast<double>(H_amp + nb_xi + 2*nb_sample) *
										    sbNumElementQ );
			    Lagrangian[nb_amp-1] = LagrangianDAmp[nb_amp-1];


				int NAmpRange = nb_amp;


				CStructBook* structBookByAmp;
				structBookByAmp = new CStructBookExp[NAmpRange];

				double* ampRangeLimit;
				ampRangeLimit = new double[NAmpRange+1];

				double* ampBar;
				ampBar = new double[NAmpRange];

				double* ampRangeNumElement;
				ampRangeNumElement = new double[NAmpRange];

				// Separate atoms among amplitude ranges and return the number of atoms in each range
				ampRangeLimit[0]=-0.001;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					ampRangeLimit[iAmpRange+1] =
						L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
					ampBar[NAmpRange-1-iAmpRange] =
						(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

					((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
																			sbNumElementQ,
																			ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																			ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

					ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

				}

				// Find the minimum Df Lagrangian and return the nb_rho and nb_phi
				double rate=0;
				for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
				{
					///////
					strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
					int sbNumElementAmp = ampRangeNumElement[iAmpRange];
					double* rhovec = new double[sbNumElementAmp];

					for (k=0;k<sbNumElementAmp;k++)
						rhovec[k] = pSBAmp[k].rho;

					bubble_srtdouble(rhovec,sbNumElementAmp );

					// Find nb_rho and nb_phi leading minimum Df Lagragian
					double chosenR_rho, chosenR_phi;
					double chosenH_rho, chosenH_phi;
					int chosenQRhoIndex, chosenQPhiIndex;
					double R_rho, R_phi;
					double H_rho, H_phi;
					double LagrangianDf;
					double minLagrangianDf=1e100;
					int nQPhiStep;
					int nQRhoStep;
					for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
					{
						if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
						nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

						H_rho = entropyRhoTable[i][j][i_qrho];
						for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
						{
							if ( (max_phase[i][j]/q_phi[i_qphi]) < 1.0) break;
							nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;

							H_phi = entropyPhiTable[i][j][i_qphi];

							LagrangianDf =  ( ampRangeNumElement[iAmpRange] *
											  (ampBar[iAmpRange]*ampBar[iAmpRange]) *
								 			  DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] ) +
											lambda *
											  ampRangeNumElement[iAmpRange] * (H_rho + H_phi);

							if (LagrangianDf<= minLagrangianDf)
							{
								minLagrangianDf = LagrangianDf;
								chosenQRhoIndex = i_qrho;
								chosenQPhiIndex = i_qphi;
								chosenH_rho = H_rho;
								chosenH_phi = H_phi;
							}

						} // nb_phase
					}  // nb_rho

					nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[chosenQRhoIndex]))*2 + 1;
#ifdef USEGGD
					R_rho = computeRateSBExpDecayGGDPdf( 	rhovec,
															sbNumElementAmp,
															rhoQStepEdge[i][j][chosenQRhoIndex],
															rhoQStepProb[i][j][chosenQRhoIndex],
															nQRhoStep);
#else
					R_rho = computeRateSBExpDecayUnifPdf(nQRhoStep,
														 ampRangeNumElement[iAmpRange]);
#endif
					chosenR_rho = R_rho;
					///////
					nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[chosenQPhiIndex])) + 1;

					R_phi = computeRateSBExpPhaseUnifPdf(	nQPhiStep,
															ampRangeNumElement[iAmpRange]);
					chosenR_phi = R_phi;

					cout << setw (10) << setfill(' ') << i+1 << " "
						 << setw (10) << setfill(' ') << j+1 << " "
						 << setw (10) << setfill(' ') << nb_amp << " "
						 << setw (10) << setfill(' ') << iAmpRange+1 << " "
						 << setw (10) << setfill(' ') << ampRangeNumElement[iAmpRange] << " "
//						 << setw (10) << setfill(' ') << nb_amp*ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << nb_xi*ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << 2*nb_sample*ampRangeNumElement[iAmpRange] << " "
						 << setw (10) << setfill(' ') << (R_amp/sbNumElementQ) * ampRangeNumElement[iAmpRange] << " "
						 << setw (14) << setfill(' ') << H_amp << " "
						 << setw (10) << setfill(' ') << chosenQRhoIndex << " "
						 << setw (10) << setfill(' ') << chosenR_rho << " "
						 << setw (10) << setfill(' ') << chosenH_rho << " "
						 << setw (10) << setfill(' ') << chosenQPhiIndex << " "
						 << setw (10) << setfill(' ') << chosenR_phi << " "
						 << setw (10) << setfill(' ') << chosenH_phi << " " << endl;


					Lagrangian[nb_amp-1] += minLagrangianDf;

					//rate+= computeRateSBExpAmp(static_cast<double>(nb_amp),ampRangeNumElement[iAmpRange]);
					rate+= (R_amp/sbNumElementQ) * ampRangeNumElement[iAmpRange];
					rate+= computeRateSBExpFreq(static_cast<double>(nb_xi),ampRangeNumElement[iAmpRange]);
					rate+= computeRateSBExpSample(static_cast<double>(nb_sample),ampRangeNumElement[iAmpRange]);
					rate+= chosenR_rho;
					rate+= chosenR_phi;

					delete [] rhovec;

				} // iAmpRange

				if (Lagrangian[nb_amp-1] < minLagrangian)
				{
					minLagrangian = Lagrangian[nb_amp-1];
					chosen_nb_amp[i][j] = static_cast<double>(nb_amp);
					chosenR_amp = R_amp;
					chosenH_amp = H_amp;
					chosenRate = rate;
				}

				delete [] (CStructBookExp*)structBookByAmp;
				delete [] ampBar;
				delete [] ampRangeLimit;
				delete [] ampRangeNumElement;
			} // nb_amp

			totalRate += static_cast<double>(chosenRate);

			delete [] pSBQ;
		} // block
	} // signal
}

void computeRDByAmpQuantArithCod(   strtSBBHeader sbbHeader, CStructBook** structBook,
									double** min_amp, double** max_amp,
									double** min_rho, double** max_rho,
									double** min_phase, double** max_phase,
									double**** DfTable,
									int N_qrho, int N_qphi,
									double* q_rho, double* q_phi,
									double nb_xi, double nb_sample,
									double** chosenNbAmp,
									strtRDByAmpQuant* rdByAmpQuant,
									double**** ampQStepEdge,
									double**** ampQStepProb,
									double**** rhoQStepEdge,
									double**** rhoQStepProb,
									double*** entropyAmpTable,
									double*** entropyRhoTable,
									double*** entropyPhiTable)
{
	int i,j,k, iAmpRange, i_qrho,i_qphi;

    int nb_amp;
    double nb_rho, nb_phi;

	double totalRate = 0.0;


	cout << setw (10) << setfill(' ') << " signal" << " "
		 << setw (10) << setfill(' ') << " block" << " "
		 << setw (10) << setfill(' ') << " nb_amp" << " "
		 << setw (10) << setfill(' ') << " iAmpRange"  << " "
		 << setw (10) << setfill(' ') << " NumEl" << " "
//		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (10) << setfill(' ') << " R_xi"  << " "
		 << setw (10) << setfill(' ') << " R_sample"  << " "
		 << setw (10) << setfill(' ') << " R_amp"  << " "
		 << setw (14) << setfill(' ') << " H_amp" << " "
		 << setw (10) << setfill(' ') << " iQRho" << " "
		 << setw (10) << setfill(' ') << " R_rho"  << " "
		 << setw (10) << setfill(' ') << " H_rho" << " "
		 << setw (10) << setfill(' ') << " iQPhi" << " "
		 << setw (10) << setfill(' ') << " R_phi" << " "
		 << setw (10) << setfill(' ') << " H_phi" << endl;

	double R_amp;
	double H_amp;
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		for (j=0;j<sbbHeader.numBlock;j++)
		{

			strtContinuousExp* pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();


			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElement];

			int sbNumElementQ;

			double L_amp = max_amp[i][j] - min_amp[i][j];
			double L_rho = max_rho[i][j] - min_rho[i][j];
			double L_phi = max_phase[i][j] - min_phase[i][j];

// 			cout << "L_amp: " << L_amp << endl;
// 			cout << "L_rho: " << L_rho << endl;
// 			cout << "L_phi: " << L_phi << endl;

			nb_amp = static_cast<int>( chosenNbAmp[i][j] );


			quantizeStructBookExpAmp(  	min_amp[i][j], max_amp[i][j], nb_amp,
										pSB, sbNumElement,
										pSBQ,sbNumElementQ);

			double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp), min_amp[i][j], max_amp[i][j]);
			int nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

#ifdef USEGGD
			double* ampvec = new double[sbNumElementQ];

			// normalize amp by max_amp (max_amp=1)
			for (k=0;k<sbNumElementQ;k++)
				ampvec[k] = ( pSBQ[k].innerProduct / max_amp[i][j]);

			bubble_srtdouble(ampvec, sbNumElementQ);

			R_amp = computeRatePdf( ampvec,
									sbNumElementQ,
									ampQStepEdge[i][j][nb_amp-1],
									ampQStepProb[i][j][nb_amp-1],
									nQAmpStep);
			delete [] ampvec;
#else
			R_amp = computeRateUnifPdf(nQAmpStep, sbNumElementQ);
#endif
			H_amp = entropyAmpTable[i][j][nb_amp-1];

			int NAmpRange = nb_amp;

			CStructBook* structBookByAmp;
			structBookByAmp = new CStructBookExp[NAmpRange];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];

			double* ampRangeNumElement;
			ampRangeNumElement = new double[NAmpRange];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));

				((CStructBookExp*)structBookByAmp)[iAmpRange].sepByAmp(	pSBQ,
																		sbNumElementQ,
																		ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																		ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

				ampRangeNumElement[iAmpRange] = ((CStructBookExp*)structBookByAmp)[iAmpRange].getNumElement();

			}

			double rate=0.0;
			int k;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				strtContinuousExp* pSBAmp = ((CStructBookExp*)structBookByAmp)[iAmpRange].getStructBook();
				int sbNumElementAmp = ampRangeNumElement[iAmpRange];
				double* rhovec = new double[sbNumElementAmp];

				for (k=0;k<sbNumElementAmp;k++)
					rhovec[k] = pSBAmp[k].rho;

				bubble_srtdouble(rhovec,sbNumElementAmp );

				// Find nb_rho and nb_phi leading minimum Df Lagragian
				int nQPhiStep;
				int nQRhoStep;
				double chosenR_rho, chosenR_phi;
				double R_rho, R_phi;
				double chosenH_rho, chosenH_phi;
				double H_rho, H_phi;
				int chosenQRhoIndex, chosenQPhiIndex;
				double LagrangianDf;
				double minLagrangianDf=1e100;
				for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
				{
					if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
					nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

					int flag=0;
					// if ( (i_qrho==20) )
// 						flag = 1;
// 					else
// 						flag = 0;


// 					R_rho = computeRateSBExpDecayGGDPdfPrint( 	rhovec,
// 															sbNumElementAmp,
// 															rhoQStepEdge[i][j][i_qrho],
// 															rhoQStepProb[i][j][i_qrho],
// 															nQRhoStep,
// 															flag);

// 					if (flag==1)
// 					{
// 						cout << "R_rho: " << R_rho << endl;
// 						cout << "iAmpRange: " << iAmpRange+1 << endl;
// 						cout << "Pressione para continuar..." << endl;
// 						int lala;
// 						cin >> lala;
// 					}
					H_rho = entropyRhoTable[i][j][i_qrho];
					for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
					{
						if ( (max_phase[i][j]/q_phi[i_qphi]) < 1.0) break;
						nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;

					    H_phi = entropyPhiTable[i][j][i_qphi];

						LagrangianDf =  ( ampRangeNumElement[iAmpRange] *
										  (ampBar[iAmpRange]*ampBar[iAmpRange]) *
										  DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] ) +
										( rdByAmpQuant->lambda *
											ampRangeNumElement[iAmpRange] * (H_rho + H_phi) );

						if (LagrangianDf<= minLagrangianDf)
						{
							minLagrangianDf = LagrangianDf;
							chosenQRhoIndex = i_qrho;
							chosenQPhiIndex = i_qphi;
							chosenH_rho = H_rho;
							chosenH_phi = H_phi;
						}
					} // nb_phase
				}  // nb_rho

				nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[chosenQRhoIndex]))*2 + 1;
#ifdef USEGGD
				R_rho = computeRateSBExpDecayGGDPdf( 	rhovec,
														sbNumElementAmp,
														rhoQStepEdge[i][j][chosenQRhoIndex],
														rhoQStepProb[i][j][chosenQRhoIndex],
														nQRhoStep);
#else
				R_rho = computeRateSBExpDecayUnifPdf(nQRhoStep,
													 ampRangeNumElement[iAmpRange]);
#endif
				chosenR_rho = R_rho;
				///////
				nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[chosenQPhiIndex])) + 1;

				R_phi = computeRateSBExpPhaseUnifPdf(	nQPhiStep,
														ampRangeNumElement[iAmpRange]);
				chosenR_phi = R_phi;


				strtQuantExp quantExp;
				quantExp.nb_amp = nb_amp;
				quantExp.nb_rho = 0;
				quantExp.nb_phase = 0;
				quantExp.nb_xi = nb_xi;
				quantExp.nb_sample = nb_sample;
				quantExp.i_qrho = chosenQRhoIndex;
				quantExp.qrho = q_rho[chosenQRhoIndex];
				quantExp.i_qphi = chosenQPhiIndex;
				quantExp.qphi = q_phi[chosenQPhiIndex];



				cout << setw (10) << setfill(' ') << i+1 << " "
					 << setw (10) << setfill(' ') << j+1 << " "
					 << setw (10) << setfill(' ') << nb_amp << " "
					 << setw (10) << setfill(' ') << iAmpRange+1 << " "
					 << setw (10) << setfill(' ') << ampRangeNumElement[iAmpRange] << " "
//					 << setw (10) << setfill(' ') << nb_amp*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') << nb_xi*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') << 2*nb_sample*ampRangeNumElement[iAmpRange] << " "
					 << setw (10) << setfill(' ') << (R_amp/sbNumElementQ) * ampRangeNumElement[iAmpRange] << " "
					 << setw (14) << setfill(' ') << H_amp << " "
					 << setw (10) << setfill(' ') << chosenQRhoIndex << " "
					 << setw (10) << setfill(' ') << chosenR_rho << " "
					 << setw (10) << setfill(' ') << chosenH_rho << " "
					 << setw (10) << setfill(' ') << chosenQPhiIndex << " "
					 << setw (10) << setfill(' ') << chosenR_phi << " "
					 << setw (10) << setfill(' ') << chosenH_phi << " " << endl;

//				rdByAmpQuant->rate[i][j][iAmpRange] =  computeRateSBExpAmp(static_cast<double>(nb_amp),ampRangeNumElement[iAmpRange]);
				rdByAmpQuant->rate[i][j][iAmpRange] =  (R_amp/sbNumElementQ) * ampRangeNumElement[iAmpRange];
				rdByAmpQuant->rate[i][j][iAmpRange] += computeRateSBExpFreq(static_cast<double>(nb_xi),ampRangeNumElement[iAmpRange]);
				rdByAmpQuant->rate[i][j][iAmpRange] += computeRateSBExpSample(static_cast<double>(nb_sample),ampRangeNumElement[iAmpRange]);
				rdByAmpQuant->rate[i][j][iAmpRange] += chosenR_rho;
				rdByAmpQuant->rate[i][j][iAmpRange] += chosenR_phi;

				rdByAmpQuant->quantExp[i][j][iAmpRange] = quantExp;
				rdByAmpQuant->nAtom[i][j][iAmpRange] = ampRangeNumElement[iAmpRange];

				rate+=rdByAmpQuant->rate[i][j][iAmpRange];

				delete [] rhovec;
			} // iAmpRange

			totalRate += static_cast<double>(rate);

			delete [] (CStructBookExp*)structBookByAmp;
			delete [] ampBar;
			delete [] ampRangeLimit;
			delete [] ampRangeNumElement;
			delete [] pSBQ;

		} // block
	} // signal
	rdByAmpQuant->totalRate = totalRate;
}


void decodeSBExpAmpRangeArithCod(char*InputFile,
								 CFileDfTable* DfTableFile,
								 CFileDictionary* dicData)
{
	////////////////
	int i, j, k, i_qrho, i_qphi;
	int iAmpRange;
	int t;
	/////
	unsigned short int dummyushint;
	unsigned int dummyuint;
	int dummyint;
	float dummyfloat;
	double dummydouble;
	////////////////
	FILE* code_file;
    char fileName[_MAX_PATH];
    char* pos;
    char aux[_MAX_PATH];
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,".mpz");
    strcpy( &pos[0], aux);
    //
    code_file = fopen (fileName,"rb");

    /// READ HEADER //////////////////
    // numSignal
    fread (&dummyint, sizeof(unsigned int), 1, code_file );
    int numSignal = dummyint;
    cout << "-numSignal: " << numSignal << endl;
    // numBlock
    //fread (&dummyint, sizeof(int), 1, code_file );
    //int numBlock = dummyint;
    //cout << "-numBlock: " << numBlock << endl;
    // signalSize
    fread (&dummyint, sizeof(unsigned int), 1, code_file );
    int signalSize = dummyint;
    cout << "-signalSize: " << signalSize << endl;
    // blockSize
    fread (&dummyint, sizeof(unsigned int), 1 , code_file );
    int blockSize = dummyint;
    cout << "-blockSize: " << blockSize << endl;
    // subBlockSize
    fread (&dummyint, sizeof(unsigned int), 1 , code_file );
    int subBlockSize = dummyint;
    cout << "-subBlockSize: " << subBlockSize << endl;
    // numBlock
    int numBlock = (int)ceil((double)signalSize/(double)blockSize);
    cout << "-numBlock: " << numBlock << endl;
    // numSubBlock
    int numSubBlock =  static_cast<int>(blockSize / subBlockSize);
     cout << "-numSubBlock: " << numSubBlock << endl;
    // Fs - sample rate
    fread (&dummyuint, sizeof(unsigned int), 1 , code_file );
    double Fs  = static_cast<double>(dummyuint);
    cout << "-Fs: " << Fs << endl;
    // lambda
    fread (&dummyfloat, sizeof(float), 1, code_file );
    double chosenLambda = static_cast<double>(dummyfloat);
    cout << "-chosenLambda: " << chosenLambda << endl;
    // signal norm
    float* signalNorm = new float[numSignal];
    fread (signalNorm, sizeof(float), numSignal, code_file );
    for (i=0;i<numSignal;i++)
	{
		cout << "Signal norm - " << i+1 << " : " << signalNorm[i] << endl;
	}

	double** min_amp= new double*[numSignal];
    double** max_amp= new double*[numSignal];
    double** min_rho= new double*[numSignal];
    double** max_rho= new double*[numSignal];
    double** min_phase= new double*[numSignal];
    double** max_phase= new double*[numSignal];
    double** nb_amp=new double*[numSignal];
    for (i=0;i<numSignal;i++)
	{
		min_amp[i]= new double[numBlock];
    	max_amp[i]= new double[numBlock];
    	min_rho[i]= new double[numBlock];
    	max_rho[i]= new double[numBlock];
    	min_phase[i]= new double[numBlock];
    	max_phase[i]= new double[numBlock];
    	nb_amp[i]= new double[numBlock];
	}

	/// READ HEADER //////////////////
	cout << "nb_amp[i][j]" << "; " << "max_amp[i][j]" << "; " << "max_rho[i][j]" << endl;
	int deltaSupMax;
	for (i=0;i<numSignal;i++)
	{
		for (j=0;j<numBlock;j++)
		{

			// nbamp/NAmpRAnge
    		fread (&dummyuint, sizeof(unsigned int), 1, code_file );
    		nb_amp[i][j] = static_cast<int>(dummyuint);
			// max_amp
    		fread (&dummyfloat, sizeof(float), 1, code_file );
    		max_amp[i][j] = static_cast<double>(dummyfloat);
			// max_rho
    		fread (&dummyfloat, sizeof(float), 1 , code_file );
    		max_rho[i][j] = static_cast<double>(dummyfloat);
    		// deltaSupMax
    		fread (&dummyuint, sizeof(unsigned int), 1, code_file );
    		deltaSupMax = static_cast<int>(dummyuint);
			// OBS:
			//   min_amp = 0
			min_amp[i][j] = 0.0;
			//   min_rho =0
			min_rho[i][j] = 0.0;
			//   min_phi = 0
			min_phase[i][j] = 0.0;
			//   max_phi = 2*pi
			max_phase[i][j] = 2*pi;

			cout << nb_amp[i][j] << "; " << max_amp[i][j]<< "; " << max_rho[i][j] << "; " << deltaSupMax << endl;

		}
	}


	int Nfreq = dicData->getNumFreq();
	int fdiscrtype = dicData->getFDiscrType(0);
	double freqi = dicData->getFreqi(0);
	double step_xi = (2*pi/Fs)*freqi;

	///////////======
	strtContinuousExp* pSB;
    int sbNumElement;

    double**** DfTable = DfTableFile->getDfTable();
    int N_qrho = DfTableFile->getNumQRho();
    int N_qphi = DfTableFile->getNumQPhi();
    double* q_rho = DfTableFile->getQRho();
    double* q_phi = DfTableFile->getQPhi();

	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
	cout << "- Compute probabilities with respect to TIME SUPPORT step quantization using" << endl;
	cout << "  Generalized Gaussian Distribution" << endl;
	//int nQDeltaSupStep = sbbHeader.blockSize;
	int nQDeltaSupStep = deltaSupMax+1;
	double* deltaSupQStepCenter = new double[nQDeltaSupStep];
	double* deltaSupQStepEdge = new double[nQDeltaSupStep+1];
	double* deltaSupQStepProb = new double[nQDeltaSupStep];
	double* deltaSupQStepCumProb = new double[nQDeltaSupStep];
	double step_deltasup = 1.0;
// 	setQStepFeatureGGD( deltaSupQStepCenter,
// 						deltaSupQStepEdge,
// 						deltaSupQStepProb,
// 						deltaSupQStepCumProb,
// 						step_deltasup,
// 						nQDeltaSupStep,
// 						0,
// 						DELTASUPGGDMEAN, DELTASUPGGDSCALE, DELTASUPGGDSHAPE);
	setQStepFeatureGGDDeltaSup( deltaSupQStepCenter,
								deltaSupQStepEdge,
								deltaSupQStepProb,
								deltaSupQStepCumProb,
								step_deltasup,
								nQDeltaSupStep,
								0);
	double H_deltasup;
#ifdef USEGGD
	H_deltasup = computeEntropyPdf(deltaSupQStepProb,nQDeltaSupStep);
#else
	H_deltasup = computeEntropyUnifPdf(nQDeltaSupStep);
#endif

   	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
	cout << "- Compute probabilities with respect to AMPLITUDE step quantization using" << endl;
	cout << "  Generalized Gaussian Distribution" << endl;
	double**** ampQStepCenter = new double***[numSignal];
	double**** ampQStepEdge = new double***[numSignal];
	double**** ampQStepProb = new double***[numSignal];
	double**** ampQStepCumProb = new double***[numSignal];
	int nQAmpStep;
	int nb_amp_;
	for (i=0;i<numSignal;i++)
	{
		ampQStepCenter[i] = new double**[numBlock];
		ampQStepEdge[i] = new double**[numBlock];
		ampQStepProb[i] = new double**[numBlock];
	    ampQStepCumProb[i] = new double**[numBlock];
		for (j=0;j<numBlock;j++)
		{
			ampQStepCenter[i][j] = new double*[16];
			ampQStepEdge[i][j] = new double*[16];
			ampQStepProb[i][j] = new double*[16];
			ampQStepCumProb[i][j] = new double*[16];
			for (nb_amp_=16; nb_amp_>0; nb_amp_--)
			{
				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp_), min_amp[i][j], max_amp[i][j]);

				nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

				// normalize step_amp
				step_amp = step_amp/max_amp[i][j];

				ampQStepCenter[i][j][nb_amp_-1] = new double[nQAmpStep];
				ampQStepEdge[i][j][nb_amp_-1] = new double[nQAmpStep+1];
				ampQStepProb[i][j][nb_amp_-1] = new double[nQAmpStep];
				ampQStepCumProb[i][j][nb_amp_-1] = new double[nQAmpStep];

				setQStepFeatureGGD( ampQStepCenter[i][j][nb_amp_-1],
									ampQStepEdge[i][j][nb_amp_-1],
									ampQStepProb[i][j][nb_amp_-1],
									ampQStepCumProb[i][j][nb_amp_-1],
									step_amp,
									nQAmpStep,
									0,
									AMPGGDMEAN, AMPGGDSCALE, AMPGGDSHAPE);
									//DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);

			}
		}



	}
	////////////////////////////////////////////////////////////////////////
	// Compute probabilities with respect to decay step quantization
// 	cout << "- Compute probabilities with respect to decay step quantization using" << endl;
// 	cout << "  Generalized Gaussian Distribution" << endl;
     int nQRhoStep;
     int nQPhiStep;
// 	double**** rhoQStepCenter = new double***[numSignal];
// 	double**** rhoQStepEdge = new double***[numSignal];
// 	double**** rhoQStepProb = new double***[numSignal];
// 	double**** rhoQStepCumProb = new double***[numSignal];
//
// 	//cout << N_qrho << endl;
// 	for (i=0;i<numSignal;i++)
// 	{
// 		rhoQStepCenter[i] = new double**[numBlock];
// 		rhoQStepEdge[i] = new double**[numBlock];
// 		rhoQStepProb[i] = new double**[numBlock];
// 		rhoQStepCumProb[i] = new double**[numBlock];
// 		for (j=0;j<numBlock;j++)
// 		{
// 			rhoQStepCenter[i][j] = new double*[N_qrho];
// 			rhoQStepEdge[i][j] = new double*[N_qrho];
// 			rhoQStepProb[i][j] = new double*[N_qrho];
// 			rhoQStepCumProb[i][j] = new double*[N_qrho];
// 			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
// 			{
// 				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
//
// 				nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;
//
// 				//cout << "- nQRhoStep: " << nQRhoStep << endl;
//
// 				rhoQStepCenter[i][j][i_qrho] = new double[nQRhoStep];
// 				rhoQStepEdge[i][j][i_qrho] = new double[nQRhoStep+1];
// 				rhoQStepProb[i][j][i_qrho] = new double[nQRhoStep];
// 				rhoQStepCumProb[i][j][i_qrho] = new double[nQRhoStep];
//
// 				setQStepFeatureGGD( rhoQStepCenter[i][j][i_qrho],
// 									rhoQStepEdge[i][j][i_qrho],
// 									rhoQStepProb[i][j][i_qrho],
// 									rhoQStepCumProb[i][j][i_qrho],
// 									q_rho[i_qrho],
// 									nQRhoStep,
// 									static_cast<int>(nQRhoStep/2),
// 									DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);
// 			}
// 		}
// 	}

	// Compute Entropy
	cout << "- Compute entropy table for Amplitude and Phase." << endl;
	double*** entropyAmpTable = new double**[numSignal];
	double*** entropyPhiTable = new double**[numSignal];
	for (i=0;i<numSignal;i++)
	{
		entropyAmpTable[i] = new double*[numBlock];
		entropyPhiTable[i] = new double*[numBlock];
		for (j=0;j<numBlock;j++)
		{
			entropyAmpTable[i][j] = new double[16];
			entropyPhiTable[i][j] = new double[N_qphi];
			for (nb_amp_=0;nb_amp_<16;nb_amp_++)
			{
				double step_amp = computeMidRiseQuantStep( static_cast<double>(nb_amp_+1), min_amp[i][j], max_amp[i][j]);
				nQAmpStep = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;

#ifdef USEGGD
				entropyAmpTable[i][j][nb_amp_] = computeEntropyPdf(	ampQStepProb[i][j][nb_amp_],
																	nQAmpStep);
#else
				entropyAmpTable[i][j][nb_amp_] = computeEntropySBExpDecayUnifPdf(nQAmpStep);
#endif
			}
			for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
			{
				if ( (max_phase[i][j]/q_phi[i_qphi]) < 1.0) break;
				nQPhiStep =  static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;

				entropyPhiTable[i][j][i_qphi] = computeEntropySBExpPhaseUnifPdf(nQPhiStep);
			}
		}
	}

	/////////////////////////////////////////////////////
	// DECODE
#ifdef DEBUG
    fstream file_decIndex("rdbyamp_decodeIndex.out",ios::out);
	file_decIndex << setw (14) << setfill(' ') << "Signal"<<" "
	              << setw (14) << setfill(' ') << "Block"<<" "
	              << setw (14) << setfill(' ') << "subBlock"<<" "
	              << setw (14) << setfill(' ') << "iAmpRange"<<" "
	              << setw (14) << setfill(' ') << "Element"<<" "
	              << setw (14) << setfill(' ') << "ampIndex"<<" "
	              << setw (14) << setfill(' ') << "rhoIndex"<<" "
	              << setw (14) << setfill(' ') << "xiIndex"<<" "
	              << setw (14) << setfill(' ') << "phiIndex"<<" "
	              << setw (14) << setfill(' ') << "initSamp"<<" "
	              << setw (14) << setfill(' ') << "deltaSamp"<<" "
	              << setw (14) << setfill(' ') << "aIndex"<<" "
	              << setw (14) << setfill(' ') << "bIndex"<<" " << endl;
#endif

	double* recBlockSignal;
    recBlockSignal = new double[blockSize];

    double** recSignal;
	recSignal = new double*[numSignal];

#ifdef	USEGGD
	double** rhoQStepCenter = new double*[N_qrho];
	double** rhoQStepEdge = new double*[N_qrho];
	double** rhoQStepProb = new double*[N_qrho];
	double** rhoQStepCumProb = new double*[N_qrho];
#endif
	double* entropyRhoTable = new double[N_qrho];


	int iSubBlock;
    for (i=0;i<numSignal;i++)
	{
		recSignal[i] = new double[signalSize];
		for (j=0;j<numBlock;j++)
		{

			////////////////////////////////////////////////////////////////////////
			// Compute probabilities with respect to decay step quantization
			cout << "- Compute probabilities and entropies with respect to decay step quantization using" << endl;
			cout << "  Generalized Gaussian Distribution" << endl;
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;

				nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;

				//cout << "- nQRhoStep: " << nQRhoStep << endl;
#ifdef	USEGGD
				rhoQStepCenter[i_qrho] = new double[nQRhoStep];
				rhoQStepEdge[i_qrho] = new double[nQRhoStep+1];
				rhoQStepProb[i_qrho] = new double[nQRhoStep];
				rhoQStepCumProb[i_qrho] = new double[nQRhoStep];

				setQStepFeatureGGD( rhoQStepCenter[i_qrho],
									rhoQStepEdge[i_qrho],
									rhoQStepProb[i_qrho],
									rhoQStepCumProb[i_qrho],
									q_rho[i_qrho],
									nQRhoStep,
									static_cast<int>(nQRhoStep/2),
									DECAYGGDMEAN, DECAYGGDSCALE, DECAYGGDSHAPE);
#endif

#ifdef USEGGD
				entropyRhoTable[i_qrho] = computeEntropySBExpDecayGGDPdf(	rhoQStepProb[i_qrho],
																			nQRhoStep);
#else
				entropyRhoTable[i_qrho] = computeEntropySBExpDecayUnifPdf(nQRhoStep);
#endif
			}


			// Define mean amplitude in range
			int NAmpRange = nb_amp[i][j];

			double* ampRangeLimit;
			ampRangeLimit = new double[NAmpRange+1];

			double* ampBar;
			ampBar = new double[NAmpRange];
			double L_amp = max_amp[i][j] - min_amp[i][j];

			ampRangeLimit[0]=-0.001;
			for (iAmpRange=0;iAmpRange<NAmpRange;iAmpRange++)
			{
				ampRangeLimit[iAmpRange+1] =
					L_amp/(pow(2.0,static_cast<double>(NAmpRange-iAmpRange-1)));
				ampBar[NAmpRange-1-iAmpRange] =
					(3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(iAmpRange)));
			}

			// Initialize block reconstruction vector (RESET)
			for (k=0;k <blockSize ;k++)
			{
				recBlockSignal[k] = 0.0;
			}

			////

			unsigned char* compressed_data=0;
			unsigned max_code_bytes = 0xFFFFFFE0U;
			////////////////////////////////////
			// Arithmetic Encoder object
			Arithmetic_Codec acd;
			acd.set_buffer(max_code_bytes, compressed_data);

			// START DECODER
			acd.read_from_file(code_file); // Alternatively acd.start_decoder();
			////////////////////////////////////

			// MODEL AMP
			Static_Data_Model ampModel;
			double step_amp = computeMidRiseQuantStep( static_cast<double>(NAmpRange), min_amp[i][j], max_amp[i][j]);
			unsigned nlevel_amp = static_cast<int>(round(max_amp[i][j]/step_amp)) + 1;
			//unsigned nlevel_amp = static_cast<unsigned>(round(pow(2.0, static_cast<double>(chosenQuantExp.nb_amp) ) ));
			cout << "- nlevel_amp: " << nlevel_amp << endl;
#ifdef USEGGD
			ampModel.set_distribution(static_cast<unsigned>(nlevel_amp),
									  ampQStepProb[i][j][NAmpRange - 1]);
#else
			ampModel.set_distribution(static_cast<unsigned>(nlevel_amp));
#endif

			unsigned* ampIndex = new unsigned[100000];
			unsigned ampIndexAux;
			int* vecNElementSubBlock =new int[numSubBlock];
			k=0;
			int nel;
			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				nel=0;
				while(1)
				{
					ampIndexAux = acd.decode(ampModel);
					if (ampIndexAux==0) break;
					ampIndex[k] = ampIndexAux;
					k++;
					nel++;
				}
				vecNElementSubBlock[iSubBlock] = nel;
			}

			// STOP DECODER
			acd.stop_decoder();

			//////////////////////////
			int sbNumElementQ = k;
			cout << "- sbNumElementQ: " << sbNumElementQ << endl;
			double* ampvec = new double[sbNumElementQ];
			dequantizeStructBookExpArithCodAmpVec(	NAmpRange, min_amp[i][j],  max_amp[i][j],
											        ampIndex, sbNumElementQ, ampvec);

			/////////
			int ** tabNumIndex =new int*[numSubBlock];
			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				tabNumIndex[iSubBlock] = new int[NAmpRange+1];
			}

			tabNumIndex[0][0] = 0;
			int count=0;
			int indAmpRange=0;

			cout<< "tabNumIndex" << endl;
			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				tabNumIndex[iSubBlock][0] = count;
				iAmpRange =0;
				for(k=0;k<vecNElementSubBlock[iSubBlock];k++)
				{
					indAmpRange = getIndexAmpRange(ampvec[count], NAmpRange, ampRangeLimit);
					count++;
					///
					while(iAmpRange!=indAmpRange)
					{
						iAmpRange++;
						tabNumIndex[iSubBlock][iAmpRange] = count-1;
						cout << tabNumIndex[iSubBlock][iAmpRange] << endl;

					}
				}
				while(iAmpRange<NAmpRange)
				{
					iAmpRange++;
					tabNumIndex[iSubBlock][iAmpRange] = count;
				}
			}


			strtContinuousExp* pSBQ = new strtContinuousExp[sbNumElementQ];

			for(k=0;k<sbNumElementQ;k++)
				pSBQ[k].innerProduct = ampvec[k];

			int* chosenQRhoIndex = new int[NAmpRange];
			int* chosenQPhiIndex = new int[NAmpRange];
			for (iAmpRange=0;iAmpRange < NAmpRange;iAmpRange++)
			{
				CStructBook* structBookByAmp = new CStructBookExp;
				((CStructBookExp*)structBookByAmp)->sepByAmp(	pSBQ,
																sbNumElementQ,
																ampRangeLimit[iAmpRange],   //lowerAmpRangeLimit
																ampRangeLimit[iAmpRange+1]); //upperAmpRangeLimit

				int numElementAmp = ((CStructBookExp*)structBookByAmp)->getNumElement();
				//////////////////////////
				getQRhoQPhiArithCodDecode(  chosenLambda,
											min_amp[i][j], max_amp[i][j],
											min_rho[i][j], max_rho[i][j],
											min_phase[i][j], max_phase[i][j],
											DfTable,
											N_qrho, N_qphi,
											q_rho, q_phi,
											numElementAmp,
											ampBar[iAmpRange],
											NAmpRange,
											iAmpRange,
											chosenQRhoIndex[iAmpRange],
											chosenQPhiIndex[iAmpRange],
											entropyRhoTable,
											entropyPhiTable[i][j]);
				cout << "- iAmpRange: " << iAmpRange << endl;
				cout << "   chosenQRhoIndex: " << chosenQRhoIndex[iAmpRange] << endl;
				cout << "   chosenQPhiIndex: " << chosenQPhiIndex[iAmpRange] << endl;
			}

			int* rhoIndex = new int[sbNumElementQ];
			unsigned* xiIndex  = new unsigned[sbNumElementQ];
			unsigned* phiIndex = new unsigned[sbNumElementQ];
			unsigned* aIndex   = new unsigned[sbNumElementQ];
			unsigned* bIndex   = new unsigned[sbNumElementQ];

			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				// START DECODER
			    acd.read_from_file(code_file); // Alternatively acd.start_decoder();

				/// MODEL FREQUENCY
				Static_Data_Model xiModel;
				cout << "- Nfreq: " << Nfreq << endl;
				xiModel.set_distribution(static_cast<unsigned>(Nfreq));

				/// MODEL TIME SUPPORT
				/// A
				Static_Data_Model initTimeSupModel;
				cout << "- subBlockSize: " << subBlockSize << endl;
				initTimeSupModel.set_distribution(static_cast<unsigned>(subBlockSize));

				/// B/ DELTA
				Static_Data_Model timeSupModel;
				//cout << "- blockSize: " << blockSize << endl;
				cout << "- deltaSupMax: " << deltaSupMax << endl;
#ifdef USEGGD
				timeSupModel.set_distribution(static_cast<unsigned>(deltaSupMax+1), deltaSupQStepProb);
#else
				timeSupModel.set_distribution(static_cast<unsigned>(deltaSupMax+1));
#endif
				for (iAmpRange=0;iAmpRange < NAmpRange;iAmpRange++)
				{

					/// MODEL DECAY
					i_qrho = chosenQRhoIndex[iAmpRange];
					nQRhoStep = static_cast<int>(ceil(max_rho[i][j]/q_rho[i_qrho]))*2 + 1;
					Static_Data_Model decayModel;
					cout << "- nQRhoStep: " << nQRhoStep << endl;
#ifdef USEGGD
					decayModel.set_distribution(static_cast<unsigned>(nQRhoStep),
												rhoQStepProb[i_qrho]);
#else
					decayModel.set_distribution(static_cast<unsigned>(nQRhoStep));
#endif
					/// MODEL PHASE
					i_qphi = chosenQPhiIndex[iAmpRange];
					nQPhiStep = static_cast<int>(ceil(max_phase[i][j]/q_phi[i_qphi])) + 1;
					Static_Data_Model phiModel;
					cout << "- nQPhiStep: " << nQPhiStep << endl;
					phiModel.set_distribution(static_cast<unsigned>(nQPhiStep));

// 					cout << "- iSubBlock: "<< iSubBlock << endl;
// 					cout << "- iAmpRange: "<< iAmpRange << endl;
// 					cout << "- tabNumIndex[iSubBlock][iAmpRange]: " << tabNumIndex[iSubBlock][iAmpRange] << endl;
// 					cout << "- tabNumIndex[iSubBlock][iAmpRange+1]: " << tabNumIndex[iSubBlock][iAmpRange+1] << endl;
					for (k=tabNumIndex[iSubBlock][iAmpRange];k<tabNumIndex[iSubBlock][iAmpRange+1];k++)
					{
						// DECAY
						rhoIndex[k] = static_cast<int>(acd.decode(decayModel)) - static_cast<int>(nQRhoStep/2);
						// PHASE
						phiIndex[k]  = acd.decode(phiModel);
						// FREQUENCY
						xiIndex[k]  = acd.decode(xiModel);
						// A
						unsigned initSamp = acd.decode(initTimeSupModel);
						aIndex[k] = initSamp + static_cast<unsigned>(subBlockSize*iSubBlock);
						// B/ DELTA
						unsigned delta = acd.decode(timeSupModel);
						bIndex[k] = aIndex[k] + delta;

#ifdef DEBUG
						file_decIndex << setw (14) << setfill(' ') << i+1<<" "
									  << setw (14) << setfill(' ') << j+1<<" "
									  << setw (14) << setfill(' ') << iSubBlock+1<<" "
									  << setw (14) << setfill(' ') << iAmpRange+1<<" "
									  << setw (14) << setfill(' ') << k+1<<" "
									  << setw (14) << setfill(' ') << ampIndex[k] <<" "
									  << setw (14) << setfill(' ') << rhoIndex[k] <<" "
									  << setw (14) << setfill(' ') << xiIndex[k]  <<" "
									  << setw (14) << setfill(' ') << phiIndex[k] <<" "
									  << setw (14) << setfill(' ') << initSamp   <<" "
									  << setw (14) << setfill(' ') << delta   <<" "
									  << setw (14) << setfill(' ') << aIndex[k]   <<" "
									  << setw (14) << setfill(' ') << bIndex[k]   <<" " << endl;
#endif
					} // k

					printf(" ****** Synthesize signal %d block %d subBlock %d amp_range %d **** \n",i+1,j+1,iSubBlock+1,iAmpRange+1);
					int numElementAmp = tabNumIndex[iSubBlock][iAmpRange+1] - tabNumIndex[iSubBlock][iAmpRange];
					int numIndex = tabNumIndex[iSubBlock][iAmpRange];
					int QRhoIndex = chosenQRhoIndex[iAmpRange];
					int QPhiIndex = chosenQPhiIndex[iAmpRange];
					dequantizeStructBookExpArithCod( nb_amp[i][j], min_amp[i][j], max_amp[i][j],
													 fdiscrtype, step_xi,
													 q_rho[QRhoIndex], q_phi[QPhiIndex],
													 &ampIndex[numIndex], &rhoIndex[numIndex], &xiIndex[numIndex],
													 &phiIndex[numIndex], &aIndex[numIndex], &bIndex[numIndex],
													 &pSBQ[numIndex], numElementAmp);

					synthSignalSBExpNoReset(recBlockSignal,
											signalNorm[i],
											blockSize,
											&pSBQ[numIndex],
											numElementAmp);

				} // iAmpRange
				// STOP DECODER
				acd.stop_decoder();

			}// iSubBlock


			int initBlockSample = j*blockSize;
			if ((initBlockSample + blockSize) >= signalSize)
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*(signalSize-initBlockSample));
			}
			else
			{
				memcpy(	&recSignal[i][initBlockSample],
						recBlockSignal,
						sizeof(double)*blockSize);
			}

#ifdef USEGGD
			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
			{
				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
				delete [] rhoQStepCenter[i_qrho];
				delete [] rhoQStepEdge[i_qrho];
				delete [] rhoQStepProb[i_qrho];
				delete [] rhoQStepCumProb[i_qrho];
			}
#endif
			delete [] ampIndex;
			delete [] rhoIndex;
			delete [] xiIndex;
			delete [] phiIndex;
			delete [] aIndex;
			delete [] bIndex;

			delete [] pSBQ;

			delete [] ampBar;
			delete [] ampRangeLimit;

			for(iSubBlock=0; iSubBlock< numSubBlock; iSubBlock++)
			{
				delete [] tabNumIndex[iSubBlock];
			}
			delete [] tabNumIndex;

		} // block
	} //signal


	delete [] entropyRhoTable;
#ifdef USEGGD
	delete [] rhoQStepCenter;
	delete [] rhoQStepEdge;
	delete [] rhoQStepProb;
	delete [] rhoQStepCumProb;
#endif

	fclose(code_file);
#ifdef DEBUG
	file_decIndex.close();
#endif
	//// Save to file the reconstructed signals
    /////
    strcpy(fileName, InputFile);
    pos = strrchr( fileName, '.');
    sprintf(aux,"_decode.wav");
    strcpy( &pos[0], aux);
    //
    CDataSignal* outSignal = new CAudioSignal;
    outSignal->setFileName(fileName);
    outSignal->setNumSignal(numSignal);
    outSignal->setSignalSize(signalSize);
    outSignal->setSignal(recSignal);
    outSignal->setNorm();
    outSignal->saveSignal();
    delete (CAudioSignal*)outSignal;

	///////////////////////////////////
    ///////// DEALLOCATE //////////////
    delete [] deltaSupQStepCenter;
	delete [] deltaSupQStepEdge;
	delete [] deltaSupQStepProb;
	delete [] deltaSupQStepCumProb;
    for (i=0;i<numSignal;i++)
	{
    	for (j=0;j<numBlock;j++)
		{
			for (nb_amp_=0;nb_amp_<16;nb_amp_++)
			{
				delete [] ampQStepCenter[i][j][nb_amp_];
				delete [] ampQStepEdge[i][j][nb_amp_];
				delete [] ampQStepProb[i][j][nb_amp_];
				delete [] ampQStepCumProb[i][j][nb_amp_];
			}
			delete [] ampQStepCenter[i][j];
			delete [] ampQStepEdge[i][j];
			delete [] ampQStepProb[i][j];
			delete [] ampQStepCumProb[i][j];
		}
		delete [] ampQStepCenter[i];
		delete [] ampQStepEdge[i];
		delete [] ampQStepProb[i];
		delete [] ampQStepCumProb[i];
	}
	delete [] ampQStepCenter;
	delete [] ampQStepEdge;
	delete [] ampQStepProb;
	delete [] ampQStepCumProb;
	////////////////////////////////
//     for (i=0;i<numSignal;i++)
// 	{
//     	for (j=0;j<numBlock;j++)
// 		{
// 			for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
// 			{
// 				if ( (max_rho[i][j]/q_rho[i_qrho]) < 1.0) break;
// 				delete [] rhoQStepCenter[i][j][i_qrho];
// 				delete [] rhoQStepEdge[i][j][i_qrho];
// 				delete [] rhoQStepProb[i][j][i_qrho];
// 				delete [] rhoQStepCumProb[i][j][i_qrho];
// 			}
// 			delete [] rhoQStepCenter[i][j];
// 			delete [] rhoQStepEdge[i][j];
// 			delete [] rhoQStepProb[i][j];
// 			delete [] rhoQStepCumProb[i][j];
// 		}
// 		delete [] rhoQStepCenter[i];
// 		delete [] rhoQStepEdge[i];
// 		delete [] rhoQStepProb[i];
// 		delete [] rhoQStepCumProb[i];
// 	}
// 	delete [] rhoQStepCenter;
// 	delete [] rhoQStepEdge;
// 	delete [] rhoQStepProb;
// 	delete [] rhoQStepCumProb;

	for (i=0;i<numSignal;i++)
	{
		delete [] min_amp[i];
    	delete [] max_amp[i];
    	delete [] min_rho[i];
    	delete [] max_rho[i];
    	delete [] min_phase[i];
    	delete [] max_phase[i];
    	delete [] nb_amp[i];
	}
	delete [] min_amp;
	delete [] max_amp;
	delete [] min_rho;
	delete [] max_rho;
	delete [] min_phase;
	delete [] max_phase;
	delete [] nb_amp;
	////
	delete [] signalNorm;
	////////////
	for (i=0;i<numSignal;i++)
	{
		for (j=0;j<numBlock;j++)
		{
			delete [] entropyAmpTable[i][j];
			delete [] entropyPhiTable[i][j];
		}
		delete [] entropyAmpTable[i];
		delete [] entropyPhiTable[i];
	}
	delete [] entropyAmpTable;
	delete [] entropyPhiTable;
	///
	for (i=0; i<numSignal; i++)
	{
		delete [] recSignal[i];
	}
	delete [] recSignal;
	delete [] recBlockSignal;

}

void getQRhoQPhiArithCodDecode( 	double lambda,
									double min_amp, double max_amp,
									double min_rho, double max_rho,
									double min_phase, double max_phase,
									double**** DfTable,
									int N_qrho, int N_qphi,
									double* q_rho, double* q_phi,
									int ampRangeNumElement,
									double ampBar,
									int nb_amp,
									int iAmpRange,
									int& chosenQRhoIndex,
									int& chosenQPhiIndex,
									double* entropyRhoTable,
									double* entropyPhiTable)
{
	int i_qrho,i_qphi;

	double L_amp = max_amp - min_amp;
	double L_rho = max_rho - min_rho;
	double L_phi = max_phase - min_phase;

// 	cout << "L_amp: " << L_amp << endl;
// 	cout << "L_rho: " << L_rho << endl;
// 	cout << "L_phi: " << L_phi << endl;

	// Find nb_rho and nb_phi leading minimum Df Lagragian
	int nQPhiStep;
	int nQRhoStep;
	double chosenH_rho, chosenH_phi;
	double H_rho, H_phi;
	double LagrangianDf;
	double minLagrangianDf=1e100;
	for (i_qrho=0;i_qrho<N_qrho;i_qrho++)
	{
		if ( (max_rho/q_rho[i_qrho]) < 1.0) break;
		nQRhoStep = static_cast<int>(ceil(max_rho/q_rho[i_qrho]))*2 + 1;

		H_rho = entropyRhoTable[i_qrho];
		for (i_qphi=0;i_qphi<N_qphi;i_qphi++)
		{
			if ( (max_phase/q_phi[i_qphi]) < 1.0) break;
			nQPhiStep =  static_cast<int>(ceil(max_phase/q_phi[i_qphi])) + 1;

			H_phi = entropyPhiTable[i_qphi];

			LagrangianDf =  ( ampRangeNumElement *
							  (ampBar*ampBar) *
							  DfTable[nb_amp-1][i_qrho][i_qphi][iAmpRange] ) +
							( lambda *
								ampRangeNumElement * (H_rho + H_phi) );

			if (LagrangianDf<= minLagrangianDf)
			{
				minLagrangianDf = LagrangianDf;
				chosenQRhoIndex = i_qrho;
				chosenQPhiIndex = i_qphi;
				chosenH_rho = H_rho;
				chosenH_phi = H_phi;
			}
		} // nb_phase
	}  // nb_rho
}

int getIndexAmpRange(double amp, int NAmpRange, double* ampRangeLimit)
{
	int i;
	for(i=0;i<NAmpRange;i++)
	{
		if ((amp>ampRangeLimit[i]) && (amp<=ampRangeLimit[i+1]))
		{
			return i;
		}
	}

	return -1;
}


// void encodeLinkedStructBook(char* InputFile, CFileGenData* genData, double rateTarget)
// {
//     printf("-> Encoding Structure Book\n");
//
//     int i,j;
//
//     //=====================================
//     // Loading Structure Books
//     printf("1) Loading the Structure Books\n");
//
//     CFileBlockRange* blockRange;
//     blockRange = new CFileBlockRange;
//     blockRange->setFileName("blockrange.dat");
//     blockRange->loadData();
//
//     strtSBBHeader sbbHeader;
//     sbbHeader = loadSBHeader("header.sbb");
//
//     CStructBook** structBook = NULL;
//     structBook = new CStructBook* [sbbHeader.numSignal];
//     for (i=0; i<sbbHeader.numSignal; i++)
//     {
//         structBook[i] = new CStructBookExp[sbbHeader.numBlock];
//     }
//
//     loadSB( InputFile,
//             blockRange,
//             sbbHeader,
//             structBook);
//
// /*    CLinkStrBook* refLinkStrBookExp;
//     refLinkStrBookExp = new CRefLinkStrBookExp[sbbHeader.numSignal];
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
//         printf(" - Set signal %d norm\n",i+1);
//         //printf("norm: %f\n",sbbHeader.norm[i]);
//         ((CRefLinkStrBookExp*)refLinkStrBookExp)[i].setNorm(sbbHeader.norm[i]);
//         printf(" - Load\n");
//         ((CRefLinkStrBookExp*)refLinkStrBookExp)[i].load(  	structBook,
// 															sbbHeader,
// 															i);
// 	 	((CRefLinkStrBookExp*)refLinkStrBookExp)[i].Print(structBook,
// 															i);
// 	}
// 	delete [] refLinkStrBookExp;
//     exit(0);*/
//
//     //=====================================
//     // Re-organizing structure book by linking atoms
//     // with continuity along the blocks
//     printf("2) Re-organizing structure book by linking atoms.\n");
//     CLinkStrBook* linkStrBookExpRef;
//     linkStrBookExpRef = new CLinkStrBookExp[sbbHeader.numSignal];
//     double** recSignal;
// 	recSignal = new double*[sbbHeader.numSignal];
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		recSignal[i] = new double[sbbHeader.signalSize];
// 	}
// 	///////////////
// 	CDictionary* expDic = new CExpDictionary;
// 	expDic->setSignalSize(sbbHeader.blockSize);
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
//         printf(" - Set signal %d norm\n",i+1);
//         //printf("norm: %f\n",sbbHeader.norm[i]);
//         ((CLinkStrBookExp*)linkStrBookExpRef)[i].setNorm(sbbHeader.norm[i]);
//         printf(" - Load\n");
//         ((CLinkStrBookExp*)linkStrBookExpRef)[i].load(  structBook,
// 									sbbHeader,
// 									i);
// 		//printf(" - Print\n");
//         //((CLinkStrBookExp*)linkStrBookExpRef)[i].Print();
//         //printf(" - Synthesize signal(s)\n");
// 		//expDic->synthSignal(recSignal[i],linkStrBookExpRef,i,sbbHeader.signalSize);
// 	}
// 	/*FILE* iorec;
// 	iorec = fopen("recsignal.out","w");
// 	for (i=0;i<sbbHeader.numSignal;i++)
//     {
//     	for (j=0;j<sbbHeader.signalSize;j++)
// 		{
// 			fprintf(iorec," %f",recSignal[i][j]);
// 		}
// 		fprintf(iorec,"\n");
// 	}
// 	fclose(iorec);
//     */
//
// 	printf("3) Rate-distortion optimization.\n");
// 	// Setting fixed quantizers
// 	strtQuantExp fixedQuantExp;
// 	if (sbbHeader.type == 2)
//     {
//         double A0 = 27.5; // Hz
//         double C8 = 4186.01;
//         double freq1= A0 * pow(2.0,-23/12);
//         double freq2= C8 * pow(2.0,23/12); // B9
//         int Nfreq = (int)(24 * ceil( log10(freq2/freq1)/log10(2.0) ) )+1;
//         fixedQuantExp.nb_xi = (int)ceil( log(Nfreq) / log(2.0) );
//     }
//     fixedQuantExp.nb_sample = (int)ceil( log(sbbHeader.signalSize) / log(2.0) );
//     fixedQuantExp.nb_block = (int)ceil( log(sbbHeader.numBlock) / log(2.0) );
// 	strtQuantExp* chosenQuantExp;
// 	chosenQuantExp = new strtQuantExp[sbbHeader.numSignal];
// 	////////////////////////////
// 	/// Load original signal
// 	CDataSignal* dataSignal;
//     if (sbbHeader.type==1)
//     {
//         dataSignal = new CComtradeSignal;
//     }
//     if (sbbHeader.type==2)
//     {
//         dataSignal = new CAudioSignal;
//     }
//     if (sbbHeader.type==3)
//     {
//         dataSignal = new CNoiseSignal;
//     }
//
//     dataSignal->setFileName(InputFile);
//     dataSignal->setBlockSize(sbbHeader.blockSize);
//     dataSignal->setBlockHop(sbbHeader.blockHop);
//     dataSignal->setSignal();
//     dataSignal->setNorm();
//
//     double** origSignal = dataSignal->getSignal();
//     /*FILE* origFile;
//     origFile = fopen("orig_signal.out","w");
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
//     	for (j=0;j<sbbHeader.signalSize;j++)
//     	{
//     		fprintf(origFile," %f",origSignal[i][j]);
//     	}
//     	fprintf(origFile,"\n");
//     }
//     fclose(origFile);
//     */
//     int rdoptFlag = genData->getFlagRDopt();
// 	optimizeRDLinkedStrBook(	origSignal,
//            		sbbHeader,
//          		linkStrBookExpRef,
//         		chosenQuantExp,
//          		rateTarget,
//          		fixedQuantExp,
//          		rdoptFlag,
//                 InputFile);
//
//     printf("4) Encoding.\n");
//
//     int nbits_signal_total=0;
//     double mse_total=0.0;
//     double mse_total_per_atom = 0.0;
//     int linkSBNumElement;
//     CLinkStrBook* linkStrBookExp;
//     linkStrBookExp = new CLinkStrBookExp[sbbHeader.numSignal];
//     double *recAtom;
//     double *recAtomQuant;
// 	recAtom = new double[sbbHeader.signalSize];
// 	recAtomQuant = new double[sbbHeader.signalSize];
//
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
//         //chosenQuantExp[i].nb_amp = 16;
// 		//chosenQuantExp[i].nb_rho = 16;
// 		//chosenQuantExp[i].nb_xi  = 8;
// 		//chosenQuantExp[i].nb_phase = 16;
// 		//chosenQuantExp[i].nb_sample = 16;
// 		//chosenQuantExp[i].nb_block = 6;
//     	// Copy the Linked Structure Book
// 		printf(" - Copy the linked structure book\n");
// 		((CLinkStrBookExp*)linkStrBookExp)[i]=((CLinkStrBookExp*)linkStrBookExpRef)[i];
// 		//((CLinkStrBookExp*)linkStrBookExpRef)[i].Print();
// 		//((CLinkStrBookExp*)linkStrBookExp)[i].Print();
// 		////////////////////////////////////////////////////////////////////////
// 		printf(" - Configure quantizer\n");
// 		((CLinkStrBookExp*)linkStrBookExp)[i].configQuantizer(  chosenQuantExp[i].nb_amp,
// 																chosenQuantExp[i].nb_rho,
// 																chosenQuantExp[i].nb_xi,
// 																chosenQuantExp[i].nb_phase,
// 																chosenQuantExp[i].nb_sample,
// 																chosenQuantExp[i].nb_block);
// 		printf(" - Quantize structure book\n");
// 		((CLinkStrBookExp*)linkStrBookExp)[i].quantize();
// 		nbits_signal_total+= ((CLinkStrBookExp*)linkStrBookExp)[i].computeNumBits();
//
// 		// Synthesize signals from structure book
// 		printf(" - Synthesize signal(s)\n");
// 		//expDic->synthSignal(recSignal[i],linkStrBookExp,i,sbbHeader.signalSize);
// 		((CLinkStrBookExp*)linkStrBookExp)[i].synthSignal(recSignal[i],sbbHeader);
//
// 		// Synthesize linked atoms from structure book
// 		printf(" - Synthesize linked atoms\n");
// 		linkSBNumElement = ((CLinkStrBookExp*)linkStrBookExp)[i].getNumElement();
// 		printf("     - Total number of linked atoms: %d\n",linkSBNumElement);
// 		for(j=0;j<linkSBNumElement;j++)
// 		{
// 		    //printf("     - linked atom: %d\n",j);
// 			//expDic->synthAtom(recAtom,linkStrBookExpRef,i,sbbHeader.signalSize,j);
// 			((CLinkStrBookExp*)linkStrBookExpRef)[i].synthAtom(recAtom,sbbHeader,j);
// 			//expDic->synthAtom(recAtomQuant,linkStrBookExp,i,sbbHeader.signalSize,j);
// 			((CLinkStrBookExp*)linkStrBookExp)[i].synthAtom(recAtomQuant,sbbHeader,j);
// 			// Compute MSE per atom
// 			mse_total_per_atom += computeMSE(	recAtom,
// 												recAtomQuant,
// 												sbbHeader.signalSize);
//
// 		}
//
// 		// Compute MSE
// 		mse_total += computeMSE(origSignal[i],
// 								recSignal[i],
// 								sbbHeader.signalSize);
//
//     }
//
//     delete [] recAtom;
//     delete [] recAtomQuant;
//
//     double bitrate = (double)nbits_signal_total/
//                      (double)(static_cast<double>(sbbHeader.signalSize)*
//                               static_cast<double>(sbbHeader.numSignal) );
//     printf("Bitrate: %f bit/sample\n",bitrate);
//     printf("Total MSE: %f \n",mse_total);
//     printf("Total MSE per Atom: %f \n",mse_total_per_atom);
//     printf("Total Square Error: %f \n",mse_total*
//                                        static_cast<double>(sbbHeader.signalSize));
//     printf("Total Square Error per Atom: %f \n",mse_total_per_atom*
//                                                 static_cast<double>(sbbHeader.signalSize));
//
//     char fileName[_MAX_PATH];
//     strcpy(fileName, InputFile);
//     char* pos;
//     pos = strrchr( fileName, '.');
//     char aux[_MAX_PATH];
//     sprintf(aux,"_rec%09.6fbps.wav",bitrate);
//     strcpy( &pos[0], aux);
//
//     CDataSignal* outSignal;
//     outSignal = new CAudioSignal;
//     outSignal->setFileName(fileName);
//     outSignal->setNumSignal(sbbHeader.numSignal);
//     outSignal->setSignalSize(sbbHeader.signalSize);
//     outSignal->setSignal(recSignal);
//     outSignal->setNorm();
//     outSignal->saveSignal();
//     delete outSignal;
//
//     // RD optimization
// //    printf("2) RD optimization\n");
// //    RDopt(  InputFile,
// //            sbbHeader,
// //            structBook);
//
//     // Quantize Structure Book
// //    printf("3) Quantize Structure Books\n");
// //    quantSB(    sbbHeader,
// //                structBook);
//
//     //////////////////////////////////////////////////////
//     // Deallocate matrix and vectors
//     delete [] linkStrBookExpRef;
//     delete [] linkStrBookExp;
//     for (i=0; i<sbbHeader.numSignal; i++)
//   	{
//     	delete [] ((CStructBookExp*)structBook[i]);
//   	}
//   	delete [] structBook;
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		delete [] recSignal[i];
// 	}
//     delete [] recSignal;
//     delete blockRange;
//     if (sbbHeader.norm!=NULL) delete [] sbbHeader.norm;
//     delete expDic;
//     delete [] chosenQuantExp;
//     delete  dataSignal;
// }

void optimizeRDLinkedStrBook(double** origSignal,
							strtSBBHeader sbbHeader,
							CLinkStrBook* linkStrBookExpAux,
							strtQuantExp* chosenQuantExp,
							double rateTarget,
							strtQuantExp fixedQuantExp,
							int rdoptFlag,
							char* InputFile)
{

	int i, j, k, t;
    //////////////////////////////////
    // Defining the quantizers range
    CFileRDBitRange* rdBitRange;
    rdBitRange = new CFileRDBitRange;
    rdBitRange->setFileName("rdbitrange.dat");
    rdBitRange->loadData();

    int init_nbit_amp, end_nbit_amp, delta_nbit_amp;
    int init_nbit_rho, end_nbit_rho, delta_nbit_rho;
    int init_nbit_phase, end_nbit_phase, delta_nbit_phase;

    init_nbit_amp   =  rdBitRange->getInitNbitAmp();
    end_nbit_amp    =  rdBitRange->getEndNbitAmp();
    delta_nbit_amp  =  rdBitRange->getDeltaNbitAmp();
    init_nbit_rho   =  rdBitRange->getInitNbitRho();
    end_nbit_rho    =  rdBitRange->getEndNbitRho();
    delta_nbit_rho  =  rdBitRange->getDeltaNbitRho();
    init_nbit_phase =  rdBitRange->getInitNbitPhase();
    end_nbit_phase  =  rdBitRange->getEndNbitPhase();
    delta_nbit_phase  = rdBitRange->getDeltaNbitPhase();

    //printf("\n amp: %d %d %d\n",init_nbit_amp,end_nbit_amp,delta_nbit_amp);
    //printf("rho: %d %d %d\n",init_nbit_rho,end_nbit_rho,delta_nbit_rho);
    //printf("phase: %d %d %d\n\n",init_nbit_phase,end_nbit_phase,delta_nbit_phase);


    /*init_nbit_amp   =  4;
    end_nbit_amp    =  8;
    delta_nbit_amp  =  2;
    init_nbit_rho   =  12;
    end_nbit_rho    =  12;
    delta_nbit_rho  =  2;
    init_nbit_phase =  12;
    end_nbit_phase  =  12;
    delta_nbit_phase  =  2;
    */
	delete rdBitRange;

    int numQuant=  (((end_nbit_amp - init_nbit_amp)/delta_nbit_amp) +1) *
                   (((end_nbit_rho - init_nbit_rho)/delta_nbit_rho) +1) *
                   (((end_nbit_phase - init_nbit_phase)/delta_nbit_phase) +1)
                   + 1;


    strtQuantExp* quantExp;
    quantExp = new strtQuantExp[numQuant];
    quantExp[0].nb_amp = 0;
    quantExp[0].nb_rho = 0;
    quantExp[0].nb_phase = 0;
    t=1;
    for (i=init_nbit_amp;i<=end_nbit_amp;i=i+delta_nbit_amp)
    {
        for (j=init_nbit_rho;j<=end_nbit_rho;j=j+delta_nbit_rho)
        {
            for (k=init_nbit_phase;k<=end_nbit_phase;k=k+delta_nbit_phase)
            {
                quantExp[t].nb_amp = i;
                quantExp[t].nb_rho = j;
                quantExp[t].nb_phase = k;
                t++;
            }
        }
    }
    strtRD** rd = new strtRD*[sbbHeader.numSignal];
    strtRDIndex* rdIndex = new strtRDIndex[sbbHeader.numSignal];
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		rd[i]= new strtRD[numQuant];
		rdIndex[i].theta_vec =NULL;
		rdIndex[i].index_vec =NULL;
	}

    int nbits_signal=0;
    ///////////////////////////
    double** recSignal;
	recSignal = new double*[sbbHeader.numSignal];
	for (i=0;i<sbbHeader.numSignal;i++)
	{
		recSignal[i] = new double[sbbHeader.signalSize];
	}
	///////////////
	CDictionary* expDic = new CExpDictionary;
	expDic->setSignalSize(sbbHeader.blockSize);
	CLinkStrBook* linkStrBookExp;
    linkStrBookExp = new CLinkStrBookExp[sbbHeader.numSignal];
    /////////////////
    FILE* rdoptFile;
    rdoptFile = fopen("rdopt.out","w");
    if (rdoptFlag==1 || rdoptFlag==2)
    {
    	loadRDFile(sbbHeader.numSignal,numQuant,rd);
    }
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		printf("Rate-distortion procedure for signal: %d\n",i+1);
		if (rdoptFlag==0)
		{
			// Performing RD calculation
			printf(" *** Performing RD calculation\n");
			printf("     Total number of quantizers: %d\n",numQuant);
            rd[i][0].rate = 0.0;
            rd[i][0].dist = (double)(sbbHeader.norm[i]*sbbHeader.norm[i])/(double)sbbHeader.signalSize;
			for (k=1; k<numQuant; k++)
			{
				printf("Quantizer: %d\n",k);
				// Copy the Linked Structure Book
				//printf(" - Copy the linked structure book\n");
				((CLinkStrBookExp*)linkStrBookExp)[i]=((CLinkStrBookExp*)linkStrBookExpAux)[i];
				//printf("Print before quant\n");
				//((CLinkStrBookExp*)linkStrBookExp)[i].Print();
				////////////////////////////////////////////////////////////////////////
				//printf(" - Configure quantizer\n");
				((CLinkStrBookExp*)linkStrBookExp)[i].configQuantizer(  quantExp[k].nb_amp,
													quantExp[k].nb_rho,
													fixedQuantExp.nb_xi,
													quantExp[k].nb_phase,
													fixedQuantExp.nb_sample,
													fixedQuantExp.nb_block);
				//printf(" - Quantize structure book\n");
				((CLinkStrBookExp*)linkStrBookExp)[i].quantize();
				nbits_signal= ((CLinkStrBookExp*)linkStrBookExp)[i].computeNumBits();
				//printf("Print after quant\n");
				//((CLinkStrBookExp*)linkStrBookExp)[i].Print();
				// Compute Rate bit/sample per signal
				rd[i][k].rate = (double)nbits_signal/
								(double)(sbbHeader.signalSize);
				// Synthesize signals from structure book
				//printf(" - Synthesize signal(s)\n");
				//expDic->synthSignal(recSignal[i],linkStrBookExp,i,sbbHeader.signalSize);
				((CLinkStrBookExp*)linkStrBookExp)[i].synthSignal(recSignal[i],sbbHeader);
				// Compute MSE
				rd[i][k].dist = computeMSE(origSignal[i],
										recSignal[i],
										sbbHeader.signalSize);
			}
		}
		///////////////////////////////////
		// Perform RD optimization
		printf(" *** Perform RD optimization\n");
		chooseBestQuantSetOptDecr(	i,
									numQuant,
									rd,
									rdIndex);
		///////////////////////////////////////////////
	    // Choose RD optimum quantizer for each signal
	    printf(" *** Choose RD optimum quantizer\n");
	    double minRateDiff=1e8;
	    int ind, chosenInd=0;
	    fprintf(rdoptFile,"Signal: %d\n",i+1);
	    for (j=0;j<rdIndex[i].numElement;j++)
		{
		    ind = rdIndex[i].index_vec[j];
		    fprintf(rdoptFile," %d",ind+1);
	    	if (fabs(rd[i][ind].rate - rateTarget)<minRateDiff)
	    	{
	    	    chosenInd = ind;
	    		chosenQuantExp[i].nb_amp = quantExp[ind].nb_amp;
	    		chosenQuantExp[i].nb_rho = quantExp[ind].nb_rho;
	    		chosenQuantExp[i].nb_phase = quantExp[ind].nb_phase;
	    		chosenQuantExp[i].nb_xi = fixedQuantExp.nb_xi;
	    		chosenQuantExp[i].nb_sample = fixedQuantExp.nb_sample;
	    		chosenQuantExp[i].nb_block = fixedQuantExp.nb_block;
	    		minRateDiff = fabs(rd[i][ind].rate - rateTarget);
	    	}
		}
		fprintf(rdoptFile,"\n");
		printf("      Chosen Quantizer:\n");
		printf("      ind: %d\n",chosenInd);
		printf("      nb_amp: %d\n",chosenQuantExp[i].nb_amp);
		printf("      nb_rho: %d\n",chosenQuantExp[i].nb_rho);
		printf("      nb_phase: %d\n",chosenQuantExp[i].nb_phase);
		printf("      nb_xi: %d\n",chosenQuantExp[i].nb_xi);
		printf("      nb_sample: %d\n",chosenQuantExp[i].nb_sample);
		printf("      nb_block: %d\n",chosenQuantExp[i].nb_block);
	}
	fclose(rdoptFile);
	/*
     * Write rate-distortion information to file
     */
	if (rdoptFlag==0)
    {
        printf("- Writing rate-distortion information to file...\n");
        FILE* rdFile;
		rdFile = fopen("rate_distortion.out","w");
		double dist_dB;
		for (k=0; k<numQuant; k++)
		{
			fprintf(rdFile," %4d %4d %4d %4d %4d %4d",	quantExp[k].nb_amp,
														quantExp[k].nb_rho,
														fixedQuantExp.nb_xi,
														quantExp[k].nb_phase,
														fixedQuantExp.nb_sample,
														fixedQuantExp.nb_block);
			for (i=0;i<sbbHeader.numSignal;i++)
			{
				//dist_dB = 20*(log(rd[i][k].dist)/log(10.0));
				fprintf(rdFile," %13.10f %13.10f",rd[i][k].rate,
												rd[i][k].dist);
			}
			fprintf(rdFile,"\n");
		}
		fclose(rdFile);
	}
    /*
     * Generate encoded files using the optimum rate distortion pairs
     */
    if (rdoptFlag==2)
    {
        int ind;
        printf("- Generating optimum R-D encoded files...\n");
        for (i=0;i<sbbHeader.numSignal;i++)
        {
            for (j=0;j<rdIndex[i].numElement;j++)
            {
                ind = rdIndex[i].index_vec[j];
                // Copy the Linked Structure Book
                ((CLinkStrBookExp*)linkStrBookExp)[i]=((CLinkStrBookExp*)linkStrBookExpAux)[i];
                ////////////////////////////////////////////////////////////////////////
                //printf(" - Configure quantizer\n");
                ((CLinkStrBookExp*)linkStrBookExp)[i].configQuantizer(  quantExp[ind].nb_amp,
                                                    quantExp[ind].nb_rho,
                                                    fixedQuantExp.nb_xi,
                                                    quantExp[ind].nb_phase,
                                                    fixedQuantExp.nb_sample,
                                                    fixedQuantExp.nb_block);
                // quantize
                ((CLinkStrBookExp*)linkStrBookExp)[i].quantize();
                // Synthesize signals from structure book
                //expDic->synthSignal(recSignal[i],linkStrBookExp,i,sbbHeader.signalSize);
                ((CLinkStrBookExp*)linkStrBookExp)[i].synthSignal(recSignal[i],sbbHeader);

                char fileName[_MAX_PATH];
                strcpy(fileName, InputFile);
                char* pos;
                pos = strrchr( fileName, '.');
                char aux[_MAX_PATH];
                sprintf(aux,"s%d_rec%09.6fbps.wav",i,rd[i][ind].rate);
                strcpy( &pos[0], aux);

                CDataSignal* outSignal;
                outSignal = new CAudioSignal;
                outSignal->setFileName(fileName);
                outSignal->setNumSignal(sbbHeader.numSignal);
                outSignal->setSignalSize(sbbHeader.signalSize);
                outSignal->setSignal(recSignal);
                outSignal->setNorm();
                outSignal->saveSignal();
                delete outSignal;
            }
        }
    }

    delete expDic;
    delete [] quantExp;
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] rd[i];
	}
    delete [] rd;
    delete [] rdIndex;
    for (i=0;i<sbbHeader.numSignal;i++)
	{
		delete [] recSignal[i];
	}
    delete [] recSignal;
    delete [] linkStrBookExp;
}

void chooseBestQuantSetOptDecr(	int nChannel,
								int numQuant,
								strtRD** rd,
								strtRDIndex* rdIndex)
{
	double* theta_vec;
	double* theta_vec_aux;
	int* index_vec;
	int* index_vec_aux;
	int numElement = 0;
	theta_vec = NULL;
	theta_vec_aux = NULL;
	index_vec = NULL;
	index_vec_aux = NULL;

	int firstPoint=0;
	int secondPoint=0;

	int i;

	double max_rate=0;
	// Find first point and max rate
	//printf("- Find first point and max rate.\n");
	for(i=0; i<numQuant; i++)
	{
		if (rd[nChannel][firstPoint].rate > rd[nChannel][i].rate)
		{
			firstPoint = i;
		}
		if (max_rate < rd[nChannel][i].rate)
		{
			max_rate = rd[nChannel][i].rate;
		}
	}

	//printf("- Building Convex hull.\n");
	double min_theta;
	double theta=0;
	double delta_rate=0;
	double delta_dist=0;
	double min_delta_rate=0;

	while(1)
	{
		min_theta = pi / 2.0;
		min_delta_rate = max_rate;
		for(i=0; i<numQuant; i++)
		{
			if ((i!=firstPoint) &&
				(rd[nChannel][i].rate > rd[nChannel][firstPoint].rate) )
			{
				delta_rate = rd[nChannel][i].rate - rd[nChannel][firstPoint].rate;
				delta_dist = rd[nChannel][i].dist - rd[nChannel][firstPoint].dist;

				if (delta_rate==0)
				{
					if(delta_dist>=0)
						theta = pi/2.0;
					else
						theta = -pi/2.0;
				}
				else
				{
					theta = atan(delta_dist/delta_rate);
				}

				if ( (theta < min_theta) && (theta <0))
				{
					min_theta = theta;
					secondPoint = i;
				}
			}
		}

		if (numElement!=0)
		{
			if (index_vec_aux!=NULL)
			{
				delete [] index_vec_aux;
			}
			index_vec_aux = new int[numElement];
			memcpy(index_vec_aux,index_vec,sizeof(int)*numElement);

			if (theta_vec_aux!=NULL)
			{
				delete [] theta_vec_aux;
			}
			theta_vec_aux = new double[numElement];
			memcpy(theta_vec_aux,theta_vec,sizeof(double)*numElement);
		}
		if (index_vec!=NULL)
		{
			delete [] index_vec;
		}
		if (theta_vec!=NULL)
		{
			delete [] theta_vec;
		}
		numElement++;
		index_vec = new int[numElement];
		theta_vec = new double[numElement];
		if ((numElement-1)!=0)
		{
			memcpy(index_vec,index_vec_aux,sizeof(int)*(numElement-1));
			memcpy(theta_vec,theta_vec_aux,sizeof(double)*(numElement-1));
		}
		index_vec[numElement-1] = firstPoint;
		theta_vec[numElement-1] = min_theta;

		if (min_theta>0)
		{
			theta_vec[numElement-1] = 0;
			break;
		}

		firstPoint = secondPoint;

		if (rd[nChannel][firstPoint].rate == max_rate)
		{
			if (index_vec_aux!=NULL)
			{
				delete [] index_vec_aux;
			}
			index_vec_aux = new int[numElement];
			memcpy(index_vec_aux,index_vec,sizeof(int)*numElement);

			if (theta_vec_aux!=NULL)
			{
				delete [] theta_vec_aux;
			}
			theta_vec_aux = new double[numElement];
			memcpy(theta_vec_aux,theta_vec,sizeof(double)*numElement);

			if (index_vec!=NULL)
			{
				delete [] index_vec;
			}
			if (theta_vec!=NULL)
			{
				delete [] theta_vec;
			}
			numElement++;

			index_vec = new int[numElement];
			memcpy(index_vec,index_vec_aux,sizeof(int)*(numElement-1));
			index_vec[numElement-1] = firstPoint;

			theta_vec = new double[numElement];
			memcpy(theta_vec,theta_vec_aux,sizeof(double)*(numElement-1));
			theta_vec[numElement-1] = 0;

			break;
		}
	}

	//printf("- Load convex hull.\n");
	//printf("- numElement %d nChannel %d\n",numElement,nChannel);
	//ofstream outFile("bestq_ind.dat",ios::out);

	if (rdIndex[nChannel].theta_vec == NULL)
	{
		rdIndex[nChannel].theta_vec = new double[numElement];
	}
	if (rdIndex[nChannel].index_vec == NULL)
	{
		rdIndex[nChannel].index_vec = new int[numElement];
	}
	rdIndex[nChannel].numElement = numElement;
	printf("      Index | Rate| Theta\n");
	for(int g=0;g<numElement;g++)
	{
		//outFile << index_vec[g] << ' ';
		rdIndex[nChannel].theta_vec[g] = theta_vec[g];
		rdIndex[nChannel].index_vec[g] = index_vec[g];

		printf("      %d %f %f\n",rdIndex[nChannel].index_vec[g],
							rd[nChannel][index_vec[g]].rate,
							rdIndex[nChannel].theta_vec[g]);
	}
	cout << endl;

	if (index_vec!=NULL) delete [] index_vec;
	if (index_vec_aux!=NULL) delete [] index_vec_aux;
	if (theta_vec!=NULL) delete [] theta_vec;
	if (theta_vec_aux!=NULL) delete [] theta_vec_aux;
}

/*
void RDopt(char* InputFile,
           strtSBBHeader sbbHeader,
           CStructBook** structBook)
{
    CExpDictionary* expDic = new CExpDictionary;
    int nbits = (int)ceil( log10( (double)(sbbHeader.blockSize) )  /  log10( (double)(2) ) );
    int sigSize = (int) pow(2.0,(double)nbits);
    expDic->setSignalSize(sigSize);

    CDataSignal* dataSignal;
    if (sbbHeader.type==1)
    {
        dataSignal = new CComtradeSignal;
    }
    if (sbbHeader.type==2)
    {
        dataSignal = new CAudioSignal;
    }
    if (sbbHeader.type==3)
    {
        dataSignal = new CNoiseSignal;
    }

    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(sbbHeader.blockSize);
    dataSignal->setBlockHop(sbbHeader.blockHop);
    dataSignal->setSignal();
    dataSignal->setNorm();

    double** origSignal = dataSignal->getSignal();
    double* pOrigBlock = new double[sbbHeader.blockSize];
    double* pRecBlock = new double[sbbHeader.blockSize];


    //////////////////////////////////
    // Defining the quantizers range
    int init_nbit_amp, end_nbit_amp;
    int init_nbit_rho, end_nbit_rho;
    int init_nbit_phase, end_nbit_phase;

    init_nbit_amp   =  1;
    end_nbit_amp    =  4;
    init_nbit_rho   =  1;
    end_nbit_rho    =  4;
    init_nbit_phase =  1;
    end_nbit_phase  =  4;

    int numQuant =  (end_nbit_amp - init_nbit_amp +1) *
                    (end_nbit_rho - init_nbit_rho +1) *
                    (end_nbit_phase - init_nbit_phase +1);

    int i, j, k, t;
    strtQuantExp* quantExp;
    quantExp = new strtQuantExp[numQuant];
    strtRD* rd = new strtRD[numQuant];

    t=0;
    for (i=init_nbit_amp;i<end_nbit_amp;i++)
    {
        for (j=init_nbit_rho;j<end_nbit_rho;j++)
        {
            for (k=init_nbit_phase;k<end_nbit_phase;k++)
            {
                quantExp[t].nb_amp = i;
                quantExp[t].nb_rho = j;
                quantExp[t].nb_phase = k;
                t++;
            }
        }
    }
    ////////////////////////////////////
    // Performing RD calculation
    CStructBook** quantSB;
    quantSB = new CStructBook* [sbbHeader.numSignal];
    strtContinuousExp* pSB;

    // frequency range
    double freq1, freq2;
    if (sbbHeader.type==2)
    {
        double A0 = 27.5; // Hz
        double C8 = 4186.01;
        freq1= A0 * pow(2.0,-23/12);
        freq2= C8 * pow(2.0,23/12); // B9
    }

    for (i=0;i<sbbHeader.numSignal;i++)
    {
        quantSB[j] = new CStructBookExp[sbbHeader.numBlock];
        for (j=0;j<sbbHeader.numBlock;j++)
        {
            pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
            if (pSB != NULL)
            {
                for (k=0; k<numQuant; k++)
                {
                    ((CStructBookExp*)quantSB[i])[j].setQuantConfig(
                                                    pSB,
                                                    ((CStructBookExp*)structBook[i])[j].getNumElement(),
                                                    ((CStructBookExp*)structBook[i])[j].getNorm(),
                                                    quantExp[k].nb_amp,   // nbits_amp,
                                                    quantExp[k].nb_rho,   // nbits_rho,
                                                    quantExp[k].nb_phase, // nbits_phase,
                                                    freq1,                // (Ffund/ Finit) (electric/audio)
                                                    freq2,                // (Fs/ Fend) (electric/audio)
                                                    sbbHeader.signalSize,
                                                    sbbHeader.type);
                    ((CStructBookExp*)quantSB[i])[j].quantStructBook(
                                                    pSB,
                                                    ((CStructBookExp*)structBook[i])[j].getNumElement());

                    // Compute Rate
                    rd[k].rate = ((CStructBookExp*)quantSB[i])[j].computeRate(sbbHeader.blockSize);

                    // Load Original Block
                    int initBlock = (j*sbbHeader.blockHop);
                    memcpy(pOrigBlock,&origSignal[i][initBlock], sbbHeader.blockSize*sizeof(double));

                    // Load Reconstructed Block
                    synthSignal(pRecBlock,
                                ((CStructBookExp*)quantSB[i])[j].getNorm(),
                                sbbHeader.blockSize,
                                ((CStructBookExp*)quantSB[i])[j].getStructBook(),
                                ((CStructBookExp*)quantSB[i])[j].getNumElement());


                    // Compute MSE
                    rd[k].dist = computeMSE(pOrigBlock,
                                            pRecBlock,
                                            sbbHeader.blockSize);
                }
            }
        }
    }

    // Deallocate vectors and matrices
    for (i=0; i<sbbHeader.numSignal; i++)
    {
        delete [] ((CStructBookExp*)quantSB[i]);
    }
    delete [] quantSB;
    delete dataSignal;
    delete expDic;
    delete [] quantExp;
    delete [] rd;
    delete [] pOrigBlock;
    delete [] pRecBlock;
}

*/
void quantSB(strtSBBHeader sbbHeader,
             CStructBook** structBook)
{
    CStructBook** quantSB;
    quantSB = new CStructBook* [sbbHeader.numSignal];
    strtContinuousExp* pSB;

    // frequency range
    double freq1, freq2;
    if (sbbHeader.type==2)
    {
        double A0 = 27.5; // Hz
        double C8 = 4186.01;
        freq1= A0 * pow(2.0,-23/12);
        freq2= C8 * pow(2.0,23/12); // B9
    }

    int i,j;
    for (i=0;i<sbbHeader.numSignal;i++)
    {
        quantSB[j] = new CStructBookExp[sbbHeader.numBlock];
        for (j=0;j<sbbHeader.numBlock;j++)
        {
            pSB = ((CStructBookExp*)structBook[i])[j].getStructBook();
            if (pSB != NULL)
            {
                ((CStructBookExp*)quantSB[i])[j].setQuantConfig(
                                                pSB,
                                                ((CStructBookExp*)structBook[i])[j].getNumElement(),
                                                ((CStructBookExp*)structBook[i])[j].getNorm(),
                                                6,// nbits_amp,
                                                6,// nbits_rho,
                                                6,// nbits_phase,
                                                freq1, // (Ffund/ Finit) (electric/audio)
                                                freq2, // (Fs/ Fend) (electric/audio)
                                                sbbHeader.signalSize,
                                                sbbHeader.type);
                ((CStructBookExp*)quantSB[i])[j].quantStructBook(
                                                pSB,
                                                ((CStructBookExp*)structBook[i])[j].getNumElement());
            }
        }
    }

    for (i=0; i<sbbHeader.numSignal; i++)
    {
        delete [] ((CStructBookExp*)quantSB[i]);
    }
    delete [] quantSB;
}


void synthSignal( double* recSignal,
				 double norm,
				 int sigSize,
				 strtContinuousExp* sb,
				 int numElement)
{
    CExpDictionary* expDic = new CExpDictionary;
    expDic->setSignalSize(sigSize);
    CParameter* parm;
    parm = new CExpParm;
    int i,j;
    double* sigAux;
    for (i=0;i<sigSize;i++)
    {
        recSignal[i] = 0;
    }
    for (i=0;i<numElement;i++)
    {
        ((CExpParm*)parm)->setParm(sb[i]);
        expDic->setRealAtom(parm);
        sigAux = expDic->getRealAtom();
        for (j=0;j<sigSize;j++)
        {
            recSignal[j] += norm*sigAux[j];
        }
    }
    delete expDic;
    delete parm;
}

void synthSignalStrBook(double* recSignal,
						CStructBook* structBook,
                        strtSBBHeader mainHeader)
{
    double norm = ((CStructBookExp*)structBook)->getNorm();
    //cout << "norm = " << norm << endl;
    CExpDictionary* expDic = new CExpDictionary;
    CParameter* parm;
    parm = new CExpParm;
    int i,j, sigSize, initBlockSample;
    double* sigAux;
    strtContinuousExp* sb;
    strtContinuousExp sb_aux;
    sb = ((CStructBookExp*)structBook)->getStructBook();
    for (i=0;i<mainHeader.signalSize;i++)
    {
        recSignal[i] = 0.0;
    }
    for (i=0;i<((CStructBookExp*)structBook)->getNumElement();i++)
    {
        sb_aux = sb[i];
        sb_aux.a = sb[i].a % mainHeader.blockSize;
        sb_aux.b = sb_aux.a + (sb[i].b - sb[i].a);
        //cout << sb[i].a << " ---- "  << sb[i].b << endl;
        ((CExpParm*)parm)->setParm(sb_aux);
        //((CExpParm*)parm)->printParm2Screen();
        sigSize = static_cast<int>(ceil(static_cast<double>(sb[i].b - sb[i].a + 1)/
                                         static_cast<double>(mainHeader.blockSize))*
                                   static_cast<double>(mainHeader.blockSize));
        expDic->setSignalSize(sigSize);
        expDic->setRealAtom(parm);
        sigAux = expDic->getRealAtom();
        initBlockSample = (int)(floor( static_cast<double>(sb[i].a) /
                                 static_cast<double>(mainHeader.blockSize) ) *
                                static_cast<double>(mainHeader.blockSize) );
        for (j=0;j<sigSize;j++)
        {
            if (initBlockSample+j >= mainHeader.signalSize) break;
            recSignal[initBlockSample+j] += norm*((CExpParm*)parm)->innerProd*sigAux[j];
        }
    }
    delete expDic;
    delete parm;
}

// void encodeRDClosedForm(char* InputFile, CFileGenData* genData, double rateTarget)
//
// {
//
// 	FILE* ioreport;
// 	ioreport = fopen("encode_report.out","w");
//
//     printf("-> Encoding Structure Book (RD Closed Form)\n");
//
//     int i,j;
//
//     //=====================================
//     // Loading Structure Books
//     printf("1) Loading the Structure Books\n");
//
//     CFileBlockRange* blockRange;
//     blockRange = new CFileBlockRange;
//     blockRange->setFileName("blockrange.dat");
//     blockRange->loadData();
//
//     strtSBBHeader sbbHeader;
//     sbbHeader = loadSBHeader("header.sbb");
//
//     CStructBook** structBook = NULL;
//     structBook = new CStructBook* [sbbHeader.numSignal];
//     for (i=0; i<sbbHeader.numSignal; i++)
//     {
//         structBook[i] = new CStructBookExp[sbbHeader.numBlock];
//     }
//
//     loadSB( InputFile,
//             blockRange,
//             sbbHeader,
//             structBook);
//
//     //=====================================
//     // Re-organizing structure book by linking atoms
//     // with continuity along the blocks
//     printf("2) Re-organizing structure book by linking atoms.\n");
//     // ----------------
//     CLinkStrBook* linkStrBookExpRef;
//     linkStrBookExpRef = new CLinkStrBookExp[sbbHeader.numSignal];
//     // -----------------
//     double** recSignal;
// 	recSignal = new double*[sbbHeader.numSignal];
// 	// -----------------
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		recSignal[i] = new double[sbbHeader.signalSize];
// 	}
// 	///////////////
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
//         printf(" - Set signal %d norm\n",i+1);
//         //printf("norm: %f\n",sbbHeader.norm[i]);
//         ((CLinkStrBookExp*)linkStrBookExpRef)[i].setNorm(sbbHeader.norm[i]);
//         printf(" - Load\n");
//         ((CLinkStrBookExp*)linkStrBookExpRef)[i].load(  structBook,
// 									sbbHeader,
// 									i);
// 		//printf(" - Print\n");
//         //((CLinkStrBookExp*)linkStrBookExpRef)[i].Print();
//         //printf(" - Synthesize signal(s)\n");
// 		//((CLinkStrBookExp*)linkStrBookExpRef)[i].synthSignal(recSignal[i],sbbHeader);
// 	}
// 	/*FILE* iorec;
// 	iorec = fopen("recsignal.out","w");
// 	for (i=0;i<sbbHeader.numSignal;i++)
//     {
//     	for (j=0;j<sbbHeader.signalSize;j++)
// 		{
// 			fprintf(iorec," %f",recSignal[i][j]);
// 		}
// 		fprintf(iorec,"\n");
// 	}
// 	fclose(iorec);
// 	*/
// 	/////////////////////
// 	// -----------------
// 	FILE* stream;
// 	stream = fopen("sblinked.out","w");
// 	CStructBook* sbLinked;
//     sbLinked = new CStructBookExp[sbbHeader.numSignal];
//     printf("3)  Convert linked structure book to common structure book.\n");
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
// 		cout << " - Convert Linked structure book to Structure book...." << endl;
// 		((CLinkStrBookExp*)linkStrBookExpRef)[i].linkStrBook2StrBook(&sbLinked[i],sbbHeader);
// 		cout << " - Set norm..." << endl;
// 		((CStructBookExp*)sbLinked)[i].setNorm(sbbHeader.norm[i]);
// 		((CStructBookExp*)sbLinked)[i].convertCoefNegativeToPositive();
// 		//cout << " - Synthesize signal" << endl;
// 		//synthSignalStrBook(recSignal[i],&sbLinked[i],sbbHeader);
// 		//cout << " - Print to screen..." << endl;
// 		//((CStructBookExp*)sbLinked)[i].printToScreen();
// 		cout << " - Print to file..." << endl;
// 		((CStructBookExp*)sbLinked)[i].saveStructBookASCII (stream );
// 	}
// 	fclose(stream);
//
// 	/*iorec = fopen("recsignal2.out","w");
// 	for (i=0;i<sbbHeader.numSignal;i++)
//     {
//     	for (j=0;j<sbbHeader.signalSize;j++)
// 		{
// 			fprintf(iorec," %f",recSignal[i][j]);
// 		}
// 		fprintf(iorec,"\n");
// 	}
// 	fclose(iorec);
// 	*/
//
//
//     /////////////////////////
//     // ----------------------
//     printf("4)  Rate distortion optimization stage.\n");
//     double max_coef, min_coef;
//     double max_rho, min_rho;
//     double max_phase, min_phase;
//     double L_amp, L_rho, L_phase;
//     int NAmpRange = 16;
//     double* ampRangeLimit;
//     ampRangeLimit = new double[NAmpRange+1];
//     double* ampBar;
//     ampBar = new double[NAmpRange];
//     CStructBook* sbLinkedSepByAmp;
//     //------------------------
//     for (i=0;i<sbbHeader.numSignal;i++)
//     {
//         cout << " - Find range of the parameters..." << endl;
//     	((CStructBookExp*)sbLinked)[i].findParmRange(min_coef, max_coef,
//     	   											 min_rho, max_rho,
//     	   											 min_phase, max_phase);
//     	cout << "#################################################################" << endl;
//     	cout << "## Range of the parameters                                     ##" << endl;
//     	cout << "=================================================================" << endl;
//     	cout <<  setw (27) << setfill(' ') << "min" <<  setw (20) << setfill(' ') << "max" << endl;
//     	cout << "Coef : " << setw (20) << setfill(' ') << min_coef << " ";
//     	cout << setw (20) << setfill(' ') << max_coef << endl;
//     	cout << "Rho  : " << setw (20) << setfill(' ') << min_rho << " ";
//     	cout << setw (20) << setfill(' ') << max_rho << endl;
//     	cout << "Phase: " << setw (20) << setfill(' ') << min_phase << " ";
//     	cout << setw (20) << setfill(' ') << max_phase << endl;
//     	cout << "#################################################################" << endl;
//     	// Parameters ranges excursion
//     	L_amp = max_coef - min_coef;
//     	L_rho  = max_rho - min_rho;
//     	L_phase= max_phase - min_phase;
//
//     	// Print the squared norm, distortion per atom
//     	printSBLinkedInfo(&sbLinked[i], sbbHeader,
//                    		L_amp, L_rho, L_phase);
//
//     	// ----------------------
// 		// Set the range limits
// 		fprintf(ioreport,"======================================\n");
// 		fprintf(ioreport,"= Amplitude Range Limits - Signal %d =\n",i+1);
// 		fprintf(ioreport,"--------------------------------------\n");
// 		ampRangeLimit[0]=-1;
// 		for (j=1;j<=NAmpRange;j++)
// 		{
// 			ampRangeLimit[j] = L_amp/(pow(2.0,static_cast<double>(NAmpRange-j)));
// 			fprintf(ioreport," Limit %4d: %12.9f - %12.9f\n",j,ampRangeLimit[j-1],ampRangeLimit[j]);
// 		}
// 		fprintf(ioreport,"===================================\n\n\n");
//
//         sbLinkedSepByAmp = new CStructBookExp[NAmpRange];
//
//     	// -----------------------
//     	cout << " - Separate atoms by amplitudes..." << endl;
//     	separateByAmp(&sbLinked[i], sbLinkedSepByAmp,NAmpRange, ampRangeLimit );
//
//     	for (j=0;j<NAmpRange;j++)
//     	{
//     	    fprintf(ioreport," ======>>>> Set of Atoms in Amplitude Range %d\n",j+1);
//     		((CStructBookExp*)sbLinkedSepByAmp)[j].saveStructBookASCII(ioreport);
//     		// Compute mean amplitude for each range
//     		ampBar[NAmpRange-1-j] = (3.0/4.0)*(L_amp/pow(2.0,static_cast<double>(j)));
// 		}
//
// 		//----------------------------------------------------------------
// 		// For a given lambda we will compute the optimum R-D solution
// 		double optR_amp,*optR_rho,*optR_phase;
// 		optR_rho = new double[NAmpRange];
// 		optR_phase = new double[NAmpRange];
// 		double *c1,*c2, *c3;
// 		c1 = new double[NAmpRange];
// 		c2 = new double[NAmpRange];
// 		c3 = new double[NAmpRange];
//
// 		cout << "- Compute RD contants..." << endl;
// 		calcRDConstant(	sbLinkedSepByAmp, ampBar, sbbHeader.Fs,
// 						c1, c2, c3, NAmpRange, sbbHeader.blockSize,
// 						L_amp,  L_rho, L_phase);
//
// 		fprintf(ioreport,"=========================================\n");
// 		fprintf(ioreport,"= RD Constants\n");
// 		fprintf(ioreport,"=========================================\n");
// 		fprintf(ioreport,"Range    C1     C2      C3        AmpBar\n");
// 		for (j=0;j<NAmpRange;j++)
//     	{
//     		fprintf(ioreport," %d %10.7f %10.7f %10.7f %10.7f \n",
//     		                  j+1,c1[j],c2[j],c3[j],ampBar[j]);
//     	}
// 		cout << "- Compute RD optimum solution..." << endl;
// 		calcRDOptimumSolution(	NAmpRange, L_amp, L_rho, L_phase,
// 								c1, c2, c3,
// 		                   		optR_amp, optR_rho, optR_phase);
//
// 		// ---------------------------
// 		delete [] ((CStructBookExp*)sbLinkedSepByAmp);
// 		delete [] optR_rho;
// 		delete [] optR_phase;
// 		delete [] c1;
// 		delete [] c2;
// 		delete [] c3;
//     }
//
//     // Close the file report
//     fclose(ioreport);
//
//     // Deallocate matrix and vectors
//     delete [] linkStrBookExpRef;
//     //--------------
//     for (i=0; i<sbbHeader.numSignal; i++)
//   	{
//     	delete [] ((CStructBookExp*)structBook[i]);
//   	}
//   	delete [] structBook;
//   	// --------------
// 	for (i=0;i<sbbHeader.numSignal;i++)
// 	{
// 		delete [] recSignal[i];
// 	}
//     delete [] recSignal;
//     // -----------------
//     delete blockRange;
//     if (sbbHeader.norm!=NULL) delete [] sbbHeader.norm;
//     // --------
//     delete [] ((CStructBookExp*)sbLinked);
//     delete [] ampRangeLimit;
//     delete [] ampBar;
// }
//
// // ==========================
// void printSBLinkedInfo(CStructBook* sbLinked,  strtSBBHeader sbbHeader,
//                       double L_amp,  double L_rho, double L_phase)
// {
//
//     strtContinuousExp* sb;
// 	int numElement;
// 	int j, initAtomSample, endAtomSample;
//
// 	CParameter* parm;
// 	parm = new CExpParm;
//
// 	CExpDictionary expDic;
// 	double R_amp = 16;
// 	double R_rho = 5;
// 	double R_phase = 5;
// 	double D;
// 	int signalSize = sbbHeader.blockSize;
// 	double Fs = sbbHeader.Fs;
// 	double norm;
//
// 	fstream file_op("sblinkedinfo.out",ios::out);
// 	numElement = ((CStructBookExp*)sbLinked)->getNumElement();
// 	sb = ((CStructBookExp*)sbLinked)->getStructBook();
// 	file_op << setw (5) << setfill(' ') << "NumEl"
// 		            <<  setw (20) << setfill(' ') << "SqrNorm"
// 		            <<  setw (20) << setfill(' ') << "SqrNormCalc"
// 					<<  setw (20) << setfill(' ') << "Func1"
// 			        <<  setw (20) << setfill(' ') << "Func2"
// 			        <<  setw (20) << setfill(' ') << "Func3"
// 			        <<  setw (20) << setfill(' ') << "AtomDistortion"
// 			        <<  setw (20) << setfill(' ') << "Coef"
// 			        <<  setw (20) << setfill(' ') << "Rho"
// 			        <<  setw (20) << setfill(' ') << "Xi"
// 			        <<  setw (20) << setfill(' ') << "Phase"
// 			        <<  setw (20) << setfill(' ') << "InitSample"
// 			        <<  setw (20) << setfill(' ') << "EndSample"
// 			        <<  setw (20) << setfill(' ') << "a"
// 			        <<  setw (20) << setfill(' ') << "b"
// 			        <<  setw (20) << setfill(' ') << "InitAtomSample"
// 			        <<  setw (20) << setfill(' ') << "EndAtomSample"
// 			        << endl;
//
// 	cout << "Delta rho: "   << L_rho/pow(2.0,R_rho)
// 	     << "; Delta phase: " << L_phase/pow(2.0,R_phase) << endl;
// 	for(j=0;j<numElement;j++)
// 	{
//
// 		((CExpParm*)parm)->setParm(sb[j]);
// 		initAtomSample = (int)( static_cast<double>(signalSize)*
// 			                  floor(static_cast<double>(sb[j].a)/static_cast<double>(signalSize)));
// 		endAtomSample = (int)( static_cast<double>(signalSize)*
// 			                  ceil(static_cast<double>(sb[j].b)/static_cast<double>(signalSize)))-1;
// 		((CExpParm*)parm)->a = sb[j].a - initAtomSample;
//
// 		((CExpParm*)parm)->b = sb[j].b - initAtomSample;
// 		expDic.setSignalSize(endAtomSample - initAtomSample + 1);
// 		expDic.setRealAtom(parm);
// 		norm = expDic.getRAtomNorm();
//
// 		D = ( ( (L_amp*L_amp)/pow(2.0,2.0*R_amp) )/12.0 ) +
// 			(((CExpParm*)parm)->innerProd*((CExpParm*)parm)->innerProd)*
// 			( ( ((L_phase*L_phase)/pow(2.0,2.0*R_phase))*
// 			    (L_rho/pow(2.0,R_rho))*
// 			    calcClosedFormDistortionFunc1(parm, Fs) ) +
// 			  ( ((L_phase*L_phase)/pow(2.0,2.0*R_phase))*
// 			    calcClosedFormDistortionFunc2(parm, Fs) ) +
// 			  ( (L_phase/pow(2.0,R_phase))*
// 			    (L_rho/pow(2.0,R_rho))*
// 			    calcClosedFormDistortionFunc3(parm, Fs) )
// 			 );
//
// 		file_op << setw (5) << setfill(' ') << j+1
// 				<<  setw (20) << setfill(' ') << calcClosedFormSqrNorm(parm,Fs)
// 				<<  setw (20) << setfill(' ') << norm*norm
// 				<<  setw (20) << setfill(' ') << calcClosedFormDistortionFunc1(parm, Fs)
// 				<<  setw (20) << setfill(' ') << calcClosedFormDistortionFunc2(parm, Fs)
// 				<<  setw (20) << setfill(' ') << calcClosedFormDistortionFunc3(parm, Fs)
// 				<<  setw (20) << setfill(' ') << D
// 				<<  setw (20) << setfill(' ') << ((CExpParm*)parm)->innerProd
// 				<<  setw (20) << setfill(' ') << ((CExpParm*)parm)->rho
// 				<<  setw (20) << setfill(' ') << ((CExpParm*)parm)->xi
// 				<<  setw (20) << setfill(' ') << ((CExpParm*)parm)->phase
// 				<<  setw (20) << setfill(' ') << sb[j].a
// 				<<  setw (20) << setfill(' ') << sb[j].b
// 				<<  setw (20) << setfill(' ') << ((CExpParm*)parm)->a
// 				<<  setw (20) << setfill(' ') << ((CExpParm*)parm)->b
// 				<<  setw (20) << setfill(' ') << initAtomSample
// 			    <<  setw (20) << setfill(' ') << endAtomSample
// 				<< endl;
// 	}
// 	file_op.close();
//
// }

//==================================================================================================
// Closed form functions
//---------------------------------
void calcRDOptimumSolution(	int NAmpRange, double L_amp,  double L_rho, double L_phase,
							double* c1, double* c2, double* c3,
		                   	double& optR_amp, double* optR_rho, double* optR_phase)
{
	double lambda=1.3449e-5;
	int i;
	int numElement;
	double a_phase, b_phase, c_phase;
	double a_rho, c_rho;
	double A,B,C;
	double x_rho,x_phase;
	CComplex root1,root2;


	optR_amp = 0.5 * (log(L_amp*L_amp*log(2.0)/(6*lambda))/log(2.0));
	cout << "optR_amp: " << optR_amp << endl;

	for (i=0;i<NAmpRange;i++)
	{
	    cout << "Range: " << i+1 << endl;
		a_phase = 2.0*log(2.0)*L_phase*L_phase*L_rho*c1[i];
		b_phase = 2.0*log(2.0)*L_phase*L_phase*c2[i];
		c_phase = log(2.0)*L_phase*L_rho*c3[i];
		a_rho = log(2.0)*L_phase*L_phase*L_rho*c1[i];
		c_rho = log(2.0)*L_phase*L_rho*c3[i];

		A = b_phase*a_rho;
		B = b_phase*c_rho;
		C = (a_phase-a_rho)*lambda;

		calcRoot2OrderBhaskara( A, B, C, root1, root2);

		/*if (A!=0.0)
		{
			x_phase = sqrt(-C/A);
		}
		else
		{
			x_phase = 1;
		}
		optR_phase[i] = - log(x_phase)/log(2.0);


		if ((a_rho*(x_phase*x_phase)+c_rho*x_phase)!=0.0)
		{
			x_rho = lambda/(a_rho*(x_phase*x_phase)+c_rho*x_phase);
		}
		else
		{
			x_rho =1;
		}
		optR_rho[i] = - log(x_rho)/log(2.0);

		cout << "X: "    << x_rho << " " << x_phase << endl;
		cout << "OptR: " << optR_rho[i] << " "<< optR_phase[i] << endl;
		*/

		cout << "A: " << A<< "; B: " << B<< "; C: " << C << endl;
		//calcRoot3OrderTartaglia(A, B, C, D,
		//						root1, root2, root3);
		cout << "Root 1 - real: " << root1.Real() << " - imag: " << root1.Imag() << endl;
		cout << "Root 2 - real: " << root2.Real() << " - imag: " << root2.Imag() << endl;
		//cout << "Root 2 - real: " << root3.Real() << " - imag: " << root3.Imag() << endl;
	}

}

void calcRDConstant(	CStructBook* sbLinkedSepByAmp,double* ampBar, double Fs,
						double* c1, double* c2, double* c3, int NAmpRange, int signalSize,
						double L_amp,  double L_rho, double L_phase)
{
	int i,j, initAtomSample, endAtomSample;
	strtContinuousExp* sb;
	int numElement;
	double c1_acum, c2_acum, c3_acum;
	CParameter* parm;
	parm = new CExpParm;

	CExpDictionary expDic;
	double R_amp = 16;
	double R_rho = 16;
	double R_phase = 16;
	double D;

	fstream file_op("rdconstant.out",ios::out);
	for (i=0;i<NAmpRange;i++)
	{
		file_op << "Range: " << i+1 << endl;
		file_op << setw (5) << setfill(' ') << "NumEl"
		            <<  setw (20) << setfill(' ') << "SqrNorm"
		            <<  setw (20) << setfill(' ') << "SqrNormCalc"
					<<  setw (20) << setfill(' ') << "Func1"
			        <<  setw (20) << setfill(' ') << "Func2"
			        <<  setw (20) << setfill(' ') << "Func3"
			        <<  setw (20) << setfill(' ') << "AtomDistortion"
			        <<  setw (20) << setfill(' ') << "Coef"
			        <<  setw (20) << setfill(' ') << "Rho"
			        <<  setw (20) << setfill(' ') << "Xi"
			        <<  setw (20) << setfill(' ') << "Phase"
			        <<  setw (20) << setfill(' ') << "InitSample"
			        <<  setw (20) << setfill(' ') << "EndSample"
			        << endl;
		sb = ((CStructBookExp*)sbLinkedSepByAmp)[i].getStructBook();
		numElement = ((CStructBookExp*)sbLinkedSepByAmp)[i].getNumElement();
		c1_acum=0.0;
		c2_acum=0.0;
		c3_acum=0.0;
		for (j=0;j<numElement;j++)
		{
			((CExpParm*)parm)->setParm(sb[j]);

			initAtomSample = (int)( static_cast<double>(signalSize)*
			                  floor(static_cast<double>(sb[j].a)/static_cast<double>(signalSize)));
			endAtomSample = (int)( static_cast<double>(signalSize)*
			                  ceil(static_cast<double>(sb[j].b)/static_cast<double>(signalSize)));
			((CExpParm*)parm)->a = ((CExpParm*)parm)->a % signalSize;

			((CExpParm*)parm)->b = sb[j].b - initAtomSample;
			expDic.setSignalSize(endAtomSample - initAtomSample);
			expDic.setRealAtom(parm);
			D = (L_amp/pow(2.0,R_phase))*(L_amp/pow(2.0,R_phase))/12 +
				(((CExpParm*)parm)->innerProd)*(((CExpParm*)parm)->innerProd)*
				( (L_phase/pow(2.0,R_phase))*
				  (L_phase/pow(2.0,R_phase))*
				  (L_phase/pow(2.0,R_rho))*
				  calcClosedFormDistortionFunc1(parm, Fs) +
				  (L_phase/pow(2.0,R_phase))*
				  (L_phase/pow(2.0,R_phase))*
				  calcClosedFormDistortionFunc2(parm, Fs) +
				  (L_phase/pow(2.0,R_phase))*
				  (L_phase/pow(2.0,R_rho))*
				  calcClosedFormDistortionFunc3(parm, Fs)
				 );

			file_op << setw (5) << setfill(' ') << j+1
			        <<  setw (20) << setfill(' ') << calcClosedFormSqrNorm(parm,Fs)
			        <<  setw (20) << setfill(' ') << (expDic.getRAtomNorm())*(expDic.getRAtomNorm())
					<<  setw (20) << setfill(' ') << calcClosedFormDistortionFunc1(parm, Fs)
			        <<  setw (20) << setfill(' ') << calcClosedFormDistortionFunc2(parm, Fs)
			        <<  setw (20) << setfill(' ') << calcClosedFormDistortionFunc3(parm, Fs)
			        <<  setw (20) << setfill(' ') << D
			        <<  setw (20) << setfill(' ') << ((CExpParm*)parm)->innerProd
			        <<  setw (20) << setfill(' ') << ((CExpParm*)parm)->rho
			        <<  setw (20) << setfill(' ') << ((CExpParm*)parm)->xi
			        <<  setw (20) << setfill(' ') << ((CExpParm*)parm)->phase
			        <<  setw (20) << setfill(' ') << ((CExpParm*)parm)->a
			        <<  setw (20) << setfill(' ') << ((CExpParm*)parm)->b
			        <<  setw (20) << setfill(' ') << initAtomSample
			        <<  setw (20) << setfill(' ') << endAtomSample
			        << endl;

			c1_acum += (((CExpParm*)parm)->innerProd)*calcClosedFormDistortionFunc1(parm, Fs);
			c2_acum += (((CExpParm*)parm)->innerProd)*calcClosedFormDistortionFunc2(parm, Fs);
			c3_acum += (((CExpParm*)parm)->innerProd)*calcClosedFormDistortionFunc3(parm, Fs);
		}
		cout << "Range: " << i+1 << endl;
		cout << ampBar[i] << "; " << c1_acum << "; " << c2_acum << "; " << c3_acum
			 << "; " << static_cast<double>(numElement) << endl;
		if (numElement!=0)
		{
			c1[i] =  c1_acum; //ampBar[i]*c1_acum / static_cast<double>(numElement);
			c2[i] =  c2_acum; //ampBar[i]*c2_acum; / static_cast<double>(numElement);
			c3[i] =  c3_acum; // ampBar[i]*c3_acum; / static_cast<double>(numElement);
		}
		else
		{
			c1[i] = 0.0;
			c2[i] = 0.0;
			c3[i] = 0.0;
		}
		cout <<  c1[i] << "; " << c2[i] << "; " << c3[i] << endl;

 	}
	delete parm;
	file_op.close();
}


double calcClosedFormSqrNorm(CParameter* parm, double Fs)
{
	/*double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho)*Fs;
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi*Fs;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    double sqrnorm=0.0;
    //int N = (int) round(xi*static_cast<double>(b-a)/(2*pi));
    int N = (int) ceil(xi*static_cast<double>(b-a)/(2*pi));
    double deltat = static_cast<double>(b-a)/Fs;
	*/
    double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho);

    double xi = ((CExpParm*)parm)->xi;
    double omega = xi;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    double sqrnorm=0.0;
    //int N = (int) round(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    int N = (int) ceil(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    double deltat = static_cast<double>((static_cast<double>(b)-static_cast<double>(a)));

    if (rho<0)
    {
    	rho = -rho;
   	}
   	else
   	{
   		//phase = xi*a + phase;
   	}



    if (deltat==0)
    {
    	return sqrnorm =1.0;
    }
    else if ( (omega!=0.0) && (rho!=0.0) )
    {
    	// Using integral interval with multiple integer of the atom frequency
    	/*return sqrnorm = (( cos(2*phase)+1 )*rho*rho - sin(2*phase)*rho*omega + omega*omega) *
							  (1 -exp(-(4*pi*static_cast<double>(N)*rho)/omega)) /
							  (4*rho*(rho*rho+omega*omega));
		*/
	    return sqrnorm = ((cos(2*phase) + 1)*rho*rho - sin(2*phase)*rho*omega + omega*omega)/(4*rho*rho*rho + 4*rho*omega*omega) -
	                      (rho*rho*(cos(2*phase + 2*deltat*omega) + 1) + omega*omega - rho*omega*sin(2*phase + 2*deltat*omega))*exp(-2*deltat*rho)/
	                      (4*rho*rho*rho + 4*rho*omega*omega);

    }
    else if ( (omega!=0.0) && (rho==0.0) )
    {
    	// Using integral interval with multiple integer of the atom frequency
    	//return sqrnorm = (pi*static_cast<double>(N))/omega;

    	return sqrnorm = (deltat/2) - (sin(2*phase) - sin(2*phase + 2*deltat*omega))/(4*omega);
    }
    else if ( (omega==0.0) && (rho!=0.0) )
    {
    	return sqrnorm = (cos(phase)*cos(phase))*(1-exp(-2*rho*deltat))/
		                      (2*rho);
    }
    else
    {
    	return sqrnorm = deltat*(cos(phase)*cos(phase));
    }
}
//-----------------------------------
double calcClosedFormDistortionFunc1(CParameter* parm, double Fs)
{
    /*double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho)*Fs;
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi*Fs;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    double f1=0.0;
    //int N = (int) round(xi*static_cast<double>(b-a)/(2*pi));
    int N = (int) ceil(xi*static_cast<double>(b-a)/(2*pi));
    double deltat = static_cast<double>(b-a)/Fs;
    */
    double f1=0.0;
    double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho);
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    //int N = (int) round(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    int N = (int) ceil(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    double deltat = static_cast<double>((static_cast<double>(b)-static_cast<double>(a)));

    double kappa2 = calcClosedFormSqrNorm(parm,Fs);


    if (deltat==0)
    {
    	return f1 = 0.0;
    }
    else if (rho!=0.0)
    {
    	return f1 = ( 1- (1+ 2*deltat*rho)*exp(-2*deltat*rho) )/ (8.0*rho*rho*kappa2);
    }
    else
    {
    	return f1 = (deltat*deltat)/(4.0*kappa2);
    }
    /*else if ( (omega!=0.0) && (rho!=0.0) )
    {
    	return f1 = (omega - exp(-(4*pi*static_cast<double>(N)*rho)/omega )*
    	                     (4*pi*static_cast<double>(N)*rho+omega) ) /
    	             (8*rho*rho*omega*kappa2);
    }
    else if ( (omega!=0.0) && (rho==0.0) )
    {
    	return f1 = (pi*static_cast<double>(N))*(pi*static_cast<double>(N)) /
    	             (omega*omega*kappa2);
    }
    else if ( (omega==0.0) && (rho!=0.0) )
    {
    	return f1 = (1-exp(-2*rho*deltat)*(2*rho*deltat+1))/
    	            (8*rho*rho*kappa2);
    }
    else
    {
    	return f1 = (deltat*deltat)/(4*kappa2);
    }*/
}
//-----------------------------------
double calcClosedFormDistortionFunc2(CParameter* parm, double Fs)
{
    /*double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho)*Fs;
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi*Fs;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    double f2=0.0;
    int N = (int) ceil(xi*static_cast<double>(b-a)/(2*pi));
    double deltat = static_cast<double>(b-a)/Fs;*/

    double f2=0.0;
    double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho);
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    //int N = (int) round(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    int N = (int) ceil(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    double deltat = static_cast<double>((static_cast<double>(b)-static_cast<double>(a)));

    double kappa2 = calcClosedFormSqrNorm(parm,Fs);

    if (deltat==0)
    {
    	return f2 = 0.0;
    }
    else if (rho!=0.0)
    {
    	return f2 = (1- exp(-2*deltat*rho))/(4*rho*rho*kappa2);
    }
    else
    {
    	return f2 = deltat/(2*kappa2);
    }
    /*else if ( (omega!=0.0) && (rho!=0.0) )
    {
    	return f2 = (1 - exp(-(4*pi*static_cast<double>(N)*rho)/omega) ) /
    	            (4*rho*kappa2);
    }
    else if ( (omega!=0.0) && (rho==0.0) )
    {
    	return f2 = (pi*static_cast<double>(N)) /
    	             (omega*kappa2);
    }
    else if ( (omega==0.0) && (rho!=0.0) )
    {
    	return f2 = (1-exp(-2*rho*deltat))/
    	            (4*rho*kappa2);
    }
    else
    {
    	return f2 = (deltat)/(2*kappa2);
    }*/
}
//-----------------------------------
double calcClosedFormDistortionFunc3(CParameter* parm, double Fs)
{
    /*double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho)*Fs;
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi*Fs;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    double f3=0.0;
    int N = (int) ceil(xi*static_cast<double>(b-a)/(2*pi));
    double deltat = static_cast<double>(b-a)/Fs;*/

    double f3=0.0;
    double innerProd = ((CExpParm*)parm)->innerProd;
    double rho = (((CExpParm*)parm)->rho);
    if (rho<0) rho = -rho;
    double xi = ((CExpParm*)parm)->xi;
    double omega = xi;
    double phase = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;
    //int N = (int) round(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    int N = (int) ceil(xi*static_cast<double>(static_cast<double>(b)-static_cast<double>(a))/(2*pi));
    double deltat = static_cast<double>((static_cast<double>(b)-static_cast<double>(a)));

    double kappa2 = calcClosedFormSqrNorm(parm,Fs);

    if (deltat==0)
    {
    	return f3 = 0.0;
    }
    else if ( (omega!=0.0) && (rho!=0.0) )
    {
    	// Using integral interval with multiple integer of the atom frequency
    	/*return f3 = (((1 - exp(-(4*pi*static_cast<double>(N)*rho)/omega) )/
    	             (4*(rho*rho+omega*omega)*(rho*rho+omega*omega)*kappa2))*
    	            (sin(2*phase)*(rho*rho-omega*omega)+2*omega*rho*cos(2*phase))) -
    	            (((pi*static_cast<double>(N)* exp(-(4*pi*static_cast<double>(N)*rho)/omega) )/
    	             (omega*(rho*rho+omega*omega)*kappa2))*
    	            (cos(2*phase)*omega + sin(2*phase)));
    	*/

    	 return f3 = -(2*omega)*(rho*rho*sin(2*phase + 2*deltat*omega)*exp(-2*deltat*rho) -
    	               omega*omega*sin(2*phase + 2*deltat*omega)*exp(-2*deltat*rho) -
    	               rho*rho*sin(2*phase) +
    	               omega*omega*sin(2*phase) +
    	               2*rho*omega*cos(2*phase + 2*deltat*omega)*exp(-2*deltat*rho) +
    	               2*deltat*omega*omega*omega*cos(2*phase + 2*deltat*omega)*exp(-2*deltat*rho) +
    	               2*deltat*rho*rho*rho*sin(2*phase + 2*deltat*omega)*exp(-2*deltat*rho) -
    	               2*rho*omega*cos(2*phase) +
    	               2*deltat*rho*rho*omega*cos(2*phase + 2*deltat*omega)*exp(-2*deltat*rho) +
    	               2*deltat*rho*omega*omega*sin(2*phase + 2*deltat*omega)*exp(-2*deltat*rho))/
    	               (4*(rho*rho + omega*omega)*(rho*rho + omega*omega)*kappa2);


    }
    else if ( (omega!=0.0) && (rho==0.0) )
    {
    	// Using integral interval with multiple integer of the atom frequency
    	/*return f3 = (-pi*static_cast<double>(N)*cos(2*phase)) /
    	             (omega*kappa2);*/

    	return f3 = -(2*omega)*(sin(2*phase) -
    	              sin(2*phase + 2*deltat*omega) +
    	              2*deltat*omega*cos(2*phase + 2*deltat*omega) ) /(4*omega*omega*kappa2);
    }
    else if ( (omega==0.0) && (rho!=0.0) )
    {
    	/*return f3 = (2*omega)*( sin(2*phase)*(1-(exp(-2*rho*deltat)*(2*rho*deltat+1)) ) )/
    	            (4*rho*rho*kappa2);
    	*/
    	return f3 =0.0;
    }
    else
    {
    	//return f3 = (2*omega)*(deltat*deltat*sin(2*phase))/(2*kappa2);
    	return f3=0.0;
    }
}
//----------------------------------
void calcRoot2OrderBhaskara(	double a, double b, double c,
								CComplex& root1, CComplex& root2)
{
	double tmp1,tmp2;

	tmp1 = b*b-4*a*c;
	if(tmp1<0)
	{
	    tmp2 = (-b/(2*a));
		root1.setReal(tmp2);
		root2.setReal(tmp2);
		//-----
		tmp2 = sqrt(fabs(tmp1))/(2*a);
		root1.setImag(tmp2);
		root2.setImag(-tmp2);
	}
	else
	{
		tmp2 = (-b/(2*a))+sqrt(tmp1)/(2*a);
		root1.setReal(tmp2);
		tmp2 = (-b/(2*a))-sqrt(tmp1)/(2*a);
		root2.setReal(tmp2);
		//-----
		tmp2 = sqrt(fabs(tmp1))/(2*a);
		root1.setImag(0.0);
		root2.setImag(0.0);
	}

}

void calcRoot3OrderTartaglia(	double a, double b, double c,double d,
								CComplex& root1, CComplex& root2, CComplex& root3)
{

	double A=b/a;
	double B=c/a;
	double C=d/a;
	double p=B-A*A/3;
	double q=C-A*B/3+2*A*A*A/27;
	double D=q*q/4+p*p*p/27;

	if(D<0.0)
	{
		double M=sqrt(-D);
		double r=sqrt(q*q/4+M*M);
		double t=acos(-(q/2.0)/r);
		double r1=2*(pow(r,(1.0/3.0)))*cos(t/3)-A/3;
		double r2=2*(pow(r,(1.0/3.0)))*cos((t+2*pi)/3)-A/3;
		double r3=2*(pow(r,(1.0/3.0)))*cos((t+4*pi)/3)-A/3;
		root1.setReal(round(r1*7E10)/7E10);
		root2.setReal(round(r2*7E10)/7E10);
		root3.setReal(round(r3*7E10)/7E10);
	}
	else
	{
		double u3=-q/2 + sqrt(D);
		double u;
		if (u3<0.0)
			u = -(pow((-u3),(1.0/3.0)));
		else
			u = (pow(u3,(1.0/3.0)));

		double v3 = -q/2 - sqrt(D);
		double v;
		if (v3<0)
			v = -(pow((-v3),(1.0/3.0)));
		else
			v = pow(v3,(1.0/3.0));


		double r1=u+v-A/3;
		r1=round(r1*10E7)/10E7;
		double Delta=(A+(r1))*(A+(r1))+4*C/r1;
		double Real=-(A+(r1))/2;
		Real=round(Real*10E7)/10E7;
		double K=fabs(Delta);
		double Imag=sqrt(K)/2;
		Imag=round(Imag*10E7)/10E7;
		if(Delta<0)
		{
			root1.setReal(r1);
			root2.setReal(Real);
			root2.setImag(Imag);
			root3.setReal(Real);
			root3.setImag(Imag);
		}
		else
		{
			double r20=Real+Imag;
			double r02=Real-Imag;
			root1.setReal(r1);
			root2.setReal(r20);
			root3.setReal(r02);
		}
	}
}

double linearQuant(double x, double step)
{
    double y;
    if (step==0.0)
        y = x;
    else
        y = floor( (x + step/2) / step) * step;

    return y;
}

int linearQuantIndex(double x, double step)
{
    int y;

    y = static_cast<int>(floor( (x + step/2) / step));

    return y;
}

double linearDeQuantIndex(int x, double step)
{
    double y;

    y =  static_cast<double>(x) * step;

    return y;
}


int geometricQuantIndex(double x, double step, int numPerOctave)
{
    int y;

    if (x==0.0)
    	y = 0;
    else
    	y =  static_cast<int>(round( static_cast<double>(numPerOctave)*log2(x/step) + 1.0));

    return y;
}

double geometricDeQuantIndex(int x, double step, int numPerOctave)
{
    double y;

	if (x==0)
		y =0.0;
	else
    	y =  step * pow(2.0, static_cast<double>(x-1)/static_cast<double>(numPerOctave)) ;

    return y;
}


