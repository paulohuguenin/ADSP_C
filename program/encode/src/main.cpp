
#include "filemgr.h"
#include "structbook.h"
#include "dictionary.h"
#include "datasignal.h"
#include "encode.h"

using namespace std;

int main(int argc, char *argv[])
{

    // Define the ".cfg" file path
    // _MAX_PATH defined in stdlib.h
    char InputFile[_MAX_PATH];
    double rateTarget;
    char auxStr[_MAX_PATH];
    if (argc > 1)
    {
        strcpy( InputFile, argv[1]);
        if (argc > 2)
        {
            strcpy( auxStr, argv[2]);
            rateTarget = atof(auxStr);
        }
        else
        {
            rateTarget = 0.0;
        }
    }
    else
    {
    //	strcpy( InputFile, "../../osc/x001.cfg");
    //	strcpy( InputFile, "MPP_PIANO_MEDIUM_A4_mono.wav");
    //	strcpy( InputFile, "../../audio/noise.");
        strcpy( InputFile, "pianoA3.wav");
        rateTarget = 0.0;
    }
    int flFixedLambdaVSUnifBlockRate = 1; // Default fixed lambda
// 	if (argc > 3)
// 	{
// 		strcpy( auxStr, argv[3]);
// 		flFixedLambdaVSUnifBlockRate = atof(auxStr);
// 	}
// 	flFixedLambdaVSUnifBlockRate = 1;

	int flLambdaBisecSearch=1;
	if (argc > 3)
	{
		strcpy( auxStr, argv[3]);
 		flLambdaBisecSearch = atoi(auxStr);
	}
	double initLambda = 0.0;
	if (argc > 4)
	{
		strcpy( auxStr, argv[4]);
 		initLambda = atof(auxStr);
	}

	/////////////////////////////

	// Loading Proceeding Switches File
    CFileEncode* proc;
    proc = new CFileEncode;
    proc->setFileName("panelEncode.dat");
    proc->loadData();

    // Loading block range input file
    CFileDecompBlockRange* blockRange;
    blockRange = new CFileDecompBlockRange;
    blockRange->setFileName("panelBlockRange.dat");
    blockRange->loadData();

    // Loading quantizers range file
    CFileRDBitRange* rdBitRange;
    rdBitRange = new CFileRDBitRange;
    rdBitRange->setFileName("panelRDBitRange.dat");
    rdBitRange->loadData();

    // Loading the dictionary file
    CFileDictionary* dicData;
    dicData = new CFileDictionary;
    dicData->setFileName("panelDictionary.dat");
    dicData->loadData();

    //=====================================
    // Loading Structure Books
    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();
	char sbbHFName[_MAX_PATH];
	char* pos;
    strcpy(sbbHFName, InputFile);
    pos = strrchr( sbbHFName, '.');
    sprintf(auxStr,"_b%d-%d_header.sbb",initBlock,finalBlock);
    strcpy( &pos[0], auxStr);
    printf("- Loading the Structure Book Header\n");
    // load header
    strtSBBHeader sbbHeader;
    sbbHeader = loadSBHeader(sbbHFName);

    //////////////////////////////////////////
    // Configuring and loading signals
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
	if (sbbHeader.type==3)
	{
		dataSignal->setNumSignal(1000);
		dataSignal->setSignalSize(sbbHeader.blockSize);
		dataSignal->setSamplingRate(44100);
	}
	dataSignal->setBlockSize(sbbHeader.blockSize);
	dataSignal->setBlockHop(sbbHeader.blockHop);
	dataSignal->setSignal();
	dataSignal->setNorm();

    /////////////////////////////
    // Procedures

    // -----------------------------------------------------------
    // Grouping blocks
    //
    if (proc->getGroupBlock()==1)
    {
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j,k;
			// Allocate structbook for exponential parameters
			CStructBook** structBook;
			structBook = new CStructBook* [sbbHeader.numSignal];
			CStructBook** structBookGroup;
			structBookGroup = new CStructBook* [sbbHeader.numSignal];
			int numBlockPerGroup = proc->getNumBlockPerGroup();
			int numBlockGroup = (int)ceil( 	static_cast<double>(sbbHeader.numBlock)/
											static_cast<double>(numBlockPerGroup) );
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				structBook[i] = new CStructBookExp[sbbHeader.numBlock];
				structBookGroup[i] = new CStructBookExp[numBlockGroup];
			}
			// load exponential structure book
			printf("- Loading the Structure Books\n");
			char sbbFName[_MAX_PATH];
			char aux[_MAX_PATH];
			char* pos;
			strcpy(sbbFName, InputFile);
			pos = strrchr( sbbFName, '.');
			int initBlock = (blockRange->getPInitBlock())[0];
			int finalBlock = (blockRange->getPFinalBlock())[0];
			sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
			strcpy( &pos[0], aux);
			loadSBExp( 	sbbFName,
						blockRange,
						sbbHeader,
						structBook);

			// Grouping blocks considering linked atoms
			printf("- Grouping blocks considering linked atoms\n");
			strtSBBHeader sbbHeaderGroup = sbbHeader;
			memcpy(sbbHeaderGroup.norm,sbbHeader.norm, sizeof(double)*sbbHeader.numSignal);
			sbbHeaderGroup.subBlockSize = sbbHeader.blockSize;
			sbbHeaderGroup.blockSize = numBlockPerGroup*sbbHeader.blockSize;
			sbbHeaderGroup.blockHop =sbbHeaderGroup.blockSize;
			sbbHeaderGroup.numBlock = numBlockGroup;

			for (i=0; i<sbbHeader.numSignal; i++)
			{
				groupBlockExp(	sbbHeader, structBook,
								sbbHeaderGroup, structBookGroup,
								numBlockPerGroup, i);
			}

			int iGrpBlock=0;
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				double norm=0.0;
				double norm_aux;
				for (j=0; j<sbbHeader.numBlock; j++)
				{
					norm_aux = ((CStructBookExp*)structBook[i])[j].getNorm();
					cout << "blocknorm: " << norm_aux << "  -- block group  " << j+1 << endl;
					norm = norm + norm_aux*norm_aux;
					k = j % numBlockPerGroup;
					if ( (k==numBlockPerGroup-1) || (j==sbbHeader.numBlock-1) )
					{
						cout << "block group norm: " << sqrt(norm) << "  -- block group  " << iGrpBlock+1 << endl;
						((CStructBookExp*)structBookGroup[i])[iGrpBlock].setNorm(sqrt(norm));
						norm =0.0;
						iGrpBlock++;
					}
				}
			}

			printf("- Save grouped structure book.\n");
			char sbbGroupFName[_MAX_PATH];
			char auxGroup[_MAX_PATH];
			char* posGroup;
			strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp_header.sbb");
			strcpy( &posGroup[0], auxGroup);

			saveSBHeaderGrp(sbbGroupFName, sbbHeaderGroup);

			strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);

			saveSBExpGrp(sbbGroupFName, sbbHeaderGroup, structBookGroup);


			printf("- Reconstruct signal based on the grouped structure book.\n");
			///////////
			double** recSignal;
			recSignal = new double*[sbbHeaderGroup.numSignal];
			double* blockSignal;
			blockSignal = new double[sbbHeaderGroup.blockSize];
			strtContinuousExp* sb;
			int signalSize = (int) ceil( static_cast<double>(sbbHeaderGroup.signalSize)/
										 static_cast<double>(sbbHeaderGroup.blockSize))  *  sbbHeaderGroup.blockSize;


			for (i=0; i<sbbHeaderGroup.numSignal; i++)
			{
				recSignal[i] = new double[signalSize];
				for (j=0; j<sbbHeaderGroup.numBlock; j++)
				{
					sb = ((CStructBookExp*)structBookGroup[i])[j].getStructBook();
					int sbNumElement = ((CStructBookExp*)structBookGroup[i])[j].getNumElement();
					synthSignalSBExp(	blockSignal,
							 			sbbHeaderGroup.norm[i],
							 			sbbHeaderGroup.blockSize,
							 			sb,
							 			sbNumElement);
					int initSample = j*sbbHeaderGroup.blockSize;
					cout << sbbHeaderGroup.norm[i] << " " <<
							sbbHeaderGroup.blockSize << " " <<
							sbNumElement << " " <<
							initSample << " " <<
							signalSize << " " << endl;
					memcpy(&recSignal[i][initSample],blockSignal, sizeof(double)*sbbHeaderGroup.blockSize);
				}
			}


			//////////////
			FILE* iorec;
			iorec = fopen("recsignal.out","w");
			for (i=0; i<sbbHeaderGroup.numSignal; i++)
			{
				for (j=0; j<signalSize; j++)
				{
					fprintf(iorec," %f",recSignal[i][j]);
				}
				fprintf(iorec,"\n");
			}
			fclose(iorec);

			CDataSignal* outSignal;
			char fileName[_MAX_PATH];
			strcpy(fileName, InputFile);
			pos = strrchr( fileName, '.');
			sprintf(aux,"_grprec.wav");
			strcpy( &pos[0],aux );
			outSignal = new CAudioSignal;
			outSignal->setFileName(fileName);
			outSignal->setNumSignal(sbbHeaderGroup.numSignal);
			outSignal->setSignalSize(sbbHeaderGroup.signalSize);
			outSignal->setSamplingRate(sbbHeaderGroup.Fs);
			outSignal->setSignal(recSignal);
			outSignal->setNorm();
			outSignal->saveSignal();
			delete outSignal;

			// Writing the Main Header
			FILE* iosba;
			strcpy(sbbGroupFName, InputFile);
			pos = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sba");
			strcpy( &pos[0], aux);
			iosba = fopen(sbbGroupFName,"w");
			fprintf(iosba,"Sign. Type :          %5i\n", sbbHeaderGroup.type);
			fprintf(iosba,"Dict. Type :          %5i\n", sbbHeaderGroup.dicType);
			fprintf(iosba,"No. Signals:          %5i\n", sbbHeaderGroup.numSignal);
			fprintf(iosba,"Signal Size:       %8i\n",    sbbHeaderGroup.signalSize);
			fprintf(iosba,"Block Hop:            %5i\n", sbbHeaderGroup.blockHop);
			fprintf(iosba,"Block Size :          %5i\n", sbbHeaderGroup.blockSize);
			fprintf(iosba,"Samp. Freq :     %10.2f\n",   sbbHeaderGroup.Fs);

			for (i=0; i<sbbHeaderGroup.numSignal; i++)
			{
				// Writing the Signal Header
				fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
				fprintf(iosba,"Signal:               %5i\n",i+1);
				fprintf(iosba,"Norm:            %10.5f\n",sbbHeaderGroup.norm[i]);
				for (j=0; j<sbbHeaderGroup.numBlock; j++)
				{

					double initBlockNorm = ((CStructBookExp*)structBookGroup[i])[j].getNorm();

					fprintf(iosba,"--------------------------------------------------------------\n");
					fprintf(iosba,"Block:                %5i\n",j+1);
					fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
					fprintf(iosba,"No.         Coef.           Decaying        Freq            Phase          Ti    Tf    PrevAtom NextAtom     Orig\n");

					sb = ((CStructBookExp*)structBookGroup[i])[j].getStructBook();
					int sbNumElement = ((CStructBookExp*)structBookGroup[i])[j].getNumElement();
					for (k=0;k<sbNumElement;k++)
					{
						fprintf ( iosba,"%10i %15.8f %15.8f %15.8f %15.8f %5i %5i    %5i    %5i    %5i \n",
								  k+1,
								  sb[k].innerProduct,
								  sb[k].rho,
								  sb[k].xi,
								  sb[k].phase,
								  sb[k].a,
								  sb[k].b,
								  sb[k].prevAtom+1,
								  sb[k].nextAtom+1,
								  sb[k].origAtomIndex+1);
					}
					fprintf(iosba,"99999\n");
				}
				fprintf(iosba,"88888\n");
			}
			fprintf(iosba,"77777\n");
			fflush(iosba);
			fclose(iosba);


			for (i=0; i<sbbHeaderGroup.numSignal; i++)
			{
				delete [] recSignal[i];
			}
			delete [] recSignal;
			delete [] blockSignal;


			// deallocate
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				delete [] (CStructBookExp*)structBook[i];
				delete [] (CStructBookExp*)structBookGroup[i];
			}
			delete [] structBook;
			delete [] structBookGroup;
		}
		if (sbbHeader.norm!=NULL) delete [] sbbHeader.norm;
	}

    // -----------------------------------------------------------
    //  Compute RD operational curves for each signal block
    //
    // if ( (genData->getProcType()==3)  ||
//          (genData->getProcType()==0) )
	if (proc->getRDOpCurve()==1)
    {
    	//=====================================
		// Loading Structure Books
		printf("- Loading the Structure Book Header\n");
    	// load header
		strtSBBHeader sbbHeader;
		char sbbGroupFName[_MAX_PATH];
        char auxGroup[_MAX_PATH];
        char* posGroup;
        strcpy(sbbGroupFName, InputFile);
        posGroup = strrchr( sbbGroupFName, '.');
        sprintf(auxGroup,"_strbookgrp_header.sbb");
        strcpy( &posGroup[0], auxGroup);
		sbbHeader = loadSBHeaderGrp(sbbGroupFName);
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;
            // Allocate structbook for exponential parameters
            CStructBook** structBook;
            structBook = new CStructBook* [sbbHeader.numSignal];

            for (i=0; i<sbbHeader.numSignal; i++)
            {
                structBook[i] = new CStructBookExp[sbbHeader.numBlock];
            }

            // load exponential structure book
            printf("- Loading the Structure Books\n");
            strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);
            loadSBExpGrp(  sbbGroupFName,
                        blockRange,
                        sbbHeader,
                        structBook);

			// RD optimization
			printf("- Rate-distortion optimization.\n");

			optimizeRDExp(	rdBitRange, dataSignal,
							sbbHeader, structBook,
					        dicData->getNumFreq(),
                            InputFile);

			// Save Operational Curves in binary format
			strtOpCurveExp* opCurve;
			int numElement;
			int numwritten;
            double dummydouble;
            strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_opcurve.bin");
			strcpy( &posGroup[0], auxGroup);
			FILE* io_opcurve;
			io_opcurve = fopen(sbbGroupFName,"wb");
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

					opCurve = ((CStructBookExp*)structBook[i])[j].getOpCurve();
					numElement = ((CStructBookExp*)structBook[i])[j].getNumElemOpCurve();
					fwrite ( &numElement, sizeof ( int ), 1, io_opcurve );
					numwritten = fwrite ( opCurve,
										  sizeof ( strtOpCurveExp ),
										  numElement,
										  io_opcurve);

				}
			}
			fclose(io_opcurve);

			// deallocate
            for (i=0; i<sbbHeader.numSignal; i++)
            {
                delete [] (CStructBookExp*)structBook[i];
            }
            delete [] structBook;
		}
		if (sbbHeader.norm!=NULL) delete [] sbbHeader.norm;
    }
    // -----------------------------------------------------------
    //  Compute RD operational curves for each signal block
    //  considering the separation of atoms among amplitude ranges
	if (proc->getRDOpCurveAmpRange()==1)
    {
    	//=====================================
		// Loading Structure Books
		printf("- Loading the Structure Book Header\n");
    	// load header
		strtSBBHeader sbbHeader;
		char sbbGroupFName[_MAX_PATH];
        char auxGroup[_MAX_PATH];
        char* posGroup;
        strcpy(sbbGroupFName, InputFile);
        posGroup = strrchr( sbbGroupFName, '.');
        sprintf(auxGroup,"_strbookgrp_header.sbb");
        strcpy( &posGroup[0], auxGroup);
		sbbHeader = loadSBHeaderGrp(sbbGroupFName);
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;
            // Allocate structbook for exponential parameters
            CStructBook** structBook;
            structBook = new CStructBook* [sbbHeader.numSignal];

            for (i=0; i<sbbHeader.numSignal; i++)
            {
                structBook[i] = new CStructBookExp[sbbHeader.numBlock];
            }

            // load exponential structure book
            printf("- Loading the Structure Books\n");
            strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);
            loadSBExpGrp(  sbbGroupFName,
                        blockRange,
                        sbbHeader,
                        structBook);

			// RD optimization
			printf("- Rate-distortion optimization.\n");

			optimizeRDExpAmpRangeOpCurve(	rdBitRange, dataSignal,
											sbbHeader, structBook,
											dicData->getNumFreq(),
                                            InputFile);

			// deallocate
            for (i=0; i<sbbHeader.numSignal; i++)
            {
                delete [] (CStructBookExp*)structBook[i];
            }
            delete [] structBook;
		}
		if (sbbHeader.norm!=NULL) delete [] sbbHeader.norm;
    }
    // -----------------------------------------------------
    //  Compute table of mean atom distortion for different
    //  number of amplitude bits
    //
	if (proc->getRDTableDf()==1)
    {
    	//=====================================
		// Loading Structure Books
		printf("- Loading the Structure Book Header\n");
    	// load header
		strtSBBHeader sbbHeader;
		char sbbGroupFName[_MAX_PATH];
        char auxGroup[_MAX_PATH];
        char* posGroup;
        strcpy(sbbGroupFName, InputFile);
        posGroup = strrchr( sbbGroupFName, '.');
        sprintf(auxGroup,"_strbookgrp_header.sbb");
        strcpy( &posGroup[0], auxGroup);
		sbbHeader = loadSBHeaderGrp(sbbGroupFName);
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;
			// Allocate structbook for exponential parameters
			CStructBook** structBook;
			structBook = new CStructBook* [sbbHeader.numSignal];
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				structBook[i] = new CStructBookExp[sbbHeader.numBlock];
			}
			// load exponential structure book
			printf("- Loading the Structure Books\n");
			strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);
			loadSBExpGrp( sbbGroupFName,
                        blockRange,
                        sbbHeader,
                        structBook);

			// RD optimization
			printf("- Rate-distortion optimization.\n");

			// optimizeRDExpByAmp(	rdBitRange, dataSignal,
// 								sbbHeaderGroup, structBookGroup,
// 						  		dicData->getNumFreq());
			calcTableRDExpByAmp(rdBitRange,
								dataSignal,
								sbbHeader, structBook,
						  		dicData->getNumFreq(),
                                InputFile);

			// deallocate
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				delete [] (CStructBookExp*)structBook[i];
			}
			delete [] structBook;
		}
		if (sbbHeader.norm!=NULL) delete [] sbbHeader.norm;
    }
    // -----------------------------------------------------
    //  Encode using RD operational curves
    //
    if (proc->getEncOpCurve()==1)
    {
    	//=====================================
		// Loading Structure Books
		printf("- Loading the Structure Book Header\n");
    	// load header
		strtSBBHeader sbbHeader;
		char sbbGroupFName[_MAX_PATH];
        char auxGroup[_MAX_PATH];
        char* posGroup;
        strcpy(sbbGroupFName, InputFile);
        posGroup = strrchr( sbbGroupFName, '.');
        sprintf(auxGroup,"_strbookgrp_header.sbb");
        strcpy( &posGroup[0], auxGroup);
		sbbHeader = loadSBHeaderGrp(sbbGroupFName);
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;
            // Allocate structbook for exponential parameters
            CStructBook** structBook;
            structBook = new CStructBook* [sbbHeader.numSignal];

            for (i=0; i<sbbHeader.numSignal; i++)
            {
                structBook[i] = new CStructBookExp[sbbHeader.numBlock];
            }
            // load exponential structure book
            printf("- Loading the Structure Books\n");
            strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);
            loadSBExpGrp(  sbbGroupFName,
                        blockRange,
                        sbbHeader,
                        structBook);

			// ENCODE
			printf("- Encode.\n");
			encodeSBExp(InputFile,
						rateTarget,
						sbbHeader,
                    	structBook,
                    	flFixedLambdaVSUnifBlockRate);
            // deallocate
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				delete [] (CStructBookExp*)structBook[i];
			}
			delete [] structBook;
    	}

    }
    // -----------------------------------------------------
    //  Encode using Df tables computed with atoms
    //  separated in amplitude ranges
    //
    if (proc->getEncAmpRange()==1)
    {
    	//=====================================
		// Loading Structure Books
		printf("- Loading the Structure Book Header\n");
    	// load header
		strtSBBHeader sbbHeader;
		char sbbGroupFName[_MAX_PATH];
        char auxGroup[_MAX_PATH];
        char* posGroup;
        strcpy(sbbGroupFName, InputFile);
        posGroup = strrchr( sbbGroupFName, '.');
        sprintf(auxGroup,"_strbookgrp_header.sbb");
        strcpy( &posGroup[0], auxGroup);
		sbbHeader = loadSBHeaderGrp(sbbGroupFName);
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;
            // Allocate structbook for exponential parameters
            CStructBook** structBook;
            structBook = new CStructBook* [sbbHeader.numSignal];

            for (i=0; i<sbbHeader.numSignal; i++)
            {
                structBook[i] = new CStructBookExp[sbbHeader.numBlock];
            }
            // load exponential structure book
            printf("- Loading the Structure Books\n");
            strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);
            loadSBExpGrp(  sbbGroupFName,
                        blockRange,
                        sbbHeader,
                        structBook);
            // Load DfTable
            printf("- Loading the Df Table\n");
            CFileDfTable* DfTableFile;
            DfTableFile = new CFileDfTable;
            DfTableFile->setFileName("./DFTABLE/4096/dftablebyamp.bin");
            DfTableFile->loadData();

			// ENCODE
			printf("- Encode.\n");
			//encodeSBExpAmpRange(InputFile,
			encodeSBExpAmpRangeBisecSearch(	InputFile,
											rateTarget,
											sbbHeader,
											structBook,
											DfTableFile,
											dicData->getNumFreq(),
											dataSignal,
											flLambdaBisecSearch,
											initLambda);
            // deallocate
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				delete [] (CStructBookExp*)structBook[i];
			}
			delete [] structBook;

			delete DfTableFile;
    	}

    }


    // -----------------------------------------------------
    //  Encode using Df tables computed with atoms
    //  separated in amplitude ranges
    //
    if (proc->getEncAmpRangeArithCod()==1)
    {
    	//=====================================
		// Loading Structure Books
		printf("- Loading the Structure Book Header\n");
    	// load header
		strtSBBHeader sbbHeader;
		char sbbGroupFName[_MAX_PATH];
        char auxGroup[_MAX_PATH];
        char* posGroup;
        strcpy(sbbGroupFName, InputFile);
        posGroup = strrchr( sbbGroupFName, '.');
        sprintf(auxGroup,"_strbookgrp_header.sbb");
        strcpy( &posGroup[0], auxGroup);
		sbbHeader = loadSBHeaderGrp(sbbGroupFName);
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;
            // Allocate structbook for exponential parameters
            CStructBook** structBook;
            structBook = new CStructBook* [sbbHeader.numSignal];

            for (i=0; i<sbbHeader.numSignal; i++)
            {
                structBook[i] = new CStructBookExp[sbbHeader.numBlock];
            }
            // load exponential structure book
            printf("- Loading the Structure Books\n");
            strcpy(sbbGroupFName, InputFile);
			posGroup = strrchr( sbbGroupFName, '.');
			sprintf(auxGroup,"_strbookgrp.sbb");
			strcpy( &posGroup[0], auxGroup);
            loadSBExpGrp(  sbbGroupFName,
                        blockRange,
                        sbbHeader,
                        structBook);
            // Load DfTable
            char* dfTableFname = rdBitRange->getDfTableFName();
            printf("- Loading the Df Table: %s \n", dfTableFname);
            CFileDfTable* DfTableFile;
            DfTableFile = new CFileDfTable;
            //DfTableFile->setFileName("./DFTABLE/4096/dftablebyamp.bin");
            DfTableFile->setFileName(dfTableFname);
            DfTableFile->loadData();

			// ENCODE
			printf("- Encode.\n");
			//encodeSBExpAmpRange(InputFile,
			// encodeSBExpAmpRangeBisecSearchArithCod(	InputFile,
// 													rateTarget,
// 													sbbHeader,
// 													structBook,
// 													DfTableFile,
// 													dicData,
// 													dataSignal,
// 													flLambdaBisecSearch,
// 											        initLambda);
			encodeSBExpAmpRangeArithCod(	InputFile,
											rateTarget,
											sbbHeader,
											structBook,
											DfTableFile,
											dicData,
											dataSignal,
											initLambda);
            // deallocate
			for (i=0; i<sbbHeader.numSignal; i++)
			{
				delete [] (CStructBookExp*)structBook[i];
			}
			delete [] structBook;

			delete DfTableFile;
    	}

    }

    // -----------------------------------------------------
    //  Decode using Df tables computed with atoms
    //  separated in amplitude ranges (arithmetic coding)
    //
    if (proc->getDecAmpRangeArithCod()==1)
    {
		//=====================================
		// Exponential dictionary
    	if (sbbHeader.dicType==1)
    	{
    		int i,j;

            // Load DfTable
            char* dfTableFname = rdBitRange->getDfTableFName();
            printf("- Loading the Df Table: %s \n", dfTableFname);
            CFileDfTable* DfTableFile;
            DfTableFile = new CFileDfTable;
            //DfTableFile->setFileName("./DFTABLE/4096/dftablebyamp.bin");
            DfTableFile->setFileName(dfTableFname);
            DfTableFile->loadData();

			// ENCODE
			printf("- Decode.\n");

			decodeSBExpAmpRangeArithCod(InputFile,
										DfTableFile,
										dicData);
            // deallocate
			delete DfTableFile;
    	}
    }

    ///////////////////
    // Delete section
    ///////////////////
    if (sbbHeader.type==1)
	{
		delete (CComtradeSignal*)dataSignal;
	}
	if (sbbHeader.type==2)
	{
		delete (CAudioSignal*)dataSignal;
	}
	if (sbbHeader.type==3)
	{
		delete (CNoiseSignal*)dataSignal;
	}
	///////
	delete proc;
    delete blockRange;
    delete rdBitRange;
    delete dicData;
  	//////
    return 0;
}
