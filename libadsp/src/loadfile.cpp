
#include "loadfile.h"

strtSBBHeader loadSBHeader(char* file)
{
    // Read the header
    strtSBBHeader sbbHeader;
    int dummyint;
    double dummydouble;
    FILE* iosbb;

    if ( (iosbb = fopen(file,"rb")) ==NULL)
    {
        printf("- File %s does not exist\n",file);
        exit(0);
    }

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.type = dummyint;
    printf("Signal type: %d \n",sbbHeader.type);

/*    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.dicType = dummyint;
    printf("Dictionary type: %d \n",sbbHeader.dicType);
*/
    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.numSignal = dummyint;
    printf("Number of signals: %d \n",sbbHeader.numSignal);

    if (sbbHeader.numSignal!=0) sbbHeader.norm = new double[sbbHeader.numSignal];

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.signalSize = dummyint;
    printf("Signal size: %d \n",sbbHeader.signalSize);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.blockHop = dummyint;
    printf("Block hop: %d \n",sbbHeader.blockHop);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.blockSize = dummyint;
    printf("Signal size: %d \n",sbbHeader.blockSize);

    fread( &dummydouble, sizeof( double ), 1, iosbb );
    sbbHeader.Fs =dummydouble;
    printf("Sampling freq.: %f \n",sbbHeader.Fs);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.initBlock = dummyint;
    printf("Initial block: %d \n",sbbHeader.initBlock);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.finalBlock = dummyint;
    printf("Final block: %d \n",sbbHeader.finalBlock);

    // Allocating memory for Structure Book Set
    sbbHeader.numBlock = (int)ceil((double)sbbHeader.signalSize/(double)sbbHeader.blockHop);
    printf("Number of blocks: %d \n",sbbHeader.numBlock);

    fclose(iosbb);

    return sbbHeader;
}

strtSBBHeader loadSBHeaderGrp(char* file)
{
    // Read the header
    strtSBBHeader sbbHeader;
    int dummyint;
    double dummydouble;
    FILE* iosbb;

    if ( (iosbb = fopen(file,"rb")) ==NULL)
    {
        printf("- File %s does not exist\n",file);
        exit(0);
    }

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.type = dummyint;
    printf("Signal type: %d \n",sbbHeader.type);

/*    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.dicType = dummyint;
    printf("Dictionary type: %d \n",sbbHeader.dicType);
*/
    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.numSignal = dummyint;
    printf("Number of signals: %d \n",sbbHeader.numSignal);

    if (sbbHeader.numSignal!=0) sbbHeader.norm = new double[sbbHeader.numSignal];

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.signalSize = dummyint;
    printf("Signal size: %d \n",sbbHeader.signalSize);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.blockHop = dummyint;
    printf("Block hop: %d \n",sbbHeader.blockHop);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.blockSize = dummyint;
    printf("Signal size: %d \n",sbbHeader.blockSize);

    fread( &dummydouble, sizeof( double ), 1, iosbb );
    sbbHeader.Fs =dummydouble;
    printf("Sampling freq.: %f \n",sbbHeader.Fs);

    fread( &dummyint, sizeof( int ), 1, iosbb );
    sbbHeader.subBlockSize = dummyint;
    printf("Sub block size: %d \n",sbbHeader.subBlockSize);

    // Allocating memory for Structure Book Set
    sbbHeader.numBlock = (int)ceil((double)sbbHeader.signalSize/(double)sbbHeader.blockHop);
    printf("Number of blocks: %d \n",sbbHeader.numBlock);

    fclose(iosbb);

    return sbbHeader;
}





void loadSBExp(	char* sbbFName,
				CFileDecompBlockRange* blockRange,
				strtSBBHeader sbbHeader,
				CStructBook** structBook)
{
    int i,j;
    FILE* iosbb;
    int dummyint;
    double dummydouble;
    strtContinuousExp dummySBExp;

    int nSignal,nBlock, nSBElement;
    CParameter* expParm;
    expParm= new CExpParm;
    //cout << "Num block range: "<< blockRange->getNumRange() << endl;
    for(i=0;i < blockRange->getNumRange();i++)
    {
        //printf("iBlockRange: %d\n",i);

        if ( (iosbb = fopen(sbbFName,"rb")) ==NULL)
        {
            printf("- File %s does not exist\n",sbbFName);
            exit(0);
        }

        printf("- Loading File %s.\n",sbbFName);

        int numread;
        //while(!feof(iosbb))
        while(1)
        {
            numread = fread( &nSignal, sizeof( int ), 1, iosbb );
            //printf("Signal: %d \n",nSignal);
            if (nSignal==77777)
            {
                printf("- Loading complete!!!\n");
                break;
            }
            numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
            sbbHeader.norm[nSignal-1] = dummydouble;
            //printf("Signal norm: %f \n",sbbHeader.norm[nSignal-1]);
            //while(!feof(iosbb))
            while(1)
            {
                numread = fread( &nBlock, sizeof( int ), 1, iosbb );
                //printf("Block: %d \n",nBlock);
                if (nBlock==88888) break;
                numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
                //printf("Block norm: %f \n",dummydouble);
                ((CStructBookExp*)structBook[nSignal-1])[nBlock-1].setNorm(dummydouble);
                //while (!feof(iosbb))
                while(1)
                {
                    numread = fread( &nSBElement, sizeof( int ), 1, iosbb );
                    //printf("Element: %d \n",nSBElement);
                    //printf("Numread: %d \n",numread);
                    if (nSBElement==99999) break;
                    numread = fread( &dummySBExp, sizeof( strtContinuousExp ), 1, iosbb );
                    ((CStructBookExp*)structBook[nSignal-1])[nBlock-1].addElement(dummySBExp);
                    //((CStructBookExp*)structBook[nSignal-1])[nBlock-1].printElementToScreen(nSBElement-1);
                    //if (feof(iosbb)) cout << "Fim de arquivo" << endl;
                    //if (ferror(iosbb)) cout << "Erro" << endl;
                }
            }
        }
        printf("- Loading has finished!!!\n");
        fclose(iosbb);
    }
    delete expParm;
}

// void loadSBExpGrp(	char* sbbFName,
// 					CFileBlockRange* blockRange,
// 					strtSBBHeader sbbHeader,
// 					CStructBook** structBook)
// {
//     int i,j;
//     FILE* iosbb;
//     int dummyint;
//     double dummydouble;
//     strtContinuousExp dummySBExp;
//
//     int nSignal,nBlock, nSBElement;
//     CParameter* expParm;
//     expParm= new CExpParm;
//     //cout << "Num block range: "<< blockRange->getNumRange() << endl;
//     for(i=0;i < blockRange->getNumRange();i++)
//     {
//         //printf("iBlockRange: %d\n",i);
//
//         if ( (iosbb = fopen(sbbFName,"rb")) ==NULL)
//         {
//             printf("- File %s does not exist\n",sbbFName);
//             exit(0);
//         }
//
//         printf("- Loading File %s.\n",sbbFName);
//
//         int numread;
//         //while(!feof(iosbb))
//         while(1)
//         {
//             numread = fread( &nSignal, sizeof( int ), 1, iosbb );
//             //printf("Signal: %d \n",nSignal);
//             if (nSignal== -3)
//             {
//                 printf("- Loading complete!!!\n");
//                 break;
//             }
//             numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
//             sbbHeader.norm[nSignal-1] = dummydouble;
//             //printf("Signal norm: %f \n",sbbHeader.norm[nSignal-1]);
//             //while(!feof(iosbb))
//             while(1)
//             {
//                 numread = fread( &nBlock, sizeof( int ), 1, iosbb );
//                 //printf("Block: %d \n",nBlock);
//                 if (nBlock== -2) break;
//                 numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
//                 //printf("Block norm: %f \n",dummydouble);
//                 ((CStructBookExp*)structBook[nSignal-1])[nBlock-1].setNorm(dummydouble);
//                 //while (!feof(iosbb))
//                 while(1)
//                 {
//                     numread = fread( &nSBElement, sizeof( int ), 1, iosbb );
//                     //printf("Element: %d \n",nSBElement);
//                     //printf("Numread: %d \n",numread);
//                     if (nSBElement== -1) break;
//                     numread = fread( &dummySBExp, sizeof( strtContinuousExp ), 1, iosbb );
//                     ((CStructBookExp*)structBook[nSignal-1])[nBlock-1].addElement(dummySBExp);
//                     //((CStructBookExp*)structBook[nSignal-1])[nBlock-1].printElementToScreen(nSBElement-1);
//                     //if (feof(iosbb)) cout << "Fim de arquivo" << endl;
//                     //if (ferror(iosbb)) cout << "Erro" << endl;
//                 }
//             }
//         }
//         printf("- Loading has finished!!!\n");
//         fclose(iosbb);
//     }
//     delete expParm;
// }

void loadSBExpGrp(	char* sbbFName,
					CFileDecompBlockRange* blockRange,
					strtSBBHeader sbbHeader,
					CStructBook** structBook)
{
    int i,j,k;
    FILE* iosbb;
    int dummyint;
    double dummydouble;

    int nSignal,nBlock, nSBElement;
    CParameter* expParm;
    expParm= new CExpParm;
    //cout << "Num block range: "<< blockRange->getNumRange() << endl;
    for(i=0;i < blockRange->getNumRange();i++)
    {
        //printf("iBlockRange: %d\n",i);

        if ( (iosbb = fopen(sbbFName,"rb")) ==NULL)
        {
            printf("- File %s does not exist\n",sbbFName);
            exit(0);
        }

        printf("- Loading File %s.\n",sbbFName);

        int numread;
        cout << sbbHeader.numSignal << "  " << sbbHeader.numBlock << endl;

        for (i=0; i < sbbHeader.numSignal; i++ )
    	{
            numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
            sbbHeader.norm[i] = dummydouble;
            printf("Signal norm: %f \n",sbbHeader.norm[i]);
            for (j=0; j < sbbHeader.numBlock; j++ )
        	{
                numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
                printf("Block norm: %f \n",dummydouble);
                ((CStructBookExp*)structBook[i])[j].setNorm(dummydouble);

                numread = fread( &nSBElement, sizeof( int ), 1, iosbb );
                printf("nSBElement: %d \n",nSBElement);

                strtContinuousExp* dummySBExp = new strtContinuousExp[nSBElement];

                numread = fread( dummySBExp, sizeof( strtContinuousExp ), nSBElement, iosbb );

				((CStructBookExp*)structBook[i])[j].addElement(dummySBExp,nSBElement);

// 				for (k=0;k<nSBElement;k++)
// 					((CStructBookExp*)structBook[i])[j].printElementToScreen(k);

				delete [] dummySBExp;
            }
        }
        printf("- Loading has finished!!!\n");
        fclose(iosbb);
    }
    delete expParm;
}


void saveSBHeader(char* file, strtSBBHeader sbbHeader)
{
    // Read the header
    int dummyint;
    double dummydouble;
    FILE* iosbb;

    if ( (iosbb = fopen(file,"wb")) ==NULL)
    {
        printf("- File %s does not exist\n","header.sbb");
        exit(0);
    }
    
	dummyint = sbbHeader.type;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Signal type: %d \n",sbbHeader.type);

    // dummyint = sbbHeader.dicType;
    // fwrite( &dummyint, sizeof( int ), 1, iosbb );
    // printf("Dictionary type: %d \n",sbbHeader.dicType);

	dummyint = sbbHeader.numSignal;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Number of signals: %d \n",sbbHeader.numSignal);

    //if (sbbHeader.numSignal!=0) sbbHeader.norm = new double[sbbHeader.numSignal];

	dummyint = sbbHeader.signalSize;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Signal size: %d \n",sbbHeader.signalSize);

	dummyint = sbbHeader.blockHop;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Block hop: %d \n",sbbHeader.blockHop);

	dummyint = sbbHeader.blockSize;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Signal size: %d \n",sbbHeader.blockSize);

    dummyint = sbbHeader.numBlock;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Num Block: %d \n",sbbHeader.numBlock);

	dummydouble = sbbHeader.Fs;
    fwrite( &dummydouble, sizeof( double ), 1, iosbb );
    printf("Sampling freq.: %f \n",sbbHeader.Fs);

    dummyint = sbbHeader.initBlock;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Init Block: %d \n",sbbHeader.initBlock);

    dummyint = sbbHeader.finalBlock;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Final Block: %d \n",sbbHeader.finalBlock);    

    // Allocating memory for Structure Book Set
    //sbbHeader.numBlock = (int)ceil((double)sbbHeader.signalSize/(double)sbbHeader.blockHop);
    //printf("Number of blocks: %d \n",sbbHeader.numBlock);

    fclose(iosbb);

}

void saveSBHeaderGrp(char* file, strtSBBHeader sbbHeader)
{
    // Read the header
    int dummyint;
    double dummydouble;
    FILE* iosbb;

    if ( (iosbb = fopen(file,"wb")) ==NULL)
    {
        printf("- File %s does not exist\n","header.sbb");
        exit(0);
    }

	dummyint = sbbHeader.type;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Signal type: %d \n",sbbHeader.type);

    // dummyint = sbbHeader.dicType;
    // fwrite( &dummyint, sizeof( int ), 1, iosbb );
    // printf("Dictionary type: %d \n",sbbHeader.dicType);

	dummyint = sbbHeader.numSignal;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Number of signals: %d \n",sbbHeader.numSignal);

    //if (sbbHeader.numSignal!=0) sbbHeader.norm = new double[sbbHeader.numSignal];

	dummyint = sbbHeader.signalSize;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Signal size: %d \n",sbbHeader.signalSize);

	dummyint = sbbHeader.blockHop;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Block hop: %d \n",sbbHeader.blockHop);

	dummyint = sbbHeader.blockSize;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Block size: %d \n",sbbHeader.blockSize);

	dummydouble = sbbHeader.Fs;
    fwrite( &dummydouble, sizeof( double ), 1, iosbb );
    printf("Sampling freq.: %f \n",sbbHeader.Fs);

    dummyint = sbbHeader.subBlockSize;
    fwrite( &dummyint, sizeof( int ), 1, iosbb );
    printf("Block size: %d \n",sbbHeader.subBlockSize);

    // Allocating memory for Structure Book Set
    //sbbHeader.numBlock = (int)ceil((double)sbbHeader.signalSize/(double)sbbHeader.blockHop);
    //printf("Number of blocks: %d \n",sbbHeader.numBlock);

    fclose(iosbb);

}

void saveSBExp(	char* InputFile,
				strtSBBHeader sbbHeader,
				CStructBook** structBook)
{
	FILE* iosbb;
	iosbb = fopen(InputFile,"wb");

   	int i,j,iNumElement;
   	int dummyint;
   	double dummydouble;
   	strtContinuousExp  *sb;
	for (i=0; i < sbbHeader.numSignal; i++ )
    {
        dummyint = i+1;
        fwrite(&dummyint, sizeof(int), 1, iosbb);
        dummydouble = sbbHeader.norm[i];
        fwrite(&dummydouble, sizeof(double), 1, iosbb);
        for (j=0; j < sbbHeader.numBlock; j++ )
        {
        	dummyint = j+1;
			fwrite(&dummyint, sizeof(int), 1, iosbb);
			dummydouble = ((CStructBookExp*)structBook[i])[j].getNorm();
			fwrite(&dummydouble, sizeof(double), 1, iosbb);

         	sb = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
			for (iNumElement=0;iNumElement<sbNumElement;iNumElement++)
			{
				((CStructBookExp*)structBook[i])[j].saveElementBin(	iosbb,
																	iNumElement+1,
																	sb[iNumElement]);
			}
			dummyint = 99999;
			fwrite( &dummyint, sizeof( int ), 1, iosbb );
        }
        dummyint = 88888;
        fwrite( &dummyint, sizeof( int ), 1, iosbb );
	}
	dummyint = 77777;
	fwrite( &dummyint, sizeof( int ), 1, iosbb );

	fclose(iosbb);
}

// void saveSBExpGrp(	char* InputFile,
// 					strtSBBHeader sbbHeader,
// 					CStructBook** structBook)
// {
// 	FILE* iosbb;
// 	iosbb = fopen(InputFile,"wb");
//
//    	int i,j,iNumElement;
//    	int dummyint;
//    	double dummydouble;
//    	strtContinuousExp  *sb;
// 	for (i=0; i < sbbHeader.numSignal; i++ )
//     {
//         dummyint = i+1;
//         fwrite(&dummyint, sizeof(int), 1, iosbb);
//         dummydouble = sbbHeader.norm[i];
//         fwrite(&dummydouble, sizeof(double), 1, iosbb);
//         for (j=0; j < sbbHeader.numBlock; j++ )
//         {
//         	dummyint = j+1;
// 			fwrite(&dummyint, sizeof(int), 1, iosbb);
// 			dummydouble = ((CStructBookExp*)structBook[i])[j].getNorm();
// 			fwrite(&dummydouble, sizeof(double), 1, iosbb);
//
//          	sb = ((CStructBookExp*)structBook[i])[j].getStructBook();
// 			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
// 			for (iNumElement=0;iNumElement<sbNumElement;iNumElement++)
// 			{
// 				((CStructBookExp*)structBook[i])[j].saveElementBin(	iosbb,
// 																	iNumElement+1,
// 																	sb[iNumElement]);
// 			}
// 			dummyint = -1;
// 			fwrite( &dummyint, sizeof( int ), 1, iosbb );
//         }
//         dummyint = -2;
//         fwrite( &dummyint, sizeof( int ), 1, iosbb );
// 	}
// 	dummyint = -3;
// 	fwrite( &dummyint, sizeof( int ), 1, iosbb );
//
// 	fclose(iosbb);
// }

void saveSBExpGrp(	char* InputFile,
					strtSBBHeader sbbHeader,
					CStructBook** structBook)
{
	FILE* iosbb;
	iosbb = fopen(InputFile,"wb");

   	int i,j,iNumElement;
   	int dummyint;
   	double dummydouble;
   	strtContinuousExp  *sb;
	for (i=0; i < sbbHeader.numSignal; i++ )
    {
        dummydouble = sbbHeader.norm[i];
        fwrite(&dummydouble, sizeof(double), 1, iosbb);
        for (j=0; j < sbbHeader.numBlock; j++ )
        {
			dummydouble = ((CStructBookExp*)structBook[i])[j].getNorm();
			fwrite(&dummydouble, sizeof(double), 1, iosbb);

         	sb = ((CStructBookExp*)structBook[i])[j].getStructBook();
			int sbNumElement = ((CStructBookExp*)structBook[i])[j].getNumElement();
			dummyint = sbNumElement;
			fwrite( &dummyint, sizeof( int ), 1, iosbb );

			fwrite(sb, sizeof(strtContinuousExp), sbNumElement, iosbb);
        }
	}

	fclose(iosbb);
}


void loadSB(	char* InputFile,
				CFileDecompBlockRange* blockRange,
				strtSBBHeader sbbHeader,
				CStructBook** structBook)
{

    int initBlock,finalBlock;
    char sbbFName[_MAX_PATH];
    char aux[_MAX_PATH];
    char* pos;

    int i,j;
    FILE* iosbb;
    int dummyint;
    double dummydouble;
    strtContinuousExp dummySBExp;

    int nSignal,nBlock, nSBElement;
    CParameter* expParm;
    expParm= new CExpParm;
    //cout << "Num block range: "<< blockRange->getNumRange() << endl;
    for(i=0;i < blockRange->getNumRange();i++)
    {
        //printf("iBlockRange: %d\n",i);
        strcpy(sbbFName, InputFile);
        pos = strrchr( sbbFName, '.');
        initBlock = (blockRange->getPInitBlock())[i];
        finalBlock = (blockRange->getPFinalBlock())[i];
        sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
        strcpy( &pos[0], aux);

        if ( (iosbb = fopen(sbbFName,"rb")) ==NULL)
        {
            printf("- File %s does not exist\n",sbbFName);
            exit(0);
        }

        printf("- Loading File %s.\n",sbbFName);

        int numread;
        //while(!feof(iosbb))
        while(1)
        {
            numread = fread( &nSignal, sizeof( int ), 1, iosbb );
            //printf("Signal: %d \n",nSignal);
            if (nSignal==77777)
            {
                printf("- Loading complete!!!\n");
                break;
            }
            numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
            sbbHeader.norm[nSignal-1] = dummydouble;
            //printf("Signal norm: %f \n",sbbHeader.norm[nSignal-1]);
            //while(!feof(iosbb))
            while(1)
            {
                numread = fread( &nBlock, sizeof( int ), 1, iosbb );
                //printf("Block: %d \n",nBlock);
                if (nBlock==88888) break;
                numread = fread( &dummydouble, sizeof( double ), 1, iosbb );
                //printf("Block norm: %f \n",dummydouble);
                ((CStructBookExp*)structBook[nSignal-1])[nBlock-1].setNorm(dummydouble);
                //while (!feof(iosbb))
                while(1)
                {
                    numread = fread( &nSBElement, sizeof( int ), 1, iosbb );
                    //printf("Element: %d \n",nSBElement);
                    //printf("Numread: %d \n",numread);
                    if (nSBElement==99999) break;
                    numread = fread( &dummySBExp, sizeof( strtContinuousExp ), 1, iosbb );
                    ((CStructBookExp*)structBook[nSignal-1])[nBlock-1].addElement(dummySBExp);
                    //((CStructBookExp*)structBook[nSignal-1])[nBlock-1].printElementToScreen(nSBElement-1);
                    //if (feof(iosbb)) cout << "Fim de arquivo" << endl;
                    //if (ferror(iosbb)) cout << "Erro" << endl;
                }
            }
        }
        printf("- Loading has finished!!!\n");
        fclose(iosbb);
    }
    delete expParm;
}

void loadRDFile(int numSignal,int numQuant,strtRD** rd)
{
    FILE* rdFile;
    rdFile = fopen("rate_distortion.out","r");

    char auxStr[_MAX_PATH];
    char auxStr2[_MAX_PATH];
    char* sNumber;
    int init_shift = 31;
    int shift = 13;
    int ind_init;
    int ind_end;
    int i,k;

    for (k=0; k<numQuant; k++)
    {
        // read line
        fgets( auxStr, _MAX_PATH, rdFile );
		for (i=0;i<numSignal;i++)
		{
			// read rate
			strcpy(auxStr2,auxStr);
			ind_init = init_shift+i;
			ind_end = init_shift+i+shift;
			auxStr2[ind_end] = NULL;
			sNumber = &auxStr[ind_init];
        	rd[i][k].rate = atof(sNumber);

        	// read distortion
			strcpy(auxStr2,auxStr);
			ind_init = ind_init+shift+1;
			ind_end = ind_end+shift+1;
			auxStr2[ind_end] = NULL;
			sNumber = &auxStr[ind_init];
        	rd[i][k].dist = atof(sNumber);
        	printf(" %13.10f %13.10f\n",rd[i][k].rate,rd[i][k].dist);
		}
	}
}
