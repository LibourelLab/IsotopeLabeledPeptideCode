#include <comutil.h>
#include <stdio.h>
#include <io.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include "mex.h"
#include <matrix.h> 
#include <string.h>
#import "C:\Program Files (x86)\Thermo\MSFileReader\Xrawfile2.dll"

using namespace std;

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    //make sure there is exactly one input
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexevalstring:nInput", "One input argument required.");
    } 
//**********************************************************************************//
//Set up MEX inout
    //Declarations 
    mxArray *yData;
    int yLength;
    char *fileName;
    
    //copy input pointer y
    yData = (mxArray*)prhs[0];
    
    //make "fileName" point to char array
    yLength = (int)(mxGetN(yData)+1);
    fileName = (char *)(mxCalloc(yLength,sizeof(char)));
    mxGetString(yData,fileName,yLength);
    
//Set up MEX output (same steps as above...)
    mxArray *numScans;
    double *output;
    numScans = plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    output = mxGetPr(numScans);
//*********************************************************************************//    
// Code that opens RAW files    
    
    //Declare for reading raw files
    //HRESULT lRet;
    TCHAR pth[MAX_PATH];
    MSFileReaderLib::IXRawfile4Ptr m_Raw;
    
    // Initialize COM interface
    HRESULT hr = CoInitialize(NULL);    
    if (FAILED(hr)) {
        mexErrMsgTxt("Cannot initialize COM interface.");
    }
 
    //CoInitialize( NULL );
    hr = m_Raw.CreateInstance("MSFileReader.XRawfile.1");
    if (FAILED(hr)) {
        mexErrMsgTxt("Cannot access XRawfile.dll.");
    }
//***********************************************************************************//
    MultiByteToWideChar(CP_ACP,0,fileName,-1,(LPWSTR)pth,MAX_PATH);
    long lr = m_Raw->Open((LPWSTR)pth);
    // Look for data that belong to the first mass spectra device in the file
     m_Raw->SetCurrentController(0, 1);
    
    long firstScanNumber = 0, lastScanNumber = 0;

    // Verifies if can get the first scan
    hr = m_Raw->GetFirstSpectrumNumber(&firstScanNumber);
    if (FAILED(hr)) {
        mexErrMsgTxt("ERROR: Unable to get first scan.");
    }

    // Ask for the last scan number to prepare memory space, for cycle 
    // and final verification
    m_Raw->GetLastSpectrumNumber(&lastScanNumber);
    long totalNumScans = (lastScanNumber - firstScanNumber) + 1;


  //output totalNumScans
    *output = totalNumScans; 
    hr=m_Raw->Close();
      if (FAILED(hr)) {
        mexErrMsgTxt("ERROR: Unable to close RAW file.");
        }

}