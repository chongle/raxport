/********************************************************/
// Raxport v3.1 by Chongle Pan, ORNL
// For exporting FT1 and FT2 files from Thermo raw files
// Dependency on MSFileReader XRawFile2.dll
/********************************************************/

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <windows.h>

#import "C:\\Program Files (x86)\\Thermo\\MSFileReader\\XRawFile2.dll"

using namespace std;
using namespace MSFileReaderLib;


// this is used only by GetPrecursorInfoFromScanNum
struct PrecursorInfo
{
    double dIsolationMass;
    double dMonoIsoMass;
    long nChargeState;
    long nParentScanNumber;
};

// this is used only by GetMassListFromScanNum
struct DataPeak
{
	double dMass;
	double dIntensity;
};

// function to convert wchar_t* to std::string.
string ConvertWCSToMBS(const wchar_t* pstr, long wslen)
{
    int len = ::WideCharToMultiByte(CP_ACP, 0, pstr, wslen, NULL, 0, NULL, NULL);

    std::string dblstr(len, '\0');
    len = ::WideCharToMultiByte(CP_ACP, 0 /* no flags */,
                                pstr, wslen /* not necessary NULL-terminated */,
                                &dblstr[0], len,
                                NULL, NULL /* no default char */);

    return dblstr;
}

// function to convert BSTR to std::string.
string ConvertBSTRToMBS(BSTR bstr)
{
    int wslen = ::SysStringLen(bstr);
    return ConvertWCSToMBS((wchar_t*)bstr, wslen);
}

// function to convert std::string to wstring.
std::wstring ConvertStringToWString(const std::string& s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
    std::wstring r(buf);
    delete[] buf;
    return r;
}

// process raw files
bool ProcessFiles( vector<wstring> vsRawFiles )
{
     MSFileReaderLib::IXRawfile5Ptr XRawfileCtrl(NULL);

     HRESULT hr = XRawfileCtrl.CreateInstance("MSFileReader.XRawfile.1", NULL, CLSCTX_INPROC_HANDLER | CLSCTX_INPROC_SERVER);
	 
     if (FAILED(hr)) {
		 cout << endl;
         cout << "Cannot access MSFileReader XRawFile2.dll" << endl;
         cout << "Please install the latest MSFileReader from Thermo Scientific freely available at" << endl;
		 cout << "http://sjsupport.thermofinnigan.com/public/detail.asp?id=703" << endl;
         return false;
     }
	 else
	 {
		 cout << "Opened MSFileReader XRawfile2.dll" << endl << endl;
	 }

	 for( int i = 0; i < (int)vsRawFiles.size(); i++)
	 {
		 wstring sRawFilename = vsRawFiles[i];

		 wstring sFT1Filename = sRawFilename.substr(0, sRawFilename.length()-4) + L".FT1";
		 wstring sFT2Filename = sRawFilename.substr(0, sRawFilename.length()-4) + L".FT2";

		 ofstream FT1stream;
		 ofstream FT2stream;

		 FT1stream.open(sFT1Filename.c_str());
		 FT2stream.open(sFT2Filename.c_str());

		 FT1stream << "H	Extractor	Raxport v3.1 Build 10/16/2012" << endl;
		 FT1stream << "H	m/z	Intensity	Resolution	Baseline	Noise	Charge" << endl;

		 FT2stream << "H	Extractor	Raxport v3.1 Build 10/16/2012" << endl;
		 FT2stream << "H	m/z	Intensity	Resolution	Baseline	Noise	Charge" << endl;

		 long fileOpenFlag = XRawfileCtrl->Open(sRawFilename.c_str());
		 if(fileOpenFlag != 0){
			wcout << "Unable to open .RAW file: " << sRawFilename << endl;
			continue;
		 }

		 XRawfileCtrl->SetCurrentController(0, 1);

		 long firstScanNumber = 0; 
		 XRawfileCtrl->GetFirstSpectrumNumber(&firstScanNumber);
		 long lastScanNumber = 0;  
		 XRawfileCtrl->GetLastSpectrumNumber(&lastScanNumber);
		 wcout << "Extracting: " << sRawFilename << endl;
		 cout << "Total Scan Range = " <<  firstScanNumber << " ... " << lastScanNumber << endl;

 	
	    for(long scanNum = firstScanNumber;  scanNum <= lastScanNumber; scanNum++) 
		{
			cout << "Scan #" <<  scanNum << "\r"; 
//			cout << "Scan #" <<  scanNum << endl;

/*
// test GetTrailerExtraForScanNum
// this block of code works
VARIANT varExtraLabels;
VariantInit(&varExtraLabels);
VARIANT varExtraValues;
VariantInit(&varExtraValues);
long nArraySize = 0;
long nRet = XRawfileCtrl->GetTrailerExtraForScanNum( scanNum, &varExtraLabels, &varExtraValues, &nArraySize);
if( nRet != 0 )
{
	cout << "ERR, GetTrailerExtraForScanNum " << endl;
	continue;
}
// Get a pointer to the SafeArray
SAFEARRAY FAR* psaLabels = varExtraLabels.parray;
varExtraLabels.parray = NULL;
SAFEARRAY FAR* psaValues = varExtraValues.parray;
varExtraValues.parray = NULL;
BSTR* pbstrLabels = NULL;
BSTR* pbstrValues = NULL;
if( FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels) ) ) )
{
	SafeArrayUnaccessData( psaLabels );
	SafeArrayDestroy( psaLabels );
	cout << "ERROR, GetTrailerExtraForScanNum " << endl;
}
if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
{
	SafeArrayUnaccessData( psaLabels );
	SafeArrayDestroy( psaLabels );
	SafeArrayUnaccessData( psaValues );
	SafeArrayDestroy( psaValues );
	cout << "ERROR GetTrailerExtraForScanNum " << endl;
}
for( long i=0; i<nArraySize; i++ )
{
	string sLabel = ConvertBSTRToMBS(pbstrLabels[i]);
	string sData = ConvertBSTRToMBS(pbstrValues[i]);

	cout << " scanNum " << scanNum << " = " << sLabel << " + " << sData << endl;
}
// Delete the SafeArray
SafeArrayUnaccessData( psaLabels );
SafeArrayDestroy( psaLabels );
SafeArrayUnaccessData( psaValues );
SafeArrayDestroy( psaValues );
// end GetTrailerExtraForScanNum
*/

			// get MS Order
			long MSOrder = 0; 
			XRawfileCtrl->GetMSOrderForScanNum(scanNum, &MSOrder); 
//			cout << "MSOrder = " << MSOrder << endl;

			// get retention time
			double retentionTime = 0; 
			XRawfileCtrl->RTFromScanNum(scanNum, &retentionTime);
//			cout << "retentionTime = " << retentionTime << endl; 

			// get the filter text
		    BSTR bstrFilter = NULL;
			XRawfileCtrl->GetFilterForScanNum( scanNum, &bstrFilter );
			string sFilterText = ConvertBSTRToMBS(bstrFilter);
//			cout << "filter text = " << sFilterText << endl; 
		
			double accurate_precursorMZ = 0;

			double isolationMZ = 0;
			int precursorCharge = 0;
			int precursorScanNumber = 0;

			if(MSOrder == 2) {

				// call GetPrecursorInfoFromScanNum
				VARIANT vPrecursorInfos;
				VariantInit(&vPrecursorInfos);
				long nPrecursorInfos = 0;
				// Get the precursor scan information
				XRawfileCtrl->GetPrecursorInfoFromScanNum(scanNum, &vPrecursorInfos, &nPrecursorInfos);
        
				// Access the safearray buffer
				BYTE* pData;
				SafeArrayAccessData(vPrecursorInfos.parray, (void**)&pData);
				for (int i=0; i < nPrecursorInfos; ++i)
				{
					// Copy the scan information from the safearray buffer
					PrecursorInfo info;
					memcpy(&info, pData + i * sizeof(MS_PrecursorInfo), sizeof(PrecursorInfo));
					precursorScanNumber = info.nParentScanNumber;
					isolationMZ = info.dIsolationMass;
					/*
					cout << " scanNum " << scanNum 
						<< " dIsolationMass " << info.dIsolationMass 
						<< " dMonoIsoMass " << info.dMonoIsoMass
						<< " nChargeState " << info.nChargeState
						<< " nParentScanNumber " << info.nParentScanNumber
						<< endl;
					*/
				}
				SafeArrayUnaccessData(vPrecursorInfos.parray);
				// end GetPrecursorInfoFromScanNum

				// two values of monoisotopic m/z retrieved by different functions
				double precursorMZ_1 = 0;
				double precursorMZ_2 = 0;

				// retrieves monoisotopic m/z. 
				VARIANT varPrecursor; 
				VariantInit(&varPrecursor);
				XRawfileCtrl->GetTrailerExtraValueForScanNum(scanNum, "Monoisotopic M/Z:" , &varPrecursor);
				if( varPrecursor.vt == VT_R4 ) 
					precursorMZ_1 = varPrecursor.fltVal;
				else if( varPrecursor.vt == VT_R8 ) 
					precursorMZ_1 = varPrecursor.dblVal;

				// retrieves parent charge state
				VARIANT varCharge; VariantInit(&varCharge);
				XRawfileCtrl->GetTrailerExtraValueForScanNum(scanNum, "Charge State:" , &varCharge);
				if( varCharge.vt == VT_I2 ) 
					precursorCharge = varCharge.iVal;

		//		cout << "precursorCharge= " << precursorCharge<< endl;
		//		cout << "precursorMZ_1 = " << fixed << setprecision(15) << precursorMZ_1 << endl;
				
				// also retrieves monoisotopic m/z.
				XRawfileCtrl->GetPrecursorMassForScanNum(scanNum, MSOrder, &precursorMZ_2);
		//		cout << "precursorMZ_2 = " << precursorMZ_2 << endl;

				// precursorMZ_1 seems to be less prone to the 1-Da error, but it's often 0.
				// use precursorMZ_2, when precursorMZ_1 is 0
				if( precursorMZ_1 > 0.1 ){
					accurate_precursorMZ = precursorMZ_1;
				}
				else{
					accurate_precursorMZ = precursorMZ_2;
				}

			}

			if( MSOrder == 1 )
			{
				FT1stream << fixed << setprecision(5);
				FT1stream << "S\t" <<  scanNum << "\t" << scanNum << endl;
				FT1stream << "I\tRetentionTime\t" <<  retentionTime << endl;
				FT1stream << "I\tFilter\t" <<  sFilterText << endl;
			}
			else
			{
				FT2stream << fixed << setprecision(5);
				FT2stream << "S\t" <<  scanNum << "\t" << scanNum << "\t"  << accurate_precursorMZ <<endl;
				FT2stream << "Z\t" << precursorCharge << "\t" << precursorCharge * accurate_precursorMZ << endl;
				FT2stream << "I\tRetentionTime\t" <<  retentionTime << endl;
		//		FT2stream << "I\tisolationMZ\t" <<  isolationMZ << endl;
				FT2stream << "I\tFilter\t" <<  sFilterText << endl;
				FT2stream << "D\tParentScanNumber\t" <<  precursorScanNumber << endl;
			}

			// get label data
			VARIANT varLabels;
			VARIANT varFlags;
			VariantInit(&varLabels);
			VariantInit(&varFlags);


			XRawfileCtrl->GetLabelData(&varLabels, &varFlags, &scanNum);
			
	        int peakNumber = varLabels.parray->rgsabound[0].cElements;
//	        cout << scanNum << "   Array size: " << peakNumber << endl;

			if (peakNumber > 0 ){
				double* pdval;
				pdval=(double*)varLabels.parray->pvData;
			
				if( MSOrder == 1 ){
					for(int n = 0; n < peakNumber; n++){
						FT1stream 
							<< fixed << setprecision(5) << pdval[n*6+0] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+1] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+2] << "\t"
							<< fixed << setprecision(0) << pdval[n*6+3] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+4] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+5] << endl;
					}
				}
				else{
					for(int n = 0; n < peakNumber; n++){
						FT2stream 
							<< fixed << setprecision(5) << pdval[n*6+0] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+1] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+2] << "\t"
							<< fixed << setprecision(0) << pdval[n*6+3] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+4] << "\t" 
							<< fixed << setprecision(0) << pdval[n*6+5] << endl;
					}
				}
			}
			else{
				VARIANT varMassList;
				VariantInit(&varMassList);
				VARIANT varPeakFlags;
				VariantInit(&varPeakFlags);
				long nArraySize = 0;
				double peakWidth = 0;
				XRawfileCtrl->GetMassListFromScanNum ( &scanNum,
					"",0, 0, 0,long(0),&peakWidth,
					&varMassList,// mass list data
					&varPeakFlags,// peak flags data
					&nArraySize ); // size of mass list array
				if( nArraySize > 0)
				{
					// Get a pointer to the SafeArray
					SAFEARRAY FAR* psa = varMassList.parray;
					DataPeak* pDataPeaks = NULL;
					SafeArrayAccessData( psa, (void**)(&pDataPeaks) );
					for( long j=0; j<nArraySize; j++ )
					{
						double dMass = pDataPeaks[j].dMass;
						double dIntensity = pDataPeaks[j].dIntensity;
						if( MSOrder == 1 ){
							FT1stream << fixed << setprecision(5) << dMass << "\t" 
									  << fixed << setprecision(0) << dIntensity << endl; 
						}
						else{
							FT2stream << fixed << setprecision(5) << dMass << "\t" 
									  << fixed << setprecision(0) << dIntensity << endl;
						}

					}
					// Release the data handle
					SafeArrayUnaccessData( psa );
				}
				if( varMassList.vt != VT_EMPTY )
				{
					SAFEARRAY FAR* psa = varMassList.parray;
					varMassList.parray = NULL;
					// Delete the SafeArray
					SafeArrayDestroy( psa );
				}
				if(varPeakFlags.vt != VT_EMPTY )
				{
					SAFEARRAY FAR* psa = varPeakFlags.parray;
					varPeakFlags.parray = NULL;
					// Delete the SafeArray
					SafeArrayDestroy( psa );
				}
			}

			VariantClear(&varLabels);
			VariantClear(&varFlags);
		}
		XRawfileCtrl->Close();
		FT1stream.close();
		FT2stream.close();
	}
	return true;
}

// recursively find matched files in sub-directories
void findMatchedFiles(wstring sPath, wstring namePattern, vector<wstring> & vsRawFilenames)
{
	if(sPath[sPath.length()-1] != L'\\'){
		sPath = sPath + L'\\';
	}

	// search for *.raw in sPath
	wstring searchPattern = sPath + namePattern;
	wcout << "Searching for: " << searchPattern << endl;

	HANDLE hFind;
	WIN32_FIND_DATA FindFileData;

	hFind = FindFirstFile(searchPattern.c_str(), &FindFileData);
	while(hFind != INVALID_HANDLE_VALUE)
	{
		wstring foundFile( FindFileData.cFileName );
		wstring pathFilename = sPath + foundFile;
		vsRawFilenames.push_back(pathFilename);
		wcout << "Found: " << pathFilename << endl;
		if(FindNextFile(hFind, &FindFileData) == FALSE)
			break;
	}
	wcout << endl;

	FindClose(hFind);

	// search for sub-directories and recursively search for raw files

	searchPattern = sPath + L'*';
//	wcout  << endl << "Searching for: " << searchPattern << endl<< endl;

	hFind = FindFirstFile(searchPattern.c_str(), &FindFileData);
	while(hFind != INVALID_HANDLE_VALUE)
	{
		wstring foundFile( FindFileData.cFileName );
		if( FindFileData.dwFileAttributes == FILE_ATTRIBUTE_DIRECTORY
            && foundFile != L"." && foundFile != L"..")
		{
			// found a subdirectory; recurse into it
			wstring pathFilename = sPath + foundFile;
//			wcout << "Go to sub-directory: " << pathFilename << endl;
			findMatchedFiles(pathFilename, namePattern, vsRawFilenames);
		}
		if(FindNextFile(hFind, &FindFileData) == FALSE)
			break;
	}

	FindClose(hFind);

	return;
}

bool initializeArguments(int argc, char* argv[], wstring & sPath)
{
    // Grab command line arguments
    vector<string> vsArguments;
    while(argc--)
		vsArguments.push_back(*argv++);

	string sWorkingDirectory  = "";
	for(int i = 1; i < (int)vsArguments.size(); i++)
	{
		if(vsArguments[i] == "-w")
		{
			i = i + 1;
			if( i < (int)vsArguments.size() )
			{
				sWorkingDirectory = vsArguments[i];
			}
			else
			{
				cout << "Usage: -w WorkingDirectory" << endl;
				return false;
			}
		}
		else if (vsArguments[i] == "-h" || vsArguments[i] == "--help")
		{
			cout << "Usage: -w WorkingDirectory" << endl;
			cout << "FT1 and FT2 files will be extracted from all raw files in the working directory and its sub-directories." << endl;
			cout << "If no option is provided, the default working directory is the current directory." << endl;
			return false;
		}
		else 
		{
			// ignore Unknown options
			cout << "Unknown option: " << vsArguments[i] << endl;
			cout << "Usage: -w WorkingDirectory" << endl;
			cout << "Use -h to get more information" << endl << endl;
		}
	}
	if( sWorkingDirectory == "" )
		sWorkingDirectory = ".";

	sPath = ConvertStringToWString(sWorkingDirectory);
	return true;
}

int main(int argc, char* argv[])
{
	cout << endl << "Raxport 3.1 Build 10/16/2012" << endl << endl;

	wstring sPath = L""; 
	if( ! initializeArguments(argc, argv, sPath) )
	{
		return 0;
	}

	wstring namePattern = L"*.raw";

	vector<wstring> vsRawFilenames;
	findMatchedFiles(sPath, namePattern, vsRawFilenames);
	if(vsRawFilenames.size() == 0)
	{
		wcout << "Cannot find any raw file in directory = " << sPath << " with Name pattern = " << namePattern << endl;
	}
	else
	{
		CoInitialize( NULL );
		ProcessFiles(vsRawFilenames);
		CoUninitialize();
	}

	return 0;

}