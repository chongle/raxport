**Export Thermo Raw files to FT1 and FT2 text files**

Chongle Pan

**Instructions:**
  1. Install the latest MSFileReader (v2.2 or higher) from Thermo Scientific. Freely available at http://sjsupport.thermofinnigan.com/public/detail.asp?id=703. Please check the boxes during installation to install both 32-bit and 64-bit versions of MSFileReader.
  1. Download raxport.exe from http://code.google.com/p/raxport/downloads/list
  1. Open a Windows console and run raxport.exe. Run "raxport.exe --help" to get usage:
```
Usage: raxport.exe -w WorkingDirectory

[-w WorkingDirectory]
Extract all raw files in the working directory and its sub-directories. 
Default: current directory

[--noFT1]
Do not extract FT1 files.
Default: false

[--noFT2]
Do not extract FT2 files.
Default: false

[--help]
Get help information

If no option is provided, extract both FT1 and FT2 files from the current directory
```

The working directory and its sub-directories contain the raw files to be exported. Raxport will export the FT1 and FT2 files to the same directory as their originating raw files. If the -w option is not provided, the working directory is defaulted to the current directory.

By default, Raxport exports both FT1 file and FT2 file from a Raw file. If FT1 files are not needed, use the --noFT1 option. If FT2 files are not needed, use the --noFT2 option.

A sample Raw file and its FT1 and FT2 files are provided at http://code.google.com/p/raxport/downloads/list

**Example 1:**
  1. Suppose that raxport.exe is saved to C:\proteomics\_software
  1. Suppose that the raw files are in D:\proteomics\_data\run1
  1. Type the command below in a Windows console
```
C:\proteomics_software\raxport.exe -w D:\proteomics_data\run1
```

This will extract FT1 and FT2 files from all raw files in D:\proteomics\_data\run1 and all its sub-directories. FT1 and FT2 files will be saved in the same directory as the raw files.

**Example 2:**
  1. Suppose that raxport.exe is saved to a directory listed in the Path variable of Windows Environment Variables. Then raxport.exe can be called without providing the full path of raxport.exe.
  1. Since the working directory is default to the current directory, you can go to the working directory and run raxport.exe without using the -w option.
  1. Type the command below in a Windows console to extract only FT2 files
```
cd D:\proteomics_data\run1
raxport.exe --noFT1
```

**FT1 File Format:**

```
S	12	12
I	RetentionTime	0.03873
I	ScanType	FT-MS1
I	ScanFilter	FTMS + p NSI Full ms [400.00-1700.00]
654.39655	628.70	12304	17.95	293.89	0
660.30737	1397.14	11801	18.18	294.56	2
660.80933	608.34	15304	18.20	294.62	2
```

  * S line: This is scan #12.
  * I line: Retention time is at 0.03873 minute
  * I line: The scan type is FT-MS1
  * I line: The filter text is "FTMS + p NSI Full ms [400.00-1700.00]", which means FTMS acquisition with positive nanospray ionization for full MS1 scan from 400 Da to 1700 Da
  * Data block: Each row represents a peak. High-resolution MS data contain six columns representing "m/z", "Intensity", "Resolution", "Baseline", "Noise", and "Charge". Low-resolution MS data contain two columns representing "m/z" and "Intensity".

**FT2 File Format:**

```
S	13	13	958.13501
Z	3	2874.40503
I	RetentionTime	0.04603
I	ScanType	FT-MS1/FT-MS2 = 958.46777 @ HCD
I	ScanFilter	FTMS + c NSI d Full ms2 958.47@hcd30.00 [110.00-2000.00]
D	ParentScanNumber	12
298.64740	331.26	13904	5.93	248.32	0
300.15582	609.44	20701	5.97	248.42	1
301.15485	326.45	14704	5.99	248.48	1
```

  * S line: This is scan #13. The monoisotopic m/z of parent ion is 958.13501
  * Z line: Parent ion charge state is 3. The charged monoisotopic mass of parent ion is 2874.40503 Da. Be aware of the 1-Da error for the monoisotopic mass.
  * I line: Retention time is at 0.04603 minute
  * I line: The scan type is FT-MS1/FT-MS2 with isolation target m/z at 958.46777. Note that the isolation target m/z is sometime different from the monoisotopic m/z of parent ion. The fragmentation method is HCD.
  * I line: The scan filter text is "FTMS + c NSI d Full ms2 958.47@hcd30.00 [110.00-2000.00]", which means centroid FTMS acquisition with positive nanospray ionization for data-dependent full MS2 scan using HCD at 30% collision energy level from 110 Da to 2000 Da
  * D line: the parent scan of this MS2 scan is scan #12
  * Data block: Each row represents a peak. High-resolution MS data contain six columns representing "m/z", "Intensity", "Resolution", "Baseline", "Noise", and "Charge". Low-resolution MS data contain two columns representing "m/z" and "Intensity".

**FAQ:**

  * What programs uses FT1 and FT2 files as input?
The Sipros program for protein identification (http://code.google.com/p/sipros), the ProRata program for protein quantification (http://code.google.com/p/prorata), the Vonode program for de novo sequencing (http://compbio.ornl.gov/Vonode/), and probably other programs that uses MS1 and MS2 file formats.

  * How are FT1 and FT2 files different from MS1 and MS2 files?
For low-resolution ion trap MS data, FT1 and FT2 files are very similar to MS1 and MS2 files, except for a few more scan header lines. For high-resolution FT MS data, FT1 and FT2 files also extract resolution, baseline, noise, and charge, provided in the labels of peaks.

  * What programs is Raxport dependent on to run?
Raxport is dependent on the freely available Thermo MSFileReader v2.2 or later. But it is not dependent on Xcalibur, which requires a license to install.

  * What operating system is supported?
We have tested Raxport on Windows XP and Windows 7. Linux is not supported, because Thermo MSFileReader is not available for Linux.


