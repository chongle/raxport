'Raxport v2.0, Chongle Pan, panc@ornl.gov

Module Module1

    Public bExtractFT1 = False 'flag for ms1 file extraction
    Public bExtractFT2 = False 'flag for ms2 file extraction
    Public bExtractMS1 = False 'flag for ms1 file extraction
    Public bExtractMS2 = False 'flag for ms2 file extraction
    Public bMS2Zline = True   'flag for printing Z line in MS2 files 


    Sub Main(ByVal sArgs() As String)

        'default PathIn to the current dir
        Dim PathIn As String = "."
        Dim PathOut As String = "NA"

        bExtractFT1 = False
        bExtractFT2 = False
        bExtractMS1 = False
        bExtractMS2 = False
        bMS2Zline = True

        Console.WriteLine("Raxport 2.0 beta; Chongle Pan; ORNL")
        Console.WriteLine("Use the --help option to get more information")

        'process arguments
        ' -w InputDir -o OutputDir -n1[Optional flag for no ms1 files] -n2[Optional flag for no ms2 files]

        Dim i As Integer = 0
        While i < sArgs.Length

            'read -w input dir
            If StrComp(sArgs(i), "-w") = 0 Then
                i = i + 1
                PathIn = sArgs(i)
            End If

            'read -o output dir
            If StrComp(sArgs(i), "-o") = 0 Then
                i = i + 1
                PathOut = sArgs(i)
            End If

            'read --FT1, to generate FT1 files
            If StrComp(sArgs(i), "--FT1") = 0 Then
                bExtractFT1 = True
                bExtractMS1 = False
            ElseIf StrComp(sArgs(i), "--MS1") = 0 Then
                bExtractMS1 = True
                bExtractFT1 = False
            End If

            'read --FT2, to generate FT2 files
            If StrComp(sArgs(i), "--FT2") = 0 Then
                bExtractFT2 = True
                bExtractMS2 = False
            ElseIf StrComp(sArgs(i), "--MS2") = 0 Then
                bExtractMS2 = True
                bExtractFT2 = False
            End If

            'read --noZ, set 
            If StrComp(sArgs(i), "--noZ") = 0 Then
                bMS2Zline = False
            End If

            'read --help
            If StrComp(sArgs(i), "--help") = 0 Then
                Console.WriteLine("Usage: -w WorkingDirectory -o OutputDiretory --FT1 [To generate FT1 files] --FT2 [To generate FT2 files] --MS1 [To generate MS1 files] --MS2 [To generate MS2 files] --noZ [Do not write Z lines in MS2 files, Z lines to be written by MS2ZAssign.exe].")
                Console.WriteLine("Default values: input directory = current working directory; output directory = input directory.")
            End If

            i = i + 1
        End While

        'if PathOut is not given, then default it to PathIn
        If PathOut = "NA" Then
            PathOut = PathIn
        End If

        If Not (bExtractFT1 Or bExtractMS1 Or bExtractFT2 Or bExtractMS2) Then
            Console.WriteLine("Please indicate which types of output files should be generated:")
            Console.WriteLine("--FT1 [To generate FT1 files] --FT2 [To generate FT2 files] --MS1 [To generate MS1 files] --MS2 [To generate MS2 files]")
        End If


        'process all RAW files in PathIn and save all FT2 files into PathOut
        Call PrecessDir(PathIn, PathOut)

    End Sub

    Sub PrecessDir(ByVal PathIn As String, ByVal PathOut As String)

        Dim fs, f, f1, fc
        Dim sFT1Filename As String
        Dim sFT2FileName As String
        Dim sRawFileName As String

        fs = CreateObject("Scripting.FileSystemObject")
        f = fs.GetFolder(PathIn)
        fc = f.Files

        'for each RAW files in the PathIN
        For Each f1 In fc
            If InStr(f1.Name, ".RAW") > 0 Then
                sRawFileName = PathIn + "\" + f1.Name
                If bExtractFT1 Then
                    sFT1Filename = PathOut + "\" + Replace(f1.Name, ".RAW", ".FT1")
                End If
                If bExtractFT2 Then
                    sFT2FileName = PathOut + "\" + Replace(f1.Name, ".RAW", ".FT2")
                End If
                If bExtractMS1 Then
                    sFT1Filename = PathOut + "\" + Replace(f1.Name, ".RAW", ".MS1")
                End If
                If bExtractMS2 Then
                    sFT2FileName = PathOut + "\" + Replace(f1.Name, ".RAW", ".MS2")
                End If

                'process this RAW file
                Call ProcessFile(sRawFileName, sFT1Filename, sFT2FileName)
            ElseIf InStr(f1.Name, ".raw") > 0 Then
                sRawFileName = PathIn + "\" + f1.Name
                If bExtractFT1 Then
                    sFT1Filename = PathOut + "\" + Replace(f1.Name, ".raw", ".FT1")
                End If
                If bExtractFT2 Then
                    sFT2FileName = PathOut + "\" + Replace(f1.Name, ".raw", ".FT2")
                End If
                If bExtractMS1 Then
                    sFT1Filename = PathOut + "\" + Replace(f1.Name, ".raw", ".MS1")
                End If
                If bExtractMS2 Then
                    sFT2FileName = PathOut + "\" + Replace(f1.Name, ".raw", ".MS2")
                End If

                'process this RAW file
                Call ProcessFile(sRawFileName, sFT1Filename, sFT2FileName)
            End If

        Next


    End Sub


    Sub ProcessFile(ByVal sRawFileName As String, ByVal sFT1FileName As String, ByVal sFT2FileName As String)

        Console.WriteLine("Extracting: " & sRawFileName)

        'Open FT1 file for output.
        If bExtractFT1 Or bExtractMS1 Then
            FileOpen(1, sFT1FileName, OpenMode.Output)
            'Write FT1 header
            PrintLine(1, "H" & vbTab & "Extractor" & vbTab & "Raxport v2.0")
            PrintLine(1, "H" & vbTab & "m/z" & vbTab & "Intensity" & vbTab & "Resolution" & vbTab & "Baseline" & vbTab & "Noise" & vbTab & "Charge")
        End If

        If bExtractFT2 Or bExtractMS2 Then
            FileOpen(2, sFT2FileName, OpenMode.Output)
            'Write FT2 header
            PrintLine(2, "H" & vbTab & "Extractor" & vbTab & "Raxport v2.0")
            PrintLine(2, "H" & vbTab & "m/z" & vbTab & "Intensity" & vbTab & "Resolution" & vbTab & "Baseline" & vbTab & "Noise" & vbTab & "Charge")
        End If

        'Open this RAW file with XCALIBURFILESLib
        Dim xRaw As XCALIBURFILESLib.XRaw
        xRaw = New XCALIBURFILESLib.XRaw
        xRaw.Open(sRawFileName)

        'Go to XDectorRead
        Dim xDetectorRead As XCALIBURFILESLib.XDetectorRead
        xDetectorRead = xRaw.Detector(0, 1)

        'Go to XFilters
        'There is a probable Bug in XFilter
        Dim xFilters As XCALIBURFILESLib.XFilters
        xFilters = xDetectorRead.Filters

        'Go to XSpectra
        Dim xSpectra As XCALIBURFILESLib.XSpectra
        xSpectra = xDetectorRead.Spectra(0)
        Dim iTotalScanCount As Integer      'Total number of spectra in this RAW file
        iTotalScanCount = xSpectra.Count()  'Get it from xSpectra

        'Cycle thru every XSpectrumRead in xSpectra and its associated XFilterRead in xFilters
        Dim iScanNumber As Integer = 1

        'Declare varibles
        Dim xSpectrumRead As XCALIBURFILESLib.XSpectrumRead
        Dim xSpectrumHeader As XCALIBURFILESLib.XSpectrumHeader

        Dim xFilter As XCALIBURFILESLib.XFilterRead
        Dim XParentScans As XCALIBURFILESLib.XParentScans
        Dim XParentScan As XCALIBURFILESLib.XParentScan
        Dim sFilterText As String

        Dim vData As Object  'This is for ITMS Data
        Dim vLabelData As Object  'This is for FTMS data
        'Dim vNoiseData As Object

        Dim dRetentionTime As Double

        Dim parentScanCount As Integer
        Dim parentScanNumber As Integer
        Dim parentMonoIsoMass As Double
        Dim parentChargeState As Double
        Dim parentIsolationMass As Double


        Dim j As Integer = 0

        Dim iMSn As Integer = 0 '1, if MS1; 2, if MS2; ...; n, if MSn
        Dim bIsFTMS As Boolean = True 'is the current scan an FTMS scan (True) or an ITMS scan (False)?


        For iScanNumber = 1 To iTotalScanCount

            'Console.WriteLine("Scan: " & iScanNumber)

            'Start to access SpectrumRead
            xSpectrumRead = xSpectra.Item(iScanNumber)

            xSpectrumHeader = xSpectrumRead.Header
            dRetentionTime = xSpectrumHeader.StartTime()

            'XFilter

            'xFilter = xFilters.ScanNumber(iScanNumber)
            'sFilterText = xFilter.Text
            'Console.WriteLine("filter: " & sFilterText)

            'determine FTMS or ITMS
            'If (InStr(sFilterText, "FTMS") > 0) Then
            '   bIsFTMS = True
            'Else
            '    bIsFTMS = False
            'End If

            'End of extracting info from xFilter

            'Determine MSn
            'This is a temperatory work-around for the sFilterText problem
            'When XParentScans is MS1, execution of xSpectrumRead.ParentScans will always throw an exception
            'If an exception is thrown, then this spectrum is MS1
            iMSn = 0
            Try
                XParentScans = xSpectrumRead.ParentScans
                parentScanCount = XParentScans.Count
                iMSn = parentScanCount + 1
            Catch ex As Exception
                iMSn = 1
            End Try


            If iMSn = 1 Then

                ' For some reason, full scan data need to be accessed to prevent memory leak

                'a work-around to determine bIsFTMS
                vLabelData = xSpectrumRead.LabelData
                vData = xSpectrumRead.Data
                ' Console.WriteLine("MS1 : " & UBound(vLabelData, 2))
                If (UBound(vLabelData, 2) <= 0) Then
                    bIsFTMS = False
                Else
                    bIsFTMS = True
                End If

                If (bExtractFT1 Or bExtractMS1) Then
                    'Write Spectrum Data
                    PrintLine(1, "S" & vbTab & iScanNumber & vbTab & iScanNumber)
                    PrintLine(1, "I" & vbTab & "RetentionTime" & vbTab & Format(dRetentionTime, "0.00"))
                    If bIsFTMS Then 'Write FTMS data
                        For j = 0 To UBound(vLabelData, 2)
                            PrintLine(1, Format(vLabelData(0, j), "0.00000") & vbTab & Format(vLabelData(1, j), "0.00") & vbTab & Format(vLabelData(2, j), "0.00") & vbTab & Format(vLabelData(3, j), "0.00") & vbTab & Format(vLabelData(4, j), "0.00") & vbTab & Format(vLabelData(5, j), "0.00"))
                        Next

                    Else 'Write ITMS data
                        For j = 0 To UBound(vData, 2)
                            PrintLine(1, Format(vData(0, j), "0.0") & vbTab & Format(vData(1, j), "0"))
                        Next

                    End If
                End If

            ElseIf iMSn = 2 Then     'Write FT2 file

                'a work-around to determine bIsFTMS
                vLabelData = xSpectrumRead.LabelData
                vData = xSpectrumRead.Data
                If (UBound(vLabelData, 2) <= 0) Then
                    bIsFTMS = False
                Else
                    bIsFTMS = True
                End If

                XParentScan = XParentScans.Item(1)

                parentMonoIsoMass = XParentScan.MonoIsoMass
                parentIsolationMass = XParentScan.IsolationMass
                parentChargeState = XParentScan.ChargeState
                parentScanNumber = XParentScan.ScanNumber


                If bExtractFT2 Or bExtractMS2 Then

                    If parentMonoIsoMass > 1 Then
                        ' FT full scans; print monoisotopic parent mass and charge state
                        PrintLine(2, "S" & vbTab & iScanNumber & vbTab & iScanNumber & vbTab & Format(parentMonoIsoMass, "0.00000"))
                        If bExtractFT2 Or bMS2Zline Then
                            PrintLine(2, "Z" & vbTab & parentChargeState & vbTab & Format((parentMonoIsoMass * parentChargeState), "0.00000"))
                        End If
                    Else
                        ' Ion trap full scans; print isolation mass and don't print Z line
                        PrintLine(2, "S" & vbTab & iScanNumber & vbTab & iScanNumber & vbTab & Format(parentIsolationMass, "0.0"))
                    End If

                    PrintLine(2, "I" & vbTab & "RetentionTime" & vbTab & Format(dRetentionTime, "0.00"))
                    PrintLine(2, "D" & vbTab & "ParentScanNumber" & vbTab & parentScanNumber)



                    'Write Spectrum Data
                    If bIsFTMS Then 'Write FTMS data
                        'PrintLine(2, "D", TAB(), "ScanType", TAB(), "FT-MS2")
                        For j = 0 To UBound(vLabelData, 2)
                            PrintLine(2, Format(vLabelData(0, j), "0.00000") & vbTab & Format(vLabelData(1, j), "0.00") & vbTab & Format(vLabelData(2, j), "0.00") & vbTab & Format(vLabelData(3, j), "0.00") & vbTab & Format(vLabelData(4, j), "0.00") & vbTab & Format(vLabelData(5, j), "0.00"))
                        Next

                    Else 'Write ITMS data
                        'PrintLine(2, "D", TAB(), "ScanType", TAB(), "IT-MS2")
                        For j = 0 To UBound(vData, 2)
                            PrintLine(2, Format(vData(0, j), "0.0") & vbTab & Format(vData(1, j), "0"))
                        Next

                    End If
                End If

            Else
                'Nothing
            End If

            'Continue to the next spectrum
        Next

        'close files
        If bExtractFT1 Or bExtractMS1 Then
            FileClose(1)
        End If

        If bExtractFT2 Or bExtractMS2 Then
            FileClose(2)
        End If

    End Sub


End Module
