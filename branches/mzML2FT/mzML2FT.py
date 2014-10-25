#!/usr/bin/python

## Import Python package modules
import sys, getopt, warnings, os, re


def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:o:",
                                    ["help",
                                     "working-dir",
				                     "output-dir",])


    # Default working dir and config file
    working_dir = "./"
    output_dir = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-w workingdirectory -o outputdirectory"
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-dir"):
            output_dir = value

    mgf_filename_list = get_file_list_with_ext(working_dir, ".mgf")
    mgf_filename_list = mgf_filename_list + get_file_list_with_ext(working_dir, ".Mgf")
    mgf_filename_list = mgf_filename_list + get_file_list_with_ext(working_dir, ".MGF")
    
    if (output_dir == "") :
        output_dir = working_dir

    return [mgf_filename_list, output_dir]
    
    
    
## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list

def outputFT2(current_spectrum_list, s_scanId, d_precursor_mz, i_precursor_z, current_FT2_file, d_retention_time):
    current_FT2_file.write("S\t"+s_scanId+"\t"+s_scanId+"\t"+str(d_precursor_mz)+"\n")
    if (i_precursor_z == 0) :
        z_mz = 0
    else :
        z_mz = i_precursor_z*d_precursor_mz
    current_FT2_file.write("Z\t"+str(i_precursor_z)+"\t"+str(z_mz)+"\n")
    current_FT2_file.write("I\tRetentionTime\t"+str(d_retention_time)+"\n")
    for each_peak in current_spectrum_list :
        current_FT2_file.write(each_peak[0]+"\t"+each_peak[1]+"\t0\t0\t0\t0\n")

def ConvertMgfFile(current_mgf_filename, output_dir) :
    MgfFileNameBase = os.path.basename(current_mgf_filename)
    (MgfFileNameRoot, MgfFileNameExt) = os.path.splitext(MgfFileNameBase)
    current_FT2_filename = output_dir + os.sep + MgfFileNameRoot + ".FT2"
    current_FT2_file = file(current_FT2_filename, "w")
    current_mgf_file = open(current_mgf_filename)
    current_FT2_file.write("H\tExtractor\tmgf2FT\n")
    current_FT2_file.write("H\tm/z\tIntensity\tResolution\tBaseline\tNoise\tCharge\n")
    current_FT2_file.write("H\tInstrument Model\tNA\n")

    current_spectrum_list = []

    for each_line in current_mgf_file :
        each_line = each_line.strip()
        if (each_line == "") :
            continue
        if each_line.startswith("BEGIN IONS") :
            current_spectrum_list = []
            s_scanId = ""
            d_precursor_mz = 0
            i_precursor_z  = 0
            continue
        if each_line.startswith("END IONS") :
            if (s_scanId == "") or (d_precursor_mz == 0) :
                print "ill scan", s_scanId
            else:
                outputFT2(current_spectrum_list, s_scanId, d_precursor_mz, i_precursor_z, current_FT2_file, d_retention_time)
            continue
        
        if each_line.startswith("TITLE=scanId=") :
            s_scanId = each_line.split("=")[2]
        elif (each_line.startswith("CHARGE=")):
            i_precursor_z = int(each_line.split("=")[1])
        elif (each_line.startswith("PEPMASS=")):
            precursor_info = each_line.split("=")[1]
            d_precursor_mz = precursor_info.split(" ")[0]
        elif (each_line.startswith("RTINSECONDS=")) :
            d_retention_time = each_line.split("=")[1]
            
        elif (each_line[0:1].isdigit()) :
            
            peak_info = each_line.split(" ")
           # print peak_info
            current_spectrum_list.append([peak_info[0], peak_info[1]])

    current_FT2_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
        # parse options
        [mgf_filename_list, output_dir] = parse_options(argv)
    for each_mgf_filename in mgf_filename_list :
        ConvertMgfFile(each_mgf_filename, output_dir)
        


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
    

    
    