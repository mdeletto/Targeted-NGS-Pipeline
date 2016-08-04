#/usr/bin/python

import optparse, sys, subprocess, datetime, os
from Variant_Detection_functions import *


#-------------------------------------------------------------------------------------------
#----------------------------Command line parser and arguments------------------------------
#-------------------------------------------------------------------------------------------

desc="""This wrapper script is intended to be used for variant detection and annotation. """

parser = optparse.OptionParser(description=desc)

parser.add_option('--variant_detection_mode', help='Enable variant detection mode and subsequent annotation', dest='mode', action='store_true',default=True)
parser.add_option('--annotation_mode', help='Enable annotation mode only', dest='mode', action='store_false',default=True)
parser.add_option('--enable_filtering', help='Enable filtering for VCFs', dest='filter', action='store_true',default=True)
parser.add_option('-s', help='Check for classes of variants in variant detection mode <All,All-HC,Somatic,Germline,LOH,Somatic_LOH', dest='status', action='store',default='Somatic_LOH')
parser.add_option("--ion_reporter_vcf", action="store", help="Input absolute path to IonReporter VCF </path/to/ionreporter.vcf>",dest="ionreporter_vcf")
parser.add_option("--ion_reporter_fusion_vcf", action="store", help="Input absolute path to IonReporter FUSION VCF (for OCP only) </path/to/ionreporter.fusion.vcf>",dest="ionreporter_fusion_vcf",default=None)
parser.add_option("--ion_reporter_tsv", action="store", help="Input absolute path to IonReporter TSV </path/to/ionreporter.tsv>",dest="ionreporter_tsv")
parser.add_option('-c', help="Base output name for files.  If a copath ID exists, please use this as the base output name.", dest='base_output', action='store')
parser.add_option('--ionreporter_url_bool', help="Flag for if VCF is coming from IonReporter url", dest='ionreporter_url_bool', action='store_true',default=False)
parser.add_option('--ionreporter_analysis_name', help="Analysis name within IonReporter", dest='ionreporter_analysis_name', action='store')
parser.add_option('--ionreporter_id', help="Analysis ID within IonReporter", dest='ionreporter_id', action='store')#,default=opts.base_output)
parser.add_option('--ionreporter_fusion_url_bool', help="Flag for if FUSION VCF is coming from IonReporter url", dest='ionreporter_fusion_url_bool', action='store_true',default=False)
parser.add_option('--ionreporter_fusion_analysis_name', help="Fusion analysis name within IonReporter", dest='ionreporter_fusion_analysis_name', action='store')
parser.add_option('--ionreporter_fusion_id', help="Fusion analysis ID within IonReporter", dest='ionreporter_fusion_id', action='store')#,default=opts.base_output)
parser.add_option('--ionreporter_fusion_vcf', help="Fusion VCF within IonReporter", dest='ionreporter_fusion_vcf', action='store')
parser.add_option('--url', help="Flag for if BAM is coming from url", dest='url', action='store_true',default=False)
parser.add_option('--regions', help="Regions to focus analysis on.  Must be in BED format", dest='regions', action='store')
parser.add_option("--ion_reporter_version", action="store", help="Please indicate IR version number (ex. 4.0, 4.4)", dest="ionreporter_version", default=4.4)

group = optparse.OptionGroup(parser, "Variant Detection Mode",
                    "Requires input of two bams--tumor and normal")
group.add_option("-t", action="store", help="Input absolute path to tumor bam </path/to/sample.bam>",dest="tumor")
group.add_option("-n", action="store", help="Input absolute path to normal bam---NOTE: If no normal bam is selected, a population normal will be used.",dest="normal",default=None)
group.add_option("--tumor-purity", action="store", help="Tumor purity <Float between 0 and 1>",dest="tumor_purity",default=1)
group.add_option("--normal-purity", action="store", help="Normal purity <Float between 0 and 1>",dest="normal_purity",default=1)
group.add_option("-p", action="store", help="Platform used for sequencing <PGM or Proton>",dest="platform",default="PGM")
parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Annotation Mode",
                    "Requires input of VCF files from both IonReporter and VarScan")
group.add_option("-v", action="store", help="Input absolute path to VarScan VCF",dest="varscan_vcf")
parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Galaxy options",
                    "ONLY TO BE USED IN A GALAXY ENVIRONMENT")
group.add_option("--galaxy_html_file", action="store", help="Path to galaxy output file---to be used only in Galaxy environment",dest="galaxy_html_file")
group.add_option("--galaxy_output_directory", action="store", help="Path to galaxy output directory---to be used only in Galaxy environment",dest="galaxy_output_directory")
parser.add_option_group(group)


(opts, args) = parser.parse_args()

#-------------------------------------------------------------------------------------------
#------------------------PATHS TO PROGRAMS/FILES USED BY PIPELINE---------------------------
#-------------------------------------------------------------------------------------------

SAMTOOLS_EXE = "/home/michael/bin/samtools-1.1/samtools"
VARSCAN_EXE = "/home/michael/bin/VarScan/VarScan.v2.3.7.jar"
SNPSIFT_EXE = "/home/michael/bin/snpEff/SnpSift.jar"
VEP_EXE = "/home/michael/bin/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl"
EDDY_EXE = "/home/michael/bin/eddy.jar"


VEP_REF_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.VEP.fasta"
REFERENCE_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.fasta"

if opts.regions=="CCP":
	REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CCPHSMV2_052013.bed"
elif opts.regions=="OCP":
	REGIONS_FILE = "/home/michael/YNHH/Reference_Files/OCP3v1.20140718.designed.bed"
elif opts.regions=="BRCA":
	REGIONS_FILE = "/home/michael/YNHH/Reference_Files/BRCA1-2/BRCA1_2.20131001.designed.bed"
elif opts.regions=="CHPv2":
	REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CHPv2/CHP2.20131001.designed.bed"
elif opts.regions=="TSC":
    REGIONS_FILE = "/home/michael/YNHH/Reference_Files/TSC1-TSC2/TSC1_2.designed.bed"
else:
	REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CCPHSMV2_052013.bed"

#POPULATION_NORMAL_PGM = 
#POPULATION_NORMAL_PROTON = 


logfile = open("%s.log.txt" % opts.base_output,"w")
logfile.write("Parameters for analysis:\n\n")
logfile.write("Date of analysis: %s\n" % datetime.datetime.now())
for key,value in opts.__dict__.iteritems():
    logfile.write("%s\t%s\n" % (key,value))
logfile.close()

if opts.mode==True:
#     mandatory_options = ['tumor','ionreporter_vcf','ionreporter_tsv','platform']
#     for m in mandatory_options:
#         # Making sure all mandatory options appeared
#         if not opts.__dict__[m]:
#             print "Mandatory option is missing!\n"
#             parser.print_help()
#             sys.exit()        

    ### COMMAND LINE CHECKS FOR REMOTE INPUT ###

    if opts.url == True:
        opts.tumor = pull_BAM_from_url(opts.tumor,opts.base_output,"tumor")
        opts.normal = pull_BAM_from_url(opts.normal,opts.base_output,"normal")
    else:
        pass

    if opts.ionreporter_url_bool == True:
        variant_link = IR_locate_variant_zip(opts.ionreporter_analysis_name,opts.ionreporter_id)
        IR_download_variant_zip(opts.ionreporter_analysis_name, variant_link)
        opts.ionreporter_vcf = "%s.ionreporter.temp.vcf" % opts.ionreporter_analysis_name
        opts.ionreporter_tsv = "%s.ionreporter.temp.tsv" % opts.ionreporter_analysis_name
    else:
        pass
    
    ### FUSION SUPPORT (added in v1.3) ####
        
    if opts.ionreporter_fusion_url_bool == True:
        print "Attempting to pull IR fusions"
        variant_link = IR_locate_variant_zip(opts.ionreporter_fusion_analysis_name,opts.ionreporter_fusion_id) # Determine download link
        IR_download_fusion_zip(opts.ionreporter_fusion_analysis_name, variant_link, opts.base_output) # Download analysis and process files
    else:
        if opts.ionreporter_fusion_vcf is None:
            if opts.regions == "OCP":
                #subprocess.call("touch %s.ionreporter.fusions.vcf" % opts.base_output,shell=True)
                subprocess.call("cat /home/michael/YNHH/Reference_Files/IR/fusions.v4.4.vcf > %s.ionreporter.fusions.vcf" % opts.base_output,shell=True)
            else:
                print "Fusion VCF is NONE - this assay did not assess gene fusions"
                pass
        else:
            print "Fusion VCF local input"
            rename_fusion_vcf(opts.ionreporter_fusion_vcf,opts.base_output) # Rename fusions.vcf file as (basename).ionreporter.fusions.vcf


    samtools_command(SAMTOOLS_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,REGIONS_FILE)
    varscan_command(VARSCAN_EXE,opts.base_output,opts.base_output,opts.tumor_purity,opts.normal_purity)
    
    if opts.status=='All':
        vcf = '%s.snp.indel.sorted.vcf' % opts.base_output
    elif opts.status=='All-HC':
        vcf = '%s.snp.indel.sorted.hc.vcf' % opts.base_output
    elif opts.status=='Somatic':
        vcf = '%s.snp.indel.sorted.somatic.vcf' % opts.base_output
    elif opts.status=='LOH':
        vcf = '%s.snp.indel.sorted.LOH.vcf' % opts.base_output
    elif opts.status=='Germline':
        vcf = '%s.snp.indel.sorted.germline.vcf' % opts.base_output
    elif opts.status=='Somatic_LOH':
        vcf = '%s.snp.indel.sorted.somatic_LOH.vcf' % opts.base_output
    else:
        vcf = '%s.snp.indel.sorted.somatic_LOH.vcf' % opts.base_output    
    
    if opts.filter==True:
        
        varscan_filter = "((SS ='1')| (SS = '2')| (SS= '3'))  & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))"
        program = 'varscan'
        SnpSift_filter(vcf,SNPSIFT_EXE,varscan_filter,opts.base_output,program)
        VEP_command(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
        
        vcf = opts.ionreporter_vcf
        vcf = IR4_4_VCF_fix(vcf,opts.base_output)
        
        
        ionreporter_no_cnv_filter = "(HRUN[*] <= 6) & (FDP[*] >= 20) & (GEN[1].FDP[*] >= 5) & !(ALT='<CNV>')"
        ionreporter_cnv_filter = "(ALT='<CNV>')"
        program = 'ionreporter.no_cnv'
        SnpSift_filter(vcf,SNPSIFT_EXE,ionreporter_no_cnv_filter,opts.base_output, program)
        VEP_command(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
        program = 'ionreporter.cnv'
        SnpSift_filter(vcf,SNPSIFT_EXE,ionreporter_cnv_filter,opts.base_output, program)
        


        
    else:
        print "Unfiltered processing is not currently supported.  Please use filtering provided or contact author for custom filtering..."
        parser.print_help()
        sys.exit()
        
            
elif opts.mode==False:
    mandatory_options = ['ionreporter_vcf','ionreporter_tsv','varscan_vcf']
    for m in mandatory_options:
        # Making sure all mandatory options appeared
        if not opts.__dict__[m]:
            print "Mandatory option is missing!\n"
            parser.print_help()
            sys.exit()

    if opts.filter==True:
        
        vcf = opts.varscan_vcf
        
        varscan_filter = "((SS ='1')| (SS = '2')| (SS= '3'))  & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))"
        program = 'varscan'
        SnpSift_filter(vcf,SNPSIFT_EXE,varscan_filter,opts.base_output,program)
        VEP_command(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
        
        vcf = opts.ionreporter_vcf
        vcf = IR4_4_VCF_fix(vcf,opts.base_output) # Fix and rename IR VCF file
        
        ionreporter_no_cnv_filter = "(HRUN[*] <= 6) & (FDP[*] >= 20) & (GEN[1].FDP[*] >= 5) & !(ALT='<CNV')"
        ionreporter_cnv_filter = "(ALT='<CNV>')"
        program = 'ionreporter.no_cnv'
        SnpSift_filter(vcf,SNPSIFT_EXE,ionreporter_no_cnv_filter,opts.base_output, program)
        VEP_command(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
        try:
        	program = 'ionreporter.cnv'
        	SnpSift_filter(vcf,SNPSIFT_EXE,ionreporter_cnv_filter,opts.base_output, program)
        except:
            print "WTF: CNV filter didn't work!"
    else:
        print "Unfiltered processing is not currently supported.  Please use filtering provided or contact author for custom filtering..."
        parser.print_help()
        sys.exit()

### FINAL POST PROCESSING ### 

edit_IR_tsv_file(opts.ionreporter_version,opts.ionreporter_tsv,opts.base_output)
Move_files_1(opts.base_output)
Eddy_command(EDDY_EXE, opts.base_output)
Move_files_2(opts.base_output)
Galaxy_special_function(opts.base_output, opts.galaxy_output_directory)
zip_files(opts.base_output,opts.galaxy_output_directory)
delete_VEP_index()


html_out = open("%s" % opts.galaxy_html_file,"w")
for x in sorted(os.listdir("%s" % opts.galaxy_output_directory)):
    html_out.write("<a href='%s'>%s</a><br>" % (x,x))


