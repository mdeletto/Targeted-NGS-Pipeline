#/usr/bin/python

import optparse, sys, subprocess, datetime, os, pysam
from Variant_Detection_functions import *


#-------------------------------------------------------------------------------------------
#----------------------------Command line parser and arguments------------------------------
#-------------------------------------------------------------------------------------------

desc="""This wrapper script is intended to be used for variant detection and annotation."""

parser = optparse.OptionParser(description=desc)

parser.add_option('--variant_detection_mode', help='Enable variant detection mode and subsequent annotation', dest='mode', action='store_true',default=True)
parser.add_option('--annotation_mode', help='Enable annotation mode only', dest='mode', action='store_false',default=True)
parser.add_option('--enable_filtering', help="""Enable consequential filtering for VCFs.  This does not apply to hard filters like read depth (these will be automatically applied to remove false positives).""", dest='filter', action='store_true',default=False)
parser.add_option('-s', help='Check for classes of variants in variant detection mode <All,All-HC,Somatic,Germline,LOH,Somatic_LOH>', dest='status', action='store',default='Somatic_LOH')
parser.add_option("--ion_reporter_vcf", action="store", help="Input absolute path to IonReporter VCF </path/to/ionreporter.vcf>",dest="ionreporter_vcf")
parser.add_option("--ion_reporter_tsv", action="store", help="Input absolute path to IonReporter TSV </path/to/ionreporter.tsv>",dest="ionreporter_tsv")
parser.add_option('-c', help="Base output name for files.  If a copath ID exists, please use this as the base output name.", dest='base_output', action='store')
parser.add_option('--ionreporter_url_bool', help="Flag for if VCF is coming from IonReporter url", dest='ionreporter_url_bool', action='store_true',default=False)
parser.add_option('--ionreporter_analysis_name', help="Analysis name within IonReporter", dest='ionreporter_analysis_name', action='store')
parser.add_option('--ionreporter_id', help="Analysis ID within IonReporter", dest='ionreporter_id', action='store')#,default=opts.base_output)

parser.add_option('--url', help="Flag for if BAM is coming from url", dest='url', action='store_true',default=False)
parser.add_option('--regions', help="Regions to focus analysis on.  Must be in BED format", dest='regions', action='store')
parser.add_option("--ion_reporter_version", action="store", help="Please indicate IR version number (ex. 4.0, 4.4)", dest="ionreporter_version", default=4.4)

group = optparse.OptionGroup(parser, "Variant Detection Mode",
                    "Requires input of two bams--tumor and normal")
group.add_option("-t", action="store", help="Input absolute path to tumor bam </path/to/sample.bam>",dest="tumor")
group.add_option("-n", action="store", help="Input absolute path to normal bam---NOTE: If no normal bam is selected, a population normal will be used.",dest="normal",default=None)
group.add_option("--tumor-purity", action="store", help="Tumor purity <Float between 0 and 1>",dest="tumor_purity",default=1)
group.add_option("--normal-purity", action="store", help="Normal purity <Float between 0 and 1>",dest="normal_purity",default=1)
group.add_option("-p", action="store", help="Platform used for sequencing <PGM or Proton>",dest="platform",default="PGM") # option not currently supported
parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Annotation Mode",
                    "Requires input of VCF files from both IonReporter and VarScan")
group.add_option("-v", action="store", help="Input absolute path to VarScan VCF",dest="varscan_vcf")
parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Fusion options (for OCP)",
                    "ONLY FOR CASES WITH FUSION DETECTION")
group.add_option('--ionreporter_fusion_url_bool', help="Flag for if FUSION VCF is coming from IonReporter url", dest='ionreporter_fusion_url_bool', action='store_true',default=False)
group.add_option('--ionreporter_fusion_analysis_name', help="Fusion analysis name within IonReporter", dest='ionreporter_fusion_analysis_name', action='store')
group.add_option('--ionreporter_fusion_id', help="Fusion analysis ID within IonReporter", dest='ionreporter_fusion_id', action='store')#,default=opts.base_output)
group.add_option('--ionreporter_fusion_vcf', help="Fusion VCF within IonReporter", dest='ionreporter_fusion_vcf', action='store')
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
global VCFLIB_DIR

SAMTOOLS_EXE = "/home/michael/bin/samtools-1.1/samtools"
VARSCAN_EXE = "/home/michael/bin/VarScan/VarScan.v2.3.9.jar"
SNPSIFT_EXE = "/home/michael/bin/snpEff/SnpSift.jar"
VEP_EXE = "/home/michael/bin/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl"
VEP_FILTER_EXE = "/home/michael/bin/ensembl-tools-release-83/scripts/variant_effect_predictor/filter_vep.pl"
EDDY_EXE = "/home/michael/bin/eddy.jar"
BAM_READCOUNT_EXE = "/home/michael/bin/bam-readcount/bin/bin/bam-readcount"
VCFLIB_DIR = "/home/michael/bin/vcflib/bin"
MUTECT_EXE = "/home/michael/bin/mutect/muTect-1.1.7.jar"
GATK_EXE = "/home/michael/bin/GATK/GenomeAnalysisTK.jar"
STRELKA_EXE = '/home/michael/bin/strelka/bin/configureStrelkaWorkflow.pl'


VEP_REF_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.VEP.fasta"
REFERENCE_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.fasta"
REF_FASTA = REFERENCE_FASTA
dbsnp_vcf = '/home/michael/YNHH/Reference_Files/Common/Common/dbsnp-common.v142.ucsc.vcf'
cosmic_vcf = '/home/michael/YNHH/Reference_Files/Common/Common/COSMIC.v74.hg19.ucsc.sorted.vcf'
STRELKA_CONFIG_TARGETED = '/home/michael/bin/strelka/etc/strelka_config.iontorrent.targeted.ini'
STRELKA_CONFIG = STRELKA_CONFIG_TARGETED

#-------------------------------------------------------------------------------------------
#--------------------------------BEGIN PREPROCESSING----------------------------------------
#-------------------------------------------------------------------------------------------
write_logfile(opts) # Write command line parameters to logfile 
REGIONS_FILE = select_target_regions(opts.regions) # Select target regions.  Defaults to CCP if none

### COMMAND LINE CHECK FOR REMOTE IR VCF INPUT ###

if opts.ionreporter_url_bool == True:
    variant_link = IR_locate_variant_zip(opts.ionreporter_analysis_name,opts.ionreporter_id)
    IR_download_variant_zip(opts.ionreporter_analysis_name, variant_link,"somatic")
    opts.ionreporter_vcf = "%s.ionreporter.somatic_temp.vcf" % opts.ionreporter_analysis_name
    opts.ionreporter_tsv = "%s.ionreporter.somatic_temp.tsv" % opts.ionreporter_analysis_name
else:
    pass


if opts.mode==True:
#     mandatory_options = ['tumor','ionreporter_vcf','ionreporter_tsv','platform']
#     for m in mandatory_options:
#         # Making sure all mandatory options appeared
#         if not opts.__dict__[m]:
#             print "Mandatory option is missing!\n"
#             parser.print_help()
#             sys.exit()        

    ### COMMAND LINE CHECKS FOR REMOTE TORRENT-SUITE BAM INPUT ###

    if opts.url == True:
        opts.tumor = pull_BAM_from_url(opts.tumor,opts.base_output,"tumor")
        opts.normal = pull_BAM_from_url(opts.normal,opts.base_output,"normal")
    else:
        pass
    
    ### FUSION SUPPORT (added in v1.2) ####
        
    if opts.ionreporter_fusion_url_bool == True:
        print "Attempting to pull IR fusions"
        variant_link = IR_locate_variant_zip(opts.ionreporter_fusion_analysis_name,opts.ionreporter_fusion_id) # Determine download link
        IR_download_fusion_zip(opts.ionreporter_fusion_analysis_name, variant_link, opts.base_output) # Download analysis and process files
    else:
        if opts.ionreporter_fusion_vcf is None:
            print "Fusion VCF is NONE"
            pass
        else:
            print "Fusion VCF local input"
            rename_fusion_vcf(opts.ionreporter_fusion_vcf,opts.base_output) # Rename fusions.vcf file as (basename).ionreporter.fusions.vcf

    ### BAM-READCOUNTS FOR TUMOR AND NORMAL BAMS ###
    
    pysam.index(opts.tumor) # create index for tumor
    pysam.index(opts.normal) # create index for normal
    #tumor_readcounts = bam_readcount_command(BAM_READCOUNT_EXE,"tumor",opts.base_output,REFERENCE_FASTA,REGIONS_FILE,opts.tumor) # calculate read counts for tumor to be used in VarScan fpfilter
    #normal_readcounts = bam_readcount_command(BAM_READCOUNT_EXE,"normal",opts.base_output,REFERENCE_FASTA,REGIONS_FILE,opts.normal) # normal(DISABLED UNTIL LARGER PANEL)
    
    ### MPILEUP FOR NORMAL/TUMOR BAMS ###
     
    samtools_mpileup(SAMTOOLS_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,REGIONS_FILE)
     
    ### VARSCAN VARIANT CALLING AND POST-PROCESSING ###
     
    varscan_base_output = str(opts.base_output)+str(".varscan.tmp") 
    varscan_command(VARSCAN_EXE,opts.base_output,varscan_base_output,opts.tumor_purity,opts.normal_purity)
     
    ### SELECT VARSCAN SUBSET VCF FOR FURTHER PROCESSING (GERMLINE,LOH,SOMATIC,ETC.) ### 
     
    varscan_vcf = select_varscan_vcf_subset(opts,varscan_base_output)
    
    ### SOMATIC VARIANT CALLING ###
        #---MUTECT VARIANT CALLING---#
        
    #mutect_vcf = muTect_caller_command(MUTECT_EXE,dbsnp_vcf,cosmic_vcf,REFERENCE_FASTA,opts.normal,opts.tumor,opts.base_output)
    
        #---STRELKA VARIANT CALLING---#
    
#     strelka_vcfs = Strelka_somatic_variant_calling_command(STRELKA_EXE,opts.normal,opts.tumor,REF_FASTA,STRELKA_CONFIG,opts.base_output)    
#     strelka_vcf = combine_vcf(VCFLIB_DIR,strelka_vcfs, opts.base_output, "strelka.snvs.indels.somatic.raw")
    
    ### SELECT IONREPORTER VCF AND APPLY IR VERSION FIX IF NECESSARY
    
    ionreporter_vcf = opts.ionreporter_vcf
    ionreporter_vcf = IR4_4_VCF_fix(ionreporter_vcf,opts.base_output)
    #ionreporter_vcf = multibreak_vcf(VCFLIB_DIR,ionreporter_vcf,opts.base_output)
    
    ### PERFORM HARD-FILTERING AND VEP ANNOTATION
    
    program_filter_vcf_list = [
                               # first index = output basename suffix
                               # second index = SnpSift filter expression
                               # third index = vcf input
                               
                                ["varscan",  
                                 """((SS ='1')| (SS = '2')| (SS= '3'))  
                                     & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))""",
                                 varscan_vcf],
                                   
                               ["ionreporter.no_cnv",
                                """(HRUN[*] <= 6)
                                    & ((FDP[*] >= 20) | (DP[*] >= 20))
                                    & ((FAO[*] >= 2) | (AO[*] >= 2)) 
                                    & ((GEN[1].FDP[*] >= 5) | (GEN[1].DP[*] >= 5)) 
                                    & !(ALT='<CNV>')""", 
                                ionreporter_vcf],
                               
                               ["ionreporter.cnv",
                                "(ALT='<CNV>')",
                                ionreporter_vcf],
                               
#                                 ["mutect",
#                                  """(GEN[0].FA <= 0.1) 
#                                      & (GEN[1].FA >= 0.05) 
#                                      & (GEN[0].DP >= 5) 
#                                      & (GEN[1].DP >= 20)""",
#                                  mutect_vcf]

#                                 ["strelka",
#                                  """(GEN[0].DP[*] >= 5)
#                                      & (GEN[1].DP[*] >= 20)""",
#                                  strelka_vcf]
                               
                               ]
    
#     if opts.platform == "PGM":
#         gatk_germline_snp_vcf = GATK_snp_caller_command(GATK_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,dbsnp_vcf)
#         gatk_germline_indel_vcf = GATK_indel_caller_command(GATK_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,dbsnp_vcf)
#         vcf_list = [gatk_germline_snp_vcf, gatk_germline_indel_vcf]
#         gatk_vcf = combine_vcf(VCFLIB_DIR,vcf_list,opts.base_output,"gatk.germline.raw.vcf")
#         program_filter_vcf_list.append(
#                                        ["gatk.germline",
#                                         """(GEN[0].DP >= 20)
#                                             & (GEN[1].DP >= 5)
#                                             & (AF[*] >= 0.30)
#                                             & (QUAL >= 30)
#                                             & (FS >= 30)
#                                         """,
#                                         gatk_vcf]
#                                        )
    
    for group in program_filter_vcf_list:
        program,filter,vcf = (i for i in group)
        SnpSift_filter(vcf,SNPSIFT_EXE,filter,opts.base_output, program)
        if opts.filter is True:
            #VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)

            VEP_filter = [""""Consequence != intron_variant or 
                              Consequence != non_coding_transcript_exon_variant or
                              Consequence != 5_prime_UTR_variant or
                              Consequence != 3_prime_UTR_variant or
                              Consequence != synonymous_variant
                              " """]
    
            VEP_command_filtered(VEP_EXE,VEP_FILTER_EXE,REF_FASTA,VEP_filter,opts.base_output,program)
        else:
            VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
        
            
elif opts.mode==False:
    
    ### OPTION NOT LONG SUPPORTED ###
    print "ERROR: TERMINATION: Annotation ONLY option no longer supported"
    sys.exit(1)
    
#     mandatory_options = ['ionreporter_vcf','ionreporter_tsv','varscan_vcf']
#     for m in mandatory_options:
#         # Making sure all mandatory options appeared
#         if not opts.__dict__[m]:
#             print "Mandatory option is missing!\n"
#             parser.print_help()
#             sys.exit()
#             
#     varscan_vcf = opts.varscan_vcf
# 
#     ionreporter_vcf = opts.ionreporter_vcf
#     ionreporter_vcf = IR4_4_VCF_fix(ionreporter_vcf,opts.base_output)
#     
#     ### PERFORM HARD-FILTERING AND VEP ANNOTATION
#     
#     program_filter_vcf_list = [
#                                # first index = output basename suffix
#                                # second index = SnpSift filter expression
#                                # third index = vcf input
#                                
#                                ["varscan",  
#                                 """((SS ='1')| (SS = '2')| (SS= '3'))  
#                                     & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))""",
#                                 varscan_vcf],
#                                    
#                                ["ionreporter.no_cnv",
#                                 """(HRUN[*] <= 6)
#                                     & ((FDP[*] >= 20) | (DP[*] >= 20))
#                                     & ((FAO[*] >= 2) | (AO[*] >= 2)) 
#                                     & ((GEN[1].FDP[*] >= 5) | (GEN[1].DP[*] >= 5)) 
#                                     & !(ALT='<CNV')""", 
#                                 ionreporter_vcf],
#                                
#                                ["ionreporter.cnv",
#                                 "(ALT='<CNV>')",
#                                 ionreporter_vcf],
#                                
#                                
#                                
#                                ]
#     
#     for group in program_filter_vcf_list:
#         program,filter,vcf = (i for i in group)
#         SnpSift_filter(vcf,SNPSIFT_EXE,filter,opts.base_output, program)
#         if opts.filter is True:
#             VEP_filter = ["""--filter "Consequence is not intron_variant or 
#                                    Consequence match not UTR or
#                                    Consequence is not synonymous_variant or
#                                    Consequence match not non_coding
#                                    " """]
#             
#             
#             VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
# 
#             #VEP_filter = 
#             #VEP_command_filtered(VEP_EXE,VEP_FILTER_EXE,REF_FASTA,VEP_filter,opts.base_output,program)
#         else:
#             VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)

### FINAL POST PROCESSING ### 

edit_IR_tsv_file(opts.ionreporter_version,opts.ionreporter_tsv,opts.base_output)
Move_files_1(opts.base_output)
#Eddy_command(EDDY_EXE, opts.base_output)
Move_files_2(opts.base_output)
Galaxy_special_function(opts.base_output, opts.galaxy_output_directory)
#zip_files(opts.base_output,opts.galaxy_output_directory)
delete_VEP_index()


html_out = open("%s" % opts.galaxy_html_file,"w")
for x in sorted(os.listdir("%s" % opts.galaxy_output_directory)):
    html_out.write("<a href='%s'>%s</a><br>" % (x,x))


