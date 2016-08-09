#!/usr/bin/python

import optparse, sys, subprocess, datetime, os, pysam, ast
from Variant_Detection_functions import *
import urllib


#-------------------------------------------------------------------------------------------
#----------------------------Command line parser and arguments------------------------------
#-------------------------------------------------------------------------------------------

desc="""A pipeline for somatic, germline, and LOH variant detection and annotation for both IonTorrent and Illumina data."""

parser = optparse.OptionParser(description=desc,version="%prog 1.2.1")

#parser.add_option('--variant_detection_mode', help='Enable variant detection mode and subsequent annotation', dest='mode', action='store_true',default=True)
#parser.add_option('--annotation_mode', help='Enable annotation mode only', dest='mode', action='store_false',default=True)
parser.add_option('--disable_filtering', help="""Disable consequential filtering for VCFs.  This does not apply to hard filters like read depth (these will be automatically applied to remove false positives).""", dest='filter_disabled', action='store',default=False)
parser.add_option('-s', help='Check for classes of variants in variant detection mode <All,All-HC,Somatic,Germline,LOH,Somatic_LOH>', dest='status', action='store',default='All-HC')
parser.add_option('-c', help="Base output name for files.  If a copath ID exists, please use this as the base output name.", dest='base_output', action='store')
parser.add_option('--url', help="Flag for if BAM is coming from url", dest='url', action='store_true',default=False)
parser.add_option('--regions', help="Regions to focus analysis on.  Must be in BED format", dest='regions', action='store')
parser.add_option("-t", action="store", help="Input absolute path to tumor bam </path/to/sample.bam>",dest="tumor",default=None)
parser.add_option("-n", action="store", help="Input absolute path to normal bam---NOTE: If no normal bam is selected, a population normal will be used.",dest="normal",default=None)
parser.add_option("-p", action="store", help="Platform used for sequencing",dest="platform",default="IonTorrent") # option not currently supported
parser.add_option("--galaxy", action="store_true", help="Apply galaxy-specific commands if this flag is supplied",dest="galaxy_flag",default=False)

group = optparse.OptionGroup(parser, "IonReporter inputs (for IonTorrent analyses)",
                             "IonReporter-specific inputs for IonTorrent data only")
group.add_option("--ionreporter_version", action="store", help="Please indicate IR version number (ex. 4.0, 4.4)", dest="ionreporter_version", default=4.4)
group.add_option('--ionreporter_somatic_url_bool', help="Flag for if VCF is coming from IonReporter url", dest='ionreporter_somatic_url_bool', action='store_true',default=False)
group.add_option('--ionreporter_somatic_analysis_name', help="Analysis name within IonReporter", dest='ionreporter_somatic_analysis_name', action='store')
group.add_option('--ionreporter_somatic_id', help="Analysis ID within IonReporter", dest='ionreporter_somatic_id', action='store')#,default=opts.base_output)
group.add_option("--ionreporter_somatic_vcf", action="store", help="Input absolute path to IonReporter VCF </path/to/ionreporter.vcf>",dest="ionreporter_somatic_vcf")
group.add_option("--ionreporter_somatic_tsv", action="store", help="Input absolute path to IonReporter TSV </path/to/ionreporter.tsv>",dest="ionreporter_somatic_tsv")
group.add_option('--ionreporter_fusion_url_bool', help="Flag for if FUSION VCF is coming from IonReporter url", dest='ionreporter_fusion_url_bool', action='store_true',default=False)
group.add_option('--ionreporter_fusion_analysis_name', help="Fusion analysis name within IonReporter", dest='ionreporter_fusion_analysis_name', action='store')
group.add_option('--ionreporter_fusion_id', help="Fusion analysis ID within IonReporter", dest='ionreporter_fusion_id', action='store')#,default=opts.base_output)
group.add_option('--ionreporter_fusion_vcf', help="Fusion VCF within IonReporter", dest='ionreporter_fusion_vcf', action='store')
group.add_option('--ionreporter_germline_url_bool', help="Flag for if VCF is coming from IonReporter url", dest='ionreporter_germline_url_bool', action='store_true',default=False)
group.add_option('--ionreporter_germline_analysis_name', help="Analysis name within IonReporter", dest='ionreporter_germline_analysis_name', action='store')
group.add_option('--ionreporter_germline_id', help="Analysis ID within IonReporter", dest='ionreporter_germline_id', action='store')#,default=opts.base_output)
group.add_option("--ionreporter_germline_vcf", action="store", help="Input absolute path to IonReporter VCF </path/to/ionreporter.vcf>",dest="ionreporter_germline_vcf")

parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Galaxy options",
                    "ONLY TO BE USED IN A GALAXY ENVIRONMENT")
group.add_option("--galaxy_html_file", action="store", help="Path to galaxy output file---to be used only in Galaxy environment",dest="galaxy_html_file")
group.add_option("--output_directory", action="store", help="Path to output directory",dest="output_directory")
parser.add_option_group(group)



(opts, args) = parser.parse_args()

#-------------------------------------------------------------------------------------------
#------------------------PATHS TO PROGRAMS/FILES USED BY PIPELINE---------------------------
#-------------------------------------------------------------------------------------------
global VCFLIB_DIR, IR_API_KEY

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
COVERAGE_ANALYSIS_EXE = '/home/michael/YNHH/Code/Github-mdeletto/TorrentTools/tumor_normal_coverage_analysis.py'

VEP_REF_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.VEP.fasta"
REFERENCE_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.fasta"
REF_FASTA = REFERENCE_FASTA
dbsnp_vcf = '/home/michael/YNHH/Reference_Files/Common/Common/dbsnp-common.v142.ucsc.vcf'
cosmic_vcf = '/home/michael/YNHH/Reference_Files/Common/Common/COSMIC.v74.hg19.ucsc.sorted.vcf'
STRELKA_CONFIG_TARGETED = '/home/michael/bin/strelka/etc/strelka_config.iontorrent.targeted.ini'
STRELKA_CONFIG = STRELKA_CONFIG_TARGETED
POPULATION_NORMAL_BAM_OCP = '/home/michael/YNHH/Reference_Files/Population_Normal/Population_Normal.OCP.v1.bam'
POPULATION_NORMAL_BAM_CCP = '/home/michael/YNHH/Reference_Files/Population_Normal/Population_Normal_PGM_v2.CCP.bam'

IR_API_KEY = "UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY"

def main():
    header = "---STARTING TPL Targeted NGS Pipeline v1.2.1---"
    print "-" * len(header)
    print header
    print "-" * len(header)

    #-------------------------------------------------------------------------------------------
    #--------------------------------BEGIN PREPROCESSING----------------------------------------
    #-------------------------------------------------------------------------------------------
    
    write_logfile(opts) # Write command line parameters to logfile 
    REGIONS_FILE = select_target_regions(opts.regions) # Select target regions.  Defaults to CCP if none
    
    # Url encode ionreporter analysis names that have excess whitespace
    
    if opts.ionreporter_somatic_url_bool is True:
        url_encoded_ionreporter_somatic_analysis_name = urllib.quote_plus(opts.ionreporter_somatic_analysis_name)
    if opts.ionreporter_fusion_url_bool is True:
        url_encoded_ionreporter_fusion_analysis_name = urllib.quote_plus(opts.ionreporter_fusion_analysis_name)

    #basename_minus_whitespace = opts.ionreporter_somatic_analysis_name.replace(" ","_") # IR converts any whitespace characters in <analysis_name> to underscores
    
    ### COMMAND LINE CHECK FOR REMOTE IR VCF INPUT ###
    
    if opts.ionreporter_somatic_url_bool == True:
        print "Remote download of SOMATIC VCF from IonReporter initiated.  Please hold...",
        try:
            variant_link = IR_locate_variant_zip(url_encoded_ionreporter_somatic_analysis_name,opts.ionreporter_somatic_id)
            IR_download_variant_zip(opts.base_output, variant_link,"somatic")
            opts.ionreporter_somatic_vcf = "%s.ionreporter.somatic_temp.vcf" % opts.base_output
            opts.ionreporter_somatic_tsv = "%s.ionreporter.somatic_temp.tsv" % opts.base_output
            print "SUCCESS"
        except:
            print "FAILED"
            sys.exit(1)
    else:
        pass

    ### FUSION SUPPORT (added in v1.2) ####
        
    if opts.ionreporter_fusion_url_bool == True:
        print "Remote download of FUSION VCF from IonReporter initiated.  Please hold...",
        try:
            variant_link = IR_locate_variant_zip(opts.ionreporter_fusion_analysis_name,opts.ionreporter_fusion_id) # Determine download link
            IR_download_fusion_zip(opts.ionreporter_fusion_analysis_name, variant_link, opts.base_output) # Download analysis and process files
            print "SUCCESS"
        except:
            print "FAILED"
            sys.exit(1)
    else:
        if opts.ionreporter_fusion_vcf is None or re.search("None", opts.ionreporter_fusion_vcf):
            if opts.regions == "OCP":
                #subprocess.call("touch %s.ionreporter.fusions.vcf" % opts.base_output,shell=True)
                subprocess.call("cat /home/michael/YNHH/Reference_Files/IR/fusions.v4.4.vcf > %s.ionreporter.fusions.vcf" % opts.base_output,shell=True)
            else:
                print "Fusion VCF is NONE - this assay did not assess gene fusions"
                pass
        else:
            print "Fusion VCF local input"
            rename_fusion_vcf(opts.ionreporter_fusion_vcf,opts.base_output) # Rename fusions.vcf file as (basename).ionreporter.fusions.vcf

#     mandatory_options = ['tumor','ionreporter_vcf','ionreporter_tsv','platform']
#     for m in mandatory_options:
#         # Making sure all mandatory options appeared
#         if not opts.__dict__[m]:
#             print "Mandatory option is missing!\n"
#             parser.print_help()
#             sys.exit()        

    ### SELECT IONREPORTER VCF AND APPLY IR VERSION FIX IF NECESSARY
    
    ionreporter_somatic_vcf = opts.ionreporter_somatic_vcf
    ionreporter_somatic_vcf = IR4_4_VCF_fix(ionreporter_somatic_vcf,opts.base_output)
    ionreporter_somatic_vcf = multibreak_vcf(VCFLIB_DIR,ionreporter_somatic_vcf,opts.base_output)
    #ionreporter_vcf = multibreak_vcf(VCFLIB_DIR,ionreporter_vcf,opts.base_output)
    
    ### PERFORM HARD-FILTERING
    
    program_filter_vcf_list = [
                               # first index = output basename suffix
                               # second index = SnpSift filter expression
                               # third index = vcf input

                                   
                               ["ionreporter.no_cnv",
                                """((FDP[*] >= 20) | (DP[*] >= 20))
                                    & ((FAO[*] >= 2) | (AO[*] >= 2)) 
                                    & ((GEN[1].FDP[*] >= 0) | (GEN[1].DP[*] >= 0)) 
                                    & !(ALT='<CNV>')""", 
                                ionreporter_somatic_vcf],
                               
                               ["ionreporter.cnv",
                                "(ALT='<CNV>')",
                                ionreporter_somatic_vcf],
                               
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


    ### COMMAND LINE CHECKS FOR REMOTE TORRENT-SUITE BAM INPUT ###

    
    if opts.tumor is None or opts.tumor == "None":
        subprocess.call("touch %s.varscan.json && touch %s.varscan.vcf" % (opts.base_output,opts.base_output), shell=True)
        varscan_vcf = "%s.varscan.vcf" % opts.base_output
    else:
        if opts.url is False:
            if re.search("Population", opts.normal):
                if opts.regions == "CCP":
                    opts.normal = POPULATION_NORMAL_BAM_CCP
                else:
                    opts.normal = POPULATION_NORMAL_BAM_OCP
        else:
            opts.tumor = pull_BAM_from_url(opts.tumor,opts.base_output,"tumor")
            if re.search("Population", opts.normal):
                if opts.regions == "CCP":
                    opts.normal = POPULATION_NORMAL_BAM_CCP
                else:
                    opts.normal = POPULATION_NORMAL_BAM_OCP
            else:
                opts.normal = pull_BAM_from_url(opts.normal,opts.base_output,"normal")
        
        ### CREATE BAM INDEXES ###
        
        pysam.index(opts.tumor) # create index for tumor
        pysam.index(opts.normal) # create index for normal
        
        ### PERFORM COVERAGE ANALYSIS WITH tumor_normal_coverage_analysis.py ###
        
        perform_coverage_analysis(COVERAGE_ANALYSIS_EXE, opts.base_output, opts.tumor, opts.normal, REGIONS_FILE)
        
        ### MPILEUP FOR NORMAL/TUMOR BAMS ###
        
        samtools_mpileup(SAMTOOLS_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,REGIONS_FILE)
        
        ### VARSCAN VARIANT CALLING AND POST-PROCESSING ###
        
        varscan_base_output = str(opts.base_output)+str(".varscan.tmp") 
        varscan_command(VARSCAN_EXE,opts.base_output,varscan_base_output)
        
        ### SELECT VARSCAN SUBSET VCF FOR FURTHER PROCESSING (GERMLINE,LOH,SOMATIC,ETC.) ### 
        
        varscan_vcf = select_varscan_vcf_subset(opts,varscan_base_output)

    



    ### BAM-READCOUNTS FOR TUMOR AND NORMAL BAMS ###

    #tumor_readcounts = bam_readcount_command(BAM_READCOUNT_EXE,"tumor",opts.base_output,REFERENCE_FASTA,REGIONS_FILE,opts.tumor) # calculate read counts for tumor to be used in VarScan fpfilter
    #normal_readcounts = bam_readcount_command(BAM_READCOUNT_EXE,"normal",opts.base_output,REFERENCE_FASTA,REGIONS_FILE,opts.normal) # normal(DISABLED UNTIL LARGER PANEL)
    
    
     
    
     

     
    
     
    
    
    ### SOMATIC VARIANT CALLING ###
        #---MUTECT VARIANT CALLING---#
        
    #mutect_vcf = muTect_caller_command(MUTECT_EXE,dbsnp_vcf,cosmic_vcf,REFERENCE_FASTA,opts.normal,opts.tumor,opts.base_output)
    
        #---STRELKA VARIANT CALLING---#
    
#     strelka_vcfs = Strelka_somatic_variant_calling_command(STRELKA_EXE,opts.normal,opts.tumor,REF_FASTA,STRELKA_CONFIG,opts.base_output)    
#     strelka_vcf = combine_vcf(VCFLIB_DIR,strelka_vcfs, opts.base_output, "strelka.snvs.indels.somatic.raw")
    

    
    if str(opts.tumor) == "None" or opts.tumor is None:
        pass
    else:
         
        varscan_filter = ["varscan",  
                          """((SS ='1') | (SS = '2') | (SS= '3') )  
                              & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))""",
                              #& !(exists INDEL)""", # filter out INDELs
                           varscan_vcf]
         
        program_filter_vcf_list.append(varscan_filter)    
    
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
    
    delete_VEP_index()
    
    for group in program_filter_vcf_list:
        program,filter,vcf = (i for i in group)
        SnpSift_filter(vcf,SNPSIFT_EXE,filter,opts.base_output, program)
        if program == "ionreporter.cnv":  # ionreporter.cnv.vcf needs to be handled via unfiltered method or will cause memory issues
            pass # for now, pass this because long CNVs add lengthy processing time by VEP
            #VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
        else:   
            
            def variant_count(vcf):
                # Counts lines without a '#' character starting the line
                # Used to determine whether to pass through VEP or not
                variant_lines = subprocess.check_output("grep -v ^# %s | wc -l" % vcf,shell=True)
                return variant_lines
            
            if variant_count(vcf) > 1:  # Test vcf to see if any variant lines exist, or create blank .json file if no lines exist
                if ast.literal_eval(str(opts.filter_disabled)) is False:
                    print "WARNING: Running consequential filters.  Only nonsynonymous coding variants will be reported!"

                    VEP_command_filtered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
                else:
                    VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
            elif variant_count(vcf) == 0:
                subprocess.call("touch %s.%s.json" % (opts.base_output, program),shell=True)


    ### FINAL POST PROCESSING ### 
    
    edit_IR_tsv_file(opts.ionreporter_version,opts.ionreporter_somatic_tsv,opts.base_output)
    
    Move_files_1(EDDY_EXE,opts.base_output,opts.galaxy_flag)
    
    extra_file_cleanup()
    
    if opts.galaxy_flag is True:
        Galaxy_special_function(opts.base_output, opts.output_directory)
        zip_files(opts.base_output,opts.output_directory)
        html_out = open("%s" % opts.galaxy_html_file,"w")
        for file in sorted(os.listdir("%s" % opts.output_directory)):
            html_out.write("<a href='%s'>%s</a><br>" % (file,file))            
    else:
        opts.output_directory = opts.base_output
        zip_files(opts.base_output,opts.output_directory)

if __name__ == "__main__":
    main()