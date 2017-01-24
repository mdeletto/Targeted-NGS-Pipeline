#!/usr/bin/python

import datetime
import optparse
import os
import pysam
import sys
import subprocess
import urllib
from Variant_Detection_functions import *


#-------------------------------------------------------------------------------------------
#----------------------------Command line parser and arguments------------------------------
#-------------------------------------------------------------------------------------------

desc="""A pipeline for somatic, germline, and LOH variant detection and annotation for both IonTorrent and Illumina data."""
vers="1.3.0"

parser = optparse.OptionParser(description=desc,version=vers)

# GENERAL OPTIONS

parser.add_option('--disable_filtering', help="""Enable consequential filtering for VCFs.  This does not apply to hard filters like read depth (these will be automatically applied to remove false positives).""", dest='filter_disabled', action='store_true',default=False)
parser.add_option('-c', help="Base output name for files.  If a copath ID exists, please use this as the base output name.", dest='base_output', action='store')
parser.add_option('--url', help="Flag for if BAM is coming from url", dest='url', action='store_true',default=False)
parser.add_option('--regions', help="Regions to focus analysis on.  Must be in BED format", dest='regions', action='store')
parser.add_option("-t", action="store", help="Input absolute path to tumor bam </path/to/sample.bam>",dest="tumor")
parser.add_option("-n", action="store", help="Input absolute path to normal bam---NOTE: If no normal bam is selected, a population normal will be used.",dest="normal",default=None)
parser.add_option("-p", action="store", help="Platform used for sequencing",dest="platform",default="IonTorrent") # option not currently supported

# IONREPORTER OPTIONS

group = optparse.OptionGroup(parser, "IonReporter inputs (for IonTorrent analyses)",
                             "IonReporter-specific inputs for IonTorrent data only")
group.add_option("--ionreporter_only", action="store_true", help="Only run IR analyses", dest="ionreporter_only", default=False)
group.add_option('--ionreporter_skip_all', help="Skip all IR analyses", dest='ionreporter_skip_all', action='store_true',default=False)
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

# MISC OPTIONS

group = optparse.OptionGroup(parser, "Miscellaneous options",
                             "Tells the pipeline how to handle other options, to allow the user finer control over output.")
group.add_option("--do-not-separate-loh", action="store_true", help="Boolean FLAG to keep GERMLINE AND LOH calls in the same VCF rather than splitting them", dest="do_not_separate_LOH", default=False)
parser.add_option_group(group)

# GALAXY OPTIONS
# Allows for the tool to be plugged into a galaxy environment (best used for testing purposes and manual analyses)

group = optparse.OptionGroup(parser, "Galaxy options",
                    "ONLY TO BE USED IN A GALAXY ENVIRONMENT")
group.add_option("--galaxy", action="store_true", help="Apply galaxy-specific commands if this flag is supplied",dest="galaxy_flag",default=False)
group.add_option("--galaxy_html_file", action="store", help="Path to galaxy output file---to be used only in Galaxy environment",dest="galaxy_html_file")
group.add_option("--output_directory", action="store", help="Path to output directory",dest="output_directory")
parser.add_option_group(group)


# LEGACY OPTIONS (UNSUPPORTED)
# Variables for VarScan variant calling.  These are disabled anyway, so there will be no effect.

group = optparse.OptionGroup(parser, "LEGACY OPTIONS",
                    "Old VarScan variant calling options.  These are disabled anyway...")
group.add_option('-s', help='Check for classes of variants in variant detection mode <All,All-HC,Somatic,Germline,LOH,Somatic_LOH>', dest='status', action='store',default='Somatic_LOH')
parser.add_option_group(group)

(opts, args) = parser.parse_args()

#-------------------------------------------------------------------------------------------
#------------------------PATHS TO PROGRAMS/FILES USED BY PIPELINE---------------------------
#-------------------------------------------------------------------------------------------
global VCFLIB_DIR, IR_API_KEY

# EXECUTABLES

SAMTOOLS_EXE = "/home/michael/bin/samtools-1.1/samtools"
BEDTOOLS_EXE = "/home/michael/bin/bedtools2/bin/bedtools"
VARSCAN_EXE = "/home/michael/bin/VarScan/VarScan.v2.3.9.jar"
SNPSIFT_EXE = "/home/michael/bin/snpEff/SnpSift.jar"
VEP_EXE = "/home/michael/bin/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl"
VEP_FILTER_EXE = "/home/michael/bin/ensembl-tools-release-83/scripts/variant_effect_predictor/filter_vep.pl"
EDDY_EXE = "/home/michael/bin/eddy.jar"
BAM_READCOUNT_EXE = "/home/michael/bin/bam-readcount/bin/bin/bam-readcount"
VCFLIB_DIR = "/home/michael/bin/vcflib/bin"
MUTECT_EXE = "/home/michael/bin/mutect/muTect-1.1.7.jar"
GATK_EXE = "/home/michael/bin/GATK/GenomeAnalysisTK.jar"
GATK_LATEST_EXE = "/home/michael/bin/GATK-3.6/GenomeAnalysisTK.jar"
STRELKA_EXE = '/home/michael/bin/strelka/bin/configureStrelkaWorkflow.pl'

# REFERENCE FILES

VEP_REF_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.VEP.fasta"
REFERENCE_FASTA = "/home/michael/YNHH/Reference_Files/FASTA/hg19.fasta"
REF_FASTA = REFERENCE_FASTA
dbsnp_vcf = '/home/michael/YNHH/Reference_Files/Common/Common/dbsnp-common.v142.ucsc.vcf'
cosmic_vcf = '/home/michael/YNHH/Reference_Files/Common/Common/COSMIC.v74.hg19.ucsc.sorted.vcf'
STRELKA_CONFIG_TARGETED = '/home/michael/bin/strelka/etc/strelka_config.iontorrent.targeted.ini'
STRELKA_CONFIG = STRELKA_CONFIG_TARGETED
MUTECT_V1_PON = "/home/michael/Development/MuTect/OCP.mutect.PON.vcf"
MUTECT2_PON_OCP = "/home/michael/YNHH/Reference_Files/tool-reference-files/MuTect2/PON/OCP/OCP.mutect2.PON.vcf"
MUTECT2_PON_CCP = "/home/michael/YNHH/Reference_Files/tool-reference-files/MuTect2/PON/CCP/CCP.mutect2.PON.vcf"

# POPULATION NORMALS

POPULATION_NORMAL_BAM_OCP = '/home/michael/YNHH/Reference_Files/Population_Normal/PopulationNormal_OCP_v2.bam'
#POPULATION_NORMAL_BAM_OCP = '/home/michael/YNHH/Reference_Files/Population_Normal/simulated/Population_Normal.OCP.v1.simulated.bam'
POPULATION_NORMAL_BAM_CCP = '/home/michael/YNHH/Reference_Files/Population_Normal/Population_Normal_PGM_v2.CCP.bam'

# OTHER GLOBAL VARIABLES

IR_API_KEY = "UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY"

def main():
    header = "---STARTING TPL Targeted NGS Pipeline v1.3---"
    print "-" * len(header)
    print header
    print "-" * len(header)    

    #-------------------------------------------------------------------------------------------
    #--------------------------------BEGIN PREPROCESSING----------------------------------------
    #-------------------------------------------------------------------------------------------
    
    # Write command line parameters to logfile 
    write_logfile(opts) 

    # Attempt to autodetect sample attributes
    sample_attribute_autodetection(opts.base_output, vers, opts.regions)

    # Select target regions.  Defaults to CCP if none
    REGIONS_FILE = select_target_regions(opts.regions)
    
    
    # Check remote TS BAM input
    if opts.tumor is None or opts.tumor == "None" and opts.ionreporter_only is True:
        subprocess.call("touch %s.varscan.json && touch %s.varscan.vcf" % (opts.base_output,opts.base_output), shell=True)
        varscan_vcf = "%s.varscan.vcf" % opts.base_output
    else:

        if opts.url is False:
            # If input is local and Population Normal exists in the name for the normal, use Population Normal BAM as normal and select appropriate PON for MuTect2
            if re.search("Population", opts.normal):
                if opts.regions == "CCP":
                    opts.normal = POPULATION_NORMAL_BAM_CCP
                else:
                    opts.normal = POPULATION_NORMAL_BAM_OCP
        else:
            opts.tumor = pull_BAM_from_url(opts.tumor,opts.base_output,"tumor")
            # If input is remote and Population Normal exists in the name for the normal, use Population Normal BAM as normal and select appropriate PON for MuTect2
            # We don't waste unnecessary time pulling in a BAM that is already stored on the server, plus the filesize for the PopNorm is very big.
            if re.search("Population", opts.normal):
                if opts.regions == "CCP":
                    opts.normal = POPULATION_NORMAL_BAM_CCP
                else:
                    opts.normal = POPULATION_NORMAL_BAM_OCP
            else:
                opts.normal = pull_BAM_from_url(opts.normal,opts.base_output,"normal")
    
    # Create index for BAMs
    pysam.index(opts.tumor) # create index for tumor
    pysam.index(opts.normal) # create index for normal
    
    # Create normal-tumor mpileup (RETIRED WITH VARSCAN in v1.3)
    #samtools_mpileup(SAMTOOLS_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,REGIONS_FILE)
    
    #---COMMAND LINE CHECK FOR REMOTE IF VCF INPUT---#
    if opts.ionreporter_skip_all is False:
        
        if opts.ionreporter_germline_url_bool == True:
            print "Remote download of GERMLINE VCF from IonReporter initiated.  Please hold..."
            url_encoded_ionreporter_germline_analysis_name = urllib.quote_plus(opts.ionreporter_germline_analysis_name)
            try:
                # Determine download link
                link_and_barcodes_list = IR_locate_germline_variant_zip(url_encoded_ionreporter_germline_analysis_name, opts.ionreporter_germline_id)
                variant_link, sample_barcode, control_barcode = (i for i in link_and_barcodes_list)
                germline_sample_dict = {"sample": sample_barcode,
                                        "control": control_barcode}
                # Download analysis and process files
                ionreporter_germline_unfiltered_vcf = IR_download_germline_variant_zip(VCFLIB_DIR, opts.base_output, variant_link, "germline", germline_sample_dict)
                print "SUCCESSFULLY DOWNLOADED GERMLINE IR VCF"
            except Exception,e:
                print str(e)
                print "FAILED TO DOWNLOAD GERMLINE IR VCF"
                sys.exit(1)
        else:
            pass
        
        if opts.ionreporter_somatic_url_bool == True:
            
            
            
            print "Remote download of SOMATIC VCF from IonReporter initiated.  Please hold..."
            # Url encode ionreporter analysis names that have excess whitespace
            url_encoded_ionreporter_somatic_analysis_name = urllib.quote_plus(opts.ionreporter_somatic_analysis_name)
            try:
                # Determine download link
                variant_link = IR_locate_variant_zip(url_encoded_ionreporter_somatic_analysis_name, opts.ionreporter_somatic_id)
                # Download analysis and process files
                ionreporter_somatic_vcf_and_tsv = IR_download_somatic_variant_zip(opts.base_output, variant_link, "somatic")
                ionreporter_somatic_unfiltered_vcf, ionreporter_somatic_tsv = (i for i in ionreporter_somatic_vcf_and_tsv)
                opts.ionreporter_somatic_vcf = ionreporter_somatic_unfiltered_vcf
                opts.ionreporter_somatic_tsv = ionreporter_somatic_tsv
                print "SUCCESSFULLY DOWNLOADED SOMATIC IR VCF"
                
                edit_IR_tsv_file(opts.ionreporter_version,opts.ionreporter_somatic_tsv,opts.base_output)
                
            except Exception, e:
                print str(e)
                print "FAILED TO DOWNLOAD SOMATIC IR VCF"
                sys.exit(1)
        else:
            pass
    
        if opts.ionreporter_fusion_url_bool == True:
            print "Remote download of FUSION VCF from IonReporter initiated.  Please hold..."
            # Url encode ionreporter analysis names that have excess whitespace
            url_encoded_ionreporter_fusion_analysis_name = urllib.quote_plus(opts.ionreporter_fusion_analysis_name)
            try:
                # Determine download link
                variant_link = IR_locate_variant_zip(url_encoded_ionreporter_fusion_analysis_name, opts.ionreporter_fusion_id) 
                # Download analysis and process files
                IR_download_fusion_zip(variant_link, opts.base_output) 
                print "SUCCESSFULLY DOWNLOADED FUSION IR VCF"
            except:
                print "FAILED TO DOWNLOAD FUSION IR VCF"
                sys.exit(1)
        else:
            if opts.ionreporter_fusion_vcf is None:
                if opts.regions == "OCP":
                    # Create pseudo fusion.vcf for OCP cases run without RNA
                    # Contains a header picked up by Downstream for insufficient RNA
                    subprocess.call("cat /home/michael/YNHH/Reference_Files/IR/fusions.v4.4.vcf > %s.ionreporter.fusions.vcf" % opts.base_output,shell=True)
                else:
                    print "Fusion VCF is NONE - this assay did not assess gene fusions"
                    pass
            else:
                print "Fusion VCF local input"
                rename_fusion_vcf(opts.ionreporter_fusion_vcf,opts.base_output) # Rename fusions.vcf file as (basename).ionreporter.fusions.vcf
    
        # Check for ionreporter.fusions.vcf
        
        ###  PROCESS FUSIONS VCF AND EXTRACT INFO ###
        
        if os.path.isfile('./%s.ionreporter.fusions.vcf' % opts.base_output):
            fusion_dict = extract_fusion_VCF_information('./%s.ionreporter.fusions.vcf' % opts.base_output)
            print fusion_dict
            # If gene expression counts exist, we will need to create a separate counts file.
            # This file will be used for differential expression analysis
            if fusion_dict['gene_expression_read_counts']:
                with open('%s.gene_expression.counts.tsv') as gene_expression_out:
                    gene_expression_out.write("Gene\t%s\n" % opts.base_output)
                    for k in fusion_dict['gene_expression_read_counts'].keys():
                        gene_expression_out.write("%s\t%s|n" % (k, fusion_dict['gene_expression_read_counts'][k]))
    
    
        #---SELECT IONREPORTER VCF AND APPLY IR VERSION FIX IF NECESSARY---#
        
        ionreporter_somatic_unfiltered_vcf = opts.ionreporter_somatic_vcf
        ionreporter_somatic_unfiltered_vcf = IR4_4_VCF_fix(ionreporter_somatic_unfiltered_vcf,opts.base_output)
        # Remove the tmp IR somatic VCF
        subprocess.call("rm %s" % opts.ionreporter_somatic_vcf, shell=True)
    

    
    #-------------------------------------------------------------------------------------------
    #----------------------------SOMATIC VARIANT CALLING----------------------------------------
    #-------------------------------------------------------------------------------------------
 
    #tumor_readcounts = bam_readcount_command(BAM_READCOUNT_EXE,"tumor",opts.base_output,REFERENCE_FASTA,REGIONS_FILE,opts.tumor) # calculate read counts for tumor to be used in VarScan fpfilter
    #normal_readcounts = bam_readcount_command(BAM_READCOUNT_EXE,"normal",opts.base_output,REFERENCE_FASTA,REGIONS_FILE,opts.normal) # normal(DISABLED UNTIL LARGER PANEL)
    
    #---MPILEUP FOR NORMAL/TUMOR BAMS (RETIRED)---#
     
    #samtools_mpileup(SAMTOOLS_EXE,opts.tumor,opts.normal,opts.base_output,REFERENCE_FASTA,REGIONS_FILE)
     
    #---VARSCAN VARIANT CALLING AND POST-PROCESSING (RETIRED)---#
     
    #varscan_base_output = str(opts.base_output)+str(".varscan.tmp") 
    #varscan_command(VARSCAN_EXE,opts.base_output,varscan_base_output)
    #varscan_vcf = select_varscan_vcf_subset(opts,varscan_base_output)

    #---MUTECT VARIANT CALLING---#
         
    #mutect_vcf = muTect_caller_command(MUTECT_EXE,REGIONS_FILE,MUTECT_V1_PON,dbsnp_vcf,cosmic_vcf,REFERENCE_FASTA,opts.normal,opts.tumor,opts.base_output)
    if opts.regions == "CCP":
        MUTECT2_PON = MUTECT2_PON_CCP
        mutect2_unfiltered_vcf = muTect2_caller_command(GATK_LATEST_EXE,REGIONS_FILE,MUTECT2_PON,dbsnp_vcf,cosmic_vcf,REFERENCE_FASTA,opts.normal,opts.tumor,opts.base_output)
    elif opts.regions in ["CHPv2","HSM","OCP","OCA"]:
        MUTECT2_PON = MUTECT2_PON_OCP
        mutect2_unfiltered_vcf = muTect2_caller_command(GATK_LATEST_EXE,REGIONS_FILE,MUTECT2_PON,dbsnp_vcf,cosmic_vcf,REFERENCE_FASTA,opts.normal,opts.tumor,opts.base_output)
    else:
        if opts.ionreporter_only is False:
            
            sys.exit("ERROR: Selected panel (%s) does not have PON for Mutect2" % (opts.regions))
             
    #---STRELKA VARIANT CALLING---#
    
    if opts.ionreporter_only is False:
        strelka_unfiltered_vcfs = Strelka_somatic_variant_calling_command(STRELKA_EXE,opts.normal,opts.tumor,REF_FASTA,STRELKA_CONFIG,opts.base_output)    
        strelka_unfiltered_vcf = combine_vcf(VCFLIB_DIR, strelka_unfiltered_vcfs, opts.base_output, "strelka.somatic.unfiltered")

    #-------------------------------------------------------------------------------------------
    #--------------------------------GERMLINE VARIANT CALLING-----------------------------------
    #-------------------------------------------------------------------------------------------

    #---PERFORM GERMLINE VARIANT CALLING IF GERMLINE CALLS DO NOT COME FROM IONREPORTER (BETA)---#
   
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



    #-------------------------------------------------------------------------------------------
    #----------------------VARIANT (VCF) FILTERING AND ANNOTATION-------------------------------
    #-------------------------------------------------------------------------------------------
    
    #---BREAK MULTI-ALLELIC SITES INTO ONE ALT ALLELE PER LINE---#
    
    if opts.ionreporter_only is True:
        vcfs_and_basenames = {"ionreporter.somatic.unfiltered" : ionreporter_somatic_unfiltered_vcf} 

    else:
        if opts.ionreporter_skip_all is True:
            vcfs_and_basenames = {"mutect2.somatic.unfiltered" : mutect2_unfiltered_vcf,
                                  "strelka.somatic.unfiltered" : strelka_unfiltered_vcf}
        else:
            vcfs_and_basenames = {"ionreporter.somatic.unfiltered" : ionreporter_somatic_unfiltered_vcf,
                                  "mutect2.somatic.unfiltered" : mutect2_unfiltered_vcf,
                                  "strelka.somatic.unfiltered" : strelka_unfiltered_vcf}
    
    if opts.ionreporter_germline_url_bool is True:
        vcfs_and_basenames.update({"ionreporter.germline.unfiltered" : ionreporter_germline_unfiltered_vcf})
    
    # If user chooses not to use IR, remove it from the list of VCFs to process
    if opts.ionreporter_skip_all is True:
        for k in vcfs_and_basenames.keys():
            if re.search("ionreporter", k):
                vcfs_and_basenames.pop(k)
    
    for description,vcf in vcfs_and_basenames.iteritems():
        if description == "ionreporter.somatic.unfiltered":
            ionreporter_somatic_unfiltered_vcf = multibreak_vcf(VCFLIB_DIR, ionreporter_somatic_unfiltered_vcf, opts.base_output, description)
        elif description == "ionreporter.germline.unfiltered":
            ionreporter_germline_unfiltered_vcf = multibreak_vcf(VCFLIB_DIR, ionreporter_germline_unfiltered_vcf, opts.base_output, description)
        elif description == "mutect2.somatic.unfiltered":
            try:
                mutect2_unfiltered_vcf = multibreak_vcf(VCFLIB_DIR, mutect2_unfiltered_vcf, opts.base_output, description)
            except:
                pass
        elif description == "strelka.somatic.unfiltered":
            strelka_unfiltered_vcf = multibreak_vcf(VCFLIB_DIR, strelka_unfiltered_vcf, opts.base_output, description)
            

    ### PERFORM HARD-FILTERING AND VEP ANNOTATION
    program_filter_vcf_list = []
    if opts.ionreporter_skip_all is False:
        program_filter_vcf_list = [
                                   # first index = output basename suffix
                                   # second index = SnpSift filter expression
                                   # third index = vcf input
    
                                   ["ionreporter.somatic",
                                    """(HRUN[*] <= 6)
                                    & ((FDP[*] >= 20) | (DP[*] >= 20))
                                    & ((FAO[*] >= 2) | (AO[*] >= 2)) 
                                    & ((GEN[1].FDP[*] >= 5) | (GEN[1].DP[*] >= 5)) 
                                    & !(ALT='<CNV>')""", 
                                    ionreporter_somatic_unfiltered_vcf],
    
                           
                                    ["ionreporter.cnv",
                                    "(ALT='<CNV>')",
                                    ionreporter_somatic_unfiltered_vcf]
                                
                                    ]
    
    if opts.ionreporter_only is False:
        if re.search("Population", opts.normal, re.IGNORECASE):
            program_filter_vcf_list = program_filter_vcf_list + [
                                       # first index = output basename suffix
                                       # second index = SnpSift filter expression
                                       # third index = vcf input
    
        #                                 ["varscan",  
        #                                  """
        #                                     ((SS ='1') | (SS = '2') | (SS= '3') )  
        #                                      & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))""",
        #                                      #& !(exists INDEL)""", # filter out INDELs
        #                                  varscan_vcf],
        #                                    
    
                               
        #                                 ["mutect.somatic",
        #                                  """(GEN[1].FA <= 0.05) 
        #                                      & (GEN[0].FA >= 0.05) 
        #                                      & (GEN[1].DP >= 5) 
        #                                      & (GEN[0].DP >= 20)""",
        #                                      #& (FILTER = 'PASS')""",
        #                                  mutect_vcf],
                                       
                                        ["mutect2.somatic",
                                         """((GEN[0].AF >= 0.05) 
                                             & (GEN[0].AD[1] >= 2))
                                             | (FILTER = 'PASS')""",
                                        mutect2_unfiltered_vcf],
         
                                        ["strelka.somatic",
                                         """(GEN[0].DP[*] >= 5)
                                             & (GEN[1].DP[*] >= 20)
                                             & (SGT !~ 'ref->ref')""",
                                         strelka_unfiltered_vcf]
                                       
                                       ]
        else:
            program_filter_vcf_list = program_filter_vcf_list + [
                                           # first index = output basename suffix
                                           # second index = SnpSift filter expression
                                           # third index = vcf input
        
            #                                 ["varscan",  
            #                                  """
            #                                     ((SS ='1') | (SS = '2') | (SS= '3') )  
            #                                      & ( (GEN[0].DP[*] >= 5) & (GEN[1].DP[*] >= 20))""",
            #                                      #& !(exists INDEL)""", # filter out INDELs
            #                                  varscan_vcf],
            #                                    
        
                                   
            #                                 ["mutect.somatic",
            #                                  """(GEN[1].FA <= 0.05) 
            #                                      & (GEN[0].FA >= 0.05) 
            #                                      & (GEN[1].DP >= 5) 
            #                                      & (GEN[0].DP >= 20)""",
            #                                      #& (FILTER = 'PASS')""",
            #                                  mutect_vcf],
                                           
                                            ["mutect2.somatic",
                                             """((GEN[1].AF <= 0.05) 
                                                 & (GEN[0].AF >= 0.05) 
                                                 & (GEN[1].AD[0] >= 5) 
                                                 & (GEN[0].AD[1] >= 2))
                                                 | (FILTER = 'PASS')""",
                                            mutect2_unfiltered_vcf],
             
                                            ["strelka.somatic",
                                             """(GEN[0].DP[*] >= 5)
                                                 & (GEN[1].DP[*] >= 20)
                                                 & (SGT !~ 'ref->ref')""",
                                             strelka_unfiltered_vcf]
                                           
                                           ]

    if opts.ionreporter_germline_url_bool is True:
        program_filter_vcf_list.append(
                                       ["ionreporter.loh",
                                         """(HRUN[*] <= 6)
                                             & ((FDP[*] >= 20) | (DP[*] >= 20))
                                             & ((FAO[*] >= 2) | (AO[*] >= 2)) 
                                             & ((GEN[1].FDP[*] >= 5) | (GEN[1].DP[*] >= 5))
                                             & ( ((GEN[0].AF[*] >= 0.65) & (GEN[1].AF[*] <= 0.65)) | (isHom(GEN[0]) & isHet(GEN[1])) )
                                             & !(ALT='<CNV>')""", 
                                         ionreporter_germline_unfiltered_vcf]
                                       )
        
        program_filter_vcf_list.append(
                                        ["ionreporter.germline",
                                        """(HRUN[*] <= 6)
                                        & ((GEN[0].AF[*] >= 0.30) & (GEN[1].AF[*] >= 0.30))
                                        & ((FDP[*] >= 20) | (DP[*] >= 20))
                                        & ((FAO[*] >= 2) | (AO[*] >= 2)) 
                                        & ((GEN[1].FDP[*] >= 5) | (GEN[1].DP[*] >= 5)) 
                                        & !(ALT='<CNV>')""", 
                                        ionreporter_germline_unfiltered_vcf]
                                       )
        



    # Delete the VEP index (i.e. force it to create a new index), so that if multiple VEP versions exist, there are no issues with the HGVS nomenclature module
    delete_VEP_index()
    # Loop over program,filter,vcf list and apply filters and annotations
    for group in program_filter_vcf_list:
        program,filter,vcf = (i for i in group)
        if os.stat(vcf).st_size == 0:
            print "WARNING: %s is empty.  Passing..." % program
        else:
            if determine_num_variants_in_vcf(vcf) < 1:
                # If no variants in VCF, create blank $program.json and pass
                print "WARNING: No variants detected in %s. \nWARNING: Final output will not have calls made by %s.  Please check error logs." % (vcf, program)        
                subprocess.call("touch %s.%s.json" % (opts.base_output, program),shell=True)
                pass
            else:
                # Apply filters via SnpSift
                filtered_vcf = SnpSift_filter(vcf,SNPSIFT_EXE,BEDTOOLS_EXE, filter, opts.base_output, program, opts.do_not_separate_LOH)
                
                # ionreporter.cnv.vcf needs to be handled via unfiltered method or will cause memory issues
                if program == "ionreporter.cnv":  
                    pass        # for now, pass this because long CNVs add lengthy processing time by VEP
                    #VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program)
                else:
                    if opts.filter_disabled is False:  
                        print "WARNING: Running consequential filters.  Only nonsynonymous coding variants will be reported!"  
                        VEP_command_filtered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program, filtered_vcf)
                    else:
                        VEP_command_unfiltered(VEP_EXE,VEP_REF_FASTA,opts.base_output,program, filtered_vcf)


    #-------------------------------------------------------------------------------------------
    #---------------------------------FINAL POST-PROCESSING-------------------------------------
    #-------------------------------------------------------------------------------------------
    
    
    
    extra_file_cleanup()
    
    move_files_to_new_subdirectory(opts.tumor, opts.normal, opts.base_output, opts.galaxy_flag)

    if opts.galaxy_flag is True:
        Galaxy_special_function(opts.base_output, opts.output_directory)
        zip_files(opts.base_output, opts.output_directory)
        html_out = open("%s" % opts.galaxy_html_file,"w")
        for file in sorted(os.listdir("%s" % opts.output_directory)):
            html_out.write("<a href='%s'>%s</a><br>" % (file,file))            
    else:
        opts.output_directory = opts.base_output
        zip_files(opts.base_output,opts.output_directory)

if __name__ == "__main__":
    main()