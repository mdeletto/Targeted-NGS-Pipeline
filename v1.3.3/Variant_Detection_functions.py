#/usr/bin/python

import subprocess
import sys
import re
import os
import json
import datetime
import openpyxl
import errno
from pprint import PrettyPrinter
import pprint
from collections import defaultdict
import vcf
import shutil
import openpyxl
import string
from openpyxl.reader.excel import load_workbook
from subprocess import CalledProcessError
import cmd
import MySQLdb


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def samtools_mpileup(Samtools,tumor_bam,normal_bam,base_output,reference_fasta,regions_file):
    print "Processing samtools mpileup..."
    try:
        print "%s mpileup -Q 7 -h 50 -o 10 -e 17 -m 4 -f %s -l %s %s %s 2>> /tmp/error > %s.mpileup" % (Samtools,reference_fasta,regions_file,normal_bam,tumor_bam,base_output)
        subprocess.call("%s mpileup -Q 7 -h 50 -o 10 -e 17 -m 4 -f %s -l %s %s %s 2>> /tmp/error > %s.mpileup" % (Samtools,reference_fasta,regions_file,normal_bam,tumor_bam,base_output),shell=True) 
    except:
        print "ERROR: Samtools mpileup failed.  Aborting..."
        sys.exit()
    
def varscan_command(VarScan,mpileup,base_output):
    """Perform VarScan variant calling and apply various filtering/grouping - RETIRED"""
    print "Processing VarScan variant calling..."
    try:
        print "java -jar VarScan.jar somatic %s.mpileup %s --mpileup 1 --somatic-p-value 0.05 --min-coverage-tumor 20 --min-coverage-normal 5 --strand-filter 1 --output-vcf 1" % (mpileup,base_output)
        subprocess.call("java -jar %s somatic %s.mpileup %s --mpileup 1 --somatic-p-value 0.05 --min-coverage-tumor 20 --min-coverage-normal 5 --strand-filter 1 --output-vcf 1 2>> /tmp/error " % (VarScan,mpileup,base_output),shell=True)
        subprocess.call("java -jar %s processSomatic %s.indel.vcf 2>> /tmp/error" % (VarScan,base_output),shell=True)
        subprocess.call("java -jar %s processSomatic %s.snp.vcf 2>> /tmp/error" % (VarScan,base_output),shell=True)
        subprocess.call("vcf-concat %s.indel.vcf %s.snp.vcf 2>> /tmp/error > %s.snp.indel.vcf" % (base_output,base_output,base_output),shell=True)
        subprocess.call("vcf-sort -c %s.snp.indel.vcf 2>> /tmp/error > %s.snp.indel.sorted.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-concat %s*.hc.vcf 2>> /tmp/error > %s.snp.indel.hc.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-sort -c %s.snp.indel.hc.vcf 2>> /tmp/error > %s.snp.indel.sorted.hc.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-concat %s*Somatic.vcf 2>> /tmp/error > %s.snp.indel.somatic.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-sort -c %s.snp.indel.somatic.vcf 2>> /tmp/error > %s.snp.indel.sorted.somatic.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-concat %s*Germline.vcf 2>> /tmp/error > %s.snp.indel.germline.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-sort -c %s.snp.indel.germline.vcf 2>> /tmp/error > %s.snp.indel.sorted.germline.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-concat %s*LOH.vcf 2>> /tmp/error > %s.snp.indel.LOH.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-sort -c %s.snp.indel.LOH.vcf 2>> /tmp/error > %s.snp.indel.sorted.LOH.vcf" % (base_output,base_output),shell=True)
        subprocess.call("vcf-concat %s.snp.Somatic.vcf %s.indel.Somatic.vcf %s.snp.LOH.vcf %s.indel.LOH.vcf 2>> /tmp/error > %s.snp.indel.somatic_LOH.vcf" % (base_output,base_output,base_output,base_output,base_output),shell=True)
        subprocess.call("vcf-sort -c %s.snp.indel.somatic_LOH.vcf 2>> /tmp/error > %s.snp.indel.sorted.somatic_LOH.vcf" % (base_output,base_output),shell=True)
    except:
        print "ERROR: VarScan somatic variant calling failed.  Aborting..."
        sys.exit()
 
def SnpSift_filter(vcf_in, SnpSift, BEDTOOLS_EXE, regex_filter, base_output, program, do_not_separate_LOH_from_germline):
    """Apply filter via SnpSift and subtract LOH calls from the GERMLINE VCF --- no point in having duplicates in there...
    ...unless...well, you want them in there...
    """
    print "Filtering VCF input for %s..." % program
    try:
        if program == "ionreporter.germline" and do_not_separate_LOH_from_germline is False:
            subprocess.call("""cat %s | \
                               java -jar %s varType - | \
                               java -jar %s filter "%s" 2>> %s.snpsift.log | \
                               %s subtract -a stdin -b %s.ionreporter.loh.vcf >> %s.%s.vcf""" % (vcf_in, SnpSift,SnpSift,regex_filter,base_output,BEDTOOLS_EXE,base_output,base_output,program),shell=True)
        else:
            subprocess.call("""cat %s | \
                                java -jar %s varType - | \
                               java -jar %s filter "%s" 2>> %s.snpsift.log >> %s.%s.vcf""" % (vcf_in,SnpSift, SnpSift,regex_filter, base_output, base_output,program),shell=True)
        return "%s.%s.vcf" % (base_output, program)
    except:
        print "ERROR: Filtering of VCF input with SnpSift failed.  Aborting..."


def VEP_command_unfiltered(VEP,REF_FASTA,base_output,program, vcf_in, DBNSFP):
    print "Annotating file..."
    try:
        vep_json_command = """perl %s --quiet \
                                    --assembly GRCh37 \
                                    --cache \
                                    --merged \
                                    --offline \
                                    --fasta %s \
                                    -i %s \
                                    --everything \
                                    --check_existing \
                                    --cache_version 90 \
                                    --json \
                                    -o %s.%s.json \
                                    --plugin LoFtool \
                                    --plugin TSSDistance \
                                    --plugin SpliceRegion \
                                    --plugin dbNSFP,%s,#chr,pos(1-coor),ref,alt,Ensembl_transcriptid,hg38_chr,hg38_pos,FATHMM_score,FATHMM_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred \
                                    -fork 16 
                                    --buffer_size 500 &> %s.vep.log""" % (VEP, REF_FASTA, vcf_in, base_output, program, DBNSFP, base_output)
        vep_json_command_process = subprocess.Popen(vep_json_command.split(" "))
        output, err = vep_json_command_process.communicate()
    except:
        print "ERROR: Could not initiate annotation on VCF file"


def VEP_command_filtered(VEP, REF_FASTA, base_output, program, vcf_in, DBNSFP, opts):
    try:
        def VEP_json_consequence_filtering(VEP_json_command_output,base_output,program):
            
            output_file = open("%s/%s.%s.json" % (os.getcwd(),base_output,program),"w")
            
            excluded_SO_terms_count = defaultdict(int)

            excluded_SO_terms = ['synonymous_variant',
                                 'start_retained_variant',
                                 'stop_retained_variant',
                                 'non_coding_transcript_exon_variant',
                                 'intron_variant',
                                 'downstream_gene_variant',
                                 'upstream_gene_variant',
                                 'UTR_variant',
                                 '3_prime_UTR_variant',
                                 '5_prime_UTR_variant',
                                 'NMD_transcript_variant',
                                 'null'
                                 '?']            
            
            for line in VEP_json_command_output.strip().split("\n"):
                loaded = json.loads(line)
                if loaded['most_severe_consequence'] in excluded_SO_terms:
                    excluded_SO_terms_count[loaded['most_severe_consequence']] += 1
                else:
                    consequences_passing_filters = []
                    try:
                        for transcript_consequence in loaded['transcript_consequences']:
                            #pprint.pprint(transcript_consequence)
                            
                            for consequence_term in transcript_consequence['consequence_terms']:
                                if consequence_term in excluded_SO_terms:
                                    
                                    try:
                                        if transcript_consequence['canonical'] == 1:
                                            consequences_passing_filters.append(transcript_consequence)
                                        else:
                                            excluded_SO_terms_count[consequence_term] += 1
                                            pass
                                    except KeyError, e:
                                        pass
                                    except Exception, e:
                                        print "ERROR: Unknown exception - %s" % str(e)
    
                                else:
                                    consequences_passing_filters.append(transcript_consequence)
                            
                
                        
                        loaded['transcript_consequences'] = consequences_passing_filters            
                        
                        output_file.write(json.dumps(loaded)+"\n")
                    
                    except Exception, e:
                        print "ERROR: Unknown exception - %s" % str(e)
                        print "ERROR: Printing line with error:"
                        print loaded
                    
            output_file.close()
            
            print "Number of variants excluded from %s due to consequential filtering..." % program
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(dict(excluded_SO_terms_count))

        with open("%s.vep.log" % base_output,"w") as err:
            vep_json_command = """perl %s --quiet \
                                        --assembly GRCh37 \
                                        --cache \
                                        --merged \
                                        --offline \
                                        --force_overwrite \
                                        --fasta %s \
                                        -i %s \
                                        --everything \
                                        --af \
                                        --check_existing \
                                        --cache_version 90 \
                                        --no_intergenic \
                                        --minimal \
                                        --allele_number \
                                        --json \
                                        --fork 16 \
                                        --buffer_size 500 \
                                        --output_file STDOUT \
                                        --plugin dbNSFP,%s,#chr,pos(1-coor),ref,alt,Ensembl_transcriptid,hg38_chr,hg38_pos,FATHMM_score,FATHMM_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred \
                                        --plugin LoFtool \
                                        --plugin SpliceRegion \
                                        --plugin TSSDistance """ % (VEP, REF_FASTA, vcf_in, DBNSFP)
                
# NEED BETTER IMPLEMENTATION - DOESN'T WORK FULLY                                     
#             if opts.filter_common_germline and re.search("germline", program):
#                 # include 1000 Genomes snps with frequency of 0.1% MAF
#                 vep_json_command += "--check_frequency --freq_pop 1KG_ALL --freq_freq 0.001 --freq_gt_lt gt --freq_filter exclude"
                
                                            
            
            vep_json_command_process = subprocess.Popen(vep_json_command.split(" "),
                                                        stdout=subprocess.PIPE,
                                                        stderr=err)

            output, err = vep_json_command_process.communicate()
        
        
        VEP_json_consequence_filtering(output, base_output, program)
        
    except Exception,e:
        print "ERROR: VEP could not annotate and filter variants"
        print "ERROR: %s" % str(e)


def move_files_to_new_subdirectory(tumor_bam, normal_bam, base_output,galaxy_flag):
    print "Moving files to new directory named %s to prepare for final spreadsheet processing..." % base_output
    try:
        #subprocess.call('/usr/lib/jvm/java-8-oracle/jre/bin/java -jar %s -f %s.varscan.json %s.ionreporter.no_cnv.json %s.ionreporter.cnv.vcf %s.ionreporter.tsv' % (EDDY,base_output,base_output,base_output,base_output),shell=True)
        if os.getcwd().split("/")[-1] == base_output:

            if galaxy_flag is False:
                mkdir_p("%s/tmp" % os.getcwd())
                subprocess.call('mv %s*tmp* tmp/ 2>> /tmp/error' % base_output, shell=True)
                subprocess.call('mv %s*temp* tmp/ 2>> /tmp/error' % base_output, shell=True)           
    
                mkdir_p("%s/BAMs" % os.getcwd())
                if re.search("Population", normal_bam, re.IGNORECASE):
                    subprocess.call('mv %s* BAMs/ 2>> /tmp/error' % tumor_bam, shell=True)
                else:
                    subprocess.call('mv %s* %s* BAMs/ 2>> /tmp/error' % (tumor_bam, normal_bam),shell=True)
            
        else:
            
            mkdir_p("%s/%s" %(os.getcwd(),base_output))
            
            if galaxy_flag is False:
                mkdir_p("%s/%s/tmp" %(os.getcwd(),base_output))
                subprocess.call('mv %s*tmp* %s/tmp/ 2>> /tmp/error' % (base_output,base_output),shell=True)
                subprocess.call('mv %s*temp* %s/tmp/ 2>> /tmp/error' % (base_output,base_output),shell=True)           
    
                mkdir_p("%s/%s/BAMs" %(os.getcwd(),base_output))
                if re.search("Population", normal_bam, re.IGNORECASE):
                    subprocess.call('mv %s* %s/BAMs/ 2>> /tmp/error' % (tumor_bam, base_output),shell=True)
                else:
                    subprocess.call('mv %s* %s* %s/BAMs/ 2>> /tmp/error' % (tumor_bam, normal_bam, base_output),shell=True)
            subprocess.call('mv %s* %s/ 2>> /tmp/error' % (base_output,base_output),shell=True)
    except:
        print "ERROR: Failed to move files"
    try:
        subprocess.call('mv %s.ionreporter.fusions.vcf %s/%s.ionreporter.fusions.vcf 2>> /tmp/error' % (base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Could not move fusion file--it may not be defined."
 
def Galaxy_special_function(base_output,galaxy_output_directory):
    print "Processing moving files to Galaxy directory..."
    try:
        subprocess.call('mkdir %s' % galaxy_output_directory,shell=True)
        subprocess.call('cp -r %s/* %s' % (base_output,galaxy_output_directory),shell=True)
    except:
        print "ERROR: Could not move to Galaxy output directory..."

def edit_IR_tsv_file(ionreporter_version,ionreporter_tsv,base_output):
    if ionreporter_version==4.0:
        subprocess.call("cp %s %s.ionreporter.cnv.tsv" % (ionreporter_tsv,base_output),shell=True)
    elif ionreporter_version==4.4:
        tsv_output = open("%s.ionreporter.cnv.tsv" % base_output, "w")
        with open("%s" % ionreporter_tsv,"r") as f:
            for line in f.readlines():
                line = re.sub("# locus","#chr\tpos",line)
                match = re.search("chr(.{1,2}):(\d*)",line)
                if match is not None:
                    line = re.sub("chr(.{1,2}):(\d*)","%s\t%s" % (match.group(1),match.group(2)),line)
                tsv_output.write(line)
    elif ionreporter_version==5.6:
        tsv_output = open("%s.tmp.ionreporter.cnv.tsv" % base_output, "w")
        with open("%s" % ionreporter_tsv,"r") as f:
            for line in f.readlines():
                line = re.sub("# locus","#chr\tpos",line)
                match = re.search("chr(.{1,2}):(\d*)",line)
                if match is not None:
                    line = re.sub("chr(.{1,2}):(\d*)","%s\t%s" % (match.group(1),match.group(2)),line)
                tsv_output.write(line)
        
        subprocess.call("cut -s --complement -f7,16,17,18 %s.tmp.ionreporter.cnv.tsv > %s.ionreporter.cnv.tsv" % (base_output, base_output), shell=True)
        
        os.remove("%s.tmp.ionreporter.cnv.tsv" % base_output)
                       
    else:
        tsv_output = open("%s.ionreporter.tsv" % base_output, "w")
        with open("%s" % ionreporter_tsv,"r") as f:
            for line in f.readlines():
                line = re.sub("# locus","#chr\tpos",line)
                match = re.search("chr(.{1,2}):(\d*)",line)
                if match is not None:
                    line = re.sub("chr(.{1,2}):(\d*)","%s\t%s" % (match.group(1),match.group(2)),line)
                tsv_output.write(line)

def IR4_4_VCF_fix(vcf,base_output):
    vcf_fix = open("%s.ionreporter.somatic.unfiltered.vcf" % base_output,"w") 
    with open("%s" % vcf,"r") as f:
        for line in f.readlines():
            if re.match("^#",line):
                if re.search('Description="+.*"+',line):
                    vcf_fix.write(line)
                elif re.search('Description=.+>',line):
                    match = re.search('Description=(.+)>',line)
                    subgroup = match.group(1)
                    line = re.sub(subgroup,'%s%s%s' % ('"',subgroup,'"'),line)
                    vcf_fix.write(line)
                else:
                    vcf_fix.write(line)
            else:
                vcf_fix.write(line)
    vcf_fix.close()
    return "%s.ionreporter.somatic.unfiltered.vcf" % base_output

def pull_BAM_from_url(url,base_output,sample_type):
    import requests
    resp = requests.get(url,auth=('downstream','downstream'))
    output = open("%s.%s.bam" % (base_output,sample_type),"w")
    output.write(resp.content)
    return "%s.%s.bam" % (base_output,sample_type)

def zip_files(base_output,galaxy_dir):
    try:
        os.chdir(galaxy_dir)
        subprocess.call('zip -qr %s.zip cnv/ %s.ionreporter.cnv.vcf %s.ionreporter.somatic.json %s.ionreporter.loh.json %s.ionreporter.germline.json %s.ionreporter.tsv %s.mutect2.somatic.json %s.strelka.somatic.json %s.ionreporter.fusions.vcf %s.specimen.json' % (base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Could not create zip file for Downstream upload. Check to make sure all output files exist."


def IR_locate_variant_zip(master_IR_dict, analysis_type):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "https://10.80.157.179/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (master_IR_dict[analysis_type]['url_encoded_analysis_name'],
                                                                                                                                                                                                                                                                                                                                                 master_IR_dict[analysis_type]['analysis_id'])] ,shell=True, stdout=subprocess.PIPE)
        output, err = proc.communicate()
    except Exception, e:
        print "Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
        print str(e)
        sys.exit(1)
    try:
        try:
            output = output.strip("[").strip("]")
            data = json.loads(output)
        except:
            print "Could not load json string"
        
        master_IR_dict[analysis_type]['unfiltered_variants_link'] = data['data_links']['unfiltered_variants']
        master_IR_dict[analysis_type]['filtered_variants_link'] = data['data_links']['filtered_variants']
        
        if analysis_type == "somatic":
            master_IR_dict[analysis_type]['control'] = data['samples']['NORMAL']
            master_IR_dict[analysis_type]['sample'] = data['samples']['TUMOR']
        elif analysis_type == "germline":
            master_IR_dict[analysis_type]['control'] = data['samples']['CONTROL']
            master_IR_dict[analysis_type]['sample'] = data['samples']['SAMPLE']
        elif analysis_type == "fusion":
            master_IR_dict[analysis_type]['control'] = None
            master_IR_dict[analysis_type]['sample'] = data['samples']['PROBAND']
        
        return master_IR_dict

    except Exception, e:
        print "Unable to process IonReporter links.  Aborting..."
        print str(e)
        sys.exit(1)

# def IR_locate_germline_variant_zip(basename,ionreporter_id):
#     try:
#         proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "https://10.80.157.179/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (basename,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
#         output, err = proc.communicate()
#         #print output
#         #print err
#     except Exception, e:
#         print "Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
#         print str(e)
#         sys.exit(1)
#     try:
#         output = output.strip("[").strip("]")
#         try:
#             data = json.loads(output)
#         except Exception, e:
#             print str(e)
#             print "Could not load json string"
#         unfiltered_variants = data['data_links']['unfiltered_variants']
#         filtered_variants = data['data_links']['filtered_variants']
#         control_barcode = data['samples']['CONTROL']
#         sample_barcode = data['samples']['SAMPLE']
#         
#         ret_vals = [unfiltered_variants,
#                     sample_barcode,
#                     control_barcode]
#         
#         return ret_vals
#         
#     except Exception, e:
#         print str(e)
#         print "Unable to process IonReporter links.  Aborting..."
#         sys.exit(1)
# 
# def IR_download_somatic_variant_zip(basename, variant_link, analysis_type, somatic_sample_dict):
# #     try:
# #         proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR.zip; \
# #                                     unzip IR.zip && unzip %s.zip; \
# #                                     cp ./Variants/*/*.vcf %s.ionreporter.%s_temp.vcf && cp ./Variants/*/*.tsv %s.ionreporter.%s_temp.tsv; \
# #                                     rm -rf %s.zip QC Variants Workflow_Settings""" % (variant_link,basename,basename,analysis_type,basename,analysis_type,basename)],shell=True,stdout=subprocess.PIPE)
# # 
# #         files = ["%s.ionreporter.%s_temp.vcf" % (basename, analysis_type),
# #                  "%s.ionreporter.%s_temp.tsv" % (basename, analysis_type)]
# #         
# #         return files
# #         
# #     except:
# #         print "Unable to download and/or unzip IonReporter files.  Aborting..."
# #         sys.exit(1)
# 
#     try:
#         proc = subprocess.Popen(["""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_somatic.zip; \
#                                   unzip -q IR_somatic.zip; \
#                                   rm -rf IR_somatic.zip; \
#                                   unzip -q %s*.zip; \
#                                   cp ./Variants/*/%s*.vcf %s.ionreporter.%s_temp.vcf && cp ./Variants/*/%s*.tsv %s.ionreporter.%s_temp.tsv; \
#                                   rm -rf IR_somatic.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (variant_link,
#                                                                                                              basename,
#                                                                                                              somatic_sample_dict['tumor'], 
#                                                                                                              basename,
#                                                                                                              analysis_type,
#                                                                                                              somatic_sample_dict['tumor'], 
#                                                                                                              basename,analysis_type,basename)], shell=True, stdout=subprocess.PIPE)
#         output, err = proc.communicate()
#         files = ["%s.ionreporter.%s_temp.vcf" % (basename, analysis_type),
#                  "%s.ionreporter.%s_temp.tsv" % (basename, analysis_type)]
#         
#         return files
#     
#     except:
#         print "Unable to download and/or unzip IonReporter files.  Aborting..."
#         sys.exit(1)
# 
# def IR_download_fusion_zip(variant_link,basename):
#     try:
#         proc = subprocess.Popen(["""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_fusion.zip; \
#                                     unzip -q IR_fusion.zip; \
#                                     rm -rf IR_fusion.zip; \
#                                     unzip -q %s*.zip; \
#                                     cp ./Variants/*/*.vcf %s.ionreporter.fusions.vcf; \
#                                     rm -rf IR_fusion.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (variant_link,basename,basename,basename)],shell=True,stdout=subprocess.PIPE)
#         proc.communicate()
#     except Exception, e:
#         print str(e)
#         print "Unable to download and/or unzip IonReporter Fusion files.  Aborting..."
#         sys.exit(1)


def IR_download_variant_zip(VCFLIB_DIR, opts, master_IR_dict, analysis_type):
    """Downloads a germline variant zip from IR.  
    
    Because IR calls germline variants separately on each BAM and reports them as such, use VCFLIB to merge the VCFs together.
    """
    try:
        if analysis_type == "germline":
        
            subprocess.call("""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_germline.zip; \
                               unzip -q IR_germline.zip; \
                               rm -rf IR_germline.zip; \
                               unzip -q %s*.zip; \
                               cp ./Variants/%s/%s*.vcf %s.ionreporter.%s_temp.sample.vcf; \
                               cp ./Variants/%s/%s*.vcf %s.ionreporter.%s_temp.control.vcf; \
                               rm -rf IR_germline.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (master_IR_dict[analysis_type]['unfiltered_variants_link'], 
                                                                                                           opts.base_output, 
                                                                                                           master_IR_dict[analysis_type]['sample'], 
                                                                                                           master_IR_dict[analysis_type]['sample'], 
                                                                                                           opts.base_output,
                                                                                                           analysis_type, 
                                                                                                           master_IR_dict[analysis_type]['control'], 
                                                                                                           master_IR_dict[analysis_type]['control'], 
                                                                                                           opts.base_output, 
                                                                                                           analysis_type, 
                                                                                                           opts.base_output), shell=True)
        
            vcf_list = ['%s.ionreporter.%s_temp.sample.vcf' % (opts.base_output, analysis_type),
                        '%s.ionreporter.%s_temp.control.vcf' % (opts.base_output, analysis_type)]
        
            master_IR_dict[analysis_type]['unfiltered_vcf'] = combine_vcf(VCFLIB_DIR, vcf_list, opts.base_output, "ionreporter.germline.unfiltered")
            
            for vcf in vcf_list:
                subprocess.call("rm %s" % vcf, shell=True)
    
        elif analysis_type == "fusion":
            subprocess.call("""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_fusion.zip; \
                                    unzip -q IR_fusion.zip; \
                                    rm -rf IR_fusion.zip; \
                                    unzip -q %s*.zip; \
                                    cp ./Variants/*/*.vcf %s.ionreporter.fusions.vcf; \
                                    rm -rf IR_fusion.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (master_IR_dict[analysis_type]['unfiltered_variants_link'],
                                                                                                              opts.base_output,
                                                                                                              opts.base_output,
                                                                                                              opts.base_output), shell=True)
            master_IR_dict[analysis_type]['unfiltered_vcf'] = "%s.ionreporter.fusions.vcf" % opts.base_output
        
        elif analysis_type == "somatic":
            
            proc = subprocess.Popen(["""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_somatic.zip; \
                                  unzip -q IR_somatic.zip; \
                                  rm -rf IR_somatic.zip; \
                                  unzip -q %s*.zip; \
                                  cp ./Variants/*/%s*.vcf %s.ionreporter.%s_temp.vcf && cp ./Variants/*/%s*.tsv %s.ionreporter.%s_temp.tsv; \
                                  rm -rf IR_somatic.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (master_IR_dict[analysis_type]['unfiltered_variants_link'],
                                                                                                             opts.base_output,
                                                                                                             master_IR_dict[analysis_type]['sample'], 
                                                                                                             opts.base_output,
                                                                                                             analysis_type,
                                                                                                             master_IR_dict[analysis_type]['sample'], 
                                                                                                             opts.base_output,
                                                                                                             analysis_type,
                                                                                                             opts.base_output)], shell=True, stdout=subprocess.PIPE)
            output, err = proc.communicate()
            
            master_IR_dict[analysis_type]['unfiltered_vcf'] = "%s.ionreporter.%s_temp.vcf" % (opts.base_output, analysis_type)
            master_IR_dict[analysis_type]['additional_tsv'] = "%s.ionreporter.%s_temp.tsv" % (opts.base_output, analysis_type)
        
        return master_IR_dict
    
    except Exception,e:
        print str(e)
        print "Unable to download and/or unzip IonReporter GERMLINE files.  Aborting..."
        sys.exit(1)

def rename_fusion_vcf(vcf,basename):
    try:
        print vcf
        subprocess.call("mv %s %s.ionreporter.fusions.vcf" % (vcf,basename),shell=True)
    except:
        print "Unable to rename IonReporter fusion.  Aborting..."
        sys.exit(1)
   
def bam_readcount_command(bam_readcount_exe,sample_type,base_output,reference_fasta,regions_file,bam_file):
    try:
        print "Producing bam readcounts for %s. This will be used for VarScan fpfilter command..." % sample_type
        subprocess.call("%s -f %s -l %s %s > %s.%s.readcounts.txt 2> /tmp/Bioinformatics_Pipeline_v1.3.bamreadcount.error.log" %(bam_readcount_exe,reference_fasta,regions_file,bam_file,base_output,sample_type),shell=True)
        readcount_filename = "%s.%s.readcounts.txt" %(base_output,sample_type)
        return readcount_filename
    except:
        print "Unable to produce bam readcounts for %s - %s." % (base_output,sample_type)
        sys.exit(1)
        
def select_target_regions(regions):
    try:
        print "Selecting target regions..."
        
        if regions=="CCP":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CCP/CCPHSMV2_052013.bed"
        elif regions=="OCPv3" or regions=="OCAv3":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/OCP/v3/OCAv3.20170110.designed.bed"
        elif regions=="OCP" or regions=="OCA" or regions=="OCPv2" or regions=="OCAv2":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/OCP/AmpliSeq_OCP/OCP.20150630.designed.bed"
        elif regions=="BRCA":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/BRCA1-2/BRCA1_2.20131001.designed.bed"
        elif regions=="CHPv2" or regions=="HSM":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CHPv2/CHP2.20131001.designed.bed"
        elif regions=="TSC":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/TSC1-TSC2/TSC1_2.designed.bed"
        elif regions=="TFNA":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/TFNA/Yale_Thyroid_DNA_WG_99191_167.1.20160607/WG_99191_167.1.20160607.TERTspike.designed.bed"
        elif regions=="AMGEN":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/AMGEN/AMGEN_WG00244_DNA.20170216.designed.bed"
        elif regions=="WhEx":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/AmpliSeqWhEx/AmpliSeqExome.20141113.designed.bed"
        else:
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CCPHSMV2_052013.bed"
            print "WARNING: No bed file was selected.  Defaulting to using CCP regions to capture as much data as possible."
        print "BED file selected as %s" %(REGIONS_FILE)
        return REGIONS_FILE
    except:
        print "ERROR: Failed to select target regions"
        sys.exit(1)

def write_logfile(opts):
    try:
        logfile = open("%s.log.txt" % opts.base_output,"w")
        logfile.write("Parameters for analysis:\n\n")
        logfile.write("Date of analysis: %s\n" % datetime.datetime.now())
        for key,value in opts.__dict__.iteritems():
            logfile.write("%s\t%s\n" % (key,value))
        logfile.close()
    except:
        print "WARNING: Failed to open log file"
        
def delete_VEP_index():
    try:
        subprocess.call("rm -rf /home/michael/YNHH/Reference_Files/FASTA/hg19.VEP.fasta.index",shell=True)
    except:
        print "ERROR: Could not delete VEP index"

def select_varscan_vcf_subset(opts,base_output):
    try:
        print "Selecting VarScan vcf subset..."
        print "%s is selected for VCF output" % opts.status
        if opts.status=='All':
            vcf = '%s.snp.indel.sorted.vcf' % base_output
        elif opts.status=='All-HC':
            vcf = '%s.snp.indel.sorted.hc.vcf' % base_output
        elif opts.status=='Somatic':
            vcf = '%s.snp.indel.sorted.somatic.vcf' % base_output
        elif opts.status=='LOH':
            vcf = '%s.snp.indel.sorted.LOH.vcf' % base_output
        elif opts.status=='Germline':
            vcf = '%s.snp.indel.sorted.germline.vcf' % base_output
        elif opts.status=='Somatic_LOH':
            vcf = '%s.snp.indel.sorted.somatic_LOH.vcf' % base_output
        else:
            vcf = '%s.snp.indel.sorted.somatic_LOH.vcf' % base_output
        return vcf
    except:
        print "ERROR: Could not select VarScan vcf subset!!!"

def multibreak_vcf(VCFLIB_DIR,vcf,base_output,description):
    try:
        subprocess.call("%s/vcfbreakmulti %s > %s.%s.multibreak.vcf" % (VCFLIB_DIR,vcf,base_output,description),shell=True)
        return "%s.%s.multibreak.vcf" % (base_output, description)
    except:
        print "ERROR: multibreak vcf failed"

def combine_vcf(VCFLIB_DIR,vcf_list,base_output,description):
    try:
        vcfs = " ".join(vcf_list)
        #print '%s/vcfcombine %s > %s.%s.vcf' % (VCFLIB_DIR,vcfs,base_output,description)
        subprocess.call('%s/vcfcombine %s > %s.%s.vcf' % (VCFLIB_DIR,vcfs,base_output,description),shell=True)
        return "%s.%s.vcf" % (base_output, description)
    except:
        print "ERROR: Could not join VCFs"
 
def muTect_caller_command(MUTECT_EXE,REGIONS_FILE,MUTECT_V1_PON,dbsnp_vcf,cosmic_vcf,REF_FASTA,normal_bam,tumor_bam,base_output):
    try:
        subprocess.call("java -jar %s --analysis_type MuTect --reference_sequence %s -L %s -normal_panel %s -cosmic %s -dbsnp %s --input_file:normal %s --input_file:tumor %s -vcf %s.mutect.somatic.unfiltered.vcf -o %s.mutect.stats.tsv --enable_extended_output --enable_qscore_output --max_alt_allele_in_normal_fraction 0.1 --max_alt_alleles_in_normal_count 10 --gap_events_threshold 30 --pir_median_threshold 4 -dcov 1000 2>> /tmp/error" % (MUTECT_EXE,REF_FASTA,REGIONS_FILE,MUTECT_V1_PON,cosmic_vcf,dbsnp_vcf,normal_bam,tumor_bam,base_output,base_output),shell=True)
        return "%s.mutect.somatic.unfiltered.vcf" % base_output
    except:
        print "ERROR: MuTect failed"

def muTect2_caller_command(GATK_LATEST_EXE,REGIONS_FILE,MUTECT2_V1_PON,dbsnp_vcf,cosmic_vcf,REF_FASTA,normal_bam,tumor_bam,base_output):
    try:
        print """#--------------------MUTECT2 SOMATIC CALLING--------------------#"""
        if re.search("Population", normal_bam, re.IGNORECASE):
            print "RUNNING MUTECT2 IN TUMOR-ONLY MODE"
            #mutect2_command = "java -jar %s --analysis_type MuTect2 --reference_sequence %s -L %s --normal_panel %s --cosmic %s --dbsnp %s --input_file:tumor %s -o %s.mutect2.somatic.unfiltered.vcf --max_alt_allele_in_normal_fraction 0.1 -nct 8 --minPruning 10 --kmerSize 60" % (GATK_LATEST_EXE,REF_FASTA,REGIONS_FILE,MUTECT2_V1_PON,cosmic_vcf,dbsnp_vcf,tumor_bam,base_output)
            mutect2_command = "java -jar %s --analysis_type MuTect2 --reference_sequence %s -L %s --normal_panel %s --cosmic %s --dbsnp %s --input_file:tumor %s -o %s.mutect2.somatic.unfiltered.vcf -nct 8 --minPruning 10 --kmerSize 60" % (GATK_LATEST_EXE,REF_FASTA,REGIONS_FILE,MUTECT2_V1_PON,cosmic_vcf,dbsnp_vcf,tumor_bam,base_output)

        else:
            mutect2_command = "java -jar %s --analysis_type MuTect2 --reference_sequence %s -L %s --normal_panel %s --cosmic %s --dbsnp %s --input_file:normal %s --input_file:tumor %s -o %s.mutect2.somatic.unfiltered.vcf --max_alt_allele_in_normal_fraction 0.1 -nct 8 --minPruning 10 --kmerSize 60" % (GATK_LATEST_EXE,REF_FASTA,REGIONS_FILE,MUTECT2_V1_PON,cosmic_vcf,dbsnp_vcf,normal_bam,tumor_bam,base_output)

        print "MUTECT2 SOMATIC CALLING COMMAND:", mutect2_command
        log = open('%s.mutect2.log' % base_output, 'w')
        subprocess.call(mutect2_command,shell=True, stdout=log, stderr=log)
        log.close()

        return "%s.mutect2.somatic.unfiltered.vcf" % base_output
    except Exception, e:
        print "ERROR: MuTect2 failed"
        print "ERROR" + str(e)

def GATK_snp_caller_command(GATK_EXE,tumor_bam,normal_bam,base_output,REF_FASTA,dbsnp_vcf):
    try:
        subprocess.call("java -jar %s --analysis_type UnifiedGenotyper -I %s -I %s -o %s.gatk.snps.germline.raw.vcf -R %s -D %s -dcov 1000 -nt 16 2>> /tmp/error" % (GATK_EXE,tumor_bam,normal_bam,base_output,REF_FASTA,dbsnp_vcf),shell=True)
        return "%s.gatk.snps.germline.raw.vcf" % base_output
    except:
        print "ERROR: GATK Germline SNP calling failed"

def GATK_indel_caller_command(GATK_EXE,tumor_bam,normal_bam,base_output,REF_FASTA,dbsnp_vcf):
    try:
        subprocess.call("java -jar %s --analysis_type UnifiedGenotyper -glm INDEL -I %s -I %s -o %s.gatk.indels.germline.raw.vcf -R %s -D %s -dcov 1000 -nt 16 2>> /tmp/error" % (GATK_EXE,tumor_bam,normal_bam,base_output,REF_FASTA,dbsnp_vcf),shell=True)
        return "%s.gatk.indels.germline.raw.vcf" % base_output
    except:
        print "ERROR: GATK Germline INDEL calling failed"

def Strelka_somatic_variant_calling_command(STRELKA_EXE,normal_bam,tumor_bam,REF_FASTA,STRELKA_CONFIG,base_output):
    try:
        strelka_analysis_dir = '%s/strelka_analysis' % os.getcwd()
        strelka_command = "perl %s --normal=%s --tumor=%s --ref=%s --config=%s --output-dir=%s; make -j 16 -C %s &> %s.strelka.log" % (STRELKA_EXE,os.path.abspath(normal_bam),os.path.abspath(tumor_bam),REF_FASTA,STRELKA_CONFIG,strelka_analysis_dir,strelka_analysis_dir, base_output)
        print """#--------------------STRELKA SOMATIC CALLING--------------------#"""
        print "STRELKA SOMATIC CALLING COMMAND:", strelka_command
        
        log = open('%s.strelka.log' % base_output, 'w')
           
        subprocess.call("perl %s --normal=%s --tumor=%s --ref=%s --config=%s --output-dir=%s" % (STRELKA_EXE,os.path.abspath(normal_bam),os.path.abspath(tumor_bam),REF_FASTA,STRELKA_CONFIG,strelka_analysis_dir),shell=True, stdout=log, stderr=log)
        subprocess.call("make -j 16 -C %s" % (strelka_analysis_dir), shell=True, stdout=log, stderr=log) 
        log.close()
        
        vcf_list = ['%s/results/passed.somatic.indels.vcf' % strelka_analysis_dir,
                    '%s/results/passed.somatic.snvs.vcf' % strelka_analysis_dir]
        print vcf_list
        
        return vcf_list
    except:
        print "ERROR: Strelka somatic variant calling failed"

def extra_file_cleanup():
    subprocess.call("rm -rf Variants Workflow_Settings QC STDOUT_summary.html VER*.log strelka_analysis",shell=True)

def determine_num_variants_in_vcf(vcf):
    """This is a quick and dirty way to count variants in a VCF (does not account for multi-allelic sites)"""
    try:
        num_variants = subprocess.check_output("grep -vc '^#' %s" % vcf,shell=True)
        num_variants = num_variants.strip("\n")
    except CalledProcessError as e:
        if int(e.output.strip("\n")) > 1:
            raise
        else:
            num_variants = 0

    return num_variants

def sample_attribute_autodetection(basename, pipeline_version, panel):
    print """#--------------------SAMPLE ATTRIBUTE AUTODETECTION--------------------#"""
    
    def redefine_Downstream_panel_name(panel):
        """Redefines panel name to be recognized by Downstream."""
        if panel=='HSM' or panel=='CHPv2':
            panel = '50 gene panel'
        elif panel=='CCP' or panel=='409':
            panel = '409 gene panel'
        elif re.search("OCP", panel) or re.search("OCA", panel):
            panel = 'Oncomine panel'
        else:
            panel = panel
        
        return panel
    
    def define_panel_version(panel):
        """Defines panel version based on abbrevation."""
        
        panel_version = re.search("v(\d)", panel)
        try:
            panel_version = str(round(float(panel_version.group(1)), 2))
        except:
            # Panel version abbreviation is not included, check to see if panel name should always be v1.0
            
            # These panels do not have anything past v1.0 yet
            v1_panels = re.compile("HSM|CCP|TFNA")
            if v1_panels.search(panel):
                panel_version = "1.0"
            else:
                if re.search("OCP", panel) or re.search("OCA", panel):
                    # All recent analyses with OCP have been 2.0
                    # panel version autodetection was not supported until v1.3.3
                    panel_version = "2.0"
                else:
                    # For unlisted panels
                    panel_version = "1.0"

        return panel_version
    
    def reformat_copath_id_format(test_string):
        """Reformats coPath ID if in X15-03344 format with 0 before numbers.
        Sometimes, the ID is entered into the spreadsheet without a preceding 0.
        """
        format = long
        if re.search("-\d+", test_string):
            match = re.search("(\w+?)-(\d+)", test_string)
            id_part_one = match.group(1)
            id_part_two = match.group(2)
            if len(id_part_two) == 5:
                print "WARNING: BASENAME IDENTIFIED AS COPATH_ID (LONG) - REFORMATTING TO COPATH_ID (SHORT)"
                reformatted_id = id_part_one + "-" + id_part_two.lstrip("0")
                reformatted_id_bool = [reformatted_id, True]
            elif len(id_part_two) > 5:
                print "ERROR: COPATH_ID not in recognize format"
                reformatted_id = None
                reformatted_id_bool = [reformatted_id, False]
            else:
                print "WARNING: BASENAME IDENTIFIED AS COPATH_ID (SHORT) - REFORMATTING TO COPATH_ID (LONG)"
                leading_zero_count = 5 - len(id_part_two)
                reformatted_id = id_part_one + "-" + ("0"*leading_zero_count) + id_part_two
                reformatted_id_bool = [reformatted_id, True]
        else:
            print "ERROR: Basename is not in COPATH_ID format."
            reformatted_id = None
            reformatted_id_bool = [reformatted_id, False]
            
        return reformatted_id_bool

    
    def pull_basename_from_TPTracker(basename, pipeline_version, panel_version):
    
        db = MySQLdb.connect("localhost","root","*23Ft198","tumor_profiling_lab")
        cursor = db.cursor()
    
        sql = "SELECT copath_id, sequencing_panel, tumor_type, tumor_source, germline_sample, germline_source, microdissection_type, percent_malignant_cells FROM sequencing_cases \
               WHERE copath_id = '%s' \
               ORDER BY req_date DESC" % basename
    
        try:
           # Execute the SQL command
            cursor.execute(sql)
           # Fetch all the rows in a list of lists.
            results = cursor.fetchall()
            if len(results) == 0:
                reformatted_id, reformatted_bool = reformat_copath_id_format(basename)
                if reformatted_bool is True:
                    print "WARNING: No matches on %s...trying reformatted ID as %s" % (basename, reformatted_id)
                    sql = "SELECT copath_id, sequencing_panel, tumor_type, tumor_source, germline_sample, germline_source, microdissection_type, percent_malignant_cells FROM sequencing_cases \
                           WHERE copath_id = '%s' \
                           ORDER BY req_date DESC" % reformatted_id
                           
                    cursor.execute(sql)
                    results = cursor.fetchall()
                    if len(results) > 0:
                        result_set = results[0]
                else:
                    print "WARNING: No matches on %s.  Basename is not in correct CoPath ID format.  Stop contacting TPTracker." % basename
            else:
                result_set = results[0]
        except Exception, e:
            print str(e)
            print "Error: unable to fetch data"
    
        # disconnect from server
        db.close()
        
        specimen_info = defaultdict(lambda: "")
        specimen_info['panel'] = panel
        specimen_info['panel_version'] = panel_version
        specimen_info['pipeline_version'] = pipeline_version
        try:
            result_set
        except NameError:
    
            specimen_info['normal'] = ""
            specimen_info['tumor'] = ""
            specimen_info['malignant_cells'] = ""
            specimen_info['dissection'] = ""  
        else:
            mysql_fields = ['copath_id',
                            'sequencing_panel',
                            'tumor_type',
                            'tumor_source',
                            'germline_sample',
                            'germline_source',
                            'microdissection_type',
                            'percent_malignant_cells']
            mysql_fields = dict(zip(mysql_fields,result_set))
            
            # Handle MySQL NULL/None values
            for k, v in mysql_fields.iteritems():
                if v is None:
                    mysql_fields[k] = ""
            
            specimen_info['tumor'] = ", ".join([mysql_fields['tumor_type'], mysql_fields['tumor_source']])
            specimen_info['normal'] = mysql_fields['germline_sample']
            specimen_info['malignant_cells'] = mysql_fields['percent_malignant_cells']
            specimen_info['dissection'] = mysql_fields['microdissection_type']
            
        if re.search("N/A", str(specimen_info['normal'])) or specimen_info['normal'] is None:
            specimen_info['normal'] = "None"
        
        with open('%s.specimen.json' % basename, 'w') as outfile:
            json.dump(dict(specimen_info), outfile)
        print json.dumps(dict(specimen_info))
    
    panel_version = define_panel_version(panel)
    panel = redefine_Downstream_panel_name(panel)
    pull_basename_from_TPTracker(basename, pipeline_version, panel_version)
    


def extract_fusion_VCF_information(vcf):
    """Extracts information from the ionreporter.fusions.vcf file."""
    
    fusion_dict = defaultdict(dict)
    fusion_dict['sum_control_count'] = 0
    
    with open(vcf) as f:
        for line in f.readlines():
            line = line.strip()
            
            if re.search("##mapd", line):
                match = re.search("##mapd=(.*)", line)
                mapd = match.group(1)
                fusion_dict['mapd'] = mapd
            elif re.search("SVTYPE=GeneExpression", line):
                split_line = line.split("\t")
                gene_match = re.search("GENE_NAME=(.+?);", line)
                gene = gene_match.group(1)
                read_count_match = re.search("READ_COUNT=(.+?);", line)
                read_count = read_count_match.group(1)
                fusion_dict['gene_expression_read_counts'][gene] = read_count
            elif re.search("SVTYPE=ExprControl", line):
                gene_match = re.search("GENE_NAME=(.+?);", line)
                gene = gene_match.group(1)
                read_count_match = re.search("READ_COUNT=(.+?);", line)
                read_count = read_count_match.group(1)
                fusion_dict['control_per_gene_read_counts'][gene] = read_count
                fusion_dict['sum_control_count'] += int(read_count)
            elif re.search("##TotalMappedFusionPanelReads", line):
                match = re.search("##TotalMappedFusionPanelReads=(.*)", line)
                total_mapped_fusion_reads = match.group(1)
                fusion_dict['total_mapped_fusion_reads'] = total_mapped_fusion_reads

    return fusion_dict

class IR_CNV_module:

    def __init__(self, remote_dir_base_path):
        self.tiles_path = remote_dir_base_path + "/outputs/CnvActor-00/diffCoverage.seg"

    def execute_CNV_module(self, REGIONS_FILE, RSCRIPT_EXE, CNV_PLOT_EXE, SNPSIFT_EXE, opts, cnv_vcf):

        def extract_CNV_segments_from_VCF(self, SnpSift, cnv_vcf, base_output):
            """Extracts CNV segment information from the IR CNV .vcf.  This is a more detailed approach than pulling the CN_Segments.seg file from IR."""
            
            try:
                if os.path.isfile("./%s" % cnv_vcf):
                    subprocess.call("java -jar %s extractFields %s CHROM POS END LEN GEN[0].CN NUMTILES CONFIDENCE PRECISION > %s.cnv.seg.detailed.tsv" % (SnpSift, cnv_vcf, base_output), shell=True)
                    cnv_seg_detailed = "%s.cnv.seg.detailed.tsv" % base_output
                    return cnv_seg_detailed
                else:
                    print "ERROR: Could not find %s." % cnv_vcf
                    print "ERROR: Not extracting fields for CNVs."
                    return None
            except Exception, e:
                print "ERROR: Unexpected error during CNV field extraction."
                print "ERROR: %s" % str(e)
                return None
        
        def IR_download_diffCoverage(self, variant_link, basename, tmp_output_name):
            """Downloads the diffCoverage.seg segment files from the IR server."""
            
            try:
                proc = subprocess.Popen(["""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o %s; \
                                            unzip -q %s; \
                                            rm -rf %s; \
                                            mv diffCoverage.seg %s.diffCoverage.seg; \
                                            rm -rf VER*.log""" % (variant_link,tmp_output_name, tmp_output_name, tmp_output_name, basename)],shell=True,stdout=subprocess.PIPE)
                output, err = proc.communicate()
                tiles_file = "%s.diffCoverage.seg" % basename
                
                return tiles_file
                
            except Exception, e:
                print str(e)
                print "Unable to download and/or unzip IonReporter CNV files.  Aborting..."
                
                return None
        
        def prepare_cnv_directory(self):
            # Create 'cnv' directory        
            mkdir_p("%s/cnv" % os.getcwd())
            
            # Change to cnv directory
            os.chdir('cnv')


        def plot_CNVs(self):
            log = open('./%s.cnv.log' % opts.base_output, 'a+')
            cmd = "%s %s %s %s %s %s" % (RSCRIPT_EXE, CNV_PLOT_EXE, opts.base_output, tiles_file, segment_file, REGIONS_FILE)
            p = subprocess.Popen(cmd.split(" "), stdout = log, stderr = log)
            p.communicate()
        
        
        
        print """#--------------------COPY NUMBER VARIANT DETECTION--------------------#"""
        
        prepare_cnv_directory(self)
        tiles_file = IR_download_diffCoverage(self, self.tiles_path, opts.base_output, "diffCoverage.seg.tmp.zip")
        segment_file = extract_CNV_segments_from_VCF(self, SNPSIFT_EXE, "../"+cnv_vcf, opts.base_output)
        # If there are no tiles or segment files, don't plot.
        if (tiles_file is not None) and (segment_file is not None):
            # plot CNVs
            plot_CNVs(self)
        os.chdir("..")
            
class QC_module:
    
    def execute_FASTQC(self, opts):
        """Execute FASTQC for BAMs in case folder."""
        
        print """#--------------------FASTQC QUALITY CONTROL--------------------#"""
        
        # Gather BAM list from current working directory
        bam_list = []
        for bam in os.listdir(os.getcwd()):
            if bam.endswith(".bam"):
                bam_list.append(bam)
        
        # make FASTQC directory if it doesn't exist
        try:
            os.makedirs("./FASTQC")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        
        
        log = open('./FASTQC/%s.fastqc.log' % opts.base_output, 'a+')
        
        # Run FASTQC command on all bams in directory using 8 threads    
        cmd = "fastqc -f bam -t 8 -o FASTQC %s" % " ".join(bam_list)
        p = subprocess.Popen(cmd.split(" "), stdout = log, stderr = log)
        p.communicate()
        
        # Close log
        log.close()
        
        # Delete .zip file created by FASTQC
        # We only need the uncompressed output directory for review
        fastqc_dir = os.path.join(os.getcwd(), "FASTQC")
        for f in os.listdir(fastqc_dir):
            if f.endswith(".zip"):
                os.remove(os.path.join(fastqc_dir, f))

def run_pipeline_parser(base_output, PP_PARSER_EXE):
    """Run pipeline-parser script to format raw pipeline files into spreadsheet."""
    try:
        subprocess.call("python %s -a --output-basename %s.unfiltered" % (PP_PARSER_EXE, base_output), shell=True)
    except Exception, e:
        print str(e)
        print "ERROR: pipeline-parser.py failed"


def filter_variant_spreadsheet(base_output, FILTER_VARIANTS_EXE):
    """Run filter_variants.py script to filter the final spreadsheet."""

    if os.path.isfile("%s.unfiltered.xlsx" % base_output):
        subprocess.call("python %s -i %s.unfiltered.xlsx -o %s.xlsx --minimum-vaf=0.05" % (FILTER_VARIANTS_EXE, base_output, base_output), shell=True)
    else:
        print "WARNING: Could not find unfiltered variant spreadsheet.  Not performing final filtering of variants"

def construct_master_IR_dict(opts):
    """Create master IR dict that holds IR analysis names and IDs."""
    
    master_IR_dict = defaultdict(lambda: None)
    
    bools_to_check = {'germline': opts.ionreporter_germline_url_bool,
                      'somatic' : opts.ionreporter_somatic_url_bool,
                      'fusion' : opts.ionreporter_fusion_url_bool
                      }
    
    for analysis_type in bools_to_check.keys():
        master_IR_dict[analysis_type] = defaultdict(lambda: None)
        if bools_to_check[analysis_type] is True:
            if analysis_type == "germline":
                analysis_id = opts.ionreporter_germline_id
                analysis_name = opts.ionreporter_germline_analysis_name
            elif analysis_type == "somatic":
                analysis_id = opts.ionreporter_somatic_id
                analysis_name = opts.ionreporter_somatic_analysis_name
            elif analysis_type == "fusion":
                analysis_id = opts.ionreporter_fusion_id
                analysis_name = opts.ionreporter_fusion_analysis_name      
          
            master_IR_dict[analysis_type]['analysis_id'] = analysis_id
            master_IR_dict[analysis_type]['analysis_name'] = analysis_name

    return master_IR_dict

def zip_files(base_output,galaxy_dir):
    try:
        os.chdir(galaxy_dir)
        subprocess.call('zip -qr %s.zip cnv/ %s.ionreporter.cnv.vcf %s.ionreporter.somatic.json %s.ionreporter.loh.json %s.ionreporter.germline.json %s.ionreporter.tsv %s.mutect2.somatic.json %s.strelka.somatic.json %s.ionreporter.fusions.vcf %s.specimen.json' % (base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Could not create zip file for Downstream upload. Check to make sure all output files exist."