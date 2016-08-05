#/usr/bin/python

import subprocess,sys,re,os,json,datetime,errno,shutil
import pprint
from collections import defaultdict
import vcf
import shutil

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
            subprocess.call("""awk 'BEGIN { FS = "\t" } ; { if ($10!="." && $11!=".") print}' %s | \
                               java -jar %s filter "%s" 2>> %s.snpsift.log | \
                               java -jar %s varType - | \
                               %s subtract -a stdin -b %s.ionreporter.loh.vcf > %s.%s.vcf""" % (vcf_in,SnpSift,regex_filter,base_output,SnpSift,BEDTOOLS_EXE,base_output,base_output,program),shell=True)
        else:
            subprocess.call("""awk 'BEGIN { FS = "\t" } ; { if ($10!="." && $11!=".") print}' %s | \
                               java -jar %s filter "%s" 2>> %s.snpsift.log | \
                               java -jar %s varType - > %s.%s.vcf""" % (vcf_in,SnpSift,regex_filter,base_output,SnpSift,base_output,program),shell=True)
        return "%s.%s.vcf" % (base_output, program)
    except:
        print "ERROR: Filtering of VCF input with SnpSift failed.  Aborting..."


def VEP_command_unfiltered(VEP,REF_FASTA,base_output,program, vcf_in):
    print "Annotating file..."
    try:
        subprocess.call('perl %s --quiet --cache --merged --offline --fasta %s -i %s --everything --check_alleles --coding_only --cache_version 83 --no_intergenic --json -o %s.%s.json -fork 16 &> %s.vep.log' % (VEP,REF_FASTA,vcf_in,base_output,program,base_output),shell=True)
    except:
        print "ERROR: Could not initiate annotation on VCF file"

def VEP_command_filtered(VEP,REF_FASTA,base_output,program,vcf_in):
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
                                 'null']            
            
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
            vep_json_command = 'perl %s --quiet --cache --merged --offline --force_overwrite --fasta %s -i %s --everything --maf_exac --check_alleles --cache_version 83 --no_intergenic --minimal --allele_number --json --fork 16 --output_file STDOUT' % (VEP, REF_FASTA, vcf_in)
            vep_json_command_process = subprocess.Popen(vep_json_command,
                                                        stdout=subprocess.PIPE,
                                                        stderr=err,
                                                        shell=True)

        output, err = vep_json_command_process.communicate()
        
        
        VEP_json_consequence_filtering(output, base_output, program)
        
    except Exception,e:
        print "ERROR: VEP could not annotate and filter variants"
        print "ERROR: %s" % str(e)


def move_files_to_new_subdirectory(EDDY,base_output,galaxy_flag):
    print "Moving files to new directory named %s to prepare for final spreadsheet processing..." % base_output
    try:
        subprocess.call('/usr/lib/jvm/java-8-oracle/jre/bin/java -jar %s -f %s.varscan.json %s.ionreporter.no_cnv.json %s.ionreporter.cnv.vcf %s.ionreporter.tsv' % (EDDY,base_output,base_output,base_output,base_output),shell=True)
        mkdir_p("%s/%s" %(os.getcwd(),base_output))
        
        if galaxy_flag is False:
            mkdir_p("%s/%s/tmp" %(os.getcwd(),base_output))
            subprocess.call('mv %s*tmp* %s/tmp/ 2>> /tmp/error' % (base_output,base_output),shell=True)
            subprocess.call('mv %s*temp* %s/tmp/ 2>> /tmp/error' % (base_output,base_output),shell=True)           

            mkdir_p("%s/%s/BAMs" %(os.getcwd(),base_output))
            subprocess.call('mv %s*bam* %s/BAMs/ 2>> /tmp/error' % (base_output,base_output),shell=True)

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
                    line = re.sub("chr(.{1,2}):(\d*)","chr%s\t%s" % (match.group(1),match.group(2)),line)
		tsv_output.write(line)
                    
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
        subprocess.call('zip -q %s.zip %s.ionreporter.cnv.vcf %s.ionreporter.somatic.json %s.ionreporter.loh.json %s.ionreporter.germline.json %s.ionreporter.tsv %s.mutect2.somatic.json %s.strelka.somatic.json %s.ionreporter.fusions.vcf' % (base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Could not create zip file for Downstream upload. Check to make sure all output files exist."


def IR_locate_variant_zip(basename,ionreporter_id):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "https://10.80.157.179/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (basename,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
        output, err = proc.communicate()
        #print output
        #print err
    except Exception, e:
        print "Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
        print str(e)
        sys.exit(1)
    try:
        output = output.strip("[").strip("]")
        #print output
        try:
            data = json.loads(output)
        except:
            print "Could not load json string"
        unfiltered_variants = data['data_links']['unfiltered_variants']
        filtered_variants = data['data_links']['filtered_variants']
        #print unfiltered_variants
        return unfiltered_variants
    except Exception, e:
        print "Unable to process IonReporter links.  Aborting..."
        print str(e)
        sys.exit(1)

def IR_locate_germline_variant_zip(basename,ionreporter_id):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "https://10.80.157.179/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (basename,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
        output, err = proc.communicate()
    except:
        print "Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
        sys.exit(1)
    try:
        output = output.strip("[").strip("]")

        try:
            data = json.loads(output)
        except:
            print "Could not load json string"
        unfiltered_variants = data['data_links']['unfiltered_variants']
        filtered_variants = data['data_links']['filtered_variants']
        control_barcode = data['samples']['CONTROL']
        sample_barcode = data['samples']['SAMPLE']
        
        ret_vals = [unfiltered_variants,
                    sample_barcode,
                    control_barcode]
        
        return ret_vals
        
    except Exception, e:
        print str(e)
        print "Unable to process IonReporter links.  Aborting..."
        sys.exit(1)

def IR_download_somatic_variant_zip(basename,variant_link,analysis_type):
#     try:
#         proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR.zip; \
#                                     unzip IR.zip && unzip %s.zip; \
#                                     cp ./Variants/*/*.vcf %s.ionreporter.%s_temp.vcf && cp ./Variants/*/*.tsv %s.ionreporter.%s_temp.tsv; \
#                                     rm -rf %s.zip QC Variants Workflow_Settings""" % (variant_link,basename,basename,analysis_type,basename,analysis_type,basename)],shell=True,stdout=subprocess.PIPE)
# 
#         files = ["%s.ionreporter.%s_temp.vcf" % (basename, analysis_type),
#                  "%s.ionreporter.%s_temp.tsv" % (basename, analysis_type)]
#         
#         return files
#         
#     except:
#         print "Unable to download and/or unzip IonReporter files.  Aborting..."
#         sys.exit(1)

    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_somatic.zip; \
                                  unzip -q IR_somatic.zip; \
                                  rm -rf IR_somatic.zip; \
                                  unzip -q %s*.zip; \
                                  cp ./Variants/*/*.vcf %s.ionreporter.%s_temp.vcf && cp ./Variants/*/*.tsv %s.ionreporter.%s_temp.tsv; \
                                  rm -rf IR_somatic.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (variant_link,basename,basename,analysis_type,basename,analysis_type,basename)], shell=True, stdout=subprocess.PIPE)

        files = ["%s.ionreporter.%s_temp.vcf" % (basename, analysis_type),
                 "%s.ionreporter.%s_temp.tsv" % (basename, analysis_type)]
        
        return files
    
    except:
        print "Unable to download and/or unzip IonReporter files.  Aborting..."
        sys.exit(1)

def IR_download_fusion_zip(variant_link,basename):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_fusion.zip; \
                                    unzip -q IR_fusion.zip; \
                                    rm -rf IR_fusion.zip; \
                                    unzip -q %s*.zip; \
                                    cp ./Variants/*/*.vcf %s.ionreporter.fusions.vcf; \
                                    rm -rf IR_fusion.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (variant_link,basename,basename,basename)],shell=True,stdout=subprocess.PIPE)
    except Exception, e:
        print str(e)
        print "Unable to download and/or unzip IonReporter Fusion files.  Aborting..."
        sys.exit(1)



def IR_download_germline_variant_zip(VCFLIB_DIR, basename, variant_link, analysis_type, germline_sample_dict):
    """Downloads a germline variant zip from IR.  
    
    Because IR calls germline variants separately on each BAM and reports them as such, use VCFLIB to merge the VCFs together.
    """
    try:
        subprocess.call("""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_germline.zip; \
                           unzip -q IR_germline.zip; \
                           rm -rf IR_germline.zip; \
                           unzip -q %s*.zip; \
                           cp ./Variants/%s/%s.vcf %s.ionreporter.%s_temp.sample.vcf; \
                           cp ./Variants/%s/%s.vcf %s.ionreporter.%s_temp.control.vcf; \
                           rm -rf IR_germline.zip %s*.zip QC Variants Workflow_Settings VER*.log""" % (variant_link, basename, germline_sample_dict['sample'], germline_sample_dict['sample'], basename, analysis_type, germline_sample_dict['control'], germline_sample_dict['control'], basename, analysis_type, basename), shell=True)
#         subprocess.call("unzip -qq IR.zip && unzip -qq %s.zip" % basename, shell=True)
#         subprocess.call("cp ./Variants/%s/%s.vcf %s.ionreporter.%s_temp.sample.vcf" % (germline_sample_dict['sample'],germline_sample_dict['sample'],basename,analysis_type), shell=True)
#         subprocess.call("cp ./Variants/%s/%s.vcf %s.ionreporter.%s_temp.control.vcf" % (germline_sample_dict['control'],germline_sample_dict['control'],basename,analysis_type), shell=True)
#         subprocess.call("rm -rf %s.zip IR.zip QC Variants Workflow_Settings" % basename, shell=True)
        
        vcf_list = ['%s.ionreporter.%s_temp.sample.vcf' % (basename, analysis_type),
                    '%s.ionreporter.%s_temp.control.vcf' % (basename, analysis_type)]
        
        germline_unfiltered_vcf = combine_vcf(VCFLIB_DIR, vcf_list, basename, "ionreporter.germline.unfiltered")
        
        for vcf in vcf_list:
            subprocess.call("rm %s" % vcf, shell=True)
            
        return germline_unfiltered_vcf
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
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CCPHSMV2_052013.bed"
        elif regions=="OCP":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/OCP3v1.20140718.designed.bed"
        elif regions=="BRCA":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/BRCA1-2/BRCA1_2.20131001.designed.bed"
        elif regions=="CHPv2" or regions=="HSM":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/CHPv2/CHP2.20131001.designed.bed"
        elif regions=="TSC":
            REGIONS_FILE = "/home/michael/YNHH/Reference_Files/TSC1-TSC2/TSC1_2.designed.bed"
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
        print '%s/vcfcombine %s > %s.%s.vcf' % (VCFLIB_DIR,vcfs,base_output,description)
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
        mutect2_command = "java -jar %s --analysis_type MuTect2 --reference_sequence %s -L %s --normal_panel %s --cosmic %s --dbsnp %s --input_file:normal %s --input_file:tumor %s -o %s.mutect2.somatic.unfiltered.vcf --max_alt_allele_in_normal_fraction 0.1 -nct 8 --minPruning 10 --kmerSize 60  2>> %s.mutect2.log" % (GATK_LATEST_EXE,REF_FASTA,REGIONS_FILE,MUTECT2_V1_PON,cosmic_vcf,dbsnp_vcf,normal_bam,tumor_bam,base_output, base_output)
        print "#" * len(mutect2_command)
        print "MUTECT2 SOMATIC CALLING COMMAND:"
        print mutect2_command
        print "#" * len(mutect2_command)   
        subprocess.call(mutect2_command,shell=True)
        return "%s.mutect2.somatic.unfiltered.vcf" % base_output
    except:
        print "ERROR: MuTect failed"

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
        print "#" * len(strelka_command)
        print "STRELKA SOMATIC CALLING COMMAND:"
        print strelka_command
        print "#" * len(strelka_command)   
        subprocess.call("perl %s --normal=%s --tumor=%s --ref=%s --config=%s --output-dir=%s 2>> %s.strelka.log" % (STRELKA_EXE,os.path.abspath(normal_bam),os.path.abspath(tumor_bam),REF_FASTA,STRELKA_CONFIG,strelka_analysis_dir,base_output),shell=True)
        subprocess.call("make -j 16 -C %s 2>> %s.strelka.log" % (strelka_analysis_dir, base_output), shell=True) 
        vcf_list = ['%s/results/passed.somatic.indels.vcf' % strelka_analysis_dir,
                    '%s/results/passed.somatic.snvs.vcf' % strelka_analysis_dir]
        print vcf_list
        return vcf_list
    except:
        print "ERROR: Strelka somatic variant calling failed"

def extra_file_cleanup():
    subprocess.call("rm -rf Variants Workflow_Settings QC STDOUT_summary.html VER*.log",shell=True)

def determine_num_variants_in_vcf(vcf):
    """This is a quick and dirty way to count variants in a VCF (does not account for multi-allelic sites)"""
    num_variants = subprocess.check_output("grep -vc '^#' %s" % vcf,shell=True)
    return num_variants