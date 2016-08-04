#/usr/bin/python

import subprocess,sys,re,os,json
import vcf

def samtools_command(Samtools,tumor_bam,normal_bam,base_output,reference_fasta,regions_file):
    print "Processing samtools mpileup..."
    try:
        print "%s mpileup -f %s -l %s %s %s 2>> /tmp/error > %s.mpileup" % (Samtools,reference_fasta,regions_file,normal_bam,tumor_bam,base_output)
        subprocess.call("%s mpileup -Q 7 -h 50 -o 10 -e 17 -m 4 -f %s -l %s %s %s 2>> /tmp/error > %s.mpileup" % (Samtools,reference_fasta,regions_file,normal_bam,tumor_bam,base_output),shell=True) 
    except:
        print "ERROR: Samtools mpileup failed.  Aborting..."
        sys.exit()
    
def varscan_command(VarScan,mpileup,base_output,tumor_purity,normal_purity):
    print "Processing VarScan variant calling..."
    try:
        print "java -jar VarScan.jar somatic %s.mpileup %s --mpileup 1 --tumor-purity %s --normal-purity %s --min-coverage 2 --strand-filter 1 --output-vcf 1" % (base_output,base_output,tumor_purity,normal_purity)
        subprocess.call("java -jar %s somatic %s.mpileup %s --mpileup 1 --tumor-purity %s --somatic-p-value 0.05 --normal-purity %s --min-coverage 2 --strand-filter 1 --output-vcf 1 2>> /tmp/error " % (VarScan,base_output,base_output,tumor_purity,normal_purity),shell=True)
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
 
def SnpSift_filter(vcf_in,SnpSift,regex_filter,base_output,program):
    print "Filtering VCF input for %s..." % program
    try:
        subprocess.call('cat %s | java -jar %s filter "%s" 2>> /tmp/error > %s.%s.vcf' % (vcf_in,SnpSift,regex_filter,base_output,program),shell=True)
    except:
        print "ERROR: Filtering of VCF input with SnpSift failed.  Aborting..."

def VEP_command(VEP,REF_FASTA,base_output,program):
    print "Annotating file..."
    print 'perl %s --quiet --cache --merged --offline --fasta %s -i %s.%s.vcf --everything --cache_version 81 --no_intergenic --filter no_stop_retained_variant --filter no_non_coding_transcript_exon_variant --filter no_intron_variant --filter no_synonymous_variant --json -o %s.%s.json -fork 16 2>> /tmp/error' % (VEP,REF_FASTA,base_output,program,base_output,program)
    try:
        #sys.stdout.write('perl %s --cache --merged --offline --fasta %s -i %s --everything --no_intergenic --filter no_stop_retained_variant --filter no_synonymous_variant --json -o %s.%s.json --fork 16' % (VEP,REF_FASTA,base_output,program,base_output,program))
        #subprocess.call('perl %s --quiet --cache --merged --offline --fasta %s -i %s.%s.vcf --everything --no_intergenic --filter no_stop_retained_variant --filter no_synonymous_variant --json -o %s.%s.json -fork 16 2>> /tmp/error' % (VEP,REF_FASTA,base_output,program,base_output,program),shell=True)
        subprocess.call('perl %s --quiet --cache --merged --offline --fasta %s -i %s.%s.vcf --everything --cache_version 81 --no_intergenic --json -o %s.%s.json --fork 16 2>> /tmp/error' % (VEP,REF_FASTA,base_output,program,base_output,program),shell=True)

    except:
        print "ERROR: Could not initiate annotation on VCF file"

def VEP_command_unfiltered(VEP,REF_FASTA,unfiltered_vcf,base_output,program):
    print "Annotating file..."
    #try:
    sys.stdout.write('perl %s --cache --merged --offline --fasta %s -i %s --everything --cache_version 81 --no_intergenic --filter no_stop_retained_variant --filter no_synonymous_variant --json -o %s.%s.json --fork 16' % (VEP,REF_FASTA,base_output,program,base_output,program))
    #subprocess.call('perl %s --quiet --cache --merged --offline --fasta %s -i %s.%s.vcf --everything --no_intergenic --filter no_intron_variant --filter no_non_coding_transcript_exon_variant --filter no_stop_retained_variant --filter no_synonymous_variant --json -o %s.%s.json --fork 16' % (VEP,REF_FASTA,base_output,program,base_output,program),shell=False)
    subprocess.call('perl %s --quiet --cache --merged --offline --fasta %s -i %s.%s.vcf --everything --cache_version 81 --no_intergenic --json -o %s.%s.json --fork 16' % (VEP,REF_FASTA,base_output,program,base_output,program),shell=False)
    #except:
     #   print "ERROR: Could not initiate annotation on VCF file"
        
def Move_files_1(base_output):
    print "Moving files to new directory named %s to prepare for final spreadsheet processing..." % base_output
    try:
        subprocess.call('mkdir %s' % base_output,shell=True)
        subprocess.call('mv %s.ionreporter.no_cnv.json %s/%s.ionreporter.no_cnv.json' % (base_output,base_output,base_output),shell=True)
        subprocess.call('mv %s.ionreporter.cnv.vcf %s/%s.ionreporter.cnv.vcf' % (base_output,base_output,base_output),shell=True)
        subprocess.call('mv %s.varscan.json %s/%s.varscan.json' % (base_output,base_output,base_output),shell=True)
        subprocess.call('mv %s.ionreporter.tsv %s/%s.ionreporter.tsv' % (base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Failed to move files"
    try:
        subprocess.call('mv %s.ionreporter.fusions.vcf %s/%s.ionreporter.fusions.vcf 2>> /tmp/error' % (base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Could not move fusion file--it may not be defined."
        
def Eddy_command(EDDY,base_output):
    print "Processing final spreadsheet for %s..." % base_output
    try:
        subprocess.call('/usr/lib/jvm/java-8-oracle/jre/bin/java -jar %s -d %s' % (EDDY,base_output),shell=True)
    except:
        print "ERROR: Could not process final spreadsheet"

def Move_files_2(base_output):
    print "Moving all final files..."
    try:
        subprocess.call('mv %s* %s/ 2>> /tmp/error' % (base_output,base_output),shell=True)
    except:
        print "ERROR: Could not move final files"
        
def Galaxy_special_function(base_output,galaxy_output_directory):
    print "Processing moving files to Galaxy directory..."
    try:
        subprocess.call('mkdir %s' % galaxy_output_directory,shell=True)
        subprocess.call('cp -r %s/* %s' % (base_output,galaxy_output_directory),shell=True)
    except:
        print "ERROR: Could not move to Galaxy output directory..."

def edit_IR_tsv_file(ionreporter_version,ionreporter_tsv,base_output):
    #subprocess.call("cp %s %s.ionreporter.tsv" % (ionreporter_tsv,base_output),shell=True)

    if ionreporter_version==4.0:
        subprocess.call("cp %s %s.ionreporter.tsv" % (ionreporter_tsv,base_output),shell=True)
    elif ionreporter_version==4.4:
        tsv_output = open("%s.ionreporter.tsv" % base_output, "w")
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
    vcf_fix = open("%s.ionreporter.vcf" % base_output,"w") 
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
    return "%s.ionreporter.vcf" % base_output

def pull_BAM_from_url(url,base_output,sample_type):
    import requests
    resp = requests.get(url,auth=('downstream','downstream'))
    output = open("%s.%s.bam" % (base_output,sample_type),"w")
    output.write(resp.content)
    return "%s.%s.bam" % (base_output,sample_type)

def zip_files(base_output,galaxy_dir):
    try:
        os.chdir(galaxy_dir)
        subprocess.call('zip -q %s.zip %s.ionreporter.cnv.vcf %s.ionreporter.no_cnv.json %s.ionreporter.tsv %s.varscan.json %s.ionreporter.fusions.vcf' % (base_output,base_output,base_output,base_output,base_output,base_output),shell=True)
    except:
        print "ERROR: Could not create zip file for Downstream upload. Check to make sure all output files exist."
        
def IR_locate_variant_zip(basename,ionreporter_id):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "https://10.80.157.179/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (basename,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
        output, err = proc.communicate()
    except:
        print "Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
        sys.exit(1)
    try:
        output = output.strip("[").strip("]")
        print output
        try:
            data = json.loads(output)
        except:
            print "Could not load json string"
        unfiltered_variants = data['data_links']['unfiltered_variants']
        filtered_variants = data['data_links']['filtered_variants']
        print unfiltered_variants
        return unfiltered_variants
    except:
        print "Unable to process IonReporter links.  Aborting..."
        sys.exit(1)
        

def IR_download_variant_zip(basename,variant_link):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR.zip; unzip IR.zip && unzip %s.zip; cp ./Variants/*/*.vcf %s.ionreporter.temp.vcf && cp ./Variants/*/*.tsv %s.ionreporter.temp.tsv; rm -rf %s.zip QC Variants Workflow_Settings""" % (variant_link,basename,basename,basename,basename)],shell=True,stdout=subprocess.PIPE)
    except:
        print "Unable to download and/or unzip IonReporter files.  Aborting..."
        sys.exit(1)

def IR_download_fusion_zip(fusion_analysis_name,variant_link,basename):
    try:
        proc = subprocess.Popen(["""curl -k -H "Authorization:UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY" "%s" 2> /dev/null -o IR_fusion.zip; unzip IR_fusion.zip && unzip %s.zip; cp ./Variants/*/*.vcf %s.ionreporter.fusions.vcf; rm -rf %s.zip""" % (variant_link,fusion_analysis_name,basename,fusion_analysis_name)],shell=True,stdout=subprocess.PIPE)
    except:
        print "Unable to download and/or unzip IonReporter Fusion files.  Aborting..."
        sys.exit(1)

def rename_fusion_vcf(vcf,basename):
    try:
        print vcf
        subprocess.call("mv %s %s.ionreporter.fusions.vcf" % (vcf,basename),shell=True)
    except:
        print "Unable to rename IonReporter fusion.  Aborting..."
        sys.exit(1)

def delete_VEP_index():
    try:
        subprocess.call("rm -rf /home/michael/YNHH/Reference_Files/FASTA/hg19.VEP.fasta.index",shell=True)
    except:
        print "ERROR: Could not delete VEP index"

def bam_readcount_command(bam_readcount_exe,sample_type,base_output,reference_fasta,regions_file,bam_file):
    try:
        print "Producing bam readcounts for %s. This will be used for VarScan fpfilter command..." % sample_type
        subprocess.call("%s -f %s -l %s %s > %s.%s.readcounts.txt " %(bam_readcount_exe,reference_fasta,regions_file,bam_file,base_output,sample_type),shell=True)
    except:
        print "Unable to produce bam readcounts for %s."
        sys.exit(1)