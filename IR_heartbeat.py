#!/usr/bin/python

import subprocess
import sys
import json
import optparse
import urllib
import re
import datetime
import pprint
import os
import pysam
import shutil
import glob
import errno
#import pymongo
import MySQLdb
from collections import defaultdict
from elasticsearch import Elasticsearch

os.chdir("/home/michael/Cases/Clinical/Pending")

global opts

desc="""IR_heartbeat.py acts as a gateway for initiating the bioinformatics analyses (variant calling, annotation, and other QC) for IonTorrent data.  The script contacts the IR server, pulls in necessary info, initiates the downstream bioinformatics analyses, and logs the activity in Elasticsearch."""
usage="./IR_heartbeat -d <analysis_date>"
version="1.0"
parser = optparse.OptionParser(description=desc,usage=usage,version=version)
 
group = optparse.OptionGroup(parser, "General options")
#group.add_option("--authorization","-a", help="Authorization API key for IonReporter",dest='auth_key',action='store')
group.add_option("--ionreporter_url","-i", help="IP address or hostname of IR server (ex: 10.80.157.179)",dest='url',action='store',default="10.80.157.179")
group.add_option("--analysis_date","-d", help="Pull IR analyses from this date.  YYYY-mm-dd format required (e.g. %default).  Default date for analysis is today!",dest='date',action='store',default=datetime.datetime.strftime(datetime.date.today(), "%Y-%m-%d"))
group.add_option("--download_bams_only","", help="Pull tumor and normal bams for a given analysis and quit",dest='bams_only',action='store_true',default=False)
group.add_option("--pipeline_version","-p", help="Pipeline version to use.  Default is to use the latest pipeline (%default)",dest='pipeline_version',action='store',default="1.2.1")#datetime.date.today())
group.add_option("-v", help="verbose <[%default]>", dest='verbose',action='store_true',default=False)
parser.add_option_group(group)

# group = optparse.OptionGroup(parser, "Single mode",
#                     "Downloads VCF for a single analysis")
#  
# group.add_option('--analysis_name', help="IR analysis name", dest='input_analysis_name', action='store')
# parser.add_option_group(group)

(opts, args) = parser.parse_args()

####################################
#---PRESET ENVIRONMENT VARIABLES---#
####################################

WORKING_DIRECTORY = "/home/michael/Cases/Clinical/Pending"
#WORKING_DIRECTORY = "/home/michael/Development/IR_heartbeat"

es = Elasticsearch()

# SSL is disabled until I get a shield license

# es = Elasticsearch(
#                    ['https://admin:tumorprofiling650@localhost:9200/'],
#                    verify_certs=True,
#                    # provide a path to CA certs on disk
#                    ca_certs='/home/michael/ca/certs/cacert.pem',
#                    # PEM formatted SSL client certificate
#                    client_cert='/etc/apache2/ssl/apache.crt',
#                    # PEM formatted SSL client key
#                    client_key='/etc/apache2/ssl/apache.key'
#                    )

###############################
#--------FUNCTIONS------------#
###############################

def set_IR_API_key_based_on_url(url):
    if url=="10.80.157.179":
        IR_API_KEY = "UmxyUXNPR3M1Q2RsbS9NYjBHQjBIaUxFTFA5RkJhRHBaMmlSSXZJTjBmUnNmQ0t1NkhOSUlrMStiNHFIQm16UjNKN2NYMzNOT2czcytqc2RveEhqK3BBSHhZNEhpNmRDVmtQaGRUZ1Z5ZXVXazJMTllQemIvV3A5c2NHOTNxRmY"
    else:
        sys.exit("ERROR: This IonReporter IP address is unsupported.")
    return IR_API_KEY

def set_pipeline_version(pipeline_version):
    if pipeline_version=="1.2.1":
        TARGETED_NGS_PIPELINE = "/home/michael/YNHH/Code/Github-mdeletto/Targeted-NGS-Pipeline/v%s/Variant_Detection_v%s.py" % (pipeline_version, pipeline_version)
    else:
        sys.exit("ERROR: Pipeline Version %s is unsupported" % pipeline_version)
    return TARGETED_NGS_PIPELINE


def command_line_parsing():
    pass

def format_date_to_string(date):
    if type(date) is str:
        date = datetime.datetime.strptime(date, "%Y-%m-%d")
    date = date.strftime("%B %d, %Y")
    return date

def IR_locate_variant_zip(analysis_name,ionreporter_id,IR_url,authorization_key):
    url_encoded_basename = urllib.quote_plus(analysis_name) # catch whitespace or other problematic characters and encode for url
    
    try:
        if ionreporter_id is None:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&name=%s" 2> /dev/null""" % (authorization_key,IR_url,url_encoded_basename)],shell=True,stdout=subprocess.PIPE)
        else:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (authorization_key,IR_url,url_encoded_basename,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
        output, err = proc.communicate()
    except:
        print "ERROR: Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
        sys.exit(1)
    try:
        try:
            data = json.loads(output)
            #print json.dumps(data,indent=4)
        except:
            print "ERROR: Could not load json string"
        unfiltered_variants = data[0]['data_links']['unfiltered_variants']
        filtered_variants = data[0]['data_links']['filtered_variants']
        return unfiltered_variants
    except Exception,e:
        #print str(e)
        print "ERROR: Unable to process IonReporter links.  Is this a valid analysis name?  Aborting..."
        return 0

def IR_download_variant_zip(basename,variant_link,authorization_key):
    if bool(variant_link) is False:
        print "ERROR: Unable to process case %s" % basename
        sys.exit(1)
    else:
        basename = basename.replace(" ","_") # IR converts any whitespace characters in <analysis_name> to underscores
        try:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "%s" 2> /dev/null -o IR.zip; unzip IR.zip && unzip %s.zip; cp ./Variants/*/*.vcf %s.vcf && cp ./Variants/*/*.tsv %s.tsv; rm -rf %s.zip QC Variants Workflow_Settings""" % (authorization_key,variant_link,basename,basename,basename,basename)],shell=True,stdout=subprocess.PIPE)
        except:
            print "Unable to download and/or unzip IonReporter files.  Aborting..."
            sys.exit(1)

def IR_analysis_summary_view():
    proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&view=summary" 2> /dev/null""" % (IR_API_KEY,"10.80.157.179")],shell=True,stdout=subprocess.PIPE)
    output, err = proc.communicate()
    all_IR_analyses = json.loads(output)
    
    return all_IR_analyses

def select_analyses(all_IR_analyses):
    names = []
    date = format_date_to_string(opts.date)
    for analysis in all_IR_analyses:
    
        if analysis['start_date'] == date: #and analysis['started_by'] == "Mike D'Eletto":
            if not re.search('-',analysis['name']) and not re.search('_',analysis['name']):
                print "ERROR: Analysis string is not in correct format.  Please consult README for using this workflow.  Trying to proceed..."
                names.append(analysis['name'].strip())
            elif not re.search('_',analysis['name']):
                name_split = analysis['name'].split("-")
                copath_id = str(name_split[0]+"-"+name_split[1])
                names.append(copath_id.strip())
            elif re.search('_',analysis['name']):
                name_split = analysis['name'].split("_")
                copath_id = name_split[0]
                names.append(copath_id.strip())
    
    
    names = list(set(names))
    
    return names

def execute_scp_bam_copy(path,bam_barcode_id):
    subprocess.call("scp -q ionadmin@10.80.157.179:%s ." % path, shell=True)
    basename_process = subprocess.check_output("""basename -s .bam %s""" % bam_barcode_id,shell=True)
    local_bam_filepath_name = basename_process.strip()  
    return local_bam_filepath_name  

def connect_to_database(database_name):

    try:
        cnx = MySQLdb.connect("localhost",
                              "root",
                              "*23Ft198",
                              database_name)
        print "SUCCESSFULLY CONNECTED TO %s" % database_name
    except MySQLdb.Error as err:
        if err.errno == errno.errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errno.errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
            
    return cnx


def download_bams_from_IR(IR_download_link,IR_download_dir,workflow_dict_key):
    # Downloads BAMs from IonReporter by downloading .rrs files, which contain samples definitions and filepaths on the IR server
    # Samtools is also called to merge BAMs (if more than one BAM exists for a sample) and index BAMs
    
    sample_name = workflow_dict_key
    bam_pipeline_flags = defaultdict(lambda: None) # to be passed back via return, to be used for pipeline automation
    bam_pipeline_flags['sample_name'] = sample_name
    
    tumor_bam_link = "%s%s%s.rrs" % (IR_download_link,IR_download_dir,workflow_dict[sample_name]['somatic_analysis']['somatic_tumor_name']) # download .rrs file that defines sample (if merged BAM)
    tumor_bam_sample_definition_file = tumor_bam_link.split("/")[-1]
    if re.search("Population",workflow_dict[sample_name]['somatic_analysis']['somatic_normal_name'],re.IGNORECASE):
        # If BAM is the Population Normal, don't waste time pulling such a large file from the IR server
        bams_to_download = {"tumor":[tumor_bam_sample_definition_file,tumor_bam_link]}
        bam_pipeline_flags['normal_bam_name'] = "Population Normal"
        bam_pipeline_flags['normal_bam_path'] = "Population Normal"
    else:
        normal_bam_link = "%s%s%s.rrs" % (IR_download_link,IR_download_dir,workflow_dict[sample_name]['somatic_analysis']['somatic_normal_name'])
        normal_bam_sample_definition_file = normal_bam_link.split("/")[-1]
        bams_to_download = {"tumor":[tumor_bam_sample_definition_file,tumor_bam_link],
                            "normal":[normal_bam_sample_definition_file,normal_bam_link]}

    # Download BAMs
    print "DOWNLOADING BAMS..."
    for sample_type,nested_list in bams_to_download.iteritems():
        sample_definition_file, bam_link = (i for i in nested_list)
        proc = subprocess.check_output("""curl -s -O -k -H "Authorization:%s" "%s" 2> /dev/null && unzip -o %s """ % (IR_API_KEY,bam_link,sample_definition_file),shell=True)
        direct_bam_filepaths = subprocess.check_output("""awk '{print $NF}' %s/%s""" % (os.getcwd(),sample_definition_file),shell=True)
        num_lines = sum(1 for line in open('%s/%s' %(os.getcwd(),sample_definition_file)))
        bam_barcode_paths_and_ids = {}
        direct_bam_filepaths = direct_bam_filepaths.strip().split("\n")
        for direct_bam_filepath in direct_bam_filepaths:
            bam_barcode_id = direct_bam_filepath.split("/")
            bam_barcode_id = bam_barcode_id[-1]
            bam_barcode_paths_and_ids[direct_bam_filepath] = bam_barcode_id
            #direct_bam_link = IR_download_link + direct_bam_filepath

        # If multiple BAMs exist in sample definition .rrs file
        if num_lines > 1 and len(bam_barcode_paths_and_ids.keys()) == 2:
            bam_pipeline_flags['merged_bams'] = True
            counter = 1
            for direct_bam_filepath, bam_barcode_id in bam_barcode_paths_and_ids.iteritems():
                local_bam_filepath_name = execute_scp_bam_copy(direct_bam_filepath, bam_barcode_id)
                shutil.move('%s/%s' % (os.getcwd(),bam_barcode_id), '%s/%s-%s-%s.bam' % (os.getcwd(),sample_name,sample_type,str(counter)))
                counter += 1
            subprocess.call("/home/michael/bin/samtools-1.1/samtools merge -f %s/%s-%s-merged.bam %s/%s-%s-%s.bam %s/%s-%s-%s.bam" % (os.getcwd(),sample_name,sample_type,os.getcwd(),sample_name,sample_type,str(1),os.getcwd(),sample_name,sample_type,str(2)),shell=True)
            subprocess.call("rm -rf %s/%s-%s-%s.bam %s/%s-%s-%s.bam" % (os.getcwd(),sample_name,sample_type,str(1),os.getcwd(),sample_name,sample_type,str(2)), shell=True)
            #pysam.index("%s/%s-%s-merged.bam" % (os.getcwd(),sample_name,sample_type))
            bam_pipeline_flags['%s_bam_name' %(sample_type)] = "%s/%s-%s-merged.bam" % (os.getcwd(),sample_name,sample_type)
        # If single BAM exists in sample definition .rrs file
        else:
            bam_pipeline_flags['merged_bams'] = False
            for direct_bam_filepath, bam_barcode_id in bam_barcode_paths_and_ids.iteritems():
                local_bam_filepath_name = execute_scp_bam_copy(direct_bam_filepath, bam_barcode_id)
                shutil.move('%s/%s' % (os.getcwd(),bam_barcode_id),'%s/%s-%s.bam' % (os.getcwd(),sample_name,sample_type))
            #pysam.index("%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type))
            bam_pipeline_flags['%s_bam_path' % (sample_type)] = "%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type)
            bam_pipeline_flags['%s_bam_name' % (sample_type)] = "%s-%s.bam" % (sample_name,sample_type)
            bam_pipeline_flags['%s_barcode' % sample_type] = local_bam_filepath_name

    return bam_pipeline_flags


def initiate_IR_download(workflow_dict,sample_name):
    if bool(workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_id']) is False:
        variant_link = IR_locate_variant_zip(workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_name'],
                                             None,
                                             opts.url,
                                             IR_API_KEY)
    else:
        variant_link = IR_locate_variant_zip(workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_name'],
                                             workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_id'],
                                             opts.url,
                                             IR_API_KEY)
    
    IR_download_URL = re.search("(http.+?download\?filepath=)(.+)",variant_link)
    if IR_download_URL:
        IR_download_link = IR_download_URL.group(1)
        IR_download_dir = IR_download_URL.group(2)
        IR_download_dir = IR_download_dir.split("/")
        IR_download_dir = "/".join(IR_download_dir[:-1])+"/"
    
    bam_dict_options = download_bams_from_IR(IR_download_link, IR_download_dir, sample_name)  # downloads BAMs and return a dictionary consisting of bam names and filepaths
    workflow_dict[sample_name].update(bam_dict_options)
    
def run_pipeline():
    
    sample_name = key
    if opts.pipeline_version=="1.2.1":
        if workflow_dict[key]['fusion_analysis']['fusion_analysis_name'] is None or re.search("None", workflow_dict[key]['fusion_analysis']['fusion_analysis_name']):
            fusion_parameters = ""
        else:
            fusion_parameters = """--ionreporter_fusion_url_bool \
                                   --ionreporter_fusion_analysis_name=%s \
                                   --ionreporter_fusion_id=%s""" % (workflow_dict[key]['fusion_analysis']['fusion_analysis_name'],
                                                                    workflow_dict[key]['fusion_analysis']['fusion_analysis_id'] )
        
        command = """%s \
            -s All-HC \
            -c %s \
            --regions=%s \
            -t %s \
            -n %s \
            -p IonTorrent \
            --ionreporter_version=4.4 \
            --ionreporter_somatic_url_bool \
            --ionreporter_somatic_analysis_name %s \
            --ionreporter_somatic_id %s \
            %s """ % (TARGETED_NGS_PIPELINE, sample_name, workflow_dict[key]['panel_name'], workflow_dict[key]['tumor_bam_name'],
                      workflow_dict[key]['normal_bam_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_id'], fusion_parameters)                               
        
        # For QC samples, we don't want any calls to be filtered based on consequence
        if re.search("QC", key) or re.search("TFNA", key):
            command += " --disable_filtering=True"
     
        command = " ".join(command.split())
        print command
        subprocess.call(command, shell=True)

#####################################
#---DYNAMIC ENVIRONMENT VARIABLES---#
#####################################

# Set pipeline executable based on selected pipeline version
TARGETED_NGS_PIPELINE = set_pipeline_version(opts.pipeline_version)

IR_API_KEY = set_IR_API_key_based_on_url(opts.url)

# Change to working directory for the script
# All files will be downloaded in this directory and moved according to their base filename (usually CoPATH ID or other identifier)

if not os.path.exists(WORKING_DIRECTORY):
    os.makedirs(WORKING_DIRECTORY)

######################################
#---CONTACT IR AND SELECT ANALYSES---#
######################################

# Connect to IR server and get basic summary view
all_IR_analyses = IR_analysis_summary_view()

# Split IR analysis names and find unique identifiers
names = select_analyses(all_IR_analyses)

# Initialize dictionary
workflow_dict = {}

# Initialize PrettyPrinter for verbose printing of dictionaries
pp = pprint.PrettyPrinter(indent=4)

# Cycle through unique identifiers and map to different analyses in the server
for name in names:
    
    header = "SEARCHING IR SERVER FOR ANALYSES RELATED TO: %s" % name
    print "-" * len(header)
    print header
    print "-" * len(header)
    
    repeat_flag = False
    workflow_nested_dict, somatic_dict, germline_dict, fusion_dict = (defaultdict(lambda: None) for i in range(4))
    workflow_dict[name] = workflow_nested_dict
    for analysis in all_IR_analyses:
        if (re.search(name,analysis['name']) or name==analysis['name'] or re.match(name,analysis['name'])) and analysis['start_date']==format_date_to_string(opts.date):
            # Pass analyses with designated flags
            pattern = re.compile('BETA|SKIP|FALSE') # Remove 'TEST|RUO'
            if pattern.search(analysis['name'], re.IGNORECASE):
                print "WARNING: Ignoring %s because it contains exclusion flag" % analysis['name']
                try:
                    workflow_dict.pop(name, None)
                except KeyError:
                    print "WARNING: %s is not in list of samples.  Did we remove it already?" % name
                except Exception, e:
                    print "ERROR: Unknown error: %s" % (str(e)) 
            else:
                
                if re.search("SOMATIC",analysis['name'],re.IGNORECASE) or re.search("tumor-normal",analysis['workflow'], re.IGNORECASE):
                    print "Case %s has somatic analysis = %s" % (name,analysis['name'])
                    
                    # Handle repeat analyses by choosing a repeat analysis if it exists.
                    if re.search("repeat", analysis['name'], re.IGNORECASE):
                        repeat_flag = True
                        
                    else:
                        
                        if repeat_flag is True:
                            pass
                        
                        else:
                            repeat_flag = False
                        
                    if repeat_flag is True and not re.search("repeat", analysis['name'], re.IGNORECASE):
                        pass
                    
                    elif repeat_flag is True and re.search("repeat", analysis['name'], re.IGNORECASE):
                        somatic_dict['somatic_analysis_status'] = analysis['status']
                        somatic_dict['somatic_analysis_name'] = analysis['name']
                        somatic_dict['somatic_analysis_id'] = analysis['id']
                        somatic_dict['somatic_tumor_name'] = analysis['samples']['TUMOR']
                        somatic_dict['somatic_normal_name'] = analysis['samples']['NORMAL']
                        somatic_dict['somatic_workflow'] = analysis['workflow']
                        somatic_dict['somatic_workflow_start_date'] = datetime.datetime.strftime(datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y"), "%Y-%m-%d")
                        try:
                            if re.search("_", analysis['name']):
                                workflow_nested_dict['panel_name'] = analysis['name'].split("_")[1]
                            else:
                                workflow_nested_dict['panel_name'] = analysis['name'].split("-")[2]
                        except Exception, e:
                            print "ERROR: Unable to define panel_name"
                            print str(e)
                        
                        try:
                            if re.search("_", analysis['name']):
                                somatic_dict['somatic_analysis_project'] = analysis['name'].split("_")[3]
                        except:
                            print "WARNING: No project has been defined for analysis %s" % analysis['name']
                        try:
                            if re.search("_", analysis['name']):
                                somatic_dict['somatic_analysis_comment'] = analysis['name'].split("_")[4]
                        except:
                            print "WARNING: No comment has been defined for analysis %s" % analysis['name']

                        workflow_nested_dict['somatic_analysis'] = somatic_dict

                    else:
                        somatic_dict['somatic_analysis_status'] = analysis['status']
                        somatic_dict['somatic_analysis_name'] = analysis['name']
                        somatic_dict['somatic_analysis_id'] = analysis['id']
                        somatic_dict['somatic_tumor_name'] = analysis['samples']['TUMOR']
                        somatic_dict['somatic_normal_name'] = analysis['samples']['NORMAL']
                        somatic_dict['somatic_workflow'] = analysis['workflow']
                        somatic_dict['somatic_workflow_start_date'] = datetime.datetime.strftime(datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y"), "%Y-%m-%d")
                        
                        try:
                            if re.search("_", analysis['name']):
                                workflow_nested_dict['panel_name'] = analysis['name'].split("_")[1]
                            else:
                                workflow_nested_dict['panel_name'] = analysis['name'].split("-")[2]
                        except Exception, e:
                            print "ERROR: Unable to define panel_name"
                            print str(e)
                        try:
                            if re.search("_", analysis['name']):
                                somatic_dict['somatic_analysis_project'] = analysis['name'].split("_")[3]
                        except:
                            print "WARNING: No project has been defined for analysis %s" % analysis['name']
                        try:
                            if re.search("_", analysis['name']):
                                somatic_dict['somatic_analysis_comment'] = analysis['name'].split("_")[4]
                        except:
                            print "WARNING: No comment has been defined for analysis %s" % analysis['name']
                            
                        workflow_nested_dict['somatic_analysis'] = somatic_dict    
                            
                else:
                    workflow_nested_dict['somatic_analysis'] = somatic_dict
                    
                if re.search("GERMLINE",analysis['name'],re.IGNORECASE):
                    print "Case %s has germline analysis = %s" % (name,analysis['name'])
                    germline_dict['germline_analysis_status'] = analysis['status']
                    germline_dict['germline_analysis_name'] = analysis['name']
                    germline_dict['germline_analysis_id'] = analysis['id']
                    germline_dict['germline_sample_name'] = analysis['samples']['SAMPLE']
                    germline_dict['germline_control_name'] = analysis['samples']['CONTROL']
                    germline_dict['germline_workflow'] = analysis['workflow']
                    germline_dict['germline_workflow_start_date'] = datetime.datetime.strftime(datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y"), "%Y-%m-%d")
                    
                    try:
                        if re.search("_", analysis['name']):
                            germline_dict['germline_analysis_project'] = analysis['name'].split("_")[3]
                    except:
                        print "WARNING: No project has been defined for analysis %s" % analysis['name']
                    try:
                        if re.search("_", analysis['name']):
                            germline_dict['germline_analysis_comment'] = analysis['name'].split("_")[4]
                    except:
                        print "WARNING: No comment has been defined for analysis %s" % analysis['name']
                    
                    workflow_nested_dict['germline_analysis'] = germline_dict
                else:
                    workflow_nested_dict['germline_analysis'] = germline_dict
                
                if re.search("FUSION",analysis['name'],re.IGNORECASE) or re.search("FUSIONS",analysis['name'],re.IGNORECASE):
                    print "Case %s has fusion analysis = %s" % (name,analysis['name'])
                    fusion_dict['fusion_analysis_status'] = analysis['status']
                    fusion_dict['fusion_analysis_name'] = analysis['name']
                    fusion_dict['fusion_analysis_id'] = analysis['id']
                    fusion_dict['fusion_workflow'] = analysis['workflow']
                    fusion_dict['fusion_workflow_start_date'] = datetime.datetime.strftime(datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y"), "%Y-%m-%d")
                    
                    try:
                        if re.search("_", analysis['name']):
                            fusion_dict['fusion_analysis_project'] = analysis['name'].split("_")[3]
                    except:
                        print "WARNING: No project has been defined for analysis %s" % analysis['name']
                    try:
                        if re.search("_", analysis['name']):
                            fusion_dict['fusion_analysis_comment'] = analysis['name'].split("_")[4]
                    except:
                        print "WARNING: No comment has been defined for analysis %s" % analysis['name']
                
                    
                    workflow_nested_dict['fusion_analysis'] = fusion_dict
                else:
                    workflow_nested_dict['fusion_analysis'] = fusion_dict
                
                workflow_nested_dict['pipeline_start_date'] = datetime.datetime.strftime(datetime.datetime.now(), "%Y-%m-%d")
                workflow_nested_dict['platform'] = 'IonTorrent'
                workflow_nested_dict['pipeline_version'] = opts.pipeline_version


# Print main dictionary for debugging purposes
if opts.verbose is True: pp.pprint(json.loads(json.dumps(workflow_dict)))

######################################
#---PING IR AND DOWNLOAD BAM FILES---#
######################################


header = "FINAL LIST OF SAMPLES FOR PIPELINE ANALYSIS:\n(NOTE: These samples will be checked across database and will fail if they already exist!  Try 'grep ERROR').\n"
print "-" * len(header)
print header
print "%s" % str("\n".join(workflow_dict.keys()))
print "-" * len(header)

counter = 1
for sample_name in workflow_dict.keys():
    
    print "(%d) PROCESSING %s..." % (counter, sample_name)
    counter += 1
    if es.exists(index='dna-seq', doc_type='pipeline-analysis-overview-test', id=sample_name):
        print "ERROR: %s exists in database...passing..." % sample_name
    else:
        if re.search("Oncomine Comprehensive Panel",workflow_dict[sample_name]['somatic_analysis']['somatic_workflow']):
            
            if re.search("SUCCESSFUL", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status']):
            
                if bool(workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_name']) is False:
                    
                    print "WARNING: Proceeding as DNA only OCP..."
                    initiate_IR_download(workflow_dict,sample_name)
                
                elif bool(workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_name']) is True:
                    
                    if re.search("SUCCESSFUL",workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_status']):
                    
                        initiate_IR_download(workflow_dict,sample_name)
                        
                    elif re.search("FAIL", workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_status'], re.IGNORECASE):
                        
                        print "WARNING: SOMATIC ANALYSIS %s FOR %s HAS %s STATUS.  AUTOMATED WORKFLOW WILL NOT BE RUN!" % (workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_name'],
                                                                                                                          sample_name,
                                                                                                                          workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_status'])
                    
                    elif re.search("PENDING", workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_status'], re.IGNORECASE) or re.search("RUNNING", workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_status'], re.IGNORECASE):
                        
                        print "ERROR: Somatic analysis is not complete.  Skipping %s until ready..." % sample_name
                        try:
                            workflow_dict.pop(sample_name, None)
                        except KeyError:
                            print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                        except Exception, e:
                            print "ERROR: Unknown error: %s" % (str(e)) 

                    else:
                        
                        print "ERROR: Unknown status (%s) for %s.  Skipping this sample..." % (workflow_dict[sample_name]['fusion_analysis']['fusion_analysis_status'],
                                                                                               workflow_dict[sample_name]['fusion_analysis']['fusion_workflow'])
                        try:
                            workflow_dict.pop(sample_name, None)
                        except KeyError:
                            print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                        except Exception, e:
                            print "ERROR: Unknown error: %s" % (str(e)) 
    
            elif re.search("FAIL", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'], re.IGNORECASE):
                
                print "WARNING: SOMATIC ANALYSIS %s FOR %s HAS %s STATUS.  AUTOMATED WORKFLOW WILL NOT BE RUN!" % (workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_name'],
                                                                                                                  sample_name,
                                                                                                                  workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'])
            
            elif re.search("PENDING", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'], re.IGNORECASE) or re.search("RUNNING", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'], re.IGNORECASE):
                
                print "ERROR: Somatic analysis is not complete.  Skipping %s until ready..." % sample_name
                try:
                    workflow_dict.pop(sample_name, None)
                except KeyError:
                    print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                except Exception, e:
                    print "ERROR: Unknown error: %s" % (str(e)) 
    
            else:
                
                print "ERROR: Unknown status (%s) for %s.  Skipping this sample..." % (workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'],
                                                                                       workflow_dict[sample_name]['somatic_analysis']['somatic_workflow'])
                try:
                    workflow_dict.pop(sample_name, None)
                except KeyError:
                    print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                except Exception, e:
                    print "ERROR: Unknown error: %s" % (str(e)) 
    
        elif re.search("Comprehensive Cancer Panel", workflow_dict[sample_name]['somatic_analysis']['somatic_workflow']) or re.search("HSM", workflow_dict[sample_name]['panel_name']) or re.search("TFNA", workflow_dict[sample_name]['panel_name']):
            
            if re.search("SUCCESSFUL", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status']):
                
                initiate_IR_download(workflow_dict,sample_name)
                
            elif re.search("FAIL", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'], re.IGNORECASE):
                
                print "WARNING: SOMATIC ANALYSIS %s FOR %s HAS %s STATUS.  AUTOMATED WORKFLOW WILL NOT BE RUN!" % (workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_name'],
                                                                                                                  sample_name,
                                                                                                                  workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'])
            
            elif re.search("PENDING", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'], re.IGNORECASE) or re.search("RUNNING", workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'], re.IGNORECASE):
                
                print "ERROR: Somatic analysis is not complete.  Skipping %s until ready..." % sample_name
                try:
                    workflow_dict.pop(sample_name, None)
                except KeyError:
                    print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                except Exception, e:
                    print "ERROR: Unknown error: %s" % (str(e)) 
    
            else:
                
                print "ERROR: Unknown status (%s) for %s.  Skipping this sample..." % (workflow_dict[sample_name]['somatic_analysis']['somatic_analysis_status'],
                                                                                       workflow_dict[sample_name]['somatic_analysis']['somatic_workflow'])
                try:
                    workflow_dict.pop(sample_name, None)
                except KeyError:
                    print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                except Exception, e:
                    print "ERROR: Unknown error: %s" % (str(e))  
        
        else:
            
            print "ERROR: At this time, this workflow %s is not supported by the IR_heartbeat." % workflow_dict[sample_name]['somatic_analysis']['somatic_workflow']
            try:
                workflow_dict.pop(sample_name, None)
            except KeyError:
                print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
            except Exception, e:
                print "ERROR: Unknown error: %s" % (str(e))  
            #print "ERROR: BUT LETS TRY IT ANYWAY!"
            #initiate_IR_download(workflow_dict, sample_name)
    

# Index BAMs in current working directory
subprocess.call("for i in *.bam; do samtools index $i; done;", shell=True)

# Remote extra log files and IonXpress barcoded BAMs
for f1 in glob.glob('%s/VER*.log' % os.getcwd()):
    os.remove(f1)
    
for f1 in glob.glob('%s/IonXpress*.bam' % os.getcwd()):
    os.remove(f1)

for f1 in glob.glob('%s/*.rrs' % os.getcwd()):
    os.remove(f1)

# If download_bams_only option is enabled, exit without executing pipeline.
# Indexed BAMs should be in the current working directory.
if opts.bams_only is True:
    sys.exit("WARNING: --download_bams_only option set to %. Exiting without executing pipeline..." % opts.bams_only)


################################################
#---DUMP TO ELASTICSEARCH AND BEGIN PIPELINE---#
################################################

# Only execute this loop if we found cases to analyze
if len(workflow_dict) > 0:
        
#     cnx = connect_to_database("tumor_profiling_lab")
#     cursor = cnx.cursor()
    
    for key in workflow_dict.keys():
        entry = json.dumps(workflow_dict[key])
        entry = json.loads(entry)
        if opts.verbose is True: pp.pprint(entry)
        try:
            es.indices.create(index='dna-seq', ignore=400)
            if es.exists(index='dna-seq', doc_type='pipeline-analysis-overview-test', id=key):
                print "ERROR: %s already exists in database...PASSING" % key
                pass
            else:
                res = es.create(index='dna-seq', doc_type='pipeline-analysis-overview-test', body=entry, id=key)
                es.indices.refresh(index="dna-seq")
                print "Was %s created in database?: %s" % (key, res['created'])
                run_pipeline()
            
            
            
#             sql = """INSERT INTO Targeted_NGS_Pipeline_Automation 
#                      (sample_name, panel_name, tumor_bam_name, tumor_barcode, normal_bam_name, normal_barcode, pipeline_start_date, pipeline_start_time, merged_bams, pipeline_status, 
#                      somatic_analysis_name, somatic_analysis_id, somatic_analysis_status, somatic_normal_name, somatic_tumor_name, somatic_workflow, somatic_workflow_start_date, 
#                      germline_analysis_name, germline_analysis_id, germline_analysis_status, germline_control_name, germline_sample_name, germline_workflow, germline_workflow_start_date, 
#                      fusion_analysis_name, fusion_analysis_id, fusion_analysis_status, fusion_workflow, fusion_workflow_start_date) 
#                      VALUES 
#                      (%s, %s, %s, %s, %s, %s, %s, %s, '%s', %s, 
#                      %s, %s, %s, %s, %s, %s, %s, 
#                      %s, %s, %s, %s, %s, %s, %s, 
#                      %s, %s, %s, %s, %s)"""
#             cursor.execute(sql, (key, workflow_dict[key]['panel_name'], workflow_dict[key]['tumor_bam_name'], workflow_dict[key]['tumor_barcode'], workflow_dict[key]['normal_bam_name'], workflow_dict[key]['normal_barcode'], workflow_dict[key]['pipeline_start_date'],workflow_dict[key]['pipeline_start_time'], workflow_dict[key]['merged_bams'], '', # Pipeline status
#                                                        workflow_dict[key]['somatic_analysis']['somatic_analysis_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_id'], workflow_dict[key]['somatic_analysis']['somatic_analysis_status'], workflow_dict[key]['somatic_analysis']['somatic_normal_name'], workflow_dict[key]['somatic_analysis']['somatic_tumor_name'], workflow_dict[key]['somatic_analysis']['somatic_workflow'], workflow_dict[key]['somatic_analysis']['somatic_workflow_start_date'],
#                                                        workflow_dict[key]['germline_analysis']['germline_analysis_name'], workflow_dict[key]['germline_analysis']['germline_analysis_id'], workflow_dict[key]['germline_analysis']['germline_analysis_status'], workflow_dict[key]['germline_analysis']['germline_control_name'], workflow_dict[key]['germline_analysis']['germline_sample_name'], workflow_dict[key]['germline_analysis']['germline_workflow'], workflow_dict[key]['germline_analysis']['germline_workflow_start_date'],
#                                                        workflow_dict[key]['fusion_analysis']['fusion_analysis_name'], workflow_dict[key]['fusion_analysis']['fusion_analysis_id'], workflow_dict[key]['fusion_analysis']['fusion_analysis_status'], workflow_dict[key]['fusion_analysis']['fusion_workflow'], workflow_dict[key]['fusion_analysis']['fusion_workflow_start_date']
#                                                        ))
#             cnx.commit()

         
        except Exception, e:
            print("Unknown error occurred")
            print(e)
            sys.exit(1)

#     cnx.close()




# connection = pymongo.MongoClient()
# tpl_db = connection['tpl']
# tpl_targeted_ngs_pipeline_log = tpl_db['tpl_targeted_ngs_pipeline_log']
# # for entry in entries:
# #     entry_id = tpl_targeted_ngs_pipeline_log.insert_one(entry).inserted_id
# 
# for document in tpl_targeted_ngs_pipeline_log.find():
#     pp.pprint(document)


# def command_line_parsing(opts):
#     
#     mandatory_options = ['auth_key','url']
#     for m in mandatory_options:
#         # Making sure all mandatory options appeared
#         if not opts.__dict__[m]:
#             print "Mandatory option is missing!\n"
#             parser.print_help()
#             sys.exit()       
# 
#     if opts.input_analysis_name is not None:
#         print "------------------------------------"
#         print "OK: Proceeding in single sample mode"
#         print "------------------------------------"
#         return "SINGLE"
#     elif opts.batch_input is not None:
#         print "------------------------------------"
#         print "OK: Proceeding in batch sample mode"
#         print "------------------------------------"
#         return "BATCH"
#     elif opts.batch_input is not None and opts.input_analysis_name is not None:
#         print "ERROR: Cannot run both single sample and batch sample modes simulataneously"    
#         parser.print_help()
#         sys.exit() 
#     elif opts.batch_input is None and opts.input_analysis_name is None:
#         print "ERROR: No input analysis name provided"
#         parser.print_help()
#         sys.exit()
#     else:
#         print "ERROR: Command line parsing failed"
#         parser.print_help()
#         sys.exit()        

      
# mode = command_line_parsing(opts)
