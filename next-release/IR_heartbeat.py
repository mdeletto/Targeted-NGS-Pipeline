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
import paramiko
import pipes
import pysam
import shutil
import glob
import errno
import time
#import pymongo
import MySQLdb
from collections import defaultdict
from elasticsearch import Elasticsearch
from globals import *  

global opts

desc="""IR_heartbeat.py acts as a gateway for initiating the bioinformatics analyses (variant calling, annotation, and other QC) for IonTorrent data.  The script contacts the IR server, pulls in necessary info, initiates the downstream bioinformatics analyses, and logs the activity in Elasticsearch."""
usage="./IR_heartbeat -d <analysis_date>"
version="1.3.3"
parser = optparse.OptionParser(description=desc,usage=usage,version=version)
 
group = optparse.OptionGroup(parser, "General options")
#group.add_option("--authorization","-a", help="Authorization API key for IonReporter",dest='auth_key',action='store')
group.add_option("--ionreporter_url","-i", help="IP address or hostname of IR server (ex: 10.80.157.179)",dest='url',action='store',default="10.80.157.179")
group.add_option("--analysis_date","-d", help="Pull IR analyses from this date.  YYYY-mm-dd format required (e.g. %default).  Default date for analysis is today!",dest='date',action='store',default=datetime.datetime.strftime(datetime.date.today(), "%Y-%m-%d"))
group.add_option("--days", help="Number of days to look back at for IR analyses.  Default is to use all IR analyses.  However, by reducing the number of days, runtime will be faster.",dest='days',action='store',default=None)
group.add_option("--download_bams_only","", help="Pull tumor and normal bams for a given analysis and quit",dest='bams_only',action='store_true',default=False)
group.add_option("--pipeline_version","-p", help="Pipeline version to use.  Default is to use the latest pipeline (%default)",dest='pipeline_version',action='store',default=version)#datetime.date.today())
group.add_option("--additional-pipeline-arguments", help="Arguments to pass to the pipeline for additional control over the pipeline.",dest='pipeline_arguments',action='store',default=None)#datetime.date.today())
group.add_option("--force","-f", help="Force unsupported panels through the pipeline anyway.  NOT RECOMMENDED.",dest='force',action='store_true',default=False)
group.add_option("--skip-es-submit", help="Skip submitting the case ID to ElasticSearch.  Default=[%default]",dest='skip_es_submit',action='store_true',default=False)
group.add_option("--enable-manual-control", help="Manual control over samples entering the pipeline.  Default=[%default]",dest='enable_manual_sample_control',action='store_true',default=False)
group.add_option("--ignore-exclusion", help="Ignore exclusion terms for analyses like BETA,SKIP,FALSE.  This will allow IR_heartbeat to detect skipped analyses.  Default=[%default]",dest='ignore_exclusion',action='store_true',default=False)
group.add_option("-v", help="verbose <[%default]>", dest='verbose',action='store_true',default=False)
parser.add_option_group(group)

(opts, args) = parser.parse_args()

####################################
#---PRESET ENVIRONMENT VARIABLES---#
####################################

try:
    if os.path.exists("/media/targetedNGSPipeline/"):
        if os.access("/media/targetedNGSPipeline/", os.W_OK):
            WORKING_DIRECTORY = "/media/targetedNGSPipeline/"
        else:
            WORKING_DIRECTORY = "/var/www/TPL/ftp/targetedNGSPipeline/"
    else:
        WORKING_DIRECTORY = "/var/www/TPL/ftp/targetedNGSPipeline/"
except:
    WORKING_DIRECTORY = "/var/www/TPL/ftp/targetedNGSPipeline/"

os.chdir(WORKING_DIRECTORY)

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
    supported_versions = ['1.2.1',
                          '1.3.0',
                          '1.3.1',
                          '1.3.2',
                          '1.3.3',
                          'next-release']
    
    if pipeline_version in supported_versions:
        if pipeline_version == '1.3.3':
            TARGETED_NGS_PIPELINE = "/home/michael/YNHH/Code/Github-mdeletto/Targeted-NGS-Pipeline/next-release/Variant_Detection.py"
        else:
            TARGETED_NGS_PIPELINE = "/home/michael/YNHH/Code/Github-mdeletto/Targeted-NGS-Pipeline/v%s/Variant_Detection_v%s.py" % (pipeline_version, pipeline_version)
    else:
        sys.exit("ERROR: Pipeline Version %s is unsupported" % pipeline_version)
        
    if pipeline_version == "next-release":
        TARGETED_NGS_PIPELINE = "/home/michael/YNHH/Code/Github-mdeletto/Targeted-NGS-Pipeline/next-release/Variant_Detection.py"
    
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
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&id=%s" 2> /dev/null""" % (authorization_key,IR_url,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
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
            proc = subprocess.Popen(["""curl -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:%s" "%s" 2> /dev/null -o IR.zip; unzip IR.zip && unzip %s.zip; cp ./Variants/*/*.vcf %s.vcf && cp ./Variants/*/*.tsv %s.tsv; rm -rf %s.zip QC Variants Workflow_Settings""" % (authorization_key,variant_link,basename,basename,basename,basename)],shell=True,stdout=subprocess.PIPE)
        except:
            print "Unable to download and/or unzip IonReporter files.  Aborting..."
            sys.exit(1)

def IR_analysis_summary_view():
    """Grabs the summary view from the IR API for analyses."""
    
    if opts.days is None:
        proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&view=summary" 2> /dev/null""" % (IR_API_KEY,"10.80.157.179")],shell=True,stdout=subprocess.PIPE)
    else:
        today = datetime.datetime.strptime(opts.date, "%Y-%m-%d")
        today_formatted = today.strftime("%Y-%m-%d")
        datediff = today - datetime.timedelta(days=int(opts.days))
        datediff_string = datediff.strftime("%Y-%m-%d")
        
        proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&view=summary&start_date=%s&end_date=%s" 2> /dev/null""" % (IR_API_KEY,"10.80.157.179", datediff_string, today_formatted)],shell=True,stdout=subprocess.PIPE)
    output, err = proc.communicate()
    all_IR_analyses = json.loads(output)
    
    return all_IR_analyses

def select_analyses(all_IR_analyses):
    
    def remove_sample_manually(name):
        input_string = raw_input("Would you like to remove %s from the list of samples for analysis? (Y\N): " % name)
        
        
        if input_string == "Y" or input_string == "Yes":
            names_to_remove.append(name)
        elif input_string == "N" or input_string == "No":
            pass
        else:
            print "ERROR: Input not recognized as a valid value.  Please try again."
            remove_sample_manually(name)
        
        return names_to_remove
    
    
    # Initialize dict for storing basename:IR_somatic_id
    basename_analysis_id_dict = {}
    
    today_date_string = format_date_to_string(opts.date)
    # Only match analyses that have SOMATIC|FUSION|GERMLINE in name
    pattern = re.compile('SOMATIC|Somatic|Fusion|FUSION|Germline|GERMLINE')
    
    for analysis in all_IR_analyses:
        if pattern.search(analysis['name'], re.IGNORECASE):
            
            # format date object into string
            today = datetime.datetime.strptime(opts.date, "%Y-%m-%d")
            today_formatted = today.strftime("%Y-%m-%d")

            if opts.days is not None: 
                # User has entered opts.days, so we look back at all analyses in this time range
                datediff = today - datetime.timedelta(days=int(opts.days))
                
                datediff_string = datediff.strftime("%Y-%m-%d")
                if datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y") <= today and datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y") >= datediff: #and analysis['started_by'] == "Mike D'Eletto":
                    if not re.search('-',analysis['name']) and not re.search('_',analysis['name']):
                        print "ERROR: Analysis string is not in correct format.  Please consult README for using this workflow.  Trying to proceed..."
                        basename_analysis_id_dict[analysis['name'].strip()] = analysis['id']
                    elif not re.search('_',analysis['name']):
                        name_split = analysis['name'].split("-")
                        copath_id = str(name_split[0]+"-"+name_split[1])
                        basename_analysis_id_dict[copath_id.strip()] = analysis['id']
                    elif re.search('_',analysis['name']):
                        name_split = analysis['name'].split("_")
                        copath_id = name_split[0]
                        basename_analysis_id_dict[copath_id.strip()] = analysis['id']
            else:
                if datetime.datetime.strptime(analysis['start_date'], "%B %d, %Y") == today: #and analysis['started_by'] == "Mike D'Eletto":
                    if not re.search('-',analysis['name']) and not re.search('_',analysis['name']):
                        print "ERROR: Analysis string is not in correct format.  Please consult README for using this workflow.  Trying to proceed..."
                        basename_analysis_id_dict[analysis['name'].strip()] = analysis['id']
                    elif not re.search('_',analysis['name']):
                        name_split = analysis['name'].split("-")
                        copath_id = str(name_split[0]+"-"+name_split[1])
                        basename_analysis_id_dict[copath_id.strip()] = analysis['id']
                    elif re.search('_',analysis['name']):
                        name_split = analysis['name'].split("_")
                        copath_id = name_split[0]
                        basename_analysis_id_dict[copath_id.strip()] = analysis['id']
    
    names = list(set(basename_analysis_id_dict.keys()))
    print "DETECTED SAMPLE LIST: %s" % ", ".join(names)
    
    if opts.enable_manual_sample_control is True:
        names_to_remove = []
        for name in names:
            names_to_remove = remove_sample_manually(name)
        
        for name in names_to_remove:
            basename_analysis_id_dict.pop(name, None)
        
        print "REVISED SAMPLE LIST: %s" % ", ".join(names)
        
    return basename_analysis_id_dict

def execute_scp_bam_copy(path,bam_barcode_id):
    subprocess.call("scp -i /home/michael/.ssh/id_rsa -q ionadmin@10.80.157.179:%s ." % path, shell=True)
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
        bam_pipeline_flags['normal_bam_name'] = "PopulationNormal"
        bam_pipeline_flags['normal_bam_path'] = "PopulationNormal"
    else:
        normal_bam_link = "%s%s%s.rrs" % (IR_download_link,IR_download_dir,workflow_dict[sample_name]['somatic_analysis']['somatic_normal_name'])
        normal_bam_sample_definition_file = normal_bam_link.split("/")[-1]
        bams_to_download = {"tumor":[tumor_bam_sample_definition_file,tumor_bam_link],
                            "normal":[normal_bam_sample_definition_file,normal_bam_link]}

    # Download BAMs
    print "DOWNLOADING BAMS..."
    
    for sample_type,nested_list in bams_to_download.iteritems():
        if not os.path.isfile("%s/%s-%s-merged.bam" % (os.getcwd(),sample_name,sample_type)) and not os.path.isfile("%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type)):
            sample_definition_file, bam_link = (i for i in nested_list)
            proc = subprocess.check_output("""curl -s -O -k -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:%s" "%s" 2> /dev/null && unzip -o %s """ % (IR_API_KEY,bam_link,sample_definition_file),shell=True)
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
                bam_pipeline_flags['%s_bam_name' %(sample_type)] = "%s-%s-merged.bam" % (sample_name,sample_type)
                bam_pipeline_flags['%s_bam_path' % (sample_type)] = "%s/%s-%s-merged.bam" % (os.getcwd(),sample_name,sample_type)
            # If single BAM exists in sample definition .rrs file
            else:
                bam_pipeline_flags['merged_bams'] = False
                for direct_bam_filepath, bam_barcode_id in bam_barcode_paths_and_ids.iteritems():
                    local_bam_filepath_name = execute_scp_bam_copy(direct_bam_filepath, bam_barcode_id)
                    shutil.move('%s/%s' % (os.getcwd(),bam_barcode_id),'%s/%s-%s.bam' % (os.getcwd(),sample_name,sample_type))
                #pysam.index("%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type))
                bam_pipeline_flags['%s_barcode' % sample_type] = local_bam_filepath_name  
                bam_pipeline_flags['%s_bam_path' % (sample_type)] = "%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type)
                bam_pipeline_flags['%s_bam_name' % (sample_type)] = "%s-%s.bam" % (sample_name,sample_type)
            
            # Remove sample definition file
            os.remove(sample_definition_file)
            
        elif os.path.isfile("%s/%s-%s-merged.bam" % (os.getcwd(),sample_name,sample_type)):
            bam_pipeline_flags['merged_bams'] = True
            bam_pipeline_flags['%s_bam_name' %(sample_type)] = "%s-%s-merged.bam" % (sample_name,sample_type)
            bam_pipeline_flags['%s_bam_path' % (sample_type)] = "%s/%s-%s-merged.bam" % (os.getcwd(),sample_name,sample_type)
        elif os.path.isfile("%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type)):
            bam_pipeline_flags['merged_bams'] = False
            bam_pipeline_flags['%s_bam_path' % (sample_type)] = "%s/%s-%s.bam" % (os.getcwd(),sample_name,sample_type)
            bam_pipeline_flags['%s_bam_name' % (sample_type)] = "%s-%s.bam" % (sample_name,sample_type)
            

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

    IR_download_URL = re.search("(http.+?download\?filePath=)(.+)",variant_link)

    if IR_download_URL:
        IR_download_link = IR_download_URL.group(1)
        IR_download_dir = IR_download_URL.group(2)
        IR_download_dir = IR_download_dir.split("/")
        IR_download_dir = "/".join(IR_download_dir[:-1])+"/"
    
    bam_dict_options = download_bams_from_IR(IR_download_link, IR_download_dir, sample_name)  # downloads BAMs and return a dictionary consisting of bam names and filepaths
    workflow_dict[sample_name].update(bam_dict_options)

def exists_remote(host, path):
    """Test if a file exists at path on a host accessible with SSH."""
    status = subprocess.call(
        ['ssh', host, 'test -f {}'.format(pipes.quote(path))])
    if status == 0:
        return True
    if status == 1:
        return False
    raise Exception('SSH failed')

def determine_IR_basename_filepath(analysis_name,ionreporter_id,IR_url,authorization_key):
    url_encoded_basename = urllib.quote_plus(analysis_name) # catch whitespace or other problematic characters and encode for url
    
    try:
        if ionreporter_id is None:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&name=%s" 2> /dev/null""" % (authorization_key,IR_url,url_encoded_basename)],shell=True,stdout=subprocess.PIPE)
        else:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&id=%s" 2> /dev/null""" % (authorization_key,IR_url,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
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

def ssh_cmd(current_host, command, sm_username, sm_key, sm_port, print_stdout):
    # make sure command ends in a new line
    if not command.endswith('\n'):
        command = command+'\n'
    # establish our dict for returning
    ssh_cmd_dict = {}
    # fill in the stuff we know
    ssh_cmd_dict['server'] = current_host
    ssh_cmd_dict['command'] = command
    # establish  our variables
    server, username, keyfile, portnumber = ( current_host, sm_username , sm_key , sm_port )
    # make the ssh object
    ssh = paramiko.SSHClient()
    # add new servers to known_hosts automatically
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # try to establish the connection to server
    try:
        ssh.connect(server, username=username, key_filename=keyfile, port=portnumber )
    except:
        # record connection failure
        ssh_cmd_dict['connection_successful'] = False
        # if it doesn't work tell 'em
        print('unable to connect to: '+server)
        return ssh_cmd_dict
    else:
        # record connection success
        ssh_cmd_dict['connection_successful'] = True
        # execute the command storing stdout, stdin, and stderr
        stdin, stdout, stderr = ssh.exec_command(command)
        # wait for exit status (that means command finished)
        exit_status = stdout.channel.recv_exit_status()
        
        # End connection if stderr does not have EOF
        timeout = 60
        endtime = time.time() + timeout
        while not stderr.channel.eof_received:
            time.sleep(1)
            if time.time() > endtime:
                stderr.channel.close()
                break

        # flush commands and cut off more writes
        stdin.flush()
        stdin.channel.shutdown_write()

        # close the connection
        ssh.close()

        # store the stdout
        ssh_cmd_dict['stdout_raw'] = stdout.readlines()

        # create an entry that is pretty to read
        if type(ssh_cmd_dict['stdout_raw']) is ListType:
            ssh_cmd_dict['stdout_print'] = "".join(ssh_cmd_dict['stdout_raw'])
        ssh_cmd_dict['stdout'] = ssh_cmd_dict['stdout_raw']

        # store stdin
        ssh_cmd_dict['stdin'] = command

        # store the stderr
        ssh_cmd_dict['stderr_raw'] = stderr.readlines()
        # create an entry for stderr that is pretty to read
        if type(ssh_cmd_dict['stderr_raw']) is ListType:
            ssh_cmd_dict['stderr_print'] = "".join(ssh_cmd_dict['stderr_raw'])
        ssh_cmd_dict['stderr'] = ssh_cmd_dict['stderr_raw']

        # store exit status
        ssh_cmd_dict['exitstatus'] = exit_status

        if exit_status > 0:
            if ssh_cmd_dict['stderr'] is not None:
                print('Unfortunately the command: '+command+' executed on server: '+server+' had a non-zero exit status. Below is the stderr:')
                print(ssh_cmd_dict['stderr'])
        else:
            if print_stdout:
                if ssh_cmd_dict['stdout_print'] is not None:
                    print(ssh_cmd_dict['stdout_print'])
        return ssh_cmd_dict


def connect_to_IR_server_and_run_command(remote_command):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect('10.80.157.179', username="ionadmin")
    stdin, stdout, stderr = client.exec_command(remote_command)
    exit_status = stdout.channel.recv_exit_status()
    #print "returning immediately"
    #return stdout

    timeout = 240
    endtime = time.time() + timeout
    while not stderr.channel.eof_received:
        time.sleep(1)
        if time.time() > endtime:
            stderr.channel.close()
            print "WARNING: Fusion VCF EOF not received.  Time limit exceeded...breaking..."
            break
     
    if stderr.read() == "":
        return stdout
    else:
        return stdout

def connect_to_IR_server_and_grab_file_contents(filename):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect('10.80.157.179', username="ionadmin")
    
    ftp = client.open_sftp()
    _file = ftp.file(filename, 'r', -1)
    data = _file.read().split("\n")
    
    return data
    


# def perform_qc_check(workflow_dict):
#     for basename in workflow_dict.keys():
#         
#         workflow_dict[basename]['qc']['comment'] = ""
#         workflow_dict[basename]['qc']['status'] = "PASS"
#         
#         if 'somatic_analysis' in workflow_dict[basename].keys():
#             try:
#                 if not workflow_dict[basename]['qc']['mapd'] < 0.5:
#                     workflow_dict[basename]['qc']['status'] = "FAIL"
#                     workflow_dict[basename]['qc']['comment'] += "MAPD >= 0.5;"
#             except KeyError:
#                 pass
# 
#         if 'fusion_analysis' in workflow_dict[basename].keys():
#             try:
#                 if workflow_dict[basename]['fusion_analysis']:
#                     if not workflow_dict[basename]['qc']['totalMappedFusionPanelReads'] >= 20000:
#                         workflow_dict[basename]['qc']['status'] = "FAIL"
#                         workflow_dict[basename]['qc']['comment'] += "totalMappedFusionPanelReads < 20,000;"
#             except KeyError:
#                 pass
#             
#         if workflow_dict[basename]['qc']['comment'] == "":
#             workflow_dict[basename]['qc']['comment'] = None
# 
#     return workflow_dict

def connect_to_mysql_and_load_qc_dict(sample_name):
    """Connect to database, choose sequencing runs that have sampleNames that match IR sampleNames, but limit sequencing runs to the most recent.
    Then perform a QC check on the IR analysis metrics, and updat the database entry.
    """
    
    cnx = connect_to_database("tumor_profiling_lab")
    cursor = cnx.cursor()
    
    # Limit 10, and Order by most recent lastModified
    sql = "SELECT resultName, runStatus, runErrorNotes FROM targetedNGSRunQualityControlMetrics \
            WHERE sampleName LIKE '%s' \
            ORDER BY datetimeCreated DESC LIMIT 10" % sample_name
    
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Fetch all the rows in a list of lists.
        results = cursor.fetchall()
        if len(results) > 0:
            row_counter = 0
            for row in results:
                row_counter += 1
                if row_counter > 1:
                    break
                else:
                    resultName = row[0]
                    runStatus = row[1]
                    runErrorNotes = row[2]
                    
            # Prepare values for MySQL update
            qc_dict = defaultdict(lambda: None, workflow_dict[sample_name]['qc'])
            
            if runErrorNotes is None or runErrorNotes == "":
               runErrorNotes = "" 
            
            if 'mapd' in workflow_dict[sample_name]['qc']: #and not re.search("Population", workflow_dict[sample_name]['somatic_analysis']['somatic_normal_name']):
                print "OK: CNV Assessment detected for %s" % sample_name
                qc_dict['cnvAssessment'] = 1
                if float(workflow_dict[sample_name]['qc']['mapd']) >= 0.5:
                    runStatus = 'FAIL'
                    runErrorNotes += 'MAPD Score >= 0.5; '
            else:
                qc_dict['cnvAssessment'] = 0
                    
            if 'totalMappedFusionPanelReads' in workflow_dict[sample_name]['qc'] and re.search("Oncomine", workflow_dict[sample_name]['somatic_analysis']['somatic_workflow']):
                print "OK: FUSION Assessment detected for %s" % sample_name
                qc_dict['geneFusionAssessment'] = 1
                if int(workflow_dict[sample_name]['qc']['totalMappedFusionPanelReads']) < 20000:
                    runStatus = 'FAIL'
                    runErrorNotes += 'totalMappedFusionPanelReads < 20,000; '
                if int(workflow_dict[sample_name]['qc']['expr_control_sum']) < 20000:
                    runStatus = 'FAIL'
                    runErrorNotes += 'sumRNAControls < 20,000; '
                if (re.search("OCPv3", workflow_dict[sample_name]['panel_name']) or re.search("OCP", workflow_dict[sample_name]['panel_name'])) and \
                    not re.search("OCPv2", workflow_dict[sample_name]['somatic_analysis']['somatic_workflow']):
                    if 'pool_counts' in workflow_dict[sample_name]['qc']:
                        
#                         poolXML = """<array _type="array">"""
#                         for pool in sorted(workflow_dict[sample_name]['qc']['pool_counts'].keys()):
#                             if workflow_dict[sample_name]['qc']['pool_counts'][pool] < 100000:
#                                 runErrorNotes += 'RNAControlPool#%s < 100,000; ' % pool
#                                 runStatus = 'FAIL'
#                             poolXML += """<XML_Serializer_Tag _originalKey="%s" _type="array">
#                                     <XML_Serializer_Tag _originalKey="Pool Number" _type="string">%s</XML_Serializer_Tag>
#                                     <XML_Serializer_Tag _originalKey="Read Count" _type="string">%s</XML_Serializer_Tag>
#                                     </XML_Serializer_Tag>""" % (str(int(pool)-1), str(pool), str(workflow_dict[sample_name]['qc']['pool_counts'][pool]))
#                         poolXML += """</array>"""
                        
                        poolXML = str(dict(workflow_dict[sample_name]['qc']['pool_counts']))
                        
                        print "HEY LOOK AT ME"
                        print poolXML
                        
                        qc_dict['pool_count_XML'] = poolXML
                    
            else:
                qc_dict['geneFusionAssessment'] = 0
        
#             # Add "FAIL" prefix to the comment, unless it already exists
#             if runErrorNotes != "":
#                 if re.search("FAIL", runErrorNotes, re.IGNORECASE) or re.search("FLAG", runErrorNotes, re.IGNORECASE):
#                     pass
#                 else:
#                     runErrorNotes = "FAIL:" + runErrorNotes
        
            # PERFORM THE MYSQL DATABASE UPDATE FOR THE ANALYSIS
        
            sql = """UPDATE targetedNGSRunQualityControlMetrics SET geneFusionAssessment = %s, 
                                                                    sumRNAControls = %s, 
                                                                    tumorRNATotalMappedFusionPanelReads = %s, 
                                                                    poolCounts = %s,
                                                                    cnvAssessment = %s, 
                                                                    mapdScore = %s, 
                                                                    runStatus = %s, 
                                                                    runErrorNotes = %s, 
                                                                    modifier = %s, 
                                                                    lastModified = %s \
                                                                WHERE resultName = %s AND sampleName = %s""" 
                    
            values = (qc_dict['geneFusionAssessment'],
                      qc_dict['expr_control_sum'],
                      qc_dict['totalMappedFusionPanelReads'],
                      qc_dict['pool_count_XML'],
                      qc_dict['cnvAssessment'],
                      qc_dict['mapd'],
                      runStatus,
                      runErrorNotes,
                      'AUTO',
                      datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                      resultName,
                      sample_name
                      )
          
            try:
                # Execute the SQL command
                cursor.execute(sql, values)
                # Commit changes
                cnx.commit()
            except Exception, e:
                print "ERROR: Could not update database for %s" % sample_name
                print str(e)
                cnx.rollback()
            
            cnx.close()

        else:
            print "WARNING: No matching entries were found in the table for %s" % sample_name
    except Exception, e:
        print "ERRROR: Unable to query database."
        print "ERRROR: %s" % str(e)
    


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
        
    elif opts.pipeline_version=="1.3.0":
        if workflow_dict[key]['fusion_analysis']['fusion_analysis_name'] is None or re.search("None", workflow_dict[key]['fusion_analysis']['fusion_analysis_name']):
            fusion_parameters = ""
        else:
            fusion_parameters = """--ionreporter_fusion_url_bool \
                                   --ionreporter_fusion_analysis_name=%s \
                                   --ionreporter_fusion_id=%s""" % (workflow_dict[key]['fusion_analysis']['fusion_analysis_name'],
                                                                    workflow_dict[key]['fusion_analysis']['fusion_analysis_id'] )
        
        if workflow_dict[key]['germline_analysis']['germline_analysis_name'] is None or re.search("None", workflow_dict[key]['germline_analysis']['germline_analysis_name']):
            germline_parameters = ""
        else:
            germline_parameters = """--ionreporter_germline_url_bool \
                                   --ionreporter_germline_analysis_name=%s \
                                   --ionreporter_germline_id=%s""" % (workflow_dict[key]['germline_analysis']['germline_analysis_name'],
                                                                    workflow_dict[key]['germline_analysis']['germline_analysis_id'] )
        


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
            %s %s""" % (TARGETED_NGS_PIPELINE, sample_name, workflow_dict[key]['panel_name'], workflow_dict[key]['tumor_bam_name'],
                      workflow_dict[key]['normal_bam_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_id'], fusion_parameters, germline_parameters)                               
        
        # For QC samples, we don't want any calls to be filtered based on consequence
        if re.search("QC", key):
            command += " --disable_filtering"
            
        if re.search("TFNA", str(workflow_dict[key]['panel_name'])):
            command += " --disable_filtering --ionreporter_only"
     
        if opts.pipeline_arguments is not None:
            command += str(opts.pipeline_arguments)
     
        command = " ".join(command.split())
        print command
        subprocess.call(command, shell=True)

    elif opts.pipeline_version=="1.3.3" or opts.pipeline_version=="next-release":
        if workflow_dict[key]['fusion_analysis']['fusion_analysis_name'] is None or re.search("None", workflow_dict[key]['fusion_analysis']['fusion_analysis_name']):
            fusion_parameters = ""
        else:
            fusion_parameters = """--ionreporter_fusion_url_bool \
                                   --ionreporter_fusion_analysis_name=%s \
                                   --ionreporter_fusion_id=%s""" % (workflow_dict[key]['fusion_analysis']['fusion_analysis_name'],
                                                                    workflow_dict[key]['fusion_analysis']['fusion_analysis_id'] )
        
        if workflow_dict[key]['germline_analysis']['germline_analysis_name'] is None or re.search("None", workflow_dict[key]['germline_analysis']['germline_analysis_name']):
            germline_parameters = ""
        else:
            germline_parameters = """--ionreporter_germline_url_bool \
                                   --ionreporter_germline_analysis_name=%s \
                                   --ionreporter_germline_id=%s""" % (workflow_dict[key]['germline_analysis']['germline_analysis_name'],
                                                                    workflow_dict[key]['germline_analysis']['germline_analysis_id'] )
        


        command = """%s \
            -s All-HC \
            -c %s \
            --regions=%s \
            -t %s \
            -n %s \
            -p IonTorrent \
            
            --ionreporter_version=5.6 \
            --ionreporter_somatic_url_bool \
            --ionreporter_somatic_analysis_name='%s' \
            --ionreporter_somatic_id=%s \
            %s %s""" % (TARGETED_NGS_PIPELINE, sample_name, workflow_dict[key]['panel_name'], workflow_dict[key]['tumor_bam_name'],
                      workflow_dict[key]['normal_bam_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_name'], workflow_dict[key]['somatic_analysis']['somatic_analysis_id'], fusion_parameters, germline_parameters)                               
        
        # For QC samples, we don't want any calls to be filtered based on consequence
        if re.search("QC", key):
            command += " --disable_filtering"
            
        if re.search("TFNA", str(workflow_dict[key]['panel_name'])):
            command += " --disable_filtering --ionreporter_only"
        
        if re.search("WhEx", str(workflow_dict[key]['panel_name'])):
            command += " --ionreporter_only"
            
        if not re.search("--disable_filtering", command):
            command += " --filter_common_germline"
     
        if opts.pipeline_arguments is not None:
            command += str(opts.pipeline_arguments)
     
        command = " ".join(command.split())
        print command
        subprocess.call(command, shell=True)

def move_bams_to_base_directory_and_change_dir(base_output):
    """Create a subdirectory for each case, move BAMs into it, cd into the dir, and index BAMs.
    This allows us to make sure files are specific to the case and standard files like "log.txt" don't get moved from another analysis.
    """
    try:
        if not os.path.exists(base_output):
            os.makedirs(base_output)
        subprocess.call("mv %s*.bam* %s/" % (base_output, base_output),shell=True)
        os.chdir(base_output)
        # Index BAMs in current working directory
        subprocess.call("for i in *.bam; do samtools index $i; done;", shell=True)
    except Exception, e:
        print "ERROR: %s" % str(e)
        print "ERROR: Could not create case directory, move BAMs, and cd into that directory"
        
######################################
#---PING IR AND DOWNLOAD BAM FILES---#
######################################

def IR_analysis_control():

        required_analyses = {
                             'Oncomine Comprehensive' : {
                                                         'somatic_analysis' : 'required',
                                                         'germline_analysis' : 'optional',
                                                         'fusion_analysis' : 'optional'
                                                         },
                             
                             'Comprehensive Cancer Panel' : {
                                                             'somatic_analysis' : 'required',
                                                             'germline_analysis' : 'required',
                                                             'fusion_analysis' : 'N/A'               
                                                             },
                             
                             'Hotspot Mutation Panel' : {
                                                        'somatic_analysis' : 'required',
                                                        'germline_analysis' : 'optional',
                                                        'fusion_analysis' : 'N/A'  
                                                         },
                             
                             'AMGEN' : {
                                        'somatic_analysis' : 'required',
                                        'germline_analysis' : 'optional',
                                        'fusion_analysis' : 'N/A'  
                                         },
                             
                             "Exome" : {
                                        'somatic_analysis' : 'required',
                                        'germline_analysis' : 'optional',
                                        'fusion_analysis' : 'N/A'  
                                         },
                             
                             'TFNA' : {
                                        'somatic_analysis' : 'required',
                                        'germline_analysis' : 'optional',
                                        'fusion_analysis' : 'optional'  
                                       }

                             }

        for panel_name in required_analyses:
            if re.search(panel_name, str(workflow_dict[sample_name]['somatic_analysis']['somatic_workflow'])):
                detected_panel_name = panel_name
            else:
                pass

        try:
            detected_panel_name = detected_panel_name
        except NameError:
            if opts.force is True:
                print "WARNING: Panel name could not be autodetected.  Attempting to force this panel through the pipeline, since --force option was invoked."
                detected_panel_name = "UNK"
                required_analyses.update({"UNK" : {
                                                    'somatic_analysis' : 'required',
                                                    'germline_analysis' : 'optional',
                                                    'fusion_analysis' : 'optional'  
                                                   }
                                           })
            else:
                print "ERROR: Unable to determine panel name.\nERROR: Somatic workflow required for the IR_heartbeat to begin"
                return True
        
        for required_analysis in required_analyses[detected_panel_name]:
            if required_analyses[detected_panel_name][required_analysis] == 'required':
                if bool(workflow_dict[sample_name][required_analysis]['%s_name' % required_analysis]) is False:
                    if re.search("QC", sample_name):
                        pass
                    else:
                        print "ERROR: %s required for %s" % (required_analysis, detected_panel_name)
                        workflow_dict.pop(sample_name, None)
                else:
                    if re.search("SUCCESSFUL", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis]):
                        
                        if required_analysis == 'somatic_analysis':
                            
                            pass
                        
                    elif re.search("FAIL", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis]):
                        
                        print "WARNING: %s FOR %s HAS %s STATUS.  AUTOMATED WORKFLOW WILL NOT BE RUN!" % (workflow_dict[sample_name][required_analysis]['%s_name' % required_analysis],
                                                                                                    sample_name,
                                                                                                    workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis])
                        
                        
                        
                    elif re.search("PENDING", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis], re.IGNORECASE) or re.search("RUNNING", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis], re.IGNORECASE):
                        
                        print "ERROR: %s is not complete.  Skipping %s until ready..." % (required_analysis, sample_name)
                        try:
                            workflow_dict.pop(sample_name, None)
                        except KeyError:
                            print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                        except Exception, e:
                            print "ERROR: Unknown error: %s" % (str(e)) 
            
                    else:
                        
                        print "ERROR: Unknown status (%s) for %s.  Skipping this sample..." % (workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis],
                                                                                               workflow_dict[sample_name][required_analysis]['%s_workflow' % required_analysis.split("_")[0]])
                        try:
                            workflow_dict.pop(sample_name, None)
                        except KeyError:
                            print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                        except Exception, e:
                            print "ERROR: Unknown error: %s" % (str(e))
                            
            elif required_analyses[detected_panel_name][required_analysis] == 'optional':
                if bool(workflow_dict[sample_name][required_analysis]['%s_name' % required_analysis]) is False:
                    print "WARNING: %s not detected for %s" % (required_analysis, detected_panel_name)
                else:
                    if re.search("SUCCESSFUL", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis]):
                        
                        if required_analysis == 'somatic_analysis':
                            
                            pass
                        
                    elif re.search("FAIL", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis]):
                        
                        "WARNING: %s FOR %s HAS %s STATUS.  AUTOMATED WORKFLOW WILL NOT BE RUN!" % (workflow_dict[sample_name][required_analysis]['%s_name' % required_analysis],
                                                                                                    sample_name,
                                                                                                    workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis])
                        
                        
                        
                    elif re.search("PENDING", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis], re.IGNORECASE) or re.search("RUNNING", workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis], re.IGNORECASE):
                        
                        print "ERROR: %s is not complete.  Skipping %s until ready..." % (required_analysis, sample_name)
                        try:
                            workflow_dict.pop(sample_name, None)
                        except KeyError:
                            print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                        except Exception, e:
                            print "ERROR: Unknown error: %s" % (str(e)) 
            
                    else:
                        
                        print "ERROR: Unknown status (%s) for %s.  Skipping this sample..." % (workflow_dict[sample_name][required_analysis]['%s_status' % required_analysis],
                                                                                               workflow_dict[sample_name][required_analysis]['%s_workflow' % required_analysis.split("_")[0]])
                        try:
                            workflow_dict.pop(sample_name, None)
                        except KeyError:
                            print "WARNING: %s is not in list of samples.  Did we remove it already?" % sample_name
                        except Exception, e:
                            print "ERROR: Unknown error: %s" % (str(e))

            else:
                
                pass

def define_user_abbreviation(analysis_dict):
    """Takes the IR analysis dict and abbreviates the 'started_by' field.
    """
    
    if re.search("Mike", str(analysis_dict['started_by'])):
        user_abbreviation = "MD"
    elif re.search("Jay", str(analysis_dict['started_by'])):
        user_abbreviation = "JW"
    elif re.search("Sand", str(analysis_dict['started_by'])):
        user_abbreviation = "SC"
    elif re.search("Jen", str(analysis_dict['started_by'])):
        user_abbreviation = "JH"
    else:
        user_abbreviation = "Default"
    
    return user_abbreviation

def query_mysql_to_check_for_basename(name):
    """Connect to database, choose sequencing runs that have sampleNames that match IR sampleNames, but limit sequencing runs to the most recent.
    Then perform a QC check on the IR analysis metrics, and updat the database entry.
    """
    
    cnx = connect_to_database("tumor_profiling_lab")
    cursor = cnx.cursor()
    
    # Limit 10, and Order by most recent lastModified
    sql = "SELECT analysisBasename, analysisID, datetimeCreated FROM targetedNGSPipelineQueue \
            WHERE analysisID LIKE '%s' AND analysisBasename LIKE '%s' \
            ORDER BY datetimeCreated DESC LIMIT 10" % (name+"_"+datetime.datetime.now().strftime('%Y-%m-%d'), name)
    
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Fetch all the rows in a list of lists.
        results = cursor.fetchall()

        if len(results) > 0:
            return False
        else:
            return True
    except Exception, e:
        print "ERROR: Unable to query database."
        print "ERROR: %s" % str(e)

def insert_case_in_queue(name):

    # PERFORM THE MYSQL DATABASE UPDATE FOR THE ANALYSIS

    cnx = connect_to_database("tumor_profiling_lab")
    cursor = cnx.cursor()

    sql = "INSERT INTO targetedNGSPipelineQueue(analysisBasename, \
           analysisID, datetimeCreated, modifier, lastModified, pipelineStage) \
           VALUES ('%s', '%s', '%s', '%s', '%s', '%s')" % \
           (name, name+"_"+datetime.datetime.now().strftime('%Y-%m-%d'), datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'AUTO', datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Initializing")

    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit changes
        cnx.commit()
    except Exception, e:
        print "ERROR: Could not update database for %s" % name
        print str(e)
        cnx.rollback()
    
    cnx.close()

def update_case_stage_in_queue(workflow_dict, name, stage):

    # PERFORM THE MYSQL DATABASE UPDATE FOR THE ANALYSIS

    cnx = connect_to_database("tumor_profiling_lab")
    cursor = cnx.cursor()
    
    sql = "UPDATE targetedNGSPipelineQueue \
           SET pipelineStage = '%s' \
           WHERE analysisID = '%s'" % (stage, name+"_"+datetime.datetime.now().strftime('%Y-%m-%d'))      
                  
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit changes
        cnx.commit()
    except Exception, e:
        print "ERROR: Could not update database for %s" % name
        print str(e)
        cnx.rollback()
    
    cnx.close()

def delete_case_from_queue(name):

    # PERFORM THE MYSQL DATABASE UPDATE FOR THE ANALYSIS

    cnx = connect_to_database("tumor_profiling_lab")
    cursor = cnx.cursor()

    sql = "DELETE FROM targetedNGSPipelineQueue \
           WHERE analysisID = '%s'" % (name+"_"+datetime.datetime.now().strftime('%Y-%m-%d'))        
                  
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit changes
        cnx.commit()
    except Exception, e:
        print "ERROR: Could not update database for %s" % name
        print str(e)
        cnx.rollback()
    
    cnx.close()

def basename_to_analyses_mapping(name):
    
    repeat_flag = False
    workflow_nested_dict, somatic_dict, germline_dict, fusion_dict, qc_dict = (defaultdict(lambda: None) for i in range(5))
    workflow_dict[name] = workflow_nested_dict
    for analysis in all_IR_analyses:
        if (re.search(name,analysis['name']) or name==analysis['name'] or re.match(name,analysis['name'])):
            # Pass analyses with designated flags
            pattern = re.compile('BETA|SKIP|FALSE') # Remove 'TEST|RUO'
            if pattern.search(analysis['name'], re.IGNORECASE) and opts.ignore_exclusion is False:
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
                    
                    # define user
                    workflow_nested_dict['user_full'] = analysis['started_by']
                    user = define_user_abbreviation(analysis)
                    workflow_nested_dict['user_abbrev'] = user
                    
                    # Get base filepath on IR server
                    somatic_base_filepath = determine_IR_basename_filepath(analysis['name'], analysis['id'], opts.url, IR_API_KEY)
                    trash, dirpath = somatic_base_filepath.split("=")
                    dirpath = "/".join(dirpath.split("/")[:-1])
                    try:
                        remote_command_output = connect_to_IR_server_and_run_command("grep '##mapd' %s/outputs/AnnotatorActor-00/annotated_variants.vcf" % dirpath)
                        match = re.search("##mapd=(.*)", remote_command_output.read().strip())
                        mapd = round(float(match.group(1)), 3)
                        qc_dict['mapd'] = mapd
                    except:
                        pass
                    
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

                    else:
                        somatic_dict['somatic_analysis_status'] = analysis['status']
                        somatic_dict['somatic_analysis_name'] = analysis['name']
                        somatic_dict['somatic_analysis_id'] = analysis['id']
                        if 'TUMOR' in analysis['samples'].keys():
                            somatic_dict['somatic_tumor_name'] = analysis['samples']['TUMOR']
                        elif 'PROBAND' in analysis['samples'].keys():
                            somatic_dict['somatic_tumor_name'] = analysis['samples']['PROBAND']

                        if 'NORMAL' in analysis['samples'].keys():
                            somatic_dict['somatic_normal_name'] = analysis['samples']['NORMAL']
                        else:
                            somatic_dict['somatic_normal_name'] = None
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
                    
                    
                    # Get base filepath on IR server
                    somatic_base_filepath = determine_IR_basename_filepath(analysis['name'], analysis['id'], opts.url, IR_API_KEY)
                    trash, dirpath = somatic_base_filepath.split("=")
                    dirpath = "/".join(dirpath.split("/")[:-1])
                    remote_command_output = connect_to_IR_server_and_grab_file_contents("%s/outputs/RNACountsActor-00/fusions.vcf" % dirpath)
                    # Extract FUSION QC info
                    qc_dict['expr_control_sum'] = 0
                    qc_dict['expr_control'] = defaultdict(dict)

                    if workflow_dict[name]['panel_name'] == 'OCPv3' or workflow_dict[name]['panel_name'] == "OCP":
                        
                        # PARSE THE BED FILE FOR POOL INFO
                        
                        control_pool_information = {}
                        #with open(PANEL_REF_DICT[workflow_dict[name]['panel_name']]['fusion_bed']) as bed:
                        with open(PANEL_REF_DICT['OCPv3']['fusion_bed']) as bed:
                            for line in bed.readlines():
                                line = line.strip()
                                if re.search("TYPE=ExpressionControl", line):
                                    gene_match = re.search("GENE_ID=(.+?);", line)
                                    gene = gene_match.group(1)
                                    pool_match = re.search("POOL=(\d+?)", line)
                                    pool = pool_match.group(1)
                                    control_pool_information[gene] = int(pool)
                    
                    # PARSE FUSION VCF and pull relevant QC metrics
                    
                    qc_dict['pool_counts'] = defaultdict(int)
                    for line in remote_command_output:
                        if re.search("SVTYPE=ExprControl", line):
                            gene_match = re.search("GENE_NAME=(.+?);", line)
                            gene = gene_match.group(1)
                            read_count_match = re.search("READ_COUNT=(\d+?);", line)
                            read_count = int(read_count_match.group(1))
                            qc_dict['expr_control'].update({gene : read_count})
                            qc_dict['expr_control_sum'] += int(read_count)
                            
                            if workflow_dict[name]['panel_name'] == 'OCPv3' or workflow_dict[name]['panel_name'] == 'OCP':
                                
                                if gene in control_pool_information.keys():
                                    #qc_dict['pool_counts'][control_pool_information[gene]] = 0
                                    qc_dict['pool_counts'][control_pool_information[gene]] += int(read_count)

                        if re.search("##TotalMappedFusionPanelReads", line):
                            match = re.search("##TotalMappedFusionPanelReads=(.*)", line)
                            total_mapped_fusion_reads = int(match.group(1))
                            qc_dict['totalMappedFusionPanelReads'] = total_mapped_fusion_reads

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
                
                workflow_nested_dict['qc'] = qc_dict
                
                workflow_nested_dict['pipeline_start_date'] = datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y-%m-%d")
                workflow_nested_dict['pipeline_start_utc_timestamp'] = str(datetime.datetime.utcnow())
                workflow_nested_dict['platform'] = 'IonTorrent'
                workflow_nested_dict['pipeline_version'] = opts.pipeline_version

    return workflow_dict

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

IR_API_KEY = set_IR_API_key_based_on_url(opts.url)

# Connect to IR server and get basic summary view
global all_IR_analyses
all_IR_analyses = IR_analysis_summary_view()

# Split IR analysis names and find unique identifiers
basenames_analysis_id_dict = select_analyses(all_IR_analyses)

# Initialize dictionary
workflow_dict = {}

# Initialize PrettyPrinter for verbose printing of dictionaries
pp = pprint.PrettyPrinter(indent=4)

# Cycle through unique identifiers and map to different analyses in the server
for name in basenames_analysis_id_dict.keys():
    header = "SEARCHING IR SERVER FOR ANALYSES RELATED TO: %s" % name
    print "-" * len(header)
    print header
    print "-" * len(header)
    
#     if es.exists(index='dna-seq', doc_type='pipeline-analysis-overview-test', id=name) and opts.skip_es_submit is False:
#         print "ERROR: %s exists in database...passing..." % name
    if query_mysql_to_check_for_basename(name) is False:
        print "ERROR: %s exists in database...passing..." % name
    else:
        # Insert case in pipeline queue
        insert_case_in_queue(name)
        
        workflow_dict = basename_to_analyses_mapping(name)
        
        # Delete case from pipeline queue if we didnt map all analyses correctly
        if not name in workflow_dict.keys():
            delete_case_from_queue(name)


# Perform QC check

#workflow_dict = perform_qc_check(workflow_dict)

# Print main dictionary for debugging purposes
if opts.verbose is True: pp.pprint(json.loads(json.dumps(workflow_dict)))


header = "FINAL LIST OF SAMPLES FOR PIPELINE ANALYSIS:\n(NOTE: These samples will be checked across database and will fail if they already exist!  Try 'grep ERROR').\n"
print "-" * len(header)
print header
print "%s" % str("\n".join(workflow_dict.keys()))
print "-" * len(header)

# print counter and sample_name
counter = 1
for sample_name in workflow_dict.keys():
    print "(%d) PROCESSING %s..." % (counter, sample_name)
    counter += 1
    workflow_dict[sample_name]['IR_heartbeat_error'] = IR_analysis_control()

for sample_name, v1 in workflow_dict.items():
    if isinstance(v1, dict):
        for k2,v2 in v1.items():
            if k2 == 'IR_heartbeat_error' and v2 is True:
                # Delete case from dictionary if there is an error or the case isn't ready yet.
                del workflow_dict[sample_name] 
                delete_case_from_queue(sample_name)

################################################
#---DUMP TO ELASTICSEARCH AND BEGIN PIPELINE---#
################################################

# Only execute this loop if we found cases to analyze
if len(workflow_dict) > 0:
    
    for key in workflow_dict.keys():
        
        # Create user subdirectory if it doesn't exist
        if not os.path.exists(WORKING_DIRECTORY+workflow_dict[key]['user_abbrev']):
            os.makedirs(workflow_dict[key]['user_abbrev'])
        
        # Change directory to user directory
        os.chdir(workflow_dict[key]['user_abbrev'])
        
        update_case_stage_in_queue(workflow_dict, key, "Downloading BAMs")
        
        # Download BAMs
        initiate_IR_download(workflow_dict, key)
        
        # If download_bams_only option is enabled, exit without executing pipeline.
        # Indexed BAMs should be in the current working directory.
        if opts.bams_only is True:
            print "WARNING: --download_bams_only option enabled.  Pipeline execution and database checks will not be performed."
            continue
        else:
    
            # JSON serialize our workflow_dict    
            entry = json.dumps(workflow_dict[key])
            entry = json.loads(entry)
            
            # Print each entry if verbose==True
            if opts.verbose is True: pp.pprint(entry)

            try:
                # Create ES Index
                es.indices.create(index='dna-seq', ignore=400)
                # If index exists for a sample, we pass the sample
                # Unless --skip_es_submit is enabled, then we don't bother checking
                if not es.exists(index='dna-seq', doc_type='pipeline-analysis-overview-test', id=key):
                    res = es.create(index='dna-seq', doc_type='pipeline-analysis-overview-test', body=entry, id=key)
                    es.indices.refresh(index="dna-seq")
                
                update_case_stage_in_queue(workflow_dict, key, "Loading QC data into TPTracker")
                
                # Load entry into mysql TPTracker interface
                connect_to_mysql_and_load_qc_dict(key)
                
                update_case_stage_in_queue(workflow_dict, key, "Moving BAMs to output directory")
                
                # Move BAMs to sample level directory, rather than user directory
                move_bams_to_base_directory_and_change_dir(key)
                
                update_case_stage_in_queue(workflow_dict, key, "Running")
                
                # Run NGSPipeline
                run_pipeline()
                
                update_case_stage_in_queue(workflow_dict, key, "Completed")
                
                # Return to higher level when pipeline has finished
                os.chdir(WORKING_DIRECTORY)

                # Remote extra log files
                files_to_delete = ['%s/VER*.log' % os.getcwd(),
                                   '%s/*.rrs' % os.getcwd()
                                   ]
             
                for filetype in files_to_delete:
                    for file in glob.glob(filetype):
                        os.remove(file)
             
            except Exception, e:
                print("Unknown error occurred")
                print(e)
                sys.exit(1)
      
# mode = command_line_parsing(opts)
