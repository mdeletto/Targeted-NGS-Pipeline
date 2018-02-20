#!/usr/bin/python

import os
from subprocess import Popen
from subprocess import PIPE
import datetime
import MySQLdb

def connect_to_mysql():
    db = MySQLdb.connect(host="localhost",
                         user="root",
                         passwd="*23Ft198",
                         db="tumor_profiling_lab"
                         )
    return db

def convert_git_timestamp_to_mysql_datetime(git_timestamp):
    """Convert a git timestamp into a mysql datetime."""

    f = "%a %b %d %H:%M:%S %Y"
    
    stripped_timezone_git_timestamp = git_timestamp.split(" ")
    stripped_timezone_git_timestamp = stripped_timezone_git_timestamp[:-1]
    stripped_timezone_git_timestamp = " ".join(stripped_timezone_git_timestamp)
    
    datetime_obj = datetime.datetime.strptime(stripped_timezone_git_timestamp, f)
    git_mysql_datetime = datetime_obj.strftime("%Y-%m-%d %H:%M:%S")
    return git_mysql_datetime
    

GIT_COMMIT_FIELDS = ['id', 'author_name', 'author_email', 'date', 'message']
GIT_LOG_FORMAT = ['%H', '%an', '%ae', '%ad', '%s']

GIT_LOG_FORMAT = '%x1f'.join(GIT_LOG_FORMAT) + '%x1e'

os.chdir("/home/michael/YNHH/Code/Github-mdeletto/Targeted-NGS-Pipeline")

p = Popen('git log --format="%s"' % GIT_LOG_FORMAT, shell=True, stdout=PIPE)
(log, _) = p.communicate()
log = log.strip('\n\x1e').split("\x1e")
log = [row.strip().split("\x1f") for row in log]
log = [dict(zip(GIT_COMMIT_FIELDS, row)) for row in log]

conn = connect_to_mysql()
c = conn.cursor()



                                      

for l in log:
    try:
        c.execute("""INSERT INTO targetedNGSPipelineCommits (gitComments, gitAuthorName, gitAuthorEmail, gitDatetimeOfCommit, gitCommitID)
                     VALUES (%s, %s, %s, %s, %s)""", (l['message'],
                                                      l['author_name'],
                                                      l['author_email'],
                                                      convert_git_timestamp_to_mysql_datetime(l['date']),
                                                      l['id'])
                  )
    except Exception, e:
        code, err_string = e
        print l['id'], code, err_string
        
        

conn.commit()
c.close()
conn.close()