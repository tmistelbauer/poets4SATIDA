import paramiko
import os
import datetime

transport = paramiko.Transport(('ftp.ipf.tuwien.ac.at', 22))
transport.connect(username='satida', password='Urofeveni468')
sftp = paramiko.SFTPClient.from_transport(transport)
basedir = '/home/tuwien/poets/DATA'
remotedir = '/_up/from_TUW'

try:
    sftp.mkdir(remotedir)
except:
    pass

dirs = ['DI', 'CDI', 'Weights']

for file in os.listdir(basedir):
    if file.endswith(".nc"):
        loc = os.path.join(basedir, file)
        rem = os.path.join(remotedir, file)
        print ('[' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
               + '] copy '+loc+' to '+rem)
        sftp.put(loc, rem)

for dir in dirs:
    sftp.chdir(remotedir)
    for file in os.listdir(os.path.join(basedir, dir)):
        try:
            sftp.mkdir(dir)
        except:
            pass
        if file.endswith(".nc"):
            loc = os.path.join(basedir, dir, file)
            rem = os.path.join(dir, file)
            print ('[' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
                   + '] copy ' + loc + ' to ' +rem)
            sftp.put(loc, rem)

sftp.close()
transport.close()
