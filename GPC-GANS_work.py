import sys
import python_client
import time
import random
import threading
import json
import os
import socket
import copy

from CallPSGFunctions import *

JE_HOST = os.getenv("JE_HOST", "http://localhost:8007/")
WORKER_THREADS = int( os.getenv("JE_THREADS", 1))
TERMINATE_ON_NO_JOBS = (os.getenv("JE_TERMINATE", "true") == "true")
JE_DEBUG = (os.getenv("JE_DEBUG", "true") == "true")


print(f"Working with JE_HOST {JE_HOST}")
print(f"TERMINATE_ON_NO_JOBS is {TERMINATE_ON_NO_JOBS}")
print(f"Number of Threads is {WORKER_THREADS}")

def myWorker(j,jec):
    payload=j['payload']
    configFile = payload['name']
    try:
        run = copy.deepcopy(configFile)
        while run.find('/') != -1:
            i = run.find('/')
            run = run[i+1:]
        run = './Data/'+run[0:-4]+'.dat'
        # print run
        # print configFile
        run = CallPSG(configFile)
#        if (os.path.exists(run) == False):
#            # print 'No file, running'
#            run = CallPSG(configFile)
#        else:
#            if True: #(os.stat(run).st_size < 100):
#                # print 'Empty file, running', configFile
#                run = CallPSG(configFile)
#            else:
#                # print 'File found', os.stat(direct+run).st_size
#                run = run
        if (os.stat(run).st_size < 100):
            raise Exception('We passed the PSG call, but the output file is still empty.')
        else:
            # PlotPSG(run)
            jec.completeJob(j['_id'], "done")
    except Exception as e:
        jec.failJob(j['_id'], f"{e}")
    # print run

def spawnme(workerid):
    ## define the API basepoint and a worker for worker example
    ## spefify the family parameter if you have more than one job family on this queue so that this client only gets jobs it can handle
    jec = python_client.GpcJobEngineClient(JE_HOST, userWorker= myWorker, terminateOnNoJobs= TERMINATE_ON_NO_JOBS, workerId=workerid, family='GPC-PSG-GANS', DEBUG = JE_DEBUG)
    jec.start()

def main():
    workers = []
    for i in range(WORKER_THREADS):
        t = threading.Thread(target=spawnme(socket.gethostname() + "-" + str(i)))
        workers.append(t)

    for w in workers:
        w.start()

    for w in workers:
        w.join()

    print(f"spawned {WORKER_THREADS} threads")

if __name__ == "__main__":
    main()
