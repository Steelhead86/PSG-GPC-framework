#### #!/usr/bin/python3

import threading
import sys
import requests
import time
import os

class GpcJobEngineClient:
    """GPC Job EngineClient"""

    def __init__(self, url, userWorker = None, noJobsSleep = 5, noJobsSignal = None, workerId = '', terminateOnNoJobs = True, DEBUG = False, family='', **kwargs):
        self.url = url
        self.processing = False
        self.workerThread = None        
        self.noJobsSleep = noJobsSleep
        self.noJobsSignal = noJobsSignal        
        self.userWorker = userWorker
        self.DEBUG = DEBUG
        self.terminateOnNoJobs = terminateOnNoJobs
        self.workerId = workerId
        self.family = family
            
    def log(self, s):
        """
        Log string
        :param s: string to log
        """
        if self.DEBUG == True:
            sys.stdout.write("%s\n" % s)

    def addJob(self, payload, friendlyName = "", jobFamily = "", priority=1):
        """
        Add a job to queue
        
        :param payload: dict of job settings
        :param priority=1: int optional priority for job
        """        
        return self.post("jobs", { 'payload': payload, "friendlyName" : friendlyName, "jobFamily" : jobFamily})

    def addJobsFromFile(self, filepath):
        if os.path.isfile(filepath):
            with open(filepath, "r") as f:
                return self.postFile("jobs", f)
        else:
            raise Exception(f"{filepath} is not a file according to os.path.isfile")

    def get(self, m, params = {}):
        if 'family' not in params:
            params['family'] = self.family
        self.log("Request GET " + self.url + m + " with parameters " + str(params))
        r = requests.get(url = self.url + m, params = params)        
        if r.status_code == 200:
            data = r.json()            
            self.log(data)
            return data
        elif r.status_code == 204:
            return False
        else:
            raise Exception("Not sure how to handle response code {}".format(r.status_code))

    def postFile(self, m, file):            
        self.log("Request POST file " + self.url + m)
        r = requests.post(url = self.url + m, files={'file' : file})        
        if r.status_code == 201:
            data = r.json()            
            self.log(data)
            return data                
        else:
            raise Exception(f"Not sure how to handle response code {r.status_code} message from api: {r.text}")    

    def post(self, m, params = {}):
        if 'family' not in params:
            params['family'] = self.family
        self.log("Request POST " + self.url + m)
        r = requests.post(url = self.url + m, data = params)        
        if r.status_code == 201:
            data = r.json()            
            self.log(data)
            return data                
        else:
            raise Exception(f"Not sure how to handle response code {r.status_code} message from api: {r.text}")

    def put(self, m, params = {}):
        self.log("Request PUT " + self.url + m)
        r = requests.put(url = self.url + m, data = params)        
        if r.status_code == 201:            
            return True       
        else:
            raise Exception(f"Not sure how to handle response code {r.status_code} message from api: {r.text}")
    
    def delete(self, m, params = {}):
        if 'family' not in params:
            params['family'] = self.family
        self.log("Request DELETE " + self.url + m)
        r = requests.delete(url = self.url + m, data = params)        
        if r.status_code == 204:
            return True        
        else:
            raise Exception(f"Not sure how to handle response code {r.status_code} message from api: {r.text}")

    def getJobById(self, _id):
        return self.get(f"jobs/{_id}")
        
    def checkoutJob(self):
        return self.get("jobs/checkout_one", {"workerid" : self.workerId})        
    
    def getJobs(self, limit=1000):
        return self.get("jobs", {"limit" : limit})
    
    def getPendingJobs(self, limit=1000):
        return self.get("jobs", {"limit" : limit, "status" : "pending"})

    def getRunningJobs(self, limit=1000):
        return self.get("jobs", {"limit" : limit, "status" : "running"})

    def getCompletedJobs(self, limit=1000):
        return self.get("jobs", {"limit" : limit, "status" : "completed"})

    def completeJob(self, _id, result):
        return self.put(f"jobs/{_id}", {"result": result })

    def progressJob(self, _id, progress):
        return self.put(f"jobs/{_id}", {"progress" : progress})
    
    def failJob(self, _id, error):
        return self.put(f"jobs/{_id}", {"error": error })

    def deleteAllJobs(self):
        return self.delete("jobs")        

    def deleteJob(self, _id):
        return self.delete(f"jobs/{_id}")
    
    def getFamilies(self):
        return self.get(f"jobs/families")
    
    def getCount(self, family):
        self.log(f"getcount(family={family})")
        return self.get(f"jobs/count", { "family" : f"{family}" })

    def getCounts(self):
        counts = {}
        families = self.getFamilies()        
        for item in families:
            family = item["_id"]
            if family == "":
                continue
            c = self.getCount(family)
            self.log(f"c={c}, family={family}")
            counts[family] = c
        return counts

### refactor, belongs in client side
    def worker(self):
        """
        Internal worker thread entry point
        """
        self.log("Worker thread active")
        
        while True:  # keep looping until signaled to stop          
            if self.processing == False:
                break
            job = None
            try:
                job = self.checkoutJob()
            except Exception as e:
                self.log(f"Exception when trying to checkout a job. Performing a sleep.\n Exception detail {e}")
                time.sleep(self.noJobsSleep)
                continue               
            
            if job == False:
                #no job
                if self.terminateOnNoJobs:
                    # self terminate worker thread since no jobs pending
                    print("Terminating because terminateOnNoJobs is true and there are no jobs pending.")
                    self.stopProcessing()
                else:
                    #long running poll mode
                    time.sleep(self.noJobsSleep)
            else:
                #got a job
                if job == "PAUSED":
                    if self.terminateOnNoJobs:
                        print("Terminating because terminateNoJobs is true and the family queue is paused")
                        self.stopProcessing()
                    else:
                        print("Job family queue is paused")
                        time.sleep(self.noJobsSleep)
                else:                    
                    self.userWorker(job,self)

    def start(self,poolid=None):
        """ starts the worker thread"""
        if self.userWorker == None:
            raise Exception("userWorker not defined")
        else:
            if poolid != None:
                self.log(f"worker pool id is {poolid}")
            self.processing = True
            self.workerThread = threading.Thread(target=self.worker)        
            self.log("Starting job engine")
            self.workerThread.start()

    def stopProcessing(self):
        """ signals that processing should be stopped """
        self.processing = False
    
    def end(self):
        """ stops processing, joins worker thread and returns"""
        self.stopProcessing()
        self.workerThread.join()
        self.log("ending job engine")
    
    def getFamily(self, family):
        self.log(f"in getFamilyConfig family={family}")
        return self.get(f"jobs/families/{family}")
