"""
Base class and sub-classes for managing parallel processing in RIOS. 

It should be emphasised at the start that it is only worth using
parallel processing in RIOS for tasks which are very computationally
intensive, as there is significant overhead in setting up the sub-jobs. 
Most image processing is I/O bound, and will not benefit from parallel 
processing. 

It should also be noted that the 'otherargs' parameter to :func:`rios.applier.apply`
works for passing data into a user function, but updated data is not passed out at present.

The base class is :class:`rios.parallel.jobmanager.JobManager`. This is an abstract base class,
and must be sub-classed before use. Any sub-class is intended to manage 
processing of the user function in a set of sub-jobs, farming out 
sub-jobs to process them all in parallel, and gathering up the 
results, and re-combining into a single set of outputs. 

Most of this work is handled in the base class, and should be generic
for different methods of parallel processing. The reason for the
sub-classes is to allow different approaches to be used, depending on 
the system configuration. In particular, one can use a cluster batch 
queue system such as PBS or SLURM to run sub-jobs as independent jobs, 
allowing it to manage scheduling of jobs and resource management. 
Alternatively, one can use MPI or Python's own multiprocessing module,
if this is more appropriate for the system configuration available. 

Sub-classes are provided for using PBS, SLURM, MPI, multiprocessing
or Python's native subprocess module. Other sub-classes can be made
as required, outside this module, and will be visible to the function::

    getJobManagerClassByName()

which is the main function used for selecting which sub-class is required.

The calling program controls the parallel processing through the
ApplierControls() object. Normal usage would be as follows::

    from rios import applier
    controls = applier.ApplierControls()
    controls.setNumThreads(5)

If a custom JobManager sub-class is used, its module should be imported 
into the calling program (in order to create the sub-class), but its use 
is selected using the same call to controls.setJobManagerType(), giving 
the jobMgrType of the custom sub-class. 

If $RIOS_DFLT_JOBMGRTYPE is set, this will be used as the default jobMgrType.
This facilitates writing of application code which can run unmodified on 
systems with different configurations. Alternatively, this can be set on
the controls object, e.g.::

    controls.setJobManagerType('pbs')
    
Environment Variables
---------------------

+---------------------------------+-----------------------------------------------------------------+
| Name                            | Description                                                     |
+=================================+=================================================================+
|RIOS_DFLT_JOBMGRTYPE             | Name string of default JobManager subclass                      |
+---------------------------------+-----------------------------------------------------------------+
|RIOS_PBSJOBMGR_QSUBOPTIONS       | String of commandline options to be used with PBS qsub.         |
|                                 | Use this for things like walltime and queue name.               |
+---------------------------------+-----------------------------------------------------------------+
|RIOS_PBSJOBMGR_INITCMDS          | String of shell command(s) which will be executed               |
|                                 | inside each PBS job, before executing the                       |
|                                 | processing commands. Not generally required, but was            |
|                                 | useful for initial testing.                                     |
+---------------------------------+-----------------------------------------------------------------+
|RIOS_SLURMJOBMGR_SBATCHOPTIONS   | String of commandline options to be used with SLURM             |
|                                 | sbatch. Use this for things like walltime and queue name.       |
|                                 | The output and error logs do not need to be specified - they    |
|                                 | are set to temporary filenames by RIOS.                         |
+---------------------------------+-----------------------------------------------------------------+
|RIOS_SLURMJOBMGR_INITCMDS        | String of shell command(s) which will be executed               |
|                                 | inside each SLURM job, before executing the                     |
|                                 | processing commands. Not generally required, but was            |
|                                 | useful for initial testing.                                     |
+---------------------------------+-----------------------------------------------------------------+

"""
# This file is part of RIOS - Raster I/O Simplification
# Copyright (C) 2012  Sam Gillingham, Neil Flood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function

import os
import sys
import math
import abc
import copy
import subprocess
import tempfile
import time
try:
    import cPickle as pickle        # For Python 2.x
except ImportError:
    import pickle

import numpy

from .. import rioserrors
# Import a pickler which can pickle functions, with their dependencies, as well
# as data. Either use the version installed with cloudpickle
# (https://github.com/cloudpipe/cloudpickle) or the bundled versin
try:
    from cloudpickle import cloudpickle
except ImportError:
    # Import from our own local copy. This is what will usually happen. 
    from . import cloudpickle

class BlockAssociations(object):
    """
    Dummy class, to mimic applier.BlockAssociations, while avoiding circular imports. 
    """

class JobInfo(object):
    """
    Abstract base class for the information that needs to be passed to a job.

    """
    __metaclass__ = abc.ABCMeta

    def prepareForPickling(self):
        """
        Returns an instance of JobInfo to be pickled.
        Normally derived classes will just return 'self'
        but in some cases more complicated processing can
        happen
        """
        return self

    @abc.abstractmethod
    def getFunctionParams(self):
        """
        Return the parameters to be passed to the actual function
        in the sub process. Should return a tuple.
        """

    @abc.abstractmethod
    def getFunctionResult(self, params):
        """
        Return the parameter(s) that were modified
        by the function so they can be returned.
        """

class JobManager(object):
    """
    Manage breaking up of RIOS processing into sub-jobs, and farming them out. 
    
    Should be sub-classed to create new ways of farming out jobs. The sub-class 
    should at least over-ride the following abstract methods::

        startOneJob()
        waitOnJobs()
        gatherAllOutputs()

    More sophisticated sub-classes might also need to over-ride::

        startAllJobs()
    
    A sub-class must also include a class attribute called jobMgrType, which has 
    string value, which is the name used to select this sub-class. 
    
    """
    __metaclass__ = abc.ABCMeta
    jobMgrType = None
    
    def __init__(self, numSubJobs):
        """
        numSubJobs is the number of sub-jobs
        """
        self.numSubJobs = numSubJobs
        self.tempdir = '.'
    
    def setTempdir(self, tempdir):
        """
        Directory to use for temporary files. This is generally set by apply(),
        using the one it has been given on the ApplierControls object. The
        default is '.'. 
        
        """
        self.tempdir = tempdir
    
    def runSubJobs(self, function, fnInputs):
        """
        Take the given list of function arguments, run the given function 
        for each one, as a separate asynchronous job. 

        Returns a list of output BlockAssociations.
        
        """
        jobIDlist = self.startAllJobs(function, fnInputs)
        self.waitOnJobs(jobIDlist)
        
        outputBlocksList = self.gatherAllOutputs(jobIDlist)
        return outputBlocksList
    
    def startAllJobs(self, function, jobInputs):
        """
        Start up all of the jobs processing blocks. Default implementation
        loops over the lists of jobs, starting each job separately. Keeps the
        first job aside, and runs it here before going off to wait for 
        the others. This means that the first element in the jobID list is not
        a jobID, but the results of the first sub-job. 

        jobInputs should be a list of JobInfo derived objects.
        
        """
        jobIDlist = [None]
        for inputs in jobInputs[1:]:
                
            jobID = self.startOneJob(function, inputs)
            jobIDlist.append(jobID)
        
        # Run the first one here
        inputs = jobInputs[0]
        params = inputs.getFunctionParams()
        function(*params)

        jobIDlist[0] = inputs.getFunctionResult(params)

        return jobIDlist
    
    @abc.abstractmethod
    def startOneJob(self, userFunc, jobInfo):
        """
        Start one job. Return a jobID object suitable for identifying the
        job, with all information required to wait for it, and 
        recover its output. This jobID is specific to the subclass. 
        
        This is an abstract method, and must be over-ridden in a sub-class.

        jobInfo should be a JobInfo derived object.        

        """
    
    @abc.abstractmethod
    def waitOnJobs(self, jobIDlist):
        """
        Wait until all the jobs in the given list have completed. This is
        an abstract method, and must be over-ridden in a sub-class. 
        
        """

    @abc.abstractmethod
    def gatherAllOutputs(self, jobIDlist):
        """
        Gather up outputs from sub-jobs, and return a list of the
        outputs objects. 
        This is an abstract method, and must be over-ridden in a sub-class. 
        
        """
    
    def __str__(self):
        """
        String representation
        """
        return "jobMgrType=%s, numSubJobs=%s" % (self.jobMgrType, self.numSubJobs)


class SubprocJobManager(JobManager):
    """
    Use Python's standard subprocess module to run individual jobs. 
    Passes input and output to/from the subprocesses using their 
    stdin/stdout. The command being executed is a simple main program
    which runs the user function on the given data, and passes back
    the resulting outputs object. 
    
    This JobManager sub-class should be used with caution, as it does not 
    involve any kind of load balancing, and all sub-processes simply run 
    concurrently. If you have enough spare cores and memory to do that, then
    no problem, but if not, you may clog the system. 
    
    """
    jobMgrType = "subproc"
    
    def startOneJob(self, userFunc, jobInfo):
        """
        Start one job. We execute the rios_subproc.py command,
        communicating via its stdin/stdout. We give it the pickled
        function and all input objects, and we get back a pickled
        outputs object. 
        
        """
        jobInfo = jobInfo.prepareForPickling()

        allInputs = (userFunc, jobInfo)
        allInputsPickled = cloudpickle.dumps(allInputs)

        proc = subprocess.Popen(['rios_subproc.py'], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE)
        proc.stdin.write(allInputsPickled)

        return proc
    
    def waitOnJobs(self, jobIDlist):
        """
        Wait until all the jobs in the given list have completed. This
        implementation doesn't wait at all, because the subprocesses may
        block on writing their output to the stdout pipe. So, we do  
        nothing here, and actually wait on the read of the stdout from
        the subprocesses. 
        
        """

    def gatherAllOutputs(self, jobIDlist):
        """
        Gather up outputs from sub-jobs, and return a list of the
        outputs objects. Note that we assume that the first element of
        jobIDlist is actually an outputs object, from running the first sub-array
        in the current process. 
        
        """
        outputBlocksList = [jobIDlist[0]]
        for proc in jobIDlist[1:]:
            pickledOutput = proc.stdout.read()
            outputObj = pickle.loads(pickledOutput)
            outputBlocksList.append(outputObj)
        return outputBlocksList


class PbsJobManager(JobManager):
    """
    Use PBS to run individual jobs
    
    """
    jobMgrType = "pbs"
    
    def startOneJob(self, userFunc, jobInfo):
        """
        Start one job. We create a shell script to submit to a PBS batch queue.
        When executed, the job will execute the rios_subproc.py command, giving
        it the names of two pickle files. The first is the pickle of all inputs
        (including the function), and the second is where it will write the 
        pickle of outputs. 
        
        Uses $RIOS_PBSJOBMGR_QSUBOPTIONS to pick up any desired options to the 
        qsub command. This should be used to control such things as requested 
        amount of memory or walltime for each job, which will otherwise be
        defaulted by PBS. 
        
        """
        jobInfo = jobInfo.prepareForPickling()

        allInputs = (userFunc, jobInfo)
        allInputsPickled = cloudpickle.dumps(allInputs)
        
        (fd, inputsfile) = tempfile.mkstemp(prefix='rios_pbsin_', dir=self.tempdir, suffix='.tmp')
        os.close(fd)
        outputsfile = inputsfile.replace('pbsin', 'pbsout')
        scriptfile = inputsfile.replace('pbsin', 'pbs').replace('.tmp', '.sh')
        logfile = outputsfile.replace('.tmp', '.log')
        
        qsubOptions = os.getenv('RIOS_PBSJOBMGR_QSUBOPTIONS')
        
        scriptCmdList = [
            "#!/bin/bash",
            "#PBS -j oe -o %s" % logfile
        ]
        if qsubOptions is not None:
            scriptCmdList.append("#PBS %s" % qsubOptions)
            
        pbsInitCmds = os.getenv('RIOS_PBSJOBMGR_INITCMDS')
        if pbsInitCmds is not None:
            scriptCmdList.append(pbsInitCmds)
            
        scriptCmdList.append("rios_subproc.py %s %s"%(inputsfile, outputsfile))
        scriptStr = '\n'.join(scriptCmdList)
        
        open(scriptfile, 'w').write(scriptStr+'\n')
        open(inputsfile, 'wb').write(allInputsPickled)
        
        submitCmdWords = ["qsub", scriptfile]
        proc = subprocess.Popen(submitCmdWords, stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, universal_newlines=True)
        # The qsub command exits almost immediately, printing the PBS job id
        # to stdout. So, we just wait for the qsub to finish, and grab the
        # jobID string.     
        (stdout, stderr) = proc.communicate()
        pbsJobID = stdout.strip()
        
        # Remove the script file, assuming that qsub took a copy of it. 
        os.remove(scriptfile)

        # If there was something in stderr from the qsub command, then probably 
        # something bad happened, so we pass it on to the user in the form of
        # an exception. 
        if len(stderr) > 0:
            msg = "Error from qsub. Message:\n"+stderr
            raise rioserrors.JobMgrError(msg)
        
        return (pbsJobID, outputsfile, logfile)
    
    def waitOnJobs(self, jobIDlist):
        """
        Wait until all jobs in the given list have completed. The jobID values
        are tuples whose first element is a PBS job id string. We poll the PBS 
        queue until none of them are left in the queue, and then return. 
        
        Note that this also assumes the technique used by the default startAllJobs()
        method, of executing the first job in the current process, and so the first
        jobID is not a jobID but the results of that. Hence we do not try to wait on
        that job, but on all the rest. 
        
        Returns only when all the listed jobID strings are no longer found in the
        PBS queue. Currently has no time-out, although perhaps it should. 
        
        """
        allFinished = False
        
        # Extract the actual PBS job ID strings, skipping the first element. 
        # Express as a set, for efficiency later on
        pbsJobIdSet = set([t[0] for t in jobIDlist[1:]])
        
        while not allFinished:
            qstatCmd = ["qstat"]
            proc = subprocess.Popen(qstatCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True)
            (stdout, stderr) = proc.communicate()
            
            stdoutLines = [line for line in stdout.split('\n') if len(line) > 0]   # No blank lines
            # Skip header lines, and grab first word on each line, which is jobID
            qstatJobIDlist = [line.split()[0] for line in stdoutLines[2:]]
            qstatJobIDset = set(qstatJobIDlist)
            
            allFinished = pbsJobIdSet.isdisjoint(qstatJobIDset)
            
            if not allFinished:
                # Sleep for a bit before checking again
                time.sleep(60)
    
    def gatherAllOutputs(self, jobIDlist):
        """
        Gather up outputs from sub-jobs, and return a list of the
        outputs objects. Note that we assume that the first element of
        jobIDlist is actually an outputs object, from running the first sub-array
        in the current process. 
        
        The jobIDlist is a list of tuples whose second element is the name of 
        the output file containing the pickled outputs object. 
        
        """
        outputBlocksList = [jobIDlist[0]]
        for (jobID, outputsfile, logfile) in jobIDlist[1:]:
            try:
                pickledOutput = open(outputsfile, 'rb').read()
                outputObj = pickle.loads(pickledOutput)
                os.remove(outputsfile)
            except Exception as e:
                logfileContents = 'No logfile found'
                if os.path.exists(logfile):
                    logfileContents = open(logfile).read()
                msg = ("Error collecting output from PBS sub-job. Exception message:\n"+str(e)+
                    "\nPBS Logfile:\n"+logfileContents)
                raise rioserrors.JobMgrError(msg)
            outputBlocksList.append(outputObj)
            os.remove(logfile)
        return outputBlocksList
        
    
class SlurmJobManager(JobManager):
    """
    Use SLURM to run individual jobs
    
    """
    jobMgrType = "slurm"

    def startOneJob(self, userFunc, jobInfo):
        """
        Start one job. We create a shell script to submit to a SLURM batch queue.
        When executed, the job will execute the rios_subproc.py command, giving
        it the names of two pickle files. The first is the pickle of all inputs
        (including the function), and the second is where it will write the 
        pickle of outputs. 
        
        Uses $RIOS_SLURMJOBMGR_SBATCHOPTIONS to pick up any desired options to the 
        sbatch command. This should be used to control such things as requested 
        amount of memory or walltime for each job, which will otherwise be
        defaulted by SLURM. 
        
        """
        jobInfo = jobInfo.prepareForPickling()

        allInputs = (userFunc, jobInfo)
        allInputsPickled = cloudpickle.dumps(allInputs)
        
        (fd, inputsfile) = tempfile.mkstemp(prefix='rios_slurmin_', dir=self.tempdir, suffix='.tmp')
        os.close(fd)
        outputsfile = inputsfile.replace('slurmin', 'slurmout')
        scriptfile = inputsfile.replace('slurmin', 'slurm').replace('.tmp', '.sl')
        logfile = outputsfile.replace('.tmp', '.log')
        
        sbatchOptions = os.getenv('RIOS_SLURMJOBMGR_SBATCHOPTIONS')
        
        scriptCmdList = [
            "#!/bin/bash",
            "#SBATCH -o %s" % logfile,
            "#SBATCH -e %s" % logfile
        ]
        if sbatchOptions is not None:
            scriptCmdList.append("#SBATCH %s" % sbatchOptions)
            
        slurmInitCmds = os.getenv('RIOS_SLURMJOBMGR_INITCMDS')
        if slurmInitCmds is not None:
            scriptCmdList.append(slurmInitCmds)
            
        scriptCmdList.append("rios_subproc.py %s %s"%(inputsfile, outputsfile))
        scriptStr = '\n'.join(scriptCmdList)
        
        open(scriptfile, 'w').write(scriptStr+'\n')
        open(inputsfile, 'wb').write(allInputsPickled)
        
        submitCmdWords = ["sbatch", scriptfile]
        proc = subprocess.Popen(submitCmdWords, stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, universal_newlines=True)
        # The sbatch command exits almost immediately, printing the SLURM job id
        # to stdout. So, we just wait for the sbatch to finish, and grab the
        # jobID string.     
        (stdout, stderr) = proc.communicate()
        slurmOutputList = stdout.strip().split()
        slurmJobID = None
        # slurm prints a sentence to the stdout:
        # 'Submitted batch job X'
        if len(slurmOutputList) >= 4:
            slurmJobID = slurmOutputList[3]
        
        # Remove the script file, assuming that sbatch took a copy of it. 
        os.remove(scriptfile)

        # If there was something in stderr from the sbatch command, then probably 
        # something bad happened, so we pass it on to the user in the form of
        # an exception. 
        if slurmJobID is None or len(stderr) > 0:
            msg = "Error from sbatch. Message:\n"+stderr
            raise rioserrors.JobMgrError(msg)
        
        return (slurmJobID, outputsfile, logfile)
    
    def waitOnJobs(self, jobIDlist):
        """
        Wait until all jobs in the given list have completed. The jobID values
        are tuples whose first element is a SLURM job id string. We poll the SLURM 
        queue until none of them are left in the queue, and then return. 
        
        Note that this also assumes the technique used by the default startAllJobs()
        method, of executing the first job in the current process, and so the first
        jobID is not a jobID but the results of that. Hence we do not try to wait on
        that job, but on all the rest. 
        
        Returns only when all the listed jobID strings are no longer found in the
        SLURM queue. Currently has no time-out, although perhaps it should. 
        
        """
        allFinished = False
        
        # Extract the actual SLURM job ID strings, skipping the first element. 
        # Express as a set, for efficiency later on
        slurmJobIdSet = set([t[0] for t in jobIDlist[1:]])
        
        while not allFinished:
            squeueCmd = ["squeue", "--noheader"]
            proc = subprocess.Popen(squeueCmd, stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, universal_newlines=True)
            (stdout, stderr) = proc.communicate()
            
            stdoutLines = [line for line in stdout.split('\n') if len(line) > 0]   # No blank lines
            # Grab first word on each line, which is jobID
            squeueJobIDlist = [line.split()[0] for line in stdoutLines]
            squeueJobIDset = set(squeueJobIDlist)
            
            allFinished = slurmJobIdSet.isdisjoint(squeueJobIDset)
            
            if not allFinished:
                # Sleep for a bit before checking again
                time.sleep(60)
    
    def gatherAllOutputs(self, jobIDlist):
        """
        Gather up outputs from sub-jobs, and return a list of the
        outputs objects. Note that we assume that the first element of
        jobIDlist is actually an outputs object, from running the first sub-array
        in the current process. 
        
        The jobIDlist is a list of tuples whose second element is the name of 
        the output file containing the pickled outputs object. 
        
        """
        outputBlocksList = [jobIDlist[0]]
        for (jobID, outputsfile, logfile) in jobIDlist[1:]:
            try:
                pickledOutput = open(outputsfile, 'rb').read()
                outputObj = pickle.loads(pickledOutput)
                os.remove(outputsfile)
            except Exception as e:
                logfileContents = 'No logfile found'
                if os.path.exists(logfile):
                    logfileContents = open(logfile).read()
                msg = ("Error collecting output from SLURM sub-job. Exception message:\n"+str(e)+
                    "\nSLURM Logfile:\n"+logfileContents)
                raise rioserrors.JobMgrError(msg)
            outputBlocksList.append(outputObj)
            os.remove(logfile)
        return outputBlocksList

def find_executable(executable):
    """
    Our own version of distutils.spawn.find_executable that finds 
    the location of a script by trying all the paths in $PATH.
    Unlike distutils.spawn.find_executable, it does not add .exe
    to the script name under Windows.

    """
    path = os.environ['PATH']
    paths = path.split(os.pathsep)

    for p in paths:
        f = os.path.join(p, executable)
        if os.path.isfile(f):
            return f
    return None
    
class MpiJobManager(JobManager):
    """
    Use MPI to run individual jobs. Requires mpi4py module. 
    
    """
    jobMgrType = "mpi"
    # returned by MPI.COMM_SELF.Spawn
    comm = None
    # current destination job, so we can spread the data around
    # to all evently
    dest = 0

    def __init__(self, numSubJobs):
        from mpi4py import MPI

        # find the path to rios_subproc_mpi.py
        subproc = find_executable('rios_subproc_mpi.py')
        if subproc is None:
            msg = 'Cannot find path to rios_subproc_mpi.py'
            raise rioserrors.FileOpenError(msg)

        # base class does one job in current process so we don't
        # need to create processes for each job
        self.comm = MPI.COMM_SELF.Spawn(sys.executable, [subproc],
                            maxprocs=(numSubJobs-1))

        # call base class implementation
        JobManager.__init__(self, numSubJobs)

    def __del__(self):
        # check constructor succeeded
        if hasattr(self, 'numSubJobs'):
            # tell all the sub jobs to exit
            for dest in range(self.numSubJobs-1):
                self.comm.send([False, 0], dest=dest)

            self.comm.Disconnect()

    def startOneJob(self, userFunc, jobInfo):
        """
        Start one job. Uses the MPI send call. 
        MPI does a certain level of pickling, but 
        we do our own here so that the function etc gets pickled.

        """
        jobInfo = jobInfo.prepareForPickling()

        allInputs = (userFunc, jobInfo)
        allInputsPickled = cloudpickle.dumps(allInputs)

        # send info off to sub process
        # we also send a flag telling the subprocess
        # not to exit and be ready for another message
        self.comm.send([True, allInputsPickled], dest=self.dest)

        # return the current one
        proc = self.dest

        # set self.dest back to zero if we have done them all
        self.dest += 1
        if self.dest >= (self.numSubJobs - 1):
            self.dest = 0

        return proc

    def waitOnJobs(self, jobIDlist):
        """
        Can't actually wait on jobs with MPI, so we do nothing here.
        The waiting happens when we gatherAllOutputs() (and call MPI.recv)
        below.

        """

    def gatherAllOutputs(self, jobIDlist):
        """
        Gather up outputs from sub-jobs, and return a list of the
        outputs objects. Note that we assume that the first element of
        jobIDlist is actually an outputs object, from running the first sub-array
        in the current process. 

        """
        outputBlocks = [jobIDlist[0]]
        for job in jobIDlist[1:]:
            pickledOutput = self.comm.recv(source=job)
            outputObj = pickle.loads(pickledOutput)

            outputBlocks.append(outputObj)

        return outputBlocks

def multiUserFunc(userFunc, jobInfo):
    """
    This function is run by the MultiJobManager to run
    one job. It runs the user function and then 
    returns the output block as multiprocessing.Pool 
    expects a function to behave.
    """

    params = jobInfo.getFunctionParams()
    userFunc(*params)

    result = jobInfo.getFunctionResult(params)
    return result
    
class MultiJobManager(JobManager):
    """
    Use Python's standard multiprocessing module to run individual jobs.
    
    This JobManager sub-class should be used with caution, as it does not 
    involve any kind of load balancing, and all sub-processes simply run 
    concurrently. If you have enough spare cores and memory to do that, then
    no problem, but if not, you may clog the system. 

    This is much faster than the SubprocJobManager presumeably due to the
    custom data pickling that Python does.
    
    """
    jobMgrType = "multiprocessing"
    # an instance of multiprocessing.Pool
    pool = None

    def __init__(self, numSubJobs):
        from multiprocessing import Pool

        # base class does one job in current process so we don't
        # need to create processes for each job
        self.pool = Pool(numSubJobs - 1)

        # call base class implementation
        JobManager.__init__(self, numSubJobs)

    def __del__(self):
        # shut down the pool object as we have finished.
        self.pool.close()
        self.pool.join()

    def startOneJob(self, userFunc, jobInfo):
        """
        Start one job. Uses the multiprocessing.Pool.apply_async call. 
        This handles all the details of getting the data into the other
        process we we don't have to do any pickling here. 

        However we do call jobInfo.prepareForPickling() since
        multiprocessing has problems with the same things that the pickler
        does so we can assume the cleanup of the function here will suffice.
        
        """
        jobInfo = jobInfo.prepareForPickling()

        proc = self.pool.apply_async(multiUserFunc, (userFunc, jobInfo))
        return proc

    def waitOnJobs(self, jobIDlist):
        """
        Wait on all the jobs with AsyncResult.wait().
        The first element of jobIDlist is actually an outputs object
        so we ignore that.

        """
        for job in jobIDlist[1:]:
            job.wait(timeout=None)
            
    def gatherAllOutputs(self, jobIDlist):
        """
        Gather up outputs from sub-jobs, and return a list of the
        outputs objects. Note that we assume that the first element of
        jobIDlist is actually an outputs object, from running the first sub-array
        in the current process. 

        """
        outputBlocks = [jobIDlist[0]]
        for job in jobIDlist[1:]:
            output = job.get(timeout=None)            
            outputBlocks.append(output)

        return outputBlocks

# This mechanism for selecting which job manager sub-class to use is important in 
# order to allow an application to run without modification on different systems.
# Our own example is that JRSRP has a system which uses PBS and another which
# uses Slurm, and we want the applications to run the same on both, which means
# that there should be a way of selecting this from the environment. 
def getJobManagerClassByType(jobMgrType):
    """
    Return a sub-class of JobManager, selected by the type name
    given. 
    
    All sub-classes of JobManager will be searched for the 
    given jobMgrType string. 
        
    """
    jobMgr = None
    subClasses = JobManager.__subclasses__()
    for c in subClasses:
        if c.jobMgrType == jobMgrType:
            jobMgr = c
    return jobMgr


def getAvailableJobManagerTypes():
    """
    Return a list of currently known job manager types
    
    """
    subClasses = JobManager.__subclasses__()
    typeList = [c.jobMgrType for c in subClasses]
    return typeList


def getJobMgrObject(controls):
    """
    Take an ApplierControls object and return a JobManager sub-class 
    object which meets the needs specified in the controls object. 
    If none is required, or none is available, then return None
    
    """
    jobmgr = None
    if controls.numThreads > 1:
        if controls.jobManagerType is None:
            raise rioserrors.JobMgrError('%d threads requested, but no jobManagerType set'%controls.numThreads)
        jobMgrTypeList = getAvailableJobManagerTypes()
        if controls.jobManagerType not in jobMgrTypeList:
            raise rioserrors.JobMgrError("JobMgrType '%s' is not known"%controls.jobManagerType)
        jobmgrClass = getJobManagerClassByType(controls.jobManagerType)
        jobmgr = jobmgrClass(controls.numThreads)
        jobmgr.setTempdir(controls.tempdir)
    return jobmgr
