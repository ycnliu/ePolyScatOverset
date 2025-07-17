---
title: Profiling Using VTune
author: Yuchen Liu
date: today
---

# Logging into Lawrencium

We need to activate an X server on our login node so add the ``-Y`` flag to you
ssh command.

This is not the ideal way to profile an application but it's a quick way to get
started. To set up a better work flow, checkout the link below...

[Set Up Remote Linux Target](https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/set-up-analysis-target/linux-targets/remote-linux-target-setup.html)

# Environment Setup

1. Load only this Intel module and do not load ``openmpi``.

  ```shell
  module load intel/2019.4.par_vtune
  ```

2. To have symbols and other information for better analysis, the ``-g`` flag
  is essential and the ``-qopenmp``, ``-debug``, and ``inline-debug-info`` flags
  are recommended by Intel for more in-depth analysis. Look at
  ``include/lrc_zsintelmkl.mk`` as an example for a setup used to run test63.

**Note**: As of 2021-08-12, this is the only module that has VTune and runs ePolyScat.
However, this module does not contain OpenMPI distribution but it does contain
the Intel MPI distribution and requires the ``mpiifort`` compiler. This must be
changed in the ``include/$MACH$COMPILER.mk`` file's ``FORTMPI`` variable.

**Note**: The default 2016 Intel module provides a lot of other modules that are
not included in this version. So you will not have access to these modules
using this Intel module

# Running VTune

1. Assuming the Intel module above is loaded, you can start the GUI

```shell
amplxe-gui
```

2. Then you need to create a **project** which just creates a directory in your home
  directory for VTune to put everything into.

3. Then point VTune to your script or application by providing the full path.

4. Finally run either *HotSpots* or *Threading* profiling.

**Note**:You should see the GUI pop up on your screen after step-1.  If not, be sure
you logged onto a node with X server running (i.e. included the ``Y`` flag).

# Test 63 Results

The major HotSpot is in ``InterpolationOverset.f90`` in the subroutine
``Determineoverlapping`` on line 869.  It was found 25.5% of CPU time is here
alone.  This seems to be because there is an all-to-all Bcast call
within a triply-nested loop. Threading profiling showed this particular call was
causing a drastic load imbalance. Pulling out the all-to-all collective communication
outside of the loops should resolve these issues.

Another part of the code found to have a long "Wait Time" is also in
``InterpolationOverset.f90`` but on line 322, which is a ``MPI_RECV`` call. Only
this particular ``MPI_RECV`` call so looks like waiting on this particular data.

# Batch Job Profiling

The above example starts the profiler by clicking the run button in the GUI.
The GUI will run an application or a script.  However, the GUI that is started
on your login node cannot run the SLURM job manager and therefore this profiling
is performed on the login node, which is not ideal.  One benefit of this is
the profiling results are obtained much faster on the login node.  However,
running multiple MPI jobs on the login node is not advised and you may get a
stern email requesting you stop.  **So know that running the profiling inside
the GUI will actually use the login node**.

A better way to profile an application is to run the profiling using VTune's
command line interface and submit a batch job to the SLURM scheduler as usual.
A tutorial on how to do this can be found below.  Though this presentation is
for KNL node on Cori, this gives you a good idea of how to use the application.
Note, there are some differences with the Cori and Lawrencium setup batch
script (e.g. on Cori there is a #SBATCH --perf=vtune option but not on
Lawrencium).

[Cori KNL node example](https://www.nersc.gov/assets/Uploads/04a-vtune-knl-20170609.pdf)

The ``profile.sh`` script found in this directory has an example of using the
``amplxe-cl`` command to run.  Following the Cori KNL node example, the
``srun`` command is used and that is essentially an ``mpiexec`` command.
Therefore, this test run is using the base ``ePolyScat.exe`` executable that
is not wrapped with the MPI meta data. This is to avoid possibly calling a MPI
job spawning command like ``mpiexec`` multiple times.

At the time of writing this, the job is still running and is far longer than
using the login node...  So long, this maybe wrong.
