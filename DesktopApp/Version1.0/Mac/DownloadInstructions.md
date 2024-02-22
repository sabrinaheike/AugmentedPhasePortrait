Mac Downloading Instructions

Download the AugmentedPhasePlane file
Dowload MATLAB Runtime
  Shortened Instructions:
  1. Download MATLAB runtime version R2023a Mac from  
      https://www.mathworks.com/products/compiler/matlab-runtime.html
  2. Unzip the MATLAB Runtime installer at the terminal using the unzip command.
      For Intel® processor-based macOS, type: unzip MATLAB_Runtime_R2023b_maci64.zip
      For Apple silicon-based macOS, type: unzip MATLAB_Runtime_R2023b_maca64.zip
  3. Double-click the DMG file to start the installer.
  4. When the MATLAB Runtime installer starts, it displays a dialog box. Read the                 information and then click Next to proceed with the installation.
  5. In the Folder Selection dialog box, specify the folder where you want to install MATLAB       Runtime. Recommendation is to install RUNTIME in the same parent folder containing the       folder of the AugmentedPhasePlane application.
  6. Confirm your choices and click Next.
Launch application from finder

More Detailed Instructions for MATLAB Runtime are below if needed:
   
AugmentedPhasePlane Executable

1. Prerequisites for Deployment 

Verify that MATLAB Runtime(R2023a) is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.
NOTE: You will need administrator rights to run the MATLAB Runtime installer. 

Alternatively, download and install the Macintosh version of the MATLAB Runtime for R2023a 
from the following link on the MathWorks website:

    https://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
"Distribute Applications" in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.

2. Files to Deploy and Package

Files to Package for Standalone 
================================
-run_AugmentedPhasePlane.sh (shell script for temporarily setting environment variables 
                             and executing the application)
   -to run the shell script, type
   
       ./run_AugmentedPhasePlane.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where MATLAB Runtime(R2023a) is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have MATLAB Runtime(R2023a) installed in 
    /mathworks/home/application/R2023a, run the shell script as:
    
       ./run_AugmentedPhasePlane.sh /mathworks/home/application/R2023a
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_AugmentedPhasePlane.sh /mathworks/devel/application/matlab
-MCRInstaller.zip 
    Note: if end users are unable to download the MATLAB Runtime using the
    instructions in the previous section, include it when building your 
    component by clicking the "Runtime included in package" link in the
    Deployment Tool.
-The Macintosh bundle directory structure AugmentedPhasePlane.app 
    Note: this can be stored in an archive file with the zip command 
    zip -r AugmentedPhasePlane.zip AugmentedPhasePlane.app
    or the tar command 
    tar -cvf AugmentedPhasePlane.tar AugmentedPhasePlane.app
-This readme file 



3. Definitions

For information on deployment terminology, go to
https://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.

4. Appendix 

A. Mac systems:
In the following directions, replace MR/R2023a by the directory on the target machine 
   where MATLAB is installed, or MR by the directory where the MATLAB Runtime is 
   installed.

If the environment variable DYLD_LIBRARY_PATH is undefined, set it to the following 
   string:

MR/R2023a/runtime/maci64:MR/R2023a/sys/os/maci64:MR/R2023a/bin/maci64

If it is defined, set it to the following:

${DYLD_LIBRARY_PATH}:MR/R2023a/runtime/maci64:MR/R2023a/sys/os/maci64:MR/R2023a/bin/maci64

    For more detailed information about setting the MATLAB Runtime paths, see Package and 
   Distribute in the MATLAB Compiler documentation in the MathWorks Documentation Center.


     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, you can
              run the shell script file run_AugmentedPhasePlane.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    



5. Launching application using Macintosh finder

If the application is purely graphical, that is, it doesn't read from standard in or 
write to standard out or standard error, it may be launched in the finder just like any 
other Macintosh application.


For more information about installing and configuring MATLAB Runtime visit: https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html#bvf6b29 
