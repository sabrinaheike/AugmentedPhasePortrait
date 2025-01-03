Mac Download Instructions:

Runtime Downloading Instructions:
1. Download MATLAB runtime version R2022a Mac from https://www.mathworks.com/products/compiler/matlab-runtime.html
2. Unzip the MATLAB Runtime installer at the terminal using the unzip command. For Intel® processor-based macOS, type: unzip MATLAB_Runtime_R2023b_maci64.zip For Apple silicon-based macOS, type: unzip MATLAB_Runtime_R2023b_maca64.zip
3. Double-click the DMG file to start the installer.
4. When the MATLAB Runtime installer starts, it displays a dialog box. Read the information and then click Next to proceed with the installation.
5. In the Folder Selection dialog box, specify the folder where you want to install MATLAB Runtime. Recommendation is to install RUNTIME in the same parent folder containing the folder of the AugmentedPhasePlane application.
6. Confirm your choices and click Next.

- Download the AugmentedPhasePlane file and the run_AugmentedPhasePlane.sh file 
- Download MATLAB Runtime version 2022a (into the same parent folder containing the AugmentedPhasePlane application and run_AugmentedPhasePlane.sh file)
- For your first time time running the app
    - In your terminal, go to the folder containing the run_AugmentedPhasePlane.sh file and AugmentedPhasePlane application file.
      For example if they are in Applications in a folder called MATLAB, run the following command:
      cd /Applications/MATLAB
    - At the command prompt, run:
      ./run_testapp.sh <mcr_directory>
      where <mcr_directory> is the directory where version 9.12 of the MATLAB Runtime is installed.
      For example, if the Runtime folder is installed in ./MATLAB_Runtime/v912 run the shell script as:
      ./run_AugmentedPhasePlane.sh ./MATLAB_Runtime/v912 
- For subsequent uses you can open the application file directly from the Applications folder or from the Launchpad.
