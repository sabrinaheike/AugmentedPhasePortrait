Mac Download Instructions:

Runtime Downloading Instructions:
1. Download MATLAB runtime version R2022a Mac from https://www.mathworks.com/products/compiler/matlab-runtime.html
2. Unzip the MATLAB Runtime installer at the terminal using the unzip command. For Intel® processor-based macOS, type: unzip MATLAB_Runtime_R2023b_maci64.zip For Apple silicon-based macOS, type: unzip MATLAB_Runtime_R2023b_maca64.zip
3. Double-click the DMG file to start the installer.
4. When the MATLAB Runtime installer starts, it displays a dialog box. Read the information and then click Next to proceed with the installation.
5. In the Folder Selection dialog box, specify the folder where you want to install MATLAB Runtime. Recommendation is to install RUNTIME in the same parent folder containing the folder of the AugmentedPhasePlane application.
6. Confirm your choices and click Next.

GUI Instructions
- Download the AugmentedPhasePlane file from this folder
- Download MATLAB Runtime version 2022a (into the same parent folder containing the AugmentedPhasePlane application)
- control + click on the AugmentedPhasePlane file and press open. This bypasses Mac's security blocks for non-verified apps. Alternatively Go to the Apple menu and select System Settings. Click Privacy & Security in the sidebar. Go to Security and click Open. Click Open Anyway. Enter your login password and click OK
- Wait for app to open in a new window.

Extra Help
- For macOS, the library path for MATLAB Runtime must be set manually. The GUI seems to work without these additional steps, but if you run into trouble here are some things you can try:
- Option 1
    1. Download the run_AugmentedPhasePlane.sh file into the same folder containing the AugmentedPhasePlane application file and MATLAB Runtime. Navigate to this folder using your terminal. For example if they are in Applications in a folder called MATLAB, run the following command: cd /Applications/MATLAB
    2. At the command prompt, run: ./run_AugmentedPhasePlane.sh <mcr_directory>, where <mcr_directory> is the directory where version 9.12 of the MATLAB Runtime is installed. For example, if the Runtime folder is installed in ./MATLAB_Runtime/v912 run the shell script using the following command: ./run_AugmentedPhasePlane.sh ./MATLAB_Runtime/v912 

- Option 2
    1. Set the environment variable DYLD_LIBRARY_PATH according to these instructions: https://www.mathworks.com/help/compiler/mcr-path-settings-for-run-time-deployment.html
    2. Run the AugmentedPhasePlane file

