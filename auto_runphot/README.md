Introduction: This script is responsible for automatically generating the required "runphot" file, which is required as an input to DOLPHOT to run the required photometry. It should be placed in /ProposalID/GalaxyID/ in order to have proper access to the required files. 

How it works: First, it ensures that DOLPHOT is installed and in the path. It reads in the files from /ProposalID/GalaxyID/raw, and then finds the common prefix between the files. It then goes through and finds the filters of each file and translates them to V or I. Then it writes the required runphot file for the specific case to /ProposalID/GalaxyID/reduced/. 

Potential Troubleshooting: Currently this has the F555W, F606W, and F814W as identifiable filters. May need to add more later on.
