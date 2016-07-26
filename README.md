This repository contains the code for SatStressGUI, a program which models stresses in icy satellites.

To use this program, download the package and run the satstressgui executable file located in Contents/MacOS.
Instructions for using the program can be found in the Help menu or by pressing f1.

Future updates to this program should either be tested in a branch before being merged into Master upon final release, or be hosted in separate repositories.

For more information, please contact Alex Patthoff at donald.a.patthoff@jpl.nasa.gov.

SatStressGUI 4.0 has been created upon efforts by
Alex Patthoff, Robert Pappalardo, Jonathan Kay, Lee Tang,
Simon Kattenhorn, C.M. Cooper, Emily S. Martin,
David Dubois, Ben J. Ayton, Jessica B. Li, 
Andre Ismailyan, Peter Sinclair.

SatStressGUI V4.0 was developed at the Jet Propulsion Laboratory, California Institute of Technology
and is based on SatStressGUI.
SatStressGUI was developed by the Planetary Geology Research group at the University of Idaho.
SatStress GUI is based on SatStress, which was designed by Zane Selvans and is available at
http://code.google.com/p/satstress and most recently at https://github.com/zaneselvans/satstress

ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any 
commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
This software may be subject to U.S. export control laws and regulations.
By accepting this document, the user agrees to comply with all applicable 
U.S. export laws and regulations. User has the responsibility to obtain export
licenses, or other export authority as may be required before exporting such
information to foreign countries or providing access to foreign persons.

Older versions of this program can be found at https://drive.google.com/folderview?id=0B_8nH6qrvb9gRXlaY1A4SkhLcGs&usp=sharing.


Getting Started
_________________

To either plot or calculate stresses, the satellite, stresses, and grid must first be defined in their respective tabs. Satellite and Grid have buttons for loading and saving their parameters (Sample .sat and .grid files are included in the sample folder).  

Information on the tabs and their parameters are given below:

Stresses -The Stresses Tab is used to select which stresses to use.

- For Diurnal and NSR stresses, the h2, k2, and l2 boxes should be left blank, unless the user wants to input their own values. 
Checking the "Input Love Numbers" box will allow you to use custom Love numbers. 
When inputting custom love numbers, you must use the format <Re> + <Im>j.  Do not use scientific notation. 
1.2 + 3e-05j would look like 1.2+0.00003j.
- The Obliquity stress must be used with Diurnal or NSR.
- The Thermal Diffusivity of the Ice Shell Thickening stress does not currently function.
- Polar Wander uses an elastic, time-independent calculation, so it should probably not be used with other stresses.
- By turning on the "Assume tidally locked satellite" option, the program will calculate the tidal axis as always perpendicular to the rotational axis.
- If you turn off the tidal locking option and the plot does not update, press 'Enter' in each of the tidal axis text boxes.
- Activating the "Despinning" box allows the user to change the initial and final rotation rate of the satellite.  
The rotational period should be input in units of hours.
- All coordinates should be input as latitude and longitude; conversion to colatitude is handled by the program.

