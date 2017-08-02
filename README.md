#About

This repository contains the code for SatStressGUI, a program which models stresses in icy satellites.

To use this program, download the package and run the satstressgui executable file located in Contents/MacOS.
Instructions for using the program can be found in the Help menu or by pressing f1.

Future updates to this program should either be tested in a branch before being merged into Master upon final release, or be hosted in separate repositories.

Please go to https://github.com/SatStressGUI/SatStressGUI/issues to suggest improvements or report bugs.

For more information, please contact Alex Patthoff at patthoff@jpl.nasa.gov.

SatStressGUI V5.0 has been created upon efforts by Zane Selvans, Mikael Beuthe, Jonathan Kay, Lee Tang,
Simon Kattenhorn, C.M. Cooper, Emily S. Martin,
David Dubois, Ben J. Ayton, Jessica B. Li, 
Andre Ismailyan, Peter Sinclair, Nhu Doan, and Chad Harper.

SatStressGUI V5.0 was developed at the Jet Propulsion Laboratory, California Institute of Technology
and is based on SatStressGUI.
SatStressGUI was developed by the Planetary Geology Research group at the University of Idaho.
SatStressGUI is based on SatStress, which was designed by Zane Selvans and is available at
http://code.google.com/p/satstress and most recently at https://github.com/zaneselvans/satstress.

ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any 
commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
This software may be subject to U.S. export control laws and regulations.
By accepting this document, the user agrees to comply with all applicable 
U.S. export laws and regulations. User has the responsibility to obtain export
licenses, or other export authority as may be required before exporting such
information to foreign countries or providing access to foreign persons.

Older versions of this program can be found here: https://drive.google.com/folderview?id=0B_8nH6qrvb9gRXlaY1A4SkhLcGs&usp=sharing.

A fork which rewrote the program in a Model-View-Controller Architecture can be found here: https://github.com/AndreI11/SatStressGui

=============================================================

#Getting Started

SatStressGUI is currently capable of calculating and plotting Diurnal, Non-Synchronous Rotation, Obliquity, Ice Shell Volume Change, and Polar Wander stresses and is able to simulate the formation of cycloids.

An explanation of these stresses and cycloids as well as additional references are given:

#Diurnal

Diurnal tidal stresses arise when a satellite is in an eccentric orbit. 
This is due to two reasons. First, the amplitude of the planet's gravitational force is greater at periapse than it is at apoapse. Secondly, the planet is rotating slightly faster (compared to its synchronous rotation rate) at periapse 
and slightly slower (again compared to its synchronous rotation rate) at apoapse. This results in a 'librational tide', where the planet appears to rock back and forth in the sky.

For more information on diurnal tides, please see:

    Wahr, J., Z. A. Selvans, M. E. Mullen, A. C. Barr, G. C. Collins, M. M. Selvans, and R. T. Pappalardo, Modeling stresses on satellites due to non-synchronous rotation and orbital eccentricity using gravitational potential theory, Icarus, Volume 200, Issue 1, March 2009, Pages 188-206.


#Obliquity

A satellite's obliquity (or axial tilt) is the angle between it rotational axis and its orbital axis. A satellite of zero obliquity will have a rotational axis perpendicular to its orbital plane. However, when the obliquity is nonzero, it causes the stresses due to diurnal tides and non-synchronous rotation to be asymmetric.

For more information on stresses due to oblique orbits, see:

    Jara-Orue, H. M., & Vermeersen, B. L. (2011). Effects of low-viscous layers and a non-zero obliquity on surface stresses induced by diurnal tides and non-synchronous rotation: The case of Europa. Icarus, 215(1), 417-438.


#Ice Shell Volume Change

As satellites age, they could become cooler or warmer. This would result in more of the liquid ocean freezing or melting, thereby increasing or decreasing the thickness of the icy crust. This process would force the ice shell to expand or shrink, putting extensional stress on the surface.

For more information on Ice Shell Volume Change as a stressing mechanism, please see:

    Nimmo, F. (2004). Stresses generated in cooling viscoelastic ice shells: Application to Europa. Journal of Geophysical Research: Planets (1991-2012), 109 (E12).

    Patthoff, D.A., et al. 2016. Viscoelastic modeling of tidal stresses on satellites with an enhanced SatStressGUI, 47th LPSC, abs. 1375.

    Wahr, J., et al., 2009. Modeling stresses on satellites due to nonsynchronous rotation and orbital eccentricity using gravitational potentialtheory, Icarus, 200, p. 186-206.


#Polar

Polar Wander is the apparent movement of a satellite's rotational pole due to nonsynchronous reorientation of the satellite's crust. If a satellite's crust is not coupled to its core, it may experience nonsynchronous rotation (NSR). 
Sometimes, this also results in a reorientation of the poles. The north pole appears to wander over the surface as the crust reorients itself. This results in stressing, due to the tidal bulge of the core and ocean moving beneath the crust, 
as well as the parent planet appearing to change its location in the sky. 
This stressing mechanism is calculated using an elastic model.

For more information on Polar Wander as a stressing mechanism, please see:

    Matsuyama, Isamu, and Francis Nimmo. "Tectonic patterns on reoriented and despun planetary bodies." Icarus 195, no. 1 (2008): 459-473.
    
    Matsuyama, Isamu, Francis Nimmo, and Jerry X. Mitrovica. "Planetary reorientation." Annual Review of Earth and Planetary Sciences 42 (2014): 605-634.


#Cycloids

Cycloids are arcuate lineaments found on the surface of Europa.  They are thought to be created when a fracture in the ice is propagated because of the stresses. In order for a cycloid to be created, the tensile stress at the location must exceed the tensile strength of the ice.Once the fracture has started, it will propagate through the ice at a certain velocity. This velocity could be constant, or could vary depending on the magnitude of the stress. During the cycloid's propagation, the satellite will continue orbiting around its primary. This causes the stress field on the satellite to change, making the cycloids curve. When the stress is no longer greater than the requisite propagation strength, the cycloid stops moving.
If the stress reaches the propagation strength again, it will continue.

For more information, please see:

    Hoppa, G.V., Tufts, B.R., Greenberg, R., Geissler, P.E., 1999b. Formation of cycloidal features on Europa. Science 285, 1899-1902
  
=============================================================  

#How To Use

To either plot or calculate stresses, the satellite, stresses, and grid must first be defined in their respective tabs. Satellite and Grid have buttons for loading and saving their parameters (sample .sat and .grid files are included in the sample folder). Alternatively, you may save or load the whole state of the program through File -> Export and File -> Load (You can use Load to load .sat, .grid, and .cycl files as well).

Information on the tabs and their parameters are given below:

#Satellite 

The Satellite Tab is used to input the physical properties of the satellite.

- SatStressGUI assumes a 4-layer (upper ice layer, lower ice layer, liquid ocean layer, and core) satellite body.
- Each entry should use the units denoted in the square brackets next to the box.
- The viscoelastic model used assumes that the satellite has two icy layers, a liquid ocean, and a solid core.
- The NSR period is usually on the order of 100,000 years.  If you are not using NSR, you can leave it as 'infinity'.
- The orbital eccentricity must be < 0.25. Otherwise, the program cannot reasonably calculate stresses.
- If you have changed a number but nothing seems to happen, try hitting 'Enter' in the box you changed.

#Stresses 

The Stresses Tab is used to select which stresses to use.

- For Diurnal and NSR stresses, the h2, k2, and l2 boxes should be left blank unless the user wants to input their own values. Checking the "Input Love Numbers" box will allow you to use custom Love numbers. 
When inputting custom love numbers, you must use the format <Re> +/- <Im>j.  Do not use scientific notation. 
For example, "1.2+3.0e-5j" should be written as "1.2+0.00003j."
- If left blank, the program will use .65 as the value for the Stefan Parameter in Ice Shell Volume Change. 
- Polar Wander uses an elastic, time-independent calculation, so it should probably not be used with other stresses.
- By turning on the "Assume tidally locked satellite" option, the program will calculate the tidal axis as always perpendicular to the rotational axis.
- If you turn off the tidal locking option and the plot does not update, press 'Enter' in each of the tidal axis text boxes.
- Activating the "Despinning" box allows the user to change the initial and final rotation rate of the satellite. 
The rotational period should be input in units of hours.
- All coordinates should be input as latitude and longitude; conversion to colatitude is handled by the program.

#Grid

The Grid Tab is used to specify what section of the satellite to look at in the plot.

- For more information about each stress, see the Information menu.
- NOTE: The number of latitude and longitude grid points must be equal.
- To examine the whole moon, use a latitude range from -90 to 90 and a longitude range of -180 to 180.
- Each row will only activate when the appropriate stress is enabled.
- The "orbital position" row is used to track diurnal stress from the satellite's orbit.  The satellite starts at the minimum position, and moves to the maximum position. 
Inputting 0 to 360 degrees will be one full orbit.  Additional orbits can be added by increasing the maximum beyond 360 degrees.
- The map will occasionally not work for certain positions.  If this happens, simply change the number of increments or the end position.
- The "amount of NSR build up" row is used to determine how long the ice shell has been rotating.
The Start Time is when the plotting starts, and the End Time is when the plotting ends.

#Cycloids

The Cycloids Tab allows the user to generate a cycloidal feature on the map.

- The cycloids are modeled and plotted on the Plot Tab.
- The Yield Threshold is how much stress must be put on the crust to break the ice and initiate a fracture.
- The Propagation Strength is how much stress must be put on the crust to make the split continue, and the split continues at the Propagation Speed.
- The Starting Latitude and Longitude determine where the cycloid begins, and the Direction determines the curvature of the cycloid.
- NOTE: The Vary Velocity option is currently untested.
- For more information on cycloids, see the Information menu.

#Plot

The Plot Tab shows a map of the stresses on the surface of the satellite.

- Tension on the map is shown as positive and compression is shown as negative.
- You can step through the plots by using the buttons to the bottom right of the graph.
- Each individual plot can be saved by using the save button to the lower left of the graph, and the series can be saved using the "save series" 
button to the lower right.
- The panel on the right allows manipulation of the map, changing the scale and type of map, as well as the stresses showed.
- The bottom panel enables and disables cycloids.
- When using Polar Wander, the initial and final locations of the rotational poles and/or sub- and anti-jove points will appear on the graph.
  - The initial north and south poles will be white circles.
  - The final north and south poles will be black circles.
  - The initial sub- and anti-jove points will be white squares.
  - The final sub- and anti-jove points will be black squares.
- The vectors created by Polar Wander do not currently appear to be generating correctly.
- When using cycloids, if the program is unable to initiate a cycloid, it will plot a black triangle at the attempted location.
  - If it creates a split, but cannot propagate it, it will plot a white triangle at the location.
- Cycloids can be saved as shape files via the appropriate button. Loading of shape files is currently not supported.
- NOTE: The cycloids cannot be saved as netcdf files currently.
- NOTE: The Lineaments feature does not function currently.
