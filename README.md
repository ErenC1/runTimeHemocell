# runTimeHemocell
Turns some of the Hemocell functions to run-time
These files are some of the changes I did in the Hemocell (https://github.com/UvaCsl/HemoCell).

Files in the "sourceFiles" turns some of the compile-time functions into run-time functions such as tauBend, and includes a scale function of RBC membrane forces. 
In hemocell.cpp only a repulsion function added. Repulsion can be activated or disabled in the config.xml (see the case files) and it can be applied every "x" iterations.

Mechanics files also includes an attempt to include Skalak model, I didn't finish it.
If you want to include that, you can use the functions as a starting point. Just uncomment the Skalak constants, and comment the current ones.

Files in the "case" is an example test case.

Step - Sinusoidal functions are NOT VALIDATED. They are achieved by changing outlet density, so larger values will probably results in crash, because its not a good way to do it.

Limiters and scale can be adjusted in RBC.xml
RDW (red cell distribution width) and PLT can be activated in config.xml, currently supports 3 different diameters, you need to create corresponding xml and pos files (such as LRBC.xml, PLT.xml).


Liked HemoCell and used these functions? Please Cite the following publications:

Závodszky, G., Van Rooij, B., Azizi, V., & Hoekstra, A. (2017). Cellular level in-silico modeling of blood rheology with an improved material model for red blood cells. Frontiers in physiology, 8, 563. https://doi.org/10.3389/fphys.2017.00563

Çolak, E., Ekici, Ö., & Erdener, Ş. E. (2025). In Silico Investigation of the RBC Velocity Fluctuations in Ex Vivo Capillaries. Applied Sciences, 15(14), 7796. https://doi.org/10.3390/app15147796
