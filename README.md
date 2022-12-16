# BaN-BoS

Basic N body simluation for Scientific computing module.

-----Info for markers------

For any markers the way you should use this is to run it, input "1" into the console. Then wait from 20s - 1.5mins for a graph of the body positions over time. After you're done with this close all matplotlib windows and enter anything into the console for the energy graph. Note precision of the simulation is traded for speed when you do this,

If you want to check out the basic case then type 3 instead of 1 in the initial prompt. This may take some time longer than the fast mode.

----------Modes---------
1 - Fast mode simulates the most important part of the solar system for a simulation at a less precise time interval. The sun, jupiter and the interloper. The interloper has a mass of about 1 solar mass by default.
2 - Normal mode is what the whole solar system interloper model runs with, a typical runtime is about 2 mins for my machine but could be significantly longer especially with more timesteps.
3 - Basic mode shows the 2-Body system test case to make sure you can check the code still works after tinkering

----Here be dragons------
If you want to adapt this code into anything feel free but I wouldnt reccomend it. There are likely far more efficient open source progeams with far less spaghetti code than this.
Be careful with the plotting routines, if you mess up the loops then it can easily crash your pc.