# AMUSE-Asteroids

Claire Olde Loohuis, Kutay Nazli, and Tünde Meijer

#### Scientific Question: Can we detect, from any point in the solar system, the light curve and orbital variation of an asteroid that is influenced by the YORP and Yarkovsky effects?

For information about the project scope/progression, see ExpectedProgression.pdf, and for an overview of all functions/classes currently being drafted, see CodeOverview.py that defines all code inputs/outputs.

The script descriptions are as follows:
1. cubesphere.py initializes the shape of the asteroid and tesselates the faces, so that each patch can be assigned its own normal direction, albedo, and emissivity.
2. asteroid.py assigns the tesselations to the chosen asteroid shape; for our purposes, a sphere with user-defined radius. Each patch is assigned an albedo and emissivity, and its normal direction is calculated. This code also contains the functions that describe the effects of YORP and Yarkovsky forces on the asteroid, as well as the framework to calculate and store the flux reflected/re-emitted from the asteroid as seen by a solar system observer.
3. system.py creates and houses the solar system, including the Sun, Earth, Jupiter, and input asteroids. It calls the flux and acceleration codes from asteroid.py, and stores them in a function that can later be called by the bridge to update the positions of all celestial bodies.
4. bridge.py evolves the system by bridging two n-body codes (the asteroids and the rest of the solar system) with the YORP/Yarkovsky field code outlined in asteroid.py. For the sake of efficiency, the YORP/Yarkovsky code affects the asteroids, but no other solar system objects as their masses are too large. Variables calculated at each timestep are stored in the system class.
5. main.py contains the plotting structure. It produces plots of the position of all solar system bodies using the chosen timestep, as well as the evolutions of the asteroid's semi-major axis and eccentricity and the flux seen by the observer. 

The simulations and plots presented can be recreated running main.py (currently the animated plotting lines are commented out).  

For a simple or less time-consuming example of what the result will look like, consider running example.py.

The Overleaf for the project proposal and report can be found at https://www.overleaf.com/2411473488tkzrntnjvjdc.
