# TheRocket

Notes to FAR-OUT Staff or any other user

1. Flight Simulator: To run a flight simulator, just download "TheRocket" subfolder and run "TheRocket" Jupyter notebook folder. "OUTPUT: NUMERICAL OUTPUTS" prints the critical parameters used in the update form. This simulator generates its own csv files for flow rate and engine thrust based on variable inputs within the file. Ensure your kernel is connected to a virtual environment with the following libraries (rocketpy, numpy).

2. Regen Solver: The propulsion team has generated a regen solver within "TheEngine" folder. This solver is being used to model pressure drop and heat transfer characteristics within the regen engine. While it does not generate an engine csv (the entire solver solves around a single thrust value), it is critical to size our cooling channels before printing. Refer to the ReadMe in that folder for specifics on how to run