# Optimization of Vehicle Speed control


*This project is based on a real life problem that is quite significant in our everyday lives. It is based on an internal combustion engine which are the most common in today’s era for transportation in light motor as well as heavy vehicles. To be precise, I will be using a high-fidelity, non-linear internal combustion engine simulator which is written in Matlab Simulink and analyze multiple methods to control the idle engine speed. </br>

*In technical terms, the primary objective is to have an optimization loop automatically “tune” the PID (closed loop) gains to achieve the best possible engine speed response in the presence of large load torque disturbances. The control input here is throttle angle (opening of the throttle). We’ll be testing the simulator in both open loop and closed loop modes. Closed loop behavior of the simulator will be investigated with the use of a PID controller. </br>

*Several optimization algorithms (constrained and unconstrained) will be implemented to find the best possible way to bring idle engine speed back up as fast as possible with the least amount of deflection.
