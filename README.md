# Control Theory Hands-on Practical

## Hovering a table tennis ball at a fixed height

Dr. Tijs W. Alleman, Daniel Illana, Dr. Gauthier Vanhaelewyn

Correspondence: tijs.alleman@ugent.be 

BIOMATH, Department of Data Analysis and Mathematical Modelling, Ghent University, Coupure Links 653, Ghent, 9000, Belgium

<img src="./tex/fig/setup.png" alt="setup" width="300"/>

### Aim

In this practical, we turn control theory into practice by hovering a table tennis ball at a fixed height in a transparent tube. We use an Arduino Uno to build a closed loop feedback controller. We measure the ball’s height using an infrared distance sensor at the top of the tube, which we the use to control a fan at the bottom of the tube. First, we study the setup’s open-loop behavior and gather response data to calibrate a mathematical model. We then use Matlab’s Simulink environment to perform an in-silico tuning of several closed-loop feedback controllers, such as the Proportional-Integral-Derivative (PID) controller and Linear-Quadratic Regulator (LQR) state feedback controller. Finally, the performance of the tuned controllers is examined and compared using the experimental setup. This setup serves as an intuitive and visual test bed for the closed-loop feedback controllers introduced in the theoretical course. Our aim is to introduce students to the use of electronics, the unfortunate reality of noisy sensors and dead time, along with their implications on closed-loop stability, the application of control theory to non-linear systems, and to the use of mathematical modeling to support the in-silico design of a suitable controller.

### Getting started

Building and assembling the setup costs roughly 50 EUR, all relevant information is included in the practical notes.