# Elliptical-Curve
Elliptical curves using GMP

## Functions: 
print_point: function to print a pair of x,y coordinates
scan_point: function to input a pair of x,y coordinates

add: adds two pairs of (x,y) coordinates according to the following conditions -> 
1.if any one point if Point of infinity, returns the other point
2.if two points are same, using point doubling formulas
3.if two points have different x coordinates, using other formulas
4.if two points have same x coordinates, returning point of infinity

negation: given (x,y) returns (x,-y). If point is infinity, returns infinity
subtract: calculates (x2,-y2) using negation function and uses add function to add (x1,y1) + (x2,-y2)
shanks_algo: performs the shank's algorithm to return a y coordinate for a given x
generate: checks if x^3+ax+b is a quadratic residue, if its not prints "No solution" else sends that to shanks_algo
point_multiplication: uses square and multiply version to calculate the point multiplication. Uses the add function
read_from_file: reads p,a,b from a file.The file name should be given as argument after compiling.
main: Takes choice from the user and performs that action


## Test Cases tried:
A. p=17, a=2, b=2

1.add: 5,1 and 6,3 = 10 gives 10,6
2.point mul: 5,1 20 times = 5,1
3.find y for x=10 : gives 6
4.Subtract: -1,-1 and 7,2 = 7,15

B. p=11, a=3, b=3
point doubling of 5,1 gives 4,5

C. p=5, a=4, b=4
given x=2, y=0
given x=1, y=2
given x=3, prints "No y"

