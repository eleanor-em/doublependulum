# doublependulum
This is a C program to generate 500 graphs of a double pendulum with arbitrary initial conditions. Bash scripts are provided to compile these images into one animated GIF.

To use the program, compile it using gcc -o pendulum pendulum.c -lm -lcpgplot, then navigate to output/ and run gen.sh. This will generate and run the Bash to resize the images, and create a GIF.
# dependencies
- Bash
- cpgplot
- ImageMagick