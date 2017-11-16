#!/bin/bash
####################################
######## COMPILING SCRIPT ##########
####################################

echo g++ --std=c++11 -o3 -Wall -Wextra -o bin/$1 src/$1.cpp -lemon -lglpk
#-Werror
g++ --std=c++14 -o3 -Wall -Wextra -o bin/$1 src/$1.cpp -lemon -lglpk
