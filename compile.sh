#!/bin/bash
####################################
######## COMPILING SCRIPT ##########
####################################

echo g++ --std=c++14 -o3 -Wall -Wextra -Werror -o bin/$1 src/$1.cpp -lemon -lglpk
g++ --std=c++14 -o3 -Wall -Wextra -Werror -o bin/$1 src/$1.cpp -lemon -lglpk
