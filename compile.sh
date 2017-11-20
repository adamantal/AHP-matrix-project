#!/bin/bash
####################################
######## COMPILING SCRIPT ##########
####################################
if [ $# -eq 0 ]; then
    declare -a arr=("_MatrixFilter" "_LPTester" "_Histogram")
    echo "Compiling all _ files..."
    for i in "${arr[@]}"
    do
       echo g++ --std=c++11 -o2 -Wall -Werror -Wextra -o bin/$i src/$i.cpp -lemon -lglpk
       g++ --std=c++11 -o2 -Wall -Werror -Wextra -o bin/$i src/$i.cpp -lemon -lglpk
    done
else
    echo g++ --std=c++11 -o2 -Wall -Werror -Wextra -o bin/$1 src/$1.cpp -lemon -lglpk
    g++ --std=c++11 -o2 -Wall -Werror -Wextra -o bin/$1 src/$1.cpp -lemon -lglpk
fi
