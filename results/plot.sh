#!/bin/bash
python3 -m http.server &
google-chrome http://localhost:8000/plot/plot.html
