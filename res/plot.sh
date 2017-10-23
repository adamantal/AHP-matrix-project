#!/bin/bash
python -m SimpleHTTPServer &
google-chrome http://localhost:8000/plot/plot.html
