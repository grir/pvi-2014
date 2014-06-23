#!/usr/bin/env bash
g++ $(gsl-config --cflags) $1.cpp $(gsl-config --libs) -lm
