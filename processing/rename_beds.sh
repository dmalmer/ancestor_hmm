#!/bin/bash

for i in $(ls ../results/ISS_*hmm-out.bed); do
    mv $i ${i%.bed}_94in.bed
done

for i in $(ls ../results/ILS_*hmm-out.bed); do
    mv $i ${i%.bed}_94in.bed
done