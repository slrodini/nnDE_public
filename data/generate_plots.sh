#!/bin/bash

./all_plots.plt
for FILE in ./*.eps
do
    epstopdf ${FILE}
    rm ${FILE}
done