#!/bin/bash

a=$(ls ~/programming/web/r_image/plots/| grep -E "^plot[0-9]{5,5}\.$1$" | sort -s | tail -1 | sed -E "s/^plot([0-9]{5,5})\.$1/\1/g")
echo $a
