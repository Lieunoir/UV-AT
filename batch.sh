#!/bin/bash

cmake --build build --target at-uv -- -j 16
search_dir=$1
for entry in "$search_dir"/*.obj
do
    echo "$entry"...
    ./build/bin/at-uv -i "$entry" --withoutGUI
done
