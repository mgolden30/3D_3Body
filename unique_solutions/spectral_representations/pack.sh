#!/bin/bash

TARBALL_NAME="solutions.tar.gz"
TXT_FILE_COUNT=$(ls -1 *.mat 2>/dev/null | wc -l)

if [ "$TXT_FILE_COUNT" -eq 0 ]; then
    # Unpack the tarball if there are no txt files
    tar -zxvf "$TARBALL_NAME"
else
    # Compress all txt files into the solutions.tar.gz tarball
    tar -czvf "$TARBALL_NAME" *.mat
    # Remove the text files after compressing
    rm -f *.mat
fi
