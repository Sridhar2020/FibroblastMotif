#!/bin/sh

name=$1
value=$2
filename=TP06_OpSplit_2D.h

if [ x$1 = x ]; then
    echo Usage: $0 parameter_name [initial_value]
    echo Set or get initial_value of parameter_name
    exit 1
fi

if [ x$2 != x ]; then
    sed -i 's/\(^#define\s\+'$name'\s\+\).*/\1'$value'/' $filename
fi

grep -E '^#define[[:space:]]+'$name'[[:space:]]*.*' $filename
