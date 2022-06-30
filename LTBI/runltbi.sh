#!/bin/bash
for CNO in {1..30}
do
    R --slave --args < GPreg.R $CNO
done

