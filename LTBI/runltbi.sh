#!/bin/bash
# for CNO in {1..30}
# do
#     R --slave --args < GPreg.R $CNO
# done


# version for machine with lots of cores:

R --slave --args < GPreg.R 1 & R --slave --args < GPreg.R 2 & R --slave --args < GPreg.R 3 & R --slave --args < GPreg.R 4 & R --slave --args < GPreg.R 5 & R --slave --args < GPreg.R 6 & R --slave --args < GPreg.R 7 & R --slave --args < GPreg.R 8 & R --slave --args < GPreg.R 9 & R --slave --args < GPreg.R 10


R --slave --args < GPreg.R 11 & R --slave --args < GPreg.R 12 & R --slave --args < GPreg.R 13 & R --slave --args < GPreg.R 14 & R --slave --args < GPreg.R 15 & R --slave --args < GPreg.R 16 & R --slave --args < GPreg.R 17 & R --slave --args < GPreg.R 18 & R --slave --args < GPreg.R 19 & R --slave --args < GPreg.R 20


R --slave --args < GPreg.R 21 & R --slave --args < GPreg.R 22 & R --slave --args < GPreg.R 23 & R --slave --args < GPreg.R 24 & R --slave --args < GPreg.R 25 & R --slave --args < GPreg.R 26 & R --slave --args < GPreg.R 27 & R --slave --args < GPreg.R 28 & R --slave --args < GPreg.R 29 & R --slave --args < GPreg.R 30

