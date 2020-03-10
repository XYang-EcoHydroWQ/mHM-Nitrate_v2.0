GFORTRAN module version '10' created from /gpfs0/home/yangx/mHM-Nitrate/github4doi/mHM-Nitrate_v2.0/src/mRM/mo_mrm_tools.f90
MD5:bc64c7468bbc6f8611c94e17f20f8b21 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

()

()

()

()

(2 'calculate_grid_properties' 'mo_mrm_tools' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 3 0 (4 5 6 7 8 9 10 11 12 13 14 15 16) () 0 () () () 0 0)
17 'get_basin_info_mrm' 'mo_mrm_tools' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 18 0 (19 20 21 22 23 24 25 26 27 28 29 30 31) () 0 () () ()
0 0)
4 'nrowsin' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
5 'ncolsin' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
6 'xllcornerin' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
7 'yllcornerin' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
8 'cellsizein' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
9 'nodata_valuein' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
10 'aimingresolution' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
11 'nrowsout' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
12 'ncolsout' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
13 'xllcornerout' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
14 'yllcornerout' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
15 'cellsizeout' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
16 'nodata_valueout' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
19 'ibasin' '' '' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
20 'ilevel' '' '' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
21 'nrows' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
22 'ncols' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
23 'ncells' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
24 'istart' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
25 'iend' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
26 'istartmask' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
27 'iendmask' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
28 'mask' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
ALLOCATABLE DIMENSION OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 ()
(2 0 DEFERRED () () () ()) 0 () () () 0 0)
29 'xllcorner' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
30 'yllcorner' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
31 'cellsize' '' '' 18 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
)

('calculate_grid_properties' 0 2 'get_basin_info_mrm' 0 17)
