GFORTRAN module version '10' created from /gpfs0/home/yangx/mHM-Nitrate/github4doi/mHM-Nitrate_v2.0/src/mHM/mo_mpr_soilmoist.f90
MD5:fa567b66b37fd663629435195b948644 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

()

()

()

()

(2 'mpr_sm' 'mo_mpr_soilmoist' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 3 0 (4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26) () 0 () () () 0 0)
4 'param' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '13')) 0 () () () 0 0)
5 'nodata' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 8 0 1 0 REAL ()) 0 0 () () 0 () () () 0 0)
6 'iflag_soil' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
7 'is_present' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
8 'nhorizons' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
9 'ntillhorizons' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
10 'sand' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
11 'clay' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
12 'dbm' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
13 'id0' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
14 'soilid0' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () (2 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0
0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
15 'lcover0' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 1 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
16 'thetas_till' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (3 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
17 'thetafc_till' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (3 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
18 'thetapw_till' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (3 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
19 'thetas' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
20 'thetafc' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
21 'thetapw' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
22 'ks' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (3 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
23 'db' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (3 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
24 'ksvar_h0' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
25 'ksvar_v0' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
26 'sms_fc0' '' '' 3 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 1 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
)

('mpr_sm' 0 2)
