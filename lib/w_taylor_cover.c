// Created by ./arrange_polyominoes.py 20 d30N20b16.dat on 13:45:35.373420
static const int Cover[4096] = {
-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
 7,21, 0, 0,22, 1, 1, 1, 1,23, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
 7,21,21,22,22,22, 1, 1,23,23,23, 2, 2, 2, 2,24,24, 3, 3, 3, 3, 3, 3,28, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,33, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,41,41, 6, 6, 6, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
 8,25,25,26,26,26,26,23,23,23,23,23,23,27,24,24,24,24, 3, 3, 3, 3,28,28,28, 4, 4, 4, 4, 4, 4, 4, 4,33,33, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,41,41, 6, 6, 6, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
 8,25,25,29,26,26,30,30,23,23,23,23,27,24,24,24,24,24,24,32,32,28,28,28,28,28, 4, 4, 4, 4, 4, 4,33,33,33,33, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,41,41,41,41, 6, 6, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
 9,31,31,29,29,30,30,30,30,23,23,27,27,24,24,24,24,24,24,32,32,28,28,28,28,28,28,28, 4, 4,33,33,33,33,33,33,33, 5, 5, 5, 5, 5, 5, 5, 5,41,41,41,41,41,41, 6, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
 9,34,34,29,29,35,30,30,36,36,36,27,27,24,24,24,24,24,24,32,28,28,28,28,28,28,28,28,28,28,33,33,33,33,33,33,33,33, 5, 5, 5, 5, 5, 5,41,41,41,41,41,41,41,41, 6, 6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,
10,34,34,37,35,35,35,38,36,36,36,36,27,27,24,24,24,24,32,32,32,28,28,28,28,28,28,28,28,33,33,33,33,33,33,33,33,33,33,33,33,41,41,41,41,41,41,41,41,41,41,41,41,41,41,56,-1,-1,-1,-1,-1,-1,-1,-1,
10,39,39,37,35,35,35,38,36,36,36,36,40,40,40,24,24,32,32,32,32,28,28,28,28,28,28,28,28,33,33,33,33,33,33,33,33,33,33,33,33,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,-1,-1,-1,-1,-1,-1,-1,-1,
11,39,39,37,37,35,38,38,38,36,36,40,40,40,40,40,40,32,32,32,32,32,28,28,28,28,28,28,47,33,33,33,33,33,33,33,33,33,33,33,33,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,-1,-1,-1,-1,-1,-1,-1,-1,
11,42,42,42,43,43,38,38,38,44,45,40,40,40,40,40,40,46,32,32,32,32,46,28,28,28,28,47,47,47,33,33,33,33,33,33,33,33,33,33,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,-1,-1,-1,-1,-1,-1,-1,-1,
12,12,42,42,43,43,43,44,44,44,44,40,40,40,40,40,40,46,46,46,46,46,46,46,46,51,47,47,47,47,33,33,33,33,33,33,33,33,33,33,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,-1,-1,-1,-1,-1,-1,-1,-1,-1,
12,12,48,48,43,43,49,44,44,44,44,45,40,40,40,40,50,46,46,46,46,46,46,46,46,47,47,47,47,47,47,33,33,33,33,33,33,33,33,57,57,41,41,41,41,41,41,41,41,41,41,41,41,41,41,-1,-1,-1,-1,-1,-1,-1,-1,-1,
13,13,48,48,48,49,49,44,44,44,44,45,45,40,40,50,50,46,46,46,46,46,46,46,46,47,47,47,47,47,47,47,47,33,33,33,33,57,57,57,57,41,41,41,41,41,41,41,41,41,41,41,41,41,41,-1,-1,-1,-1,-1,-1,-1,-1,-1,
13,13,48,48,52,49,49,49,44,44,45,45,45,50,50,50,50,50,46,46,46,46,46,46,51,47,47,47,47,47,47,47,47,47,47,47,57,57,57,57,57,57,41,41,41,41,41,41,41,41,41,41,41,41,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,53,52,52,52,49,49,49,49,54,54,54,55,50,50,50,50,50,50,46,46,46,46,51,51,47,47,47,47,47,47,47,47,47,47,47,57,57,57,57,57,57,57,41,41,41,41,41,41,41,41,41,41,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,14,52,52,52,52,49,49,54,54,54,54,54,50,50,50,50,50,50,50,50,51,51,51,51,47,47,47,47,47,47,47,47,47,47,47,57,57,57,57,57,57,57,56,41,41,41,41,41,41,41,41,56,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,14,53,52,52,52,58,54,54,54,54,54,54,54,50,50,50,50,50,50,59,51,51,51,51,51,47,47,47,47,47,47,47,47,47,57,57,57,57,57,57,57,57,57,56,56,41,41,41,41,56,56,56,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,53,53,53,58,58,58,54,54,54,54,54,54,54,55,50,50,50,50,59,59,59,51,51,51,51,51,47,47,47,47,47,47,47,57,57,57,57,57,57,57,57,57,57,56,56,56,56,56,56,56,56,56,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,60,60,58,58,58,58,54,54,54,54,54,55,55,55,55,61,59,59,59,59,51,51,51,51,51,51,47,47,47,47,47,64,57,57,57,57,57,57,57,57,57,57,57,56,56,56,56,56,56,56,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,60,60,58,58,58,58,58,54,54,54,55,55,55,55,55,59,59,59,59,59,59,59,51,51,51,63,63,63,63,63,63,63,57,57,57,57,57,57,57,57,57,57,57,57,56,56,56,56,56,56,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,60,60,60,58,58,58,62,62,62,62,55,55,55,55,61,59,59,59,59,59,59,59,59,59,63,63,63,63,63,63,63,63,63,57,57,57,57,57,57,57,57,57,57,57,57,69,56,56,56,56,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,60,60,60,60,65,62,62,62,62,62,62,61,61,61,61,59,59,59,59,59,59,59,59,59,63,63,63,63,63,63,63,63,63,57,57,57,57,57,57,57,57,57,57,57,57,69,69,69,69,72,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,66,60,60,60,60,65,62,62,62,62,62,62,61,61,61,61,61,59,59,59,59,59,59,59,63,63,63,63,63,63,63,63,63,63,63,57,57,57,57,57,57,57,57,57,57,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,66,66,65,65,65,62,62,62,62,62,62,61,61,61,61,61,59,59,59,59,59,59,59,63,63,63,63,63,63,63,63,63,63,63,64,57,57,57,57,57,57,57,57,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,16,65,65,65,65,62,62,62,62,62,62,67,61,61,61,61,61,68,59,59,59,68,68,63,63,63,63,63,63,63,63,63,63,63,64,64,64,57,57,57,57,69,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,16,66,65,65,65,65,62,62,62,62,67,67,67,61,61,61,68,68,68,68,68,68,68,63,63,63,63,63,63,63,63,63,63,63,64,64,64,64,64,64,64,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,66,66,65,65,65,65,65,70,67,67,67,67,67,67,67,67,68,68,68,68,68,68,68,68,63,63,63,63,63,63,63,63,63,64,64,64,64,64,64,64,69,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,66,66,66,71,70,65,70,70,70,67,67,67,67,67,67,67,67,68,68,68,68,68,68,68,68,63,63,63,63,63,63,63,63,63,64,64,64,64,64,64,64,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,73,71,70,70,70,70,70,70,67,67,67,67,67,67,67,67,68,68,68,68,68,68,68,68,68,63,63,63,63,63,63,63,64,64,64,64,64,64,64,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,71,70,70,70,70,70,70,67,67,67,67,67,67,67,67,68,68,68,68,68,68,68,68,68,74,74,63,63,63,64,64,64,64,64,64,64,64,69,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,71,70,70,70,70,70,70,67,67,67,67,67,67,67,67,68,68,68,68,68,68,68,68,74,74,74,74,74,74,74,64,64,64,64,64,64,69,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,71,71,70,70,70,70,70,70,67,67,67,67,67,67,75,68,68,68,68,68,68,68,68,74,74,74,74,74,74,74,74,74,74,74,80,80,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,71,71,71,70,70,70,70,76,77,77,67,67,75,75,75,75,68,68,68,68,68,68,74,74,74,74,74,74,74,74,74,74,74,74,80,80,69,69,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,73,73,71,71,71,71,76,76,76,76,76,77,75,75,75,75,75,75,75,75,75,79,74,74,74,74,74,74,74,74,74,74,74,74,74,74,80,80,69,69,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,73,73,78,78,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,79,74,74,74,74,74,74,74,74,74,74,74,74,74,74,80,80,80,69,69,69,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,78,78,78,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,75,75,74,74,74,74,74,74,74,74,74,74,74,74,74,74,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,78,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,75,75,79,74,74,74,74,74,74,74,74,74,74,74,74,80,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,78,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,75,75,79,74,74,74,74,74,74,74,74,74,74,74,74,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,78,78,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,79,79,79,74,74,74,74,74,74,74,74,74,74,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,78,78,76,76,76,76,76,76,76,76,75,75,75,75,75,75,75,75,75,79,79,79,79,74,74,74,74,74,74,74,74,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,78,78,78,78,76,76,76,76,76,76,77,77,75,75,75,75,75,75,75,79,79,79,79,79,79,74,74,74,74,74,74,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,78,78,78,78,78,78,81,76,76,77,77,82,82,82,82,75,75,75,79,79,79,79,79,79,79,79,79,79,79,80,80,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,83,78,78,83,81,81,81,81,81,81,81,82,82,82,82,82,82,79,79,79,79,79,79,79,79,79,79,80,80,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,83,83,81,81,81,81,81,81,81,81,81,82,82,82,82,82,79,79,79,79,79,79,79,79,79,79,84,80,80,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,81,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,79,79,79,79,79,79,79,79,84,84,84,80,80,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,81,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,79,79,79,79,84,84,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,19,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,82,84,84,84,84,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,19,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,82,84,84,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,81,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,84,84,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,81,81,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,20,20,81,81,81,81,81,81,81,81,81,82,82,82,82,82,82,82,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,20,20,20,83,81,81,81,81,81,81,81,82,82,82,82,82,82,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,20,83,81,81,81,81,81,85,85,85,85,85,85,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,20,85,85,85,85,85,85,85,85,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,20,85,85,85,85,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
};
