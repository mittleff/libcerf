// Created by ./arrange_polyominoes.py 20 d30N20b16.dat on 13:22:22.247172
static const int Cover[4096] = {
 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
 8, 8, 1, 1,21, 2, 2, 2, 2,22, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
 8, 8,21,21,21,23, 2, 2,22,22,22, 3, 3, 3, 3,24,24, 4, 4, 4, 4, 4, 4,27, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,32, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,39,39, 7, 7, 7, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
 9, 9,25,25,23,23,23,22,22,22,22,22,22,26,24,24,24,24, 4, 4, 4, 4,27,27,27, 5, 5, 5, 5, 5, 5, 5, 5,32,32, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,39,39, 7, 7, 7, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
10,28,25,25,23,23,23,29,22,22,22,22,26,24,24,24,24,24,24,31,31,27,27,27,27,27, 5, 5, 5, 5, 5, 5,32,32,32,32, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,39,39,39,39, 7, 7, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
10,28,28,30,30,23,29,29,29,22,22,26,26,24,24,24,24,24,24,31,31,27,27,27,27,27,27,27, 5, 5,32,32,32,32,32,32,32, 6, 6, 6, 6, 6, 6, 6, 6,39,39,39,39,39,39, 7, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
11,11,33,30,30,34,29,29,29,29,29,26,26,24,24,24,24,24,24,31,27,27,27,27,27,27,27,27,27,27,32,32,32,32,32,32,32,32, 6, 6, 6, 6, 6, 6,39,39,39,39,39,39,39,39, 7, 7, 7, 7,-1,-1,-1,-1,-1,-1,-1,-1,
11,11,33,30,30,34,34,29,29,29,35,35,26,26,24,24,24,24,31,31,31,27,27,27,27,27,27,27,27,32,32,32,32,32,32,32,32,32,32,32,32,39,39,39,39,39,39,39,39,39,39,39,39,39,39,52,-1,-1,-1,-1,-1,-1,-1,-1,
12,36,33,33,34,34,34,37,37,35,35,35,35,35,38,24,24,31,31,31,31,27,27,27,27,27,27,27,27,32,32,32,32,32,32,32,32,32,32,32,32,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,-1,-1,-1,-1,-1,-1,-1,-1,
12,36,36,36,40,34,37,37,37,35,35,35,35,35,38,38,38,31,31,31,31,31,27,27,27,27,27,27,43,32,32,32,32,32,32,32,32,32,32,32,32,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,-1,-1,-1,-1,-1,-1,-1,-1,
13,41,41,40,40,40,37,37,37,35,35,35,35,35,38,38,38,38,31,31,31,31,44,27,27,27,27,43,43,43,32,32,32,32,32,32,32,32,32,32,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,-1,-1,-1,-1,-1,-1,-1,-1,
13,13,41,40,40,40,42,37,37,42,35,35,35,38,38,38,38,38,38,44,44,44,44,44,44,48,43,43,43,43,32,32,32,32,32,32,32,32,32,32,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,-1,-1,-1,-1,-1,-1,-1,-1,-1,
13,13,41,45,40,46,42,42,42,42,47,47,47,38,38,38,38,38,38,44,44,44,44,44,44,43,43,43,43,43,43,32,32,32,32,32,32,32,32,53,53,39,39,39,39,39,39,39,39,39,39,39,39,39,39,-1,-1,-1,-1,-1,-1,-1,-1,-1,
13,45,45,45,45,46,42,42,42,42,47,47,47,47,38,38,38,38,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,32,32,32,32,53,53,53,53,39,39,39,39,39,39,39,39,39,39,39,39,39,39,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,49,45,45,46,46,46,42,42,47,47,47,47,47,50,50,50,50,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,43,43,53,53,53,53,53,53,39,39,39,39,39,39,39,39,39,39,39,39,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,14,49,49,51,46,46,51,51,47,47,47,47,47,50,50,50,50,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,43,43,53,53,53,53,53,53,53,39,39,39,39,39,39,39,39,39,39,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,14,49,49,51,51,51,51,51,47,47,47,47,47,50,50,50,50,44,44,44,44,44,44,44,43,43,43,43,43,43,43,43,43,43,43,53,53,53,53,53,53,53,52,39,39,39,39,39,39,39,39,52,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
14,49,49,49,54,51,51,51,54,55,47,47,47,50,50,50,50,50,50,44,44,44,44,44,44,48,43,43,43,43,43,43,43,43,43,53,53,53,53,53,53,53,53,53,52,52,39,39,39,39,52,52,52,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,56,56,54,54,54,54,54,55,55,55,55,50,50,50,50,50,50,57,57,44,44,48,48,48,48,43,43,43,43,43,43,43,53,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,52,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,15,56,54,54,54,54,54,55,55,55,55,55,50,50,50,50,57,57,57,57,48,48,48,48,48,48,43,43,43,43,43,63,53,53,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,15,56,54,54,54,54,54,55,55,55,55,55,58,58,59,57,57,57,57,57,57,57,48,48,48,62,62,62,62,62,62,62,53,53,53,53,53,53,53,53,53,53,53,53,52,52,52,52,52,52,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
15,15,56,56,56,60,54,60,55,55,55,55,55,55,58,58,57,57,57,57,57,57,57,57,57,57,62,62,62,62,62,62,62,62,62,53,53,53,53,53,53,53,53,53,53,53,53,67,52,52,52,52,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,61,61,56,60,60,60,60,60,55,55,55,55,58,58,58,57,57,57,57,57,57,57,57,57,57,62,62,62,62,62,62,62,62,62,53,53,53,53,53,53,53,53,53,53,53,53,67,67,67,67,71,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,61,60,60,60,60,60,60,64,64,58,58,58,58,58,57,57,57,57,57,57,57,57,57,57,62,62,62,62,62,62,62,62,62,62,53,53,53,53,53,53,53,53,53,53,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,16,60,60,60,60,60,60,64,64,58,58,58,58,58,58,57,57,57,57,57,57,57,57,62,62,62,62,62,62,62,62,62,62,62,63,53,53,53,53,53,53,53,53,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,16,61,60,60,60,60,64,64,64,64,58,58,58,58,58,58,57,57,57,57,57,57,66,62,62,62,62,62,62,62,62,62,62,62,63,63,63,53,53,53,53,67,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,16,61,61,61,65,65,64,64,64,64,64,58,58,58,58,58,58,59,57,57,57,57,66,66,62,62,62,62,62,62,62,62,62,62,62,63,63,63,63,63,63,63,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
16,68,68,68,65,65,65,64,64,64,64,64,64,64,58,58,59,59,59,59,66,66,66,66,66,66,62,62,62,62,62,62,62,62,62,63,63,63,63,63,63,63,67,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,68,68,65,65,65,65,64,64,64,64,64,64,69,70,59,59,70,66,66,66,66,66,66,66,62,62,62,62,62,62,62,62,62,63,63,63,63,63,63,63,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,65,65,65,65,65,65,64,64,64,64,69,69,69,69,70,70,66,66,66,66,66,66,66,66,62,62,62,62,62,62,62,63,63,63,63,63,63,63,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,68,65,65,65,65,65,65,72,69,69,69,69,69,69,69,66,66,66,66,66,66,66,66,66,66,66,62,62,62,63,63,63,63,63,63,63,63,67,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,68,65,65,65,65,65,65,69,69,69,69,69,69,69,69,69,66,66,66,66,66,66,66,66,66,66,75,75,75,75,63,63,63,63,63,63,67,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,17,68,73,65,65,65,65,72,69,69,69,69,69,69,69,69,69,66,66,66,66,66,66,66,66,66,66,75,75,75,75,75,75,76,76,76,76,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
17,17,74,73,73,73,72,72,72,72,69,69,69,69,69,69,69,69,69,70,66,66,66,66,66,66,66,66,75,75,75,75,75,75,75,75,76,76,76,67,67,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,73,73,73,73,72,72,72,72,69,69,69,69,69,69,69,69,69,70,66,66,66,66,66,66,66,66,75,75,75,75,75,75,75,75,75,76,76,76,67,67,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,73,73,73,72,72,72,72,72,69,69,69,69,69,69,69,70,70,70,70,66,66,66,66,75,75,75,75,75,75,75,75,75,75,75,76,76,76,76,67,67,67,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,73,73,72,72,72,72,72,72,69,69,69,69,69,70,70,70,77,77,77,77,75,75,75,75,75,75,75,75,75,75,75,75,75,75,76,76,76,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,73,73,73,72,72,72,72,72,72,78,69,79,79,79,77,77,77,77,77,77,75,75,75,75,75,75,75,75,75,75,75,75,75,75,76,76,76,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,73,73,73,73,73,72,72,78,78,78,78,78,78,77,77,77,77,77,77,77,75,75,75,75,75,75,75,75,75,75,75,75,75,75,76,76,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,18,73,73,73,73,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,75,75,75,75,75,75,75,75,75,75,75,75,75,75,76,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,18,74,74,80,80,80,80,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,75,75,75,75,75,75,75,75,75,75,75,75,76,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
18,18,81,81,80,80,80,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,75,75,75,75,75,75,75,75,75,75,75,75,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,81,80,80,80,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,75,75,75,75,75,75,75,75,75,75,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,80,80,80,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,77,75,75,75,75,75,75,75,75,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,80,80,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,77,77,82,75,75,75,75,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,80,80,80,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,82,82,83,83,83,76,76,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,19,80,80,80,78,78,78,78,78,78,78,78,77,77,77,77,77,77,77,77,77,77,82,82,82,83,83,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,19,80,80,80,80,78,78,78,78,78,78,78,79,77,77,77,77,77,77,77,77,82,82,82,82,83,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,80,80,80,80,80,80,80,84,78,84,79,79,79,82,82,77,77,77,77,82,82,82,82,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,19,81,81,80,80,80,84,84,84,84,84,84,82,82,82,82,82,82,82,82,82,82,82,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,19,20,81,81,84,84,84,84,84,84,84,84,84,82,82,82,82,82,82,82,82,82,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
19,19,19,20,20,20,84,84,84,84,84,84,84,84,84,84,84,82,82,82,82,82,82,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,84,84,84,84,84,84,84,84,84,84,84,82,82,82,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,20,84,84,84,84,84,84,84,84,84,84,84,82,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,20,84,84,84,84,84,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
20,20,20,20,20,20,20,84,84,84,84,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
};
