using BrainWave, Base.Test

x = [1 2 3 4;
     5 6 7 8;
     9 10 11 12;
     13 14 15 16
     ]

@test isequal(deleteSlice2D(x, 1, 1), x[2:end,:])
@test isequal(deleteSlice2D(x, [1,4], 1), x[2:3,:])
@test isequal(deleteSlice2D(x, [4,1], 1), x[2:3,:])

@test isequal(deleteSlice2D(x, 1, 2), x[:,2:end])
@test isequal(deleteSlice2D(x, [1,4], 2), x[:,2:3])
@test isequal(deleteSlice2D(x, [4,1], 2), x[:,2:3])
