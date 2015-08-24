using BrainWave, Base.Test

x = reshape(collect(1:27), 3,3,3)
@test isequal(deleteSlice3D(x, 2, 1), x[[1,3],:,:])
@test isequal(deleteSlice3D(x, [1,3], 2), x[:,[2], :])
@test isequal(deleteSlice3D(x, [3,1], 2), x[:,[2], :])
@test isequal(deleteSlice3D(x, [2,3], 3), x[:,:, [1]])

