t1 = time_ns() #tic()

include("test_baseline_correct.jl")
include("test_workflow.jl")
include("test_delete_slice_2D.jl")
include("test_delete_slice_3D.jl")
include("test_filterContinuous_DArray.jl")
include("test_filterContinuous_SharedArray.jl")

include("run_examples.jl")
include("run_spectrum.jl")
t2 = time_ns() #toc()
println(float(t2-t1)*1e9)
