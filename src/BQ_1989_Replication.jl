module BQ_1989_Replication

using XLSX
using DataFrames
using Dates
using LinearAlgebra
using Distributions
using Plots
using Statistics
using Pkg
using Test
using TestItemRunner



export run, runtests

#Import the main functions 
function run()
    include("src/files_path.jl")
    include("src/IRFs.jl") #Generates the IRFs (figures 1-2-3-4-5-6)
    #include("src/VarianceDecomp.jl") #Generates the Variance Decomposition tables
end

function runtests()
    include("test/runtests.jl")
end



end # module 
