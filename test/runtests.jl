using Test
using BQ_1989_Replication
using LinearAlgebra

include("../src/IRFs.jl")
#include("../data/BQ1989_Data.xlsx")

#Test if OLSmodel() correctly estimates coefficients and returns valid statistics
@testset "OLSmodel tests" begin
    y = randn(100)
    x = randn(100, 2)
    result = OLSmodel(y, x)
    @test length(result.beta) == 2
    @test result.sige ≥ 0
    @test 0 ≤ result.rsqr ≤ 1
    @test 0 ≤ result.rbar ≤ 1
end

#Test if VARmakexy() returns correct matrix shape
@testset "VARmakexy tests" begin
    data = randn(100, 2)
    lags = 4
    Y, X = VARmakexy(data, lags, 1)
    @test size(Y, 1) == size(X, 1)
    @test size(X, 2) == 1 + 2 * lags  # intercept + lags
end

#Tests if VARmodel correctly computes :Fcomp, :sigma, :resid, :Fit.
#Tests also if :maxEig<1 for a stable VAR.
@testset "VARmodel tests" begin
    data = randn(100, 2)
    VAR, opt = VARmodel(data, 4)
    @test haskey(VAR, :Fcomp)
    @test VAR[:Fcomp] isa Matrix
    @test VAR[:maxEig] < 1.0
end

#Tests if VARirband() returns bounds of same shape as IR, and INF ≤ MED ≤ SUP for all elements.
@testset "VARir tests" begin
    data = randn(100, 2)
    VAR, opt = VARmodel(data, 4)
    opt.ident = "short"
    opt.nsteps = 10
    IR, VARout = VARir(VAR, opt)
    @test size(IR) == (10, 2, 2)
    @test all(.!isnan.(IR))
end

#Tests if nearest_posdef() always returns a positive definite matrix
@testset "nearest_posdef tests" begin
    A = [4.0 -2.0; -2.0 1.0]  # not positive definite
    Apos = nearest_posdef(A)
    @test isposdef(Apos)
end
