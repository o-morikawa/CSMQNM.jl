using Test
using CSMQNM

@testset "configuration" begin
    @test CSMQNM.canonical_potential_label(Dict(:ipot=>1)) == :schwarzschild_rw
    @test CSMQNM.canonical_potential_label(Dict(:ipot=>:sds_dim_rw)) == :sds_dim_rw
    @test CSMQNM.basis_dimension(BasisConfig(imax=4, nmax=2, range=:real)) == 12
    @test CSMQNM.basis_dimension(BasisConfig(imax=4, nmax=2, range=:complex)) == 24
end

@testset "Schwarzschild potential map" begin
    physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2)
    @test isfinite(CSMQNM.potential_value(physics, 0.0))
end

@testset "full small Schwarzschild solve" begin
    cfg = RunConfig(
        basis = BasisConfig(imax=4, nmax=1, r0=0.2, rmax=4.0, range=:real),
        csm = CSMConfig(the0=25.0),
        integration = IntegrationConfig(xmin=-5.0, xmax=5.0, dx=0.5),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm(cfg)
    @test length(result.energy) == result.basis_size_after
    @test all(imag.(result.omega) .<= 1e-12)
end
