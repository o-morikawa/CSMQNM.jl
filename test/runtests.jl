using Test
using CSMQNM

@testset "configuration" begin
    @test CSMQNM.canonical_potential_label(Dict(:ipot=>1)) == :schwarzschild_rw
    @test CSMQNM.canonical_potential_label(Dict(:ipot=>:sds_dim_rw)) == :sds_dim_rw
    @test CSMQNM.basis_dimension(BasisConfig(imax=4, nmax=2, range=:real)) == 12
    @test CSMQNM.basis_dimension(BasisConfig(imax=4, nmax=2, range=:complex)) == 12
end

@testset "complex range validation" begin
    valid = RunConfig(
        basis = BasisConfig(imax=4, nmax=1, range=:complex, beta=π/2),
        csm = CSMConfig(the0=42.0),
        physics = Dict(:ipot=>:schwarzschild_rw),
    )
    @test CSMQNM.validate_config(valid) === nothing

    bad_beta = RunConfig(
        basis = BasisConfig(imax=4, nmax=1, range=:complex, beta=0.1),
        csm = CSMConfig(the0=80.0),
        physics = Dict(:ipot=>:schwarzschild_rw),
    )
    @test_throws ErrorException CSMQNM.validate_config(bad_beta)
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

@testset "full small complex-range solve" begin
    cfg = RunConfig(
        basis = BasisConfig(imax=4, nmax=3, r0=0.2, rmax=4.0, range=:complex, beta=π/2),
        csm = CSMConfig(the0=42.0),
        integration = IntegrationConfig(xmin=-5.0, xmax=5.0, dx=0.5),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm(cfg)
    @test length(result.energy) == cfg.basis.imax * (cfg.basis.nmax + 1)
    @test result.dropped_basis_vectors == 0
    @test all(isfinite, real.(result.energy))
    @test all(isfinite, imag.(result.energy))
end

@testset "ECS configuration and contour" begin
    ecs = ECSConfig(theta_deg=60.0, xmin=-10.0, xmax=10.0, nx=101, x0_left=3.0, x0_right=3.0)
    @test CSMQNM.ecs_contour(0.0, ecs) == 0.0 + 0.0im
    @test CSMQNM.ecs_jacobian(0.0, ecs) == 1.0 + 0.0im
    @test CSMQNM.ecs_jacobian(5.0, ecs) ≈ exp(im * 60.0 * π / 180)
end


@testset "ECS kinetic matrix" begin
    ecs = ECSConfig(
        theta_deg=1e-12,
        xmin=-1.0,
        xmax=1.0,
        nx=5,
        x0_left=0.5,
        x0_right=0.5,
    )
    x = collect(range(ecs.xmin, ecs.xmax; length=ecs.nx))
    T = CSMQNM.ecs_kinetic_matrix(x, ecs)
    dx = x[2] - x[1]
    expected = [
        2 / dx^2  -1 / dx^2  0;
        -1 / dx^2  2 / dx^2  -1 / dx^2;
        0  -1 / dx^2  2 / dx^2
    ]
    @test T ≈ ComplexF64.(expected)
    @test count(!iszero, T) == 7
end

@testset "small Schwarzschild ECS solve" begin
    cfg = ECSRunConfig(
        ecs = ECSConfig(theta_deg=50.0, xmin=-8.0, xmax=8.0, nx=81, x0_left=3.0, x0_right=3.0),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm_ecs(cfg)
    @test length(result.energy) == cfg.ecs.nx - 2
    @test length(result.x) == cfg.ecs.nx - 2
    @test all(isfinite, real.(result.energy))
    @test all(isfinite, imag.(result.energy))
end

@testset "FEDVR utilities" begin
    ξ, w = CSMQNM.legendre_gll(6)
    @test length(ξ) == 6
    @test length(w) == 6
    @test ξ[1] ≈ -1.0
    @test ξ[end] ≈ 1.0
    @test sum(w) ≈ 2.0
    D = CSMQNM.lagrange_derivative_matrix(ξ)
    @test size(D) == (6, 6)
    @test isapprox(D * ones(6), zeros(6); atol=1e-10)
end


@testset "FEDVR bridge connectivity" begin
    e = FEDVRECSConfig(xmin=-2.0, xmax=2.0, nelements=4, order=6)
    breaks = CSMQNM.fedvr_element_breaks(e)
    ξ, w = CSMQNM.legendre_gll(e.order)
    x, l2g, scale = CSMQNM.fedvr_bridge_connectivity(e, breaks, ξ, w)
    @test length(x) == e.nelements * (e.order - 1) - 1
    @test l2g[1, 1] == 0
    @test l2g[end, end] == 0
    @test l2g[1, end] == l2g[2, 1]
    @test scale[1, end] > 0
    @test scale[2, 1] > 0
end

@testset "small Schwarzschild FEDVR-ECS solve" begin
    cfg = FEDVRECSRunConfig(
        ecs = FEDVRECSConfig(theta_deg=50.0, xmin=-8.0, xmax=8.0, nelements=8, order=6, x0_left=3.0, x0_right=3.0),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm_ecs_fedvr(cfg)
    expected_size = cfg.ecs.nelements * (cfg.ecs.order - 1) - 1
    @test length(result.energy) == expected_size
    @test size(result.overlap) == (expected_size, expected_size)
    @test all(isfinite, real.(result.energy))
    @test all(isfinite, imag.(result.energy))
end

#=
@testset "Takagi-ECS Gaussian matrix utilities" begin
    cfg = TakagiECSRunConfig(
        ecs = TakagiECSConfig(
            theta_deg=50.0,
            xmin=-8.0,
            xmax=8.0,
            nbasis=12,
            sigma_scale=1.2,
            x0_left=3.0,
            x0_right=3.0,
            quadrature_panels=16,
            quadrature_order=6,
            orthogonalization=:hermitian,
        ),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    H, N, centers, sigmas = CSMQNM.gaussian_ecs_matrices(cfg)
    @test size(H) == (cfg.ecs.nbasis, cfg.ecs.nbasis)
    @test size(N) == (cfg.ecs.nbasis, cfg.ecs.nbasis)
    @test length(centers) == cfg.ecs.nbasis
    @test length(sigmas) == cfg.ecs.nbasis
    @test norm(N - transpose(N)) / max(1, norm(N)) < 1e-10
    @test all(isfinite, real.(H))
    @test all(isfinite, imag.(H))
end

@testset "small Schwarzschild Gaussian-ECS Hermitian solve" begin
    # CI should not depend on experimental Takagi convergence.
    # This test exercises the Gaussian-ECS backend with deterministic
    # Hermitian auxiliary-metric conditioning only.
    cfg = TakagiECSRunConfig(
        ecs = TakagiECSConfig(
            theta_deg=50.0,
            xmin=-8.0,
            xmax=8.0,
            nbasis=12,
            sigma_scale=1.2,
            x0_left=3.0,
            x0_right=3.0,
            quadrature_panels=16,
            quadrature_order=6,
            overlap_cutoff=1e-9,
            orthogonalization=:hermitian,
        ),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm_ecs_takagi(cfg)
    @test result.orthogonalization_method == :hermitian
    @test result.basis_size_after <= result.basis_size_before
    @test length(result.energy) == result.basis_size_after
    @test all(isfinite, real.(result.energy))
    @test all(isfinite, imag.(result.energy))
    @test isfinite(result.orthogonalization_error)
end
=#
