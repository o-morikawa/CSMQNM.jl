using CSMQNM

cfg = RunConfig(
    basis = BasisConfig(imax=30, nmax=3, r0=0.1, rmax=60.0, range=:real),
    csm = CSMConfig(the0=42.0),
    integration = IntegrationConfig(xmin=-60.0, xmax=80.0, dx=0.02),
    physics = Dict(
        :ipot => :sds_dim_rw,
        :M => 1.0,
        :ell => 2,
        :s => 0,
        :lam => 0.4,
        :dim => 5,
        :c => :tensor,
    ),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm(cfg)

println("basis size before orthogonalization: ", result.basis_size_before)
println("basis size after orthogonalization:  ", result.basis_size_after)
println("first 10 omega values:")
println(result.omega[1:min(10, length(result.omega))])
