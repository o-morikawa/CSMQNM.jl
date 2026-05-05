using CSMQNM

result = solve_qnm_from_dict(
    basis = Dict(:imax=>30, :nmax=>3, :r0=>0.1, :rmax=>60.0, :range=>:complex),
    csm = Dict(:the0=>42.0),
    integration = Dict(:xmin=>-60.0, :xmax=>80.0, :dx=>0.02),
    physics = Dict(:ipot=>:sds4_rw, :M=>1.0, :ell=>2, :s=>2, :lam=>0.08),
    output = Dict(:write_spectrum=>true),
)

println(result.omega[1:min(10, length(result.omega))])
