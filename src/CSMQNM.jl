module CSMQNM

using LinearAlgebra
using Printf
using LambertW
using AutoTortoise

export BasisConfig, CSMConfig, IntegrationConfig, OutputConfig, RunConfig, QNMResult
export solve_qnm, solve_qnm_from_dict, write_potential, potential_value, potential_labels
export ECSConfig, ECSRunConfig, ECSResult, solve_qnm_ecs, solve_qnm_ecs_from_dict
export ecs_contour, ecs_jacobian, write_ecs_spectrum, write_ecs_potential
export FEDVRECSConfig, FEDVRECSRunConfig, FEDVRECSResult
export solve_qnm_ecs_fedvr, solve_qnm_ecs_fedvr_from_dict
export write_fedvr_ecs_spectrum, write_fedvr_ecs_potential

const DEG = π / 180
const ORTHOGONALIZATION_CUTOFF = 1e-5

include("CSM.jl")
include("ECS.jl")
include("FEDVR_ECS.jl")

end
