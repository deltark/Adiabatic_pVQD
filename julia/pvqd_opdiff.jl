using StatsBase
using Yao.AD
using Yao.AD: Rotor
using YaoExtensions

using Yao, PastaQ

import YaoExtensions.faithful_grad
import Yao.AbstractBlock

include("YaoPastaQ.jl")

const CPhaseGate{N, T} = ControlBlock{N,<:ShiftGate{T},<:Any}

############# manipulate quantum differential node ############
Yao.AD.generator(c::CPhaseGate{N}) where N = ControlBlock{N}(c.ctrl_locs, c.ctrl_config, Z, c.locs)

function get_diffblocks(circuit::AbstractBlock)
    diffblocks = collect_blocks(Union{RotationGate, CPhaseGate}, circuit)
    length(diffblocks) == nparameters(circuit) || throw(ArgumentError("Input circuit contains non-differentiable/unsupported parameters."))
    return diffblocks
end

#### interface #####
export faithful_grad

@inline function _perturb(func, gate::AbstractBlock, δ::Real)
    dispatch!(-, gate, (δ,))
    r1 = func()
    dispatch!(+, gate, (2δ,))
    r2 = func()
    dispatch!(-, gate, (δ,))
    r1, r2
end

@inline function _expect(circuit::AbstractBlock, reg::ArrayReg; nshots=nothing)
    if nshots === nothing
        expect(op, reg)
    else
        mean(measure(op, reg; nshots=nshots))
    end
end

@inline function _expect_pasta(circuit::AbstractBlock, noise::Union{Nothing, Tuple, NamedTuple}; nshots=100)
    if nshots===nothing
        throw(ArgumentError("Noisy simulations require shots, nshots cannot be nothing"))
    end

    N = YaoBase.nqubits(circuit)

    pasta_circ = YaoPastaQ.genlist(circuit)
    bases = randombases(N, nshots; local_basis = ["Z"])
    ψ = runcircuit(N, pasta_circ, noise = noise)
    pasta_data = getsamples(ψ, bases)

    processed_pasta = prod(last.(pasta_data); dims=2)
    #last. is to get the second member of each pair
    #The product does the opposite of projector_zero, then we flip it
    processed_pasta = 1 .- processed_pasta

    mean(processed_pasta)
end

## same but forcing the diffblocks (added by user)
## now with noise (expectation value only on projector_zero for nown, so only usable for overlap)
@inline function faithful_grad(op::AbstractBlock, pair::Pair{<:ArrayReg, <:AbstractBlock}, diffblock_forced::AbstractBlock, noise::Union{Nothing, Tuple, NamedTuple}=nothing; nshots=nothing)
    map(get_diffblocks(diffblock_forced)) do diffblock
        if noise===nothing
            r1, r2 = _perturb(()->_expect(op, copy(pair.first) |> pair.second; nshots=nshots) |> real, diffblock, π/2)
        else
            r1, r2 = _perturb(()->_expect_pasta(pair.second, noise; nshots) |> real, diffblock, π/2)
        end
        (r2 - r1)/2
    end
end