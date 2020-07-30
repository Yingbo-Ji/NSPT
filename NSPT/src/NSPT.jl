module NSPT
export Langevin, φ4
import DSP, Statistics

function Interaction(a::Array, n::Integer)
    a = transpose(a)
    c = DSP.conv(a, DSP.conv(a, a))
    return c[n]
end
function Langevin(n, m, δ)
    #=where n:Langevin time
    m :perturbative order
    δ:time step
    =#
    container = zeros(n, m - 1)
    noise = [randn() for i = 1:n-1] .* sqrt(2δ)
    for i ∈ 1:m-1
        if i == 1
            for j ∈ 1:n-1
                container[j+1, i] = (1 - δ)container[j, i] + noise[j]
            end
        else
            for j ∈ 1:n-1
                container[j+1, i] =
                    (1 - δ)container[j, i] +
                    δ * Interaction(container[j, 1:i-1], i - 1)
            end
        end
    end
    return container
end

function φ4(a, n, m)
    a = transpose(a)
    container = Array{Float64}(undef, m, n)
    for i ∈ 1:n
        tem = DPS.conv(
            DPS.conv(DPS.conv(a[1:end, i], a[1:end, i]), a[1:end, i]),
            a[1:end, i],
        )
        container[1:m, i] = tem[1:m,]
    end
    container = Statistics.mean(container, dims = 2)
end


end # module
