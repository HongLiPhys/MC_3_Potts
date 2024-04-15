module MC
using Plots
using DataStructures
using DelimitedFiles

struct NematicOrders
    size::Int64 #Sample size: unit cell numbers along one axis
    steps::Int64 #MC simulation sweep numbers
    T::Float64 #Temperature in K units
    kB::Float64 #Boltzmann constant
    h::Float64 #Strain field strength
    grid::Matrix{Int64} #matrix of Vector order map
    Field::Matrix{Int64} #Field distribution matrix, must be the same size as grid
    M::Dict{Int64,Vector{Int64}} #Magnetization for n = -1,0,1; total number of a same n 
end

function initialize_NematicOrders(h)
    size = 80
    steps = 10^8
    T = 5 #K
    kB = 8.6 * 10^(-2) # meV/T
    board = rand((-1, 0, 1), size, size)
    field = zeros(Float64, size, size)
    for i in 1:size
        for j in 1:size
            if (i - size / 3)^2 + (j - 3 * size / 5)^2 <= (size / 5)^2 || (i - 3 * size / 5)^2 + (j - size / 3)^2 <= (size / 5)^2
                field[i, j] = -1
            else
                field[i, j] = 1
            end
        end
    end
    Ms = Dict(-1 => zeros(Int64, steps + 1), 0 => zeros(Int64, steps + 1), 1 => zeros(Int64, steps + 1))
    m0 = counter(board)
    Ms[-1][1] = m0[-1]
    Ms[0][1] = m0[0]
    Ms[1][1] = m0[1]
    return NematicOrders(
        size,
        steps,
        T,
        kB,
        h,
        board,
        field,
        Ms
    )
end



function δE(i, j, s, n::NematicOrders)
    size = n.size
    function PBC(site)
        if site == 0
            size
        elseif site == size + 1
            1
        else
            site
        end
    end
    neighbours = [n.grid[PBC(i - 1), j],
        n.grid[PBC(i + 1), j],
        n.grid[i, PBC(j - 1)],
        n.grid[i, PBC(j + 1)],
        n.grid[PBC(i - 1), PBC(j + 1)],
        n.grid[PBC(i + 1), PBC(j - 1)]] #This is a Triangular lattice model, 6 nearest neighbour hopping term counted with Periodic Boundary Condition imposed.
    c = counter(neighbours)
    if n.grid[i, j] == s
        0
    else
        -(c[s] - c[n.grid[i, j]]) + n.h * (
            if n.grid[i, j] == n.Field[i, j]
                1
            else
                0
            end
        ) - n.h * (
            if s == n.Field[i, j]
                1
            else
                0
            end
        )
    end
end



function simulate(n)
    for step in 1:n.steps
        n.M[-1][step+1] = n.M[-1][step]
        n.M[0][step+1] = n.M[0][step]
        n.M[1][step+1] = n.M[1][step]
        i = rand(1:n.size)
        j = rand(1:n.size)
        old_s = n.grid[i, j]
        ss = filter!(e -> e ≠ old_s, [-1, 0, 1])
        s = rand(ss)
        δe = δE(i, j, s, n)
        if δe <= 0
            n.grid[i, j] = s
        elseif rand() < exp(-δe / (n.kB * n.T))
            n.grid[i, j] = s
        end

        if n.grid[i, j] != old_s
            n.M[old_s][step+1] = n.M[old_s][step] - 1
            n.M[s][step+1] = n.M[s][step] + 1
        end
    end
end

function show_Ordermap_M(n)
    steps = n.steps
    Mn1 = zeros(Int64, 100)
    M0 = zeros(Int64, 100)
    Mp1 = zeros(Int64, 100)

    for i in 1:100
        step = floor(Int, 1 + (i - 1) * steps / 100)
        Mn1[i] = n.M[-1][step]
        M0[i] = n.M[0][step]
        Mp1[i] = n.M[1][step]
    end
    ordermap = heatmap(1:n.size,
        1:n.size, n.grid,
        c=cgrad([:blue, :white, :red, :yellow]),
        xlabel="x axis", ylabel="y axis",
        title="Nematic vector map",
        aspect_ratio=1)
    magnetization = plot(range(0, 100, length=100), [Mn1, M0, Mp1])
    display(plot(ordermap, magnetization, layout=(1, 2)))
end

function generate_data(h_min, δh, h_max)
    for h in h_min:δh:h_max
        n = initialize_NematicOrders(h)
        simulate(n)
        steps = n.steps
        points_number = 200
        M_short = zeros(Int64, 3, points_number)
        for i in 1:points_number
            step = floor(Int, 1 + (i - 1) * steps / points_number)
            M_short[1, i] = n.M[-1][step]
            M_short[2, i] = n.M[0][step]
            M_short[3, i] = n.M[1][step]
        end
        writedlm("n_map_" * string(h) * ".dat", n.grid)
        writedlm("M_curve_" * string(h) * ".dat", M_short)
    end

end

n = initialize_NematicOrders(0.2)
@time simulate(n)
@time show_Ordermap_M(n)
@time generate_data(0.3, 0.1, 1)


end # module MC
