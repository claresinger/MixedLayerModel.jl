using OrdinaryDiffEq
using StochasticDiffEq
using Plots

function CF(x)
    # u[5] < 0.6 ? 0.2 : 1.0
    return tanh(8*(x-0.6))*0.5*0.75 + 0.6
end

function fd(du, u, p, t)
    du[:] .= 0.
    τCF = 2*24*3600.
    du[5] = (CF(u[5]) - u[5]) / τCF
end

function fs(du, u, p, t)
    du[:] .= 0.
    τCF = 2*24*3600.
    τsyn = 5*24*3600.
    CFnew = CF(u[5]) + (CF(u[5]) - u[5]) / τsyn
    du[5] = (CFnew - u[5]) / τCF
end

function gu(du, u, p, t)
    du[:] .= 0.
    σ² = 1e-7;
    du[5] = sqrt(σ²)
end

function gb(du, u, p, t)
    du[:] .= 0.
    σ = sqrt(5 / (5*24*3600.))
    du[5] = σ * (u[5] - 0) * (1 - u[5])
end

function gd(du, u, p, t)
    noise = 3e-3
    du[:] .= 0.
    du[5] = noise * (u[5]-0) * (1 - u[5])
end

u0 = [0, 0, 0, 0, 1.0]
params = 0.
dtmax, tspan = 1*3600., (0., 1000*24*3600.);
steptol = 1e-3;

prob = ODEProblem(
        ODEFunction(fd, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
        u0, 
        tspan, 
        params
);
sol = solve(prob, Euler(), dt=dtmax);
plot(sol.t/(24*3600.), getindex.(sol.u,5), lw=2,
    xlabel="time (days)", ylabel="cloud fraction", label="deterministic")

prob = SDEProblem(
    SDEFunction(fs, gb, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
    gb,
    u0,
    tspan,
    seed = 18
);
sol = solve(prob, EM(), dt=dtmax);
plot!(sol.t/(24*3600.), getindex.(sol.u,5), label="stochastic")

τsyn = 5*24*3600. # 5 days
W = OrnsteinUhlenbeckProcess(1/τsyn, 0., 1., 0., 0.)
prob = SDEProblem(
    fd,
    gd,
    u0,
    tspan,
    seed = 18,
    noise = W
);
sol = solve(prob, EM(), dt=dtmax);
plot!(sol.t/(24*3600.), getindex.(sol.u,5), label="O-U")

savefig("SDEtest-down.png")

############
# ensemble
############

prob = ODEProblem(
        ODEFunction(fd, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
        u0, 
        tspan, 
        params
);
sol = solve(prob, Euler(), dt=dtmax);
plot(sol, lw=2, xlabel="time (s)", ylabel="cloud fraction", label="deterministic")

prob = SDEProblem(
    SDEFunction(fs, gb, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
    gb,
    u0,
    tspan,
);
ensembleprob = EnsembleProblem(prob);
sol = solve(ensembleprob, EM(), EnsembleThreads(), trajectories=100, dt=dtmax);
plot!(EnsembleSummary(sol), label="stochastic")

τsyn = 5*24*3600. # 5 days
W = OrnsteinUhlenbeckProcess(1/τsyn, 0., 1., 0., 0.)
prob = SDEProblem(
    fd,
    gd,
    u0,
    tspan,
    noise = W
);
ensembleprob = EnsembleProblem(prob);
sol = solve(ensembleprob, EM(), EnsembleThreads(), trajectories=100, dt=dtmax);
plot!(EnsembleSummary(sol), ls=:dot, label="OU")

savefig("SDEtest-ensemble-down.png")