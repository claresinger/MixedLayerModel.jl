using OrdinaryDiffEq
using StochasticDiffEq
using Plots

function f_deterministic(du, u, p, t)
    du[:] .= 0.
    
    CF0 = u[5] < 0.6 ? 0.2 : 1.0
    τCF = 2*24*3600.
    du[5] = (CF0 - u[5]) / τCF
end

function f(du, u, p, t)
    du[:] .= 0.
    
    CF0 = u[5] < 0.6 ? 0.2 : 1.0
    τCF = 2*24*3600.
    τsyn = 5*24*3600.
    CFnew = CF0 + (CF0 - u[5]) / τsyn
    du[5] = (CFnew - u[5]) / τCF
end

function g(du, u, p, t)
    du[:] .= 0.
    σ² = 1e-7;
    du[5] = sqrt(σ²)
end

function g_bounded(du, u, p, t)
    du[:] .= 0.
    σ = sqrt(5 / (5*24*3600.))
    du[5] = σ * (u[5] - 0.1) * (1 - u[5])
end

u0 = [0, 0, 0, 0, 0.2]
params = 0.
dtmax, tspan = 5*3600., (0., 100*24*3600.);
steptol = 1e-3;

prob = ODEProblem(
        ODEFunction(f_deterministic, tgrad=(du, u, p, t) -> fill!(du, 0.0)), 
        u0, 
        tspan, 
        params);
# sol = solve(prob, Rodas5(autodiff=false), abstol=0.0, reltol=steptol, dtmax=dtmax);
sol = solve(prob, Euler(), dt=dtmax);
plot(sol.t/(24*3600.), getindex.(sol.u,5), xlabel="time (days)", ylabel="cloud fraction", lw=2)

prob = SDEProblem(
    SDEFunction(f, g_bounded, tgrad=(du, u, p, t) -> fill!(du, 0.0)),
    g_bounded,
    u0,
    tspan,
    seed = 18);
sol = solve(prob, EM(), dt=dtmax);
# sol = solve(prob, SOSRI(), abstol=0.0, reltol=steptol, dtmax=dtmax);
plot!(sol.t/(24*3600.), getindex.(sol.u,5))
