
@everywhere module RegionalCF

using Distributed
using NCDatasets
using MixedLayerModel
include("mlm_solve_funcs.jl")

function one_day(dayi, ds, skip1, skip2, par)
    println(dayi)
    println(ds)
    lon = ds["lon"][1:skip1:end]
    lat = ds["lat"][1:skip2:end]
    N, M = length(lon), length(lat)

    CFday = fill(NaN, (N, M))
    for (i1, _) in enumerate(lon)
        println(i1)
        for (i2, _) in enumerate(lat)
            println(i2)
            local j1 = 1 + (i1-1)*skip1
            local j2 = 1 + (i2-1)*skip2
            
            if typeof(ds["sst"][j1,j2,dayi]) == Missing
                println("skip")
            else
                par.SST0 = ds["sst"][j1,j2,dayi] # (K)
                par.V = ds["WS"][j1,j2,dayi] # m/s
                par.D = ds["D500"][j1,j2,dayi] # (1/s)
                par.RHft = ds["RH500"][j1,j2,dayi] # (-)
                par.EIS0 = ds["EIS"][j1,j2,dayi] # (K)
                par.CO2 = 400 # (ppm)
                try
                    _, sol = run_mlm(par, dt=3600.0*24.0*5, tspan=(0.0,3600.0*24.0*100), quiet=true);
                    CFday[i1,i2] = sol.u[end][5];
                catch
                    println("fail")
                end
            end
        end
    end
    return CFday
end
end