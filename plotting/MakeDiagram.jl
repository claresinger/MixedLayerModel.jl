using PyPlots
using LaTeXStrings

function make_MLM_diagram()
    fig = figure(figsize=(8,6))
    rcParams.update({'font.size': 15})

    plot([0,80],[0,0],"k-")
    fill_between([0,80],[50,50],[100,100],color="grey",alpha=0.5)
    plot([20,20,10,10],[0,100,102,150],"k-")
    plot([65,65,50,45],[0,100,102,150],"k-")
    
    RH = concatenate((58 + linspace(0,7,num=50), 65 + zeros(50), 35 - linspace(0,2,num=5)))
    z = concatenate((linspace(0,100,num=100), linspace(102,150,num=5)))
    plot(RH,z,"k:")

    text(-4,50,"$z_b$")
    text(-4,100,"$z_i$")

    text(66,10,"$q_{tM}$")
    text(21,10,"$h_M$")
    text(51,10,"RH(z)")
    
    text(11,120,"$h^{ft}(z)$")
    text(49,120,"$q_t^{ft}(z)$")
    text(35,120,"RH$^{ft}$")

    #title("Mixed-layer model schematic",fontsize=15)
    axis("off")
    box(False)

    tight_layout()
    savefig("./figs/mlm-diagram.png",dpi=300)
    #show()

end