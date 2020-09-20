import numpy as np
import matplotlib.pyplot as plt

def main():
    make_MLM_diagram()

def make_MLM_diagram():
    fig = plt.figure(figsize=(8,6))
    plt.rcParams.update({'font.size': 15})

    plt.plot([0,80],[0,0],"k-")
    plt.fill_between([0,80],[50,50],[100,100],color="grey",alpha=0.5)
    plt.plot([20,20,10,10],[0,100,102,150],"k-")
    plt.plot([65,65,50,45],[0,100,102,150],"k-")
    
    RH = np.concatenate((58+np.linspace(0,7,num=50), 65+np.zeros(50), 35-np.linspace(0,2,num=5)))
    z = np.concatenate((np.linspace(0,100,num=100),np.linspace(102,150,num=5)))
    plt.plot(RH,z,"k:")

    plt.text(-4,50,"$z_b$")
    plt.text(-4,100,"$z_i$")

    plt.text(66,10,"$q_{tM}$")
    plt.text(21,10,"$h_M$")
    plt.text(51,10,"RH(z)")
    
    plt.text(11,120,"$h^{ft}(z)$")
    plt.text(49,120,"$q_t^{ft}(z)$")
    plt.text(35,120,"RH$^{ft}$")

    #plt.title("Mixed-layer model schematic",fontsize=15)
    plt.axis("off")
    plt.box(False)

    plt.tight_layout()
    plt.savefig("./figs/mlm-diagram.png",dpi=300)
    #plt.show()

if __name__ == '__main__':
    main()