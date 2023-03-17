using Statistics

CO2updn_list = [200,300,400,800,1000,1200,1300,1400,1600,1400,1300,1200,1000,800,400,300,200];
nd = length(CO2updn_list);
SSTupdn_list = [287.7,289.1,290.0,292.2,293.2,294.3,304.5,305.8,308.0,306.2,304.7,303.7,302.0,300.9,297.9,296.8,287.6];
LHFupdn_list = [97.8,103.9,107.1,112.8,115.3,120.5,208.7,213.5,220.9,214.4,209.4,206.0,199.7,195.2,183.7,179.1,96.1];
ziupdn_list = [1442,1349,1266,1078,1011,972,834,781,703,766,826,866,949,1013,1230,1340,1416];
μSST, σSST = mean(SSTupdn_list), std(SSTupdn_list);
μLHF, σLHF = mean(LHFupdn_list), std(LHFupdn_list);

function normalize_data(x, name="SST")
    if name == "SST"
        x = (x .- μSST) ./ σSST;
    end
    if name == "LHF"
        x = (x .- μLHF) ./ σLHF;
    end
    return x
end

function unnormalize_data(x, name="SST")
    if name == "SST"
        x = x .* σSST .+ μSST;
    end
    if name == "LHF"
        x = x .* σLHF .+ μLHF;
    end
    return x
end