import cube_sphere
import asteroid
import system as sys
import bridge as b
import main
%matplotlib inline


config = main.get_config()
system = sys.System(config[0])

experiment_config = {
    'time step': 0.1,
    'end time' : 15
}

b.Evolve(system, 0.1, 15)

main.plots(system, experiment_config)
