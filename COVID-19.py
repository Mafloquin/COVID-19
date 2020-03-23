import numpy as np
from matplotlib import pyplot as plt

def seir_model_with_soc_dist(init_vals, params, t):
    S_0, E_0, I_0, R_0 = init_vals
    S, E, I, R = [S_0], [E_0], [I_0], [R_0]
    alpha, beta, gamma, rho = params
    dt = t[1] - t[0]
    for _ in t[1:]:
        next_S = S[-1] - (rho*beta*S[-1]*I[-1])*dt
        next_E = E[-1] + (rho*beta*S[-1]*I[-1] - alpha*E[-1])*dt
        next_I = I[-1] + (alpha*E[-1] - gamma*I[-1])*dt
        next_R = R[-1] + (gamma*I[-1])*dt
        S.append(next_S)
        E.append(next_E)
        I.append(next_I)
        R.append(next_R)
    return np.stack([S, E, I, R]).T

# Define parameters
t_max = 200 #jours
periode_incubation = 5 #jours
beta = 1.75 #taux de contact de la population
periode_infection_mean = 2 #jours
population_france = 66990000

dt = .1
t = np.linspace(0, t_max, int(t_max/dt) + 1)
N = 10000
init_vals = 1 - 1/N, 1/N, 0, 0
alpha = 1/periode_incubation #Inverse de la période d'incubation
gamma = 1/periode_infection_mean
rho = [1.0,0.9,0.8,0.7,0.6,0.5,0.4]

plt.figure(0)

#SIMULATION
for counter, value in enumerate(rho):
    plt.cla()
    plt.ylim(0,30)
    plt.grid()
    params = alpha, beta, gamma, value
    results = seir_model_with_soc_dist(init_vals, params, t)
    
    text_legend_exposed = r'Exposé ($\rho$=' + str(value) + ')'
    text_legend_infected = r'Infecté ($\rho$=' + str(value) + ')'
    couleur = np.random.rand(3,)
    plt.plot(np.arange(0,t_max+dt,dt), results[:,1] * 100,c=couleur,label=text_legend_exposed)
    plt.plot(np.arange(0,t_max+dt,dt), results[:,2] * 100,c=couleur,linestyle='--',label=text_legend_infected)
    plt.xlabel('Jours [-]')
    plt.ylabel('Pourcentage de la population [%]')
    plt.title(r'COVID-19, Modèle SEIR avec confinement ($\alpha$=' + str(alpha) + r', $\beta$=' + str(beta) + r', $\gamma$=' + str(gamma) + ')')
    plt.legend()
    plt.savefig('rho=' + str(value) + '.jpg')