import subprocess
import os
import shlex
import csv
import torch
import numpy as np
import matplotlib.pyplot as plt
from botorch.models import SingleTaskGP
from gpytorch.likelihoods import FixedNoiseGaussianLikelihood
from gpytorch.kernels import MaternKernel, ScaleKernel, RBFKernel
from gpytorch.means import ConstantMean, ZeroMean
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.fit import fit_gpytorch_mll
from botorch.acquisition.monte_carlo import qExpectedImprovement
from botorch.acquisition import ExpectedImprovement
from botorch.optim import optimize_acqf
from botorch.models.transforms.input import InputStandardize
from botorch.models.transforms.outcome import Standardize
from gpytorch.priors import GammaPrior
#import variables

#print(variables.end_time)
#variables.end_time = 5
#print(variables.end_time)

########################################################################################################################
#Customize code here

#Major changes:
nu = 1.5
matern = True
rbf = False

const = False
zero = True

num_consecutive_trials = 4
EI = True
max_stddev = False
freq = 3

#Minor changes:
fixed_Yvar = 1e-6
lower_bound = 0.
upper_bound = 10.
sobol_on = True
improvement_threshold = 1e-2
num_initial_trials = 2 #this needs to be >=2
########################################################################################################################


def simulation(force):
    force = force.numpy()[0]
    command = shlex.split(f"./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin {force}")
    subprocess.run(command)

    f = open("muscle_length_prestretch.csv")
    reader = csv.reader(f, delimiter=",")
    for row in reader:
        prestretch = float(row[1]) - float(row[0])
        print("The muscle was stretched ", prestretch)
    f.close()

    f = open("muscle_length_contraction.csv")
    reader = csv.reader(f, delimiter=",")
    muscle_length_process = []
    for row in reader:
        for j in row:
            muscle_length_process.append(j)
        
    contraction = float(muscle_length_process[0]) - float(muscle_length_process[-2])
    print("The muscle contracted ", contraction)
    f.close()

    return contraction


class CustomSingleTaskGP(SingleTaskGP):
    def __init__(self, train_X, train_Y):
        train_Yvar = torch.full_like(train_Y, fixed_Yvar, dtype=torch.double)
        likelihood = FixedNoiseGaussianLikelihood(noise=train_Yvar)
        if matern:
            kernel = ScaleKernel(MaternKernel(nu=nu))
        elif rbf:
            kernel = ScaleKernel(RBFKernel())
        else:
            kernel = ScaleKernel(MaternKernel(nu=nu))
        if const:
            mean = ConstantMean()
        elif zero:
            mean = ZeroMean()
        else:
            mean = ConstantMean()
        input_transform = InputStandardize(d=1)
        output_tansform = Standardize(m=1)

        super().__init__(train_X,
                         train_Y,
                         likelihood=likelihood,
                         covar_module=kernel,
                         mean_module=mean,
                         #input_transform=input_transform,
                         #outcome_transform=output_tansform,
                        )


os.chdir("..")
os.chdir("cuboid_muscle/build_release")


def scale_samples(samples, lower_bound, upper_bound):
    return lower_bound + (upper_bound - lower_bound) * samples

sobol = torch.quasirandom.SobolEngine(dimension=1, scramble=True)
if sobol_on:
    initial_x = scale_samples(sobol.draw(num_initial_trials, dtype=torch.double), lower_bound, upper_bound)
else:
    initial_x = torch.tensor([[0.],[10.],[5.]], dtype=torch.double)

initial_y = torch.tensor([])
for force in initial_x:
    y = torch.tensor([[simulation(force)]], dtype=torch.double)
    initial_y = torch.cat([initial_y, y])
#initial_y = torch.tensor(simulation(initial_x),dtype=torch.double)
initial_yvar = torch.full_like(initial_y, fixed_Yvar, dtype=torch.double)


gp = CustomSingleTaskGP(initial_x, initial_y)
mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
fit_gpytorch_mll(mll)

print("Lengthscale:", gp.covar_module.base_kernel.lengthscale.item())
print("Outputscale:", gp.covar_module.outputscale.item())
print("Noise:", gp.likelihood.noise.mean().item())

num_iterations = 100
best_value = -float('inf')
no_improvement_trials = 0
counter = num_initial_trials

for i in range(num_iterations):
    if EI:
        acq_fct = ExpectedImprovement(model=gp, best_f=initial_y.max())
    else:
        acq_fct = ExpectedImprovement(model=gp, best_f=initial_y.max())

    if not max_stddev or counter%freq!=0:
        candidate, acq_value = optimize_acqf(
            acq_function=acq_fct,
            bounds=torch.tensor([[lower_bound], [upper_bound]], dtype=torch.double),
            q=1,
            num_restarts=200,
            raw_samples=512,
        )
    else:
        x_query = torch.linspace(lower_bound, upper_bound, 1000).unsqueeze(-1)
        posterior = gp.posterior(x_query)
        variance = posterior.variance.squeeze(-1)
        stddev = torch.sqrt(variance).detach().numpy()
        argmax_stddev = stddev.argmax()
        candidate = torch.tensor([x_query[argmax_stddev].numpy()])
        print("max stddev used")

    new_y = torch.tensor([[simulation(candidate[0])]], dtype=torch.double)
    #new_y = simulation(candidate).clone().detach()
    new_yvar = torch.full_like(new_y, fixed_Yvar, dtype=torch.double)

    initial_x = torch.cat([initial_x, candidate])
    initial_y = torch.cat([initial_y, new_y])
    initial_yvar = torch.cat([initial_yvar, new_yvar])
    gp = CustomSingleTaskGP(initial_x, initial_y)
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    x_query = torch.linspace(lower_bound, upper_bound, 1000).unsqueeze(-1)
    posterior = gp.posterior(x_query)

    mean = posterior.mean.squeeze(-1).detach().numpy()
    variance = posterior.variance.squeeze(-1)
    stddev = torch.sqrt(variance).detach().numpy()

    print("Lengthscale:", gp.covar_module.base_kernel.lengthscale.item())
    print("Outputscale:", gp.covar_module.outputscale.item())
    print("Noise:", gp.likelihood.noise.mean().item())

    current_value = new_y.item()
    if current_value > best_value + improvement_threshold:
        best_value = current_value
        no_improvement_trials = 0
    elif len(initial_x) > num_initial_trials:
        no_improvement_trials += 1

    counter += 1

    print(f"Trial {i + 1 + num_initial_trials}: x = {candidate.item()}, Value = {current_value}, Best Value = {best_value}")

    if no_improvement_trials >= num_consecutive_trials:
        print("Stopping criterion met. No significant improvement for consecutive trials.")
        print("Number of total trials: ", i+1+num_initial_trials)
        break


    plt.scatter(initial_x.numpy(), initial_y.numpy(), color="red", label="Trials", zorder=3)
    plt.plot(initial_x.numpy(), initial_y.numpy(), color="red", linestyle="", markersize=3)
    plt.plot(x_query, mean)
    plt.scatter(candidate.numpy(), new_y.numpy(), color="green", s=30, zorder=5, label="New query point")
    plt.fill_between(x_query.numpy().squeeze(), mean - 2 * stddev, mean + 2 * stddev, alpha=0.3, label="GP 95% CI")
    plt.xlabel("x")
    plt.ylabel("Objective Value")
    plt.title("Optimization Process")
    plt.legend()
    plt.show()

x_query = torch.linspace(0, 10, 1000).unsqueeze(-1)
posterior = gp.posterior(x_query)

mean = posterior.mean.squeeze(-1).detach().numpy()
variance = posterior.variance.squeeze(-1)
stddev = torch.sqrt(variance).detach().numpy()

plt.scatter(initial_x.numpy(), initial_y.numpy(), color="red", label="Trials", zorder=3)
plt.plot(initial_x.numpy(), initial_y.numpy(), color="red", linestyle="", markersize=3)
plt.plot(x_query, mean)
plt.fill_between(x_query.numpy().squeeze(), mean - 2 * stddev, mean + 2 * stddev, alpha=0.3, label="GP 95% CI")
plt.xlabel("x")
plt.ylabel("Objective Value")
plt.title("Optimization Results")
plt.legend()
plt.show()

continuing = input("Do you want to add another query point? (y/n)")

while continuing == "y":
    candidate = input("Which point do you want to add?")
    candidate = torch.tensor([[float(candidate)]], dtype=torch.double)

    new_y = torch.tensor([[simulation(candidate[0])]], dtype=torch.double)
    new_yvar = torch.full_like(new_y, fixed_Yvar, dtype=torch.double)

    initial_x = torch.cat([initial_x, candidate])
    initial_y = torch.cat([initial_y, new_y])
    initial_yvar = torch.cat([initial_yvar, new_yvar])
    gp = CustomSingleTaskGP(initial_x, initial_y)
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    x_query = torch.linspace(lower_bound, upper_bound, 1000).unsqueeze(-1)
    posterior = gp.posterior(x_query)

    mean = posterior.mean.squeeze(-1).detach().numpy()
    variance = posterior.variance.squeeze(-1)
    stddev = torch.sqrt(variance).detach().numpy()

    plt.scatter(initial_x.numpy(), initial_y.numpy(), color="red", label="Trials", zorder=3)
    plt.plot(initial_x.numpy(), initial_y.numpy(), color="red", linestyle="", markersize=3)
    plt.plot(x_query, mean)
    plt.scatter(candidate.numpy(), new_y.numpy(), color="green", s=30, zorder=5, label="New query point")
    plt.fill_between(x_query.numpy().squeeze(), mean - 2 * stddev, mean + 2 * stddev, alpha=0.3, label="GP 95% CI")
    plt.xlabel("x")
    plt.ylabel("Objective Value")
    plt.title("Optimization Process")
    plt.legend()
    plt.show()

    continuing = input("Do you want to add another query point? (y/n)")