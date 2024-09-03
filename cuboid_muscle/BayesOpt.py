import subprocess
import sys
import os
import shlex
import csv
import torch
import numpy as np
import matplotlib.pyplot as plt
from botorch.models import SingleTaskGP
from gpytorch.likelihoods import GaussianLikelihood
from gpytorch.kernels import MaternKernel, ScaleKernel, RBFKernel
from gpytorch.means import ConstantMean, ZeroMean
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.fit import fit_gpytorch_mll
from botorch.acquisition import ExpectedImprovement, ProbabilityOfImprovement, PosteriorMean
from botorch.acquisition.knowledge_gradient import qKnowledgeGradient
from botorch.acquisition.max_value_entropy_search import qMaxValueEntropy
from botorch.optim import optimize_acqf
from botorch.models.transforms.input import Normalize
from botorch.models.transforms.outcome import Standardize
import time


script_path = os.path.dirname(os.path.abspath(__file__))
var_path = os.path.join(script_path, "variables")
sys.path.insert(0, var_path)

import variables

########################################################################################################################
#Customize code here

#Major changes:
nu = 1.5
matern = True
rbf = False

const = False
zero = True

fixed_noise = False
variable_noise = True

EI = False
PI = False
KG = False
ES = True

stopping_y = False
improvement_threshold = 1e-4
stopping_xy = True
x_range = 5e-2
num_consecutive_trials = 3

#Minor changes:
fixed_Yvar = 1e-6
lower_bound = 0.
#upper_bound = 20. ########### what should this be?
sobol_on = True
num_initial_trials = 2 #this needs to be >=2
visualize = True
add_points = False
relative_prestretch_min = 1.5
relative_prestretch_max = 1.6
########################################################################################################################

global_individuality_parameter = ""
if matern:
    global_individuality_parameter = global_individuality_parameter + "_matern_" + str(nu)
elif rbf:
    global_individuality_parameter = global_individuality_parameter + "_rbf"
if const:
    global_individuality_parameter = global_individuality_parameter + "_const"
elif zero:
    global_individuality_parameter = global_individuality_parameter + "_zero"
if fixed_noise:
    global_individuality_parameter = global_individuality_parameter + "_fixed_noise"
elif variable_noise:
    global_individuality_parameter = global_individuality_parameter + "_variable_noise"
if EI:
    global_individuality_parameter = global_individuality_parameter + "_EI"
elif PI:
    global_individuality_parameter = global_individuality_parameter + "_PI"
elif KG:
    global_individuality_parameter = global_individuality_parameter + "_KG"
elif ES:
    global_individuality_parameter = global_individuality_parameter + "_ES"
if stopping_y:
    global_individuality_parameter = global_individuality_parameter + "_stopping_y"
elif stopping_xy:
    global_individuality_parameter = global_individuality_parameter + "_stopping_xy"


def simulation(force):
    x = force.numpy()[0]
    f = -0.001678*x**2 + 0.05034*x
    return f
"""

def simulation(force):
    force = force.numpy()[0]
    print("start simulation with force", force)
    individuality_parameter = str(int(time.time()))+str(force)
    command = shlex.split(f"./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin {force} {individuality_parameter}")
    subprocess.run(command)

    print("end simulation")

    f = open("muscle_length_prestretch"+individuality_parameter+".csv")
    reader = csv.reader(f, delimiter=",")
    for row in reader:
        prestretch = float(row[1]) - float(row[0])
        print("The muscle was stretched ", prestretch)
    f.close()

    f = open("muscle_length_contraction"+individuality_parameter+".csv")
    reader = csv.reader(f, delimiter=",")
    muscle_length_process = []
    for row in reader:
        for j in row:
            muscle_length_process.append(j)
        
    contraction = float(muscle_length_process[0]) - float(muscle_length_process[-2])
    print("The muscle contracted ", contraction)
    f.close()

    return contraction
"""

class CustomSingleTaskGP(SingleTaskGP):
    def __init__(self, train_X, train_Y):
        train_Yvar = torch.full_like(train_Y, fixed_Yvar, dtype=torch.double)
        if fixed_noise:
            likelihood = GaussianLikelihood(noise=train_Yvar)
        elif variable_noise:
            likelihood = GaussianLikelihood()
        else:
            print("Wrong input, used variable noise instead")
            likelihood = GaussianLikelihood()
        if matern:
            kernel = ScaleKernel(MaternKernel(nu=nu))
        elif rbf:
            kernel = ScaleKernel(RBFKernel())
        else:
            print("Wrong input, used Matern Kernel with nu=1.5 instead")
            kernel = ScaleKernel(MaternKernel(nu=1.5))

        if const:
            mean = ConstantMean()
        elif zero:
            mean = ZeroMean()
        else:
            print("Wrong input, used Constant Mean instead")
            mean = ConstantMean()

        input_transform = Normalize(d=train_X.shape[-1])
        output_transform = Standardize(m=1)

        super().__init__(train_X,
                         train_Y,
                         likelihood=likelihood,
                         covar_module=kernel,
                         mean_module=mean,
                         input_transform=input_transform,
                         outcome_transform=output_transform,
                        )

def find_relative_prestretch(force):
    print(force)
    individuality_parameter = str(int(time.time()))+"_"+str(force)
    command = shlex.split(f"./incompressible_mooney_rivlin_2 ../prestretch_tensile_test.py incompressible_mooney_rivlin_2 {force} {individuality_parameter}")
    subprocess.run(command)

    f = open("muscle_length_prestretch"+individuality_parameter+".csv")
    reader = csv.reader(f, delimiter=",")
    for row in reader:
        relative_prestretch = float(row[1]) / float(row[0])
    f.close()

    command2 = shlex.split("rm muscle_length_prestretch"+individuality_parameter+".csv")
    subprocess.run(command2)

    return relative_prestretch

def find_upper_bound():
    lower_guess = 0
    upper_guess = 10
    relative_prestretch = find_relative_prestretch(upper_guess)

    while relative_prestretch < relative_prestretch_min or relative_prestretch > relative_prestretch_max:
        if relative_prestretch < relative_prestretch_min:
            not_relevant = upper_guess
            upper_guess = 2*upper_guess - lower_guess
            lower_guess = not_relevant
        if relative_prestretch > relative_prestretch_max:
            upper_guess = (upper_guess + lower_guess)/2
        relative_prestretch = find_relative_prestretch(upper_guess)
        
    return upper_guess

os.chdir("build_release")

upper_bound = find_upper_bound()

print(upper_bound)

starting_time = time.time()

sobol = torch.quasirandom.SobolEngine(dimension=1, scramble=True)
if sobol_on:
    initial_x = sobol.draw(num_initial_trials, dtype=torch.double)
else:
    initial_x = torch.linspace(0, 1, num_initial_trials)

with open("BayesOpt_outputs"+global_individuality_parameter+".csv", "w"):
    pass

initial_y = torch.tensor([])
for force in initial_x:
    y = torch.tensor([[simulation(force*(upper_bound-lower_bound)+lower_bound)]], dtype=torch.double)

    with open("BayesOpt_outputs"+global_individuality_parameter+".csv", "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([force.numpy()[0]*(upper_bound-lower_bound)+lower_bound, y.numpy()[0,0]])

    initial_y = torch.cat([initial_y, y])
initial_yvar = torch.full_like(initial_y, fixed_Yvar, dtype=torch.double)


gp = CustomSingleTaskGP(initial_x, initial_y)
mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
fit_gpytorch_mll(mll)

print("Lengthscale:", gp.covar_module.base_kernel.lengthscale.item())
#print("Outputscale:", gp.covar_module.outputscale.item())
print("Noise:", gp.likelihood.noise.mean().item())

num_iterations = 100
best_value = -float('inf')
no_improvement_trials = 0
counter = num_initial_trials

for i in range(num_iterations):
    if EI:
        acq_fct = ExpectedImprovement(model=gp, best_f=initial_y.max())
    elif PI:
        acq_fct = ProbabilityOfImprovement(model=gp, best_f=initial_y.max())
    elif KG:
        pass
    elif ES:
        bounds=torch.tensor([[0], [1]], dtype=torch.double)
        candidate_set = torch.rand(1000, bounds.size(1))
        candidate_set = bounds[0] + (bounds[1] - bounds[0]) * candidate_set
        acq_fct = qMaxValueEntropy(model=gp, candidate_set=candidate_set)
    else:
        print("Wrong input, used Expected Improvement instead.")
        acq_fct = ExpectedImprovement(model=gp, best_f=initial_y.max())

    if KG:
        SMOKE_TEST = os.environ.get("SMOKE_TEST")
        NUM_FANTASIES = 128 if not SMOKE_TEST else 4
        NUM_RESTARTS = 10 if not SMOKE_TEST else 2
        RAW_SAMPLES = 128
        bounds = torch.stack([torch.zeros(1, dtype=torch.double), torch.ones(1, dtype=torch.double)])
        acq_fct = qKnowledgeGradient(model=gp, num_fantasies=NUM_FANTASIES)
        candidates, acq_value = optimize_acqf(
            acq_function=acq_fct,
            bounds=bounds,
            q=1,
            num_restarts=NUM_RESTARTS,
            raw_samples=RAW_SAMPLES,
        )

        argmax_pmean, max_pmean = optimize_acqf(
            acq_function=PosteriorMean(gp),
            bounds=bounds,
            q=1,
            num_restarts=NUM_RESTARTS,
            raw_samples=RAW_SAMPLES,
        )
        qKG_proper = qKnowledgeGradient(
            gp,
            num_fantasies=NUM_FANTASIES,
            sampler=acq_fct.sampler,
            current_value=max_pmean,
        )

        candidate, acq_value_proper = optimize_acqf(
            acq_function=qKG_proper,
            bounds=bounds,
            q=1,
            num_restarts=NUM_RESTARTS,
            raw_samples=RAW_SAMPLES,
        )
    else:
        candidate, acq_value = optimize_acqf(
            acq_function=acq_fct,
            bounds=torch.tensor([[0], [1]], dtype=torch.double),
            q=1,
            num_restarts=20,
            raw_samples=256,
        )

    new_y = torch.tensor([[simulation(candidate[0]*(upper_bound-lower_bound)+lower_bound)]], dtype=torch.double)
    new_yvar = torch.full_like(new_y, fixed_Yvar, dtype=torch.double)

    with open("BayesOpt_outputs"+global_individuality_parameter+".csv", "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([candidate.numpy()[0,0]*(upper_bound-lower_bound)+lower_bound, new_y.numpy()[0,0]])

    initial_x = torch.cat([initial_x, candidate])
    initial_y = torch.cat([initial_y, new_y])
    initial_yvar = torch.cat([initial_yvar, new_yvar])
    gp = CustomSingleTaskGP(initial_x, initial_y)
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    x_query = torch.linspace(0, 1, 1000).unsqueeze(-1)
    posterior = gp.posterior(x_query)

    mean = posterior.mean.squeeze(-1).detach().numpy()
    variance = posterior.variance.squeeze(-1)
    stddev = torch.sqrt(variance).detach().numpy()

    print("Lengthscale:", gp.covar_module.base_kernel.lengthscale.item())
    #print("Outputscale:", gp.covar_module.outputscale.item())
    print("Noise:", gp.likelihood.noise.mean().item())

    if visualize:
        plt.scatter(initial_x.numpy()*(upper_bound-lower_bound)+lower_bound, initial_y.numpy(), color="red", label="Trials", zorder=3)
        plt.plot(initial_x.numpy()*(upper_bound-lower_bound)+lower_bound, initial_y.numpy(), color="red", linestyle="", markersize=3)
        plt.plot(x_query*(upper_bound-lower_bound)+lower_bound, mean)
        plt.scatter(candidate.numpy()*(upper_bound-lower_bound)+lower_bound, new_y.numpy(), color="green", s=30, zorder=5, label="New query point")
        plt.fill_between(x_query.numpy().squeeze()*(upper_bound-lower_bound)+lower_bound, mean - 2 * stddev, mean + 2 * stddev, alpha=0.3, label="GP 95% CI")
        plt.xlabel("prestretch force")
        plt.ylabel("contraction of muscle")
        plt.title("Optimization Process")
        plt.legend()
        plt.show()

    counter += 1

    if stopping_y:
        current_value = new_y.item()
        if current_value > best_value + improvement_threshold:
            best_value = current_value
            no_improvement_trials = 0
        elif len(initial_x) > num_initial_trials:
            no_improvement_trials += 1
        if no_improvement_trials >= num_consecutive_trials:
            print(f"Trial {i + 1 + num_initial_trials}: x = {candidate.item()*(upper_bound-lower_bound)+lower_bound}, Value = {current_value}, Best Value = {best_value}")
            print("Stopping criterion met. No significant improvement for consecutive trials.")
            print("Number of total trials: ", i+1+num_initial_trials)
            break
    elif stopping_xy:
        max_index = torch.argmax(initial_y)
        for k in range(len(initial_x)):
            number_x_in_epsilon_neighborhood = 0
            breaking = False
            max_y_in_range = False
            for j in range(len(initial_x)):
                if np.abs(initial_x[k,0].numpy() - initial_x[j,0].numpy()) < x_range:
                    number_x_in_epsilon_neighborhood += 1
                    if initial_x[max_index,0].numpy() == initial_x[k,0].numpy() or initial_x[max_index,0].numpy() == initial_x[j,0].numpy():
                        max_y_in_range = True
            if number_x_in_epsilon_neighborhood >= num_consecutive_trials and max_y_in_range:
                print("Stopping criterion met. No significant improvement for consecutive trials.")
                print("Number of total trials: ", i+1+num_initial_trials)
                breaking = True
                break
        if breaking:
            break
    else:
        print("Wrong input, used stopping_y instead.")
        stopping_y = True

    current_value = new_y.item()
    if current_value > best_value + improvement_threshold:
        best_value = current_value

    print(f"Trial {i + 1 + num_initial_trials}: x = {candidate.item()*(upper_bound-lower_bound)+lower_bound}, Value = {current_value}, Best Value = {best_value}")

    
x_query = torch.linspace(0, 1, 1000).unsqueeze(-1)
posterior = gp.posterior(x_query)

mean = posterior.mean.squeeze(-1).detach().numpy()
variance = posterior.variance.squeeze(-1)
stddev = torch.sqrt(variance).detach().numpy()

if visualize:
    plt.scatter(initial_x.numpy()*(upper_bound-lower_bound)+lower_bound, initial_y.numpy(), color="red", label="Trials", zorder=3)
    plt.plot(initial_x.numpy()*(upper_bound-lower_bound)+lower_bound, initial_y.numpy(), color="red", linestyle="", markersize=3)
    plt.plot(x_query*(upper_bound-lower_bound)+lower_bound, mean)
    plt.fill_between(x_query.numpy().squeeze()*(upper_bound-lower_bound)+lower_bound, mean - 2 * stddev, mean + 2 * stddev, alpha=0.3, label="GP 95% CI")
    plt.xlabel("prestretch force")
    plt.ylabel("contraction of muscle")
    plt.title("Optimization Results")
    plt.legend()
    plt.show()

if add_points:
    continuing = input("Do you want to add another query point? (y/n)")
else:
    continuing = "n"

while continuing == "y":
    candidate = input("Which point do you want to add?")
    candidate = torch.tensor([[float(candidate)]], dtype=torch.double)

    new_y = torch.tensor([[simulation(candidate[0]*(upper_bound-lower_bound)+lower_bound)]], dtype=torch.double)
    new_yvar = torch.full_like(new_y, fixed_Yvar, dtype=torch.double)

    with open("BayesOpt_outputs"+global_individuality_parameter+".csv", "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([candidate.numpy()[0,0]*(upper_bound-lower_bound)+lower_bound, new_y.numpy()[0,0]])

    initial_x = torch.cat([initial_x, candidate])
    initial_y = torch.cat([initial_y, new_y])
    initial_yvar = torch.cat([initial_yvar, new_yvar])
    gp = CustomSingleTaskGP(initial_x, initial_y)
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    counter += 1

    print(f"Trial {i + 1 + num_initial_trials}: x = {candidate.item()}, Value = {current_value}, Best Value = {best_value}")

    x_query = torch.linspace(0, 1, 1000).unsqueeze(-1)
    posterior = gp.posterior(x_query)

    mean = posterior.mean.squeeze(-1).detach().numpy()
    variance = posterior.variance.squeeze(-1)
    stddev = torch.sqrt(variance).detach().numpy()

    if visualize:
        plt.scatter(initial_x.numpy()*(upper_bound-lower_bound)+lower_bound, initial_y.numpy(), color="red", label="Trials", zorder=3)
        plt.plot(initial_x.numpy()*(upper_bound-lower_bound)+lower_bound, initial_y.numpy(), color="red", linestyle="", markersize=3)
        plt.plot(x_query*(upper_bound-lower_bound)+lower_bound, mean)
        plt.scatter(candidate.numpy()*(upper_bound-lower_bound)+lower_bound, new_y.numpy(), color="green", s=30, zorder=5, label="New query point")
        plt.fill_between(x_query.numpy().squeeze()*(upper_bound-lower_bound)+lower_bound, mean - 2 * stddev, mean + 2 * stddev, alpha=0.3, label="GP 95% CI")
        plt.xlabel("prestretch force")
        plt.ylabel("contraction of muscle")
        plt.title("Optimization Process")
        plt.legend()
        plt.show()

    continuing = input("Do you want to add another query point? (y/n)")

max_index = torch.argmax(initial_y)
maximizer = initial_x[max_index]
best_y = initial_y[max_index]

with open("BayesOpt_outputs"+global_individuality_parameter+".csv", "a", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(np.linspace(lower_bound, upper_bound, 1000))
    writer.writerow(mean)
    writer.writerow(stddev)
    writer.writerow([counter])
    writer.writerow([maximizer.numpy()[0]*(upper_bound-lower_bound)+lower_bound, best_y.numpy()[0]])
    writer.writerow([time.time()-starting_time])

print(global_individuality_parameter)

with open("BayesOpt_global_individuality_parameters.csv", "a") as f:
    writer = csv.writer(f)
    writer.writerow([global_individuality_parameter])