# This settings file can be used for two different equations:
# - Isotropic hyperelastic material
# - Linear elasticity
#
# arguments: <scenario_name> <force>


import numpy as np
import sys, os

script_path = os.path.dirname(os.path.abspath(__file__))
var_path = os.path.join(script_path, "variables")
sys.path.insert(0, var_path)

import variables

n_ranks = (int)(sys.argv[-1])

# parameters
force = 0.0                       # [N] load on top
material_parameters = [3.176e-10, 1.813, 1.075e-2, 1.0]     # [c1, c2, b, d]
physical_extent = [3.0, 3.0, 12.0]
constant_body_force = None                                                                      
scenario_name = "tensile_test"
dirichlet_bc_mode = "fix_floating"                                                              
 
if len(sys.argv) > 3:                                                                           
  scenario_name = sys.argv[0]
  force = float(sys.argv[1])
  print("scenario_name: {}".format(scenario_name))
  print("force: {}".format(force))
    
  # set material parameters depending on scenario name
  if scenario_name == "compressible_mooney_rivlin":
    material_parameters = [3.176e-10, 1.813, 10]      # c1, c2, c
    
  elif scenario_name == "compressible_mooney_rivlin_decoupled":
    material_parameters = [3.176e-10, 1.813, 10.0]      # c1, c2, kappa
    
  elif scenario_name == "incompressible_mooney_rivlin":
    material_parameters = [3.176e-10, 1.813]      # c1, c2
    
  elif scenario_name == "nearly_incompressible_mooney_rivlin":
    material_parameters = [3.176e-10, 1.813, 1e3]      # c1, c2, kappa

  elif scenario_name == "nearly_incompressible_mooney_rivlin_decoupled":
    material_parameters = [3.176e-10, 1.813, 1e3]      # c1, c2, kappa

  elif scenario_name == "linear":
    pass

  elif scenario_name == "nearly_incompressible_mooney_rivlin_febio":
    material_parameters = [3.176e-10, 1.813, 1e3]      # c1, c2, kappa

  else:
    print("Error! Please specify the correct scenario, see settings.py for allowed values.\n")
    quit()

nx, ny, nz = 3, 3, 12                     # number of elements
mx, my, mz = 2*nx+1, 2*ny+1, 2*nz+1 # quadratic basis functions

fb_x, fb_y = 10, 10         # number of fibers
fb_points = 100             # number of points per fiber
fiber_direction = [0, 0, 1] # direction of fiber in element

def get_fiber_no(fiber_x, fiber_y):
    return fiber_x + fiber_y*fb_x

meshes = { # create 3D mechanics mesh
    "mesh3D": {
        "nElements":            [nx, ny, nz],
        "physicalExtent":       [nx, ny, nz],
        "physicalOffset":       [0, 0, 0],
        "logKey":               "mesh3D",
        "inputMeshIsGlobal":    True,
        "nRanks":               n_ranks
    },
    "3Dmesh_quadratic": { 
      "inputMeshIsGlobal":          True,                       # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      "nElements":                  [nx, ny, nz],               # number of quadratic elements in x, y and z direction
      "physicalExtent":             physical_extent,            # physical size of the box
      "physicalOffset":             [0, 0, 0],                  # offset/translation where the whole mesh begins
    },
    "3Dmesh_febio": { 
      "inputMeshIsGlobal":          True,                       # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      "nElements":                  [2*nx, 2*ny, 2*nz],               # number of quadratic elements in x, y and z direction
      "physicalExtent":             physical_extent,            # physical size of the box
      "physicalOffset":             [0, 0, 0],                  # offset/translation where the whole mesh begins
    }
}

for fiber_x in range(fb_x):
    for fiber_y in range(fb_y):
        fiber_no = get_fiber_no(fiber_x, fiber_y)
        x = nx * fiber_x / (fb_x - 1)
        y = ny * fiber_y / (fb_y - 1)
        nodePositions = [[x, y, nz * i / (fb_points - 1)] for i in range(fb_points)]
        meshName = "fiber{}".format(fiber_no)
        meshes[meshName] = { # create fiber meshes
            "nElements":            [fb_points - 1],
            "nodePositions":        nodePositions,
            "inputMeshIsGlobal":    True,
            "nRanks":               n_ranks
        }

# set Dirichlet BC, fix bottom
elasticity_dirichlet_bc = {}
k = 0

# fix z value on the whole x-y-plane
for j in range(my):
  for i in range(mx):
    elasticity_dirichlet_bc[k*mx*my + j*mx + i] = [None,None,0.0,None,None,None]

# fix left edge 
for j in range(my):
  elasticity_dirichlet_bc[k*mx*my + j*mx + 0][0] = 0.0
  
# fix front edge 
for i in range(mx):
  elasticity_dirichlet_bc[k*mx*my + 0*mx + i][1] = 0.0
       
# set Neumann BC, set traction at the top
k = nz-1
traction_vector = [0, 0, force]     # the traction force in specified in the reference configuration

elasticity_neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": traction_vector, "face": "2+"} for j in range(ny) for i in range(nx)]

# callback for result
def handle_result_hyperelasticity(result):
  data = result[0]

  number_of_nodes = mx * my
  average_z_start = 0
  average_z_end = 0

  z_data = data["data"][0]["components"][2]["values"]

  for i in range(number_of_nodes):
    average_z_start += z_data[i]
    average_z_end += z_data[number_of_nodes*(mz -1) + i]

  average_z_start /= number_of_nodes
  average_z_end /= number_of_nodes

  length_of_muscle = np.abs(average_z_end - average_z_start)
  print("length of muscle: ", length_of_muscle)

  if data["timeStepNo"] == 0:
    f = open("muscle_length_prestretch.csv", "w")
    f.write(str(length_of_muscle))
    f.write(",")
    f.close()
  else:
    f = open("muscle_length_prestretch.csv", "a")
    f.write(str(length_of_muscle))
    f.write(",")
    f.close()
  
  if data["timeStepNo"] == 1:
    field_variables = data["data"]
    
    strain = max(field_variables[1]["components"][2]["values"])
    stress = max(field_variables[5]["components"][2]["values"])
    
    print("strain: {}, stress: {}".format(strain, stress))
    
    with open("result.csv","a") as f:
      f.write("{},{},{}\n".format(scenario_name,strain,stress))

# callback for result
def handle_result_febio(result):
  data = result[0]
  
  if data["timeStepNo"] == 1:
    field_variables = data["data"]
    
    strain = max(field_variables[2]["components"][2]["values"])
    stress = max(field_variables[5]["components"][2]["values"])
    
    print("strain: {}, stress: {}".format(strain, stress))
    
    with open("result.csv","a") as f:
      f.write("{},{},{}\n".format(scenario_name,strain,stress))

# callback for result
def handle_result_linear_elasticity(result):
  data = result[0]
  
  if data["timeStepNo"] == -1:
    field_variables = data["data"]
    
    strain = max(field_variables[1]["components"][2]["values"])
    K = 50    # parameters as given in config
    mu = 100
    stress = (K + 4./3*mu) * strain
    
    print("strain: {}, stress: {}".format(strain, stress))
    
    with open("result.csv","a") as f:
      f.write("{},{},{}\n".format(scenario_name,strain,stress))

def callback_function_contraction(raw_data):
  t = raw_data[0]["currentTime"]
  if True:
    number_of_nodes = variables.bs_x * variables.bs_y
    average_z_start = 0
    average_z_end = 0

    z_data = raw_data[0]["data"][0]["components"][2]["values"]

    for i in range(number_of_nodes):
      average_z_start += z_data[i]
      average_z_end += z_data[number_of_nodes*(variables.bs_z -1) + i]

    average_z_start /= number_of_nodes
    average_z_end /= number_of_nodes

    length_of_muscle = np.abs(average_z_end - average_z_start)
    print("length of muscle: ", length_of_muscle)

    if t == variables.dt_3D:
      f = open("muscle_length_contraction.csv", "w")
      f.write(str(length_of_muscle))
      f.write(",")
      f.close()
    else:
      f = open("muscle_length_contraction.csv", "a")
      f.write(str(length_of_muscle))
      f.write(",")
      f.close()


config = {
  "scenarioName":                 scenario_name,                # scenario name to identify the simulation runs in the log file
  "logFormat":                    "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":   "solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile": "mappings_between_meshes_log.txt",    # log file for mappings 

  "Meshes": meshes,
  "MappingsBetweenMeshes": { 
    "mesh3D" : ["fiber{}".format(variables.get_fiber_no(fiber_x, fiber_y)) for fiber_x in range(variables.fb_x) for fiber_y in range(variables.fb_y)]
  },

  "Solvers": {
    "linearElasticitySolver": {           # solver for linear elasticity
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    ,
      "maxIterations":      1e4,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    }, 
    "diffusionSolver": {
      "solverType":                     "cg",
      "preconditionerType":             "none",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    },
    "mechanicsSolver": {
      "solverType":                     "preonly",
      "preconditionerType":             "lu",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "snesLineSearchType":             "l2",
      "snesRelativeTolerance":          1e-5,
      "snesAbsoluteTolerance":          1e-5,
      "snesMaxIterations":              10,
      "snesMaxFunctionEvaluations":     1e8,
      "snesRebuildJacobianFrequency":   5,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    }
  },

  "Coupling": {

    #I dont know what to put here

    "timeStepWidth": variables.end_time,
    "endTime": variables.end_time,
    "connectedSlotsTerm1To2": {},

    "Term1": {
      "HyperelasticitySolver": {
        "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
        
        "materialParameters":         material_parameters,          # material parameters of the Mooney-Rivlin material
        "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
        "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
        "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
        "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
          
        "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
        # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
        
        # mesh
        "meshName":                   "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
        "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
        
        #"fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
        #"fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
        
        # nonlinear solver
        "relativeTolerance":          1e-5,                         # 1e-10 relative tolerance of the linear solver
        "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
        "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
        "preconditionerType":         "lu",                         # type of the preconditioner
        "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
        "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
        "snesMaxIterations":          100,                           # maximum number of iterations in the nonlinear solver
        "snesRelativeTolerance":      1e-5,                         # relative tolerance of the nonlinear solver
        "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
        "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
        "snesRebuildJacobianFrequency": 1,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
        
        #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
        "dumpFilename":               "",                           # dump disabled
        "dumpFormat":                 "default",                     # default, ascii, matlab
        
        #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
        #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
        "loadFactors":                [],                           # no load factors, solve problem directly
        "loadFactorGiveUpThreshold":    0.1,                        # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
        "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
        
        # boundary and initial conditions
        "dirichletBoundaryConditions": elasticity_dirichlet_bc,             # the initial Dirichlet boundary conditions that define values for displacements u
        "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        "neumannBoundaryConditions":   elasticity_neumann_bc,               # Neumann boundary conditions that define traction forces on surfaces of elements
        "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
        "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
        "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
        
        "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
        "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
        "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
        "constantBodyForce":           constant_body_force,                 # a constant force that acts on the whole body, e.g. for gravity
        
        "dirichletOutputFilename":      "out/"+scenario_name+"/dirichlet_boundary_conditions",           # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        
        # define which file formats should be written
        # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
        "OutputWriter" : [
          
          # Paraview files
          {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/prestretch", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          
          # Python files and callback
          {"format": "PythonFile", "outputInterval": 1, "filename": "out/all/"+scenario_name, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          {"format": "PythonCallback", "outputInterval": 1, "filename": "out/all/"+scenario_name, "callback": handle_result_hyperelasticity, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ],
        # 2. additional output writer that writes also the hydrostatic pressure
        "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
          "OutputWriter" : [
            #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
        },
        # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
        "LoadIncrements": {   
          "OutputWriter" : [
            #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
        },
      },
    },
    "Term2": {
      "Coupling": {
        "timeStepWidth":            variables.dt_3D,
        "logTimeStepWidthAsKey":    "dt_3D",
        "durationLogKey":           "duration_3D",
        "endTime":                  variables.end_time,
        "timeStepOutputInterval":   1,
        "connectedSlotsTerm1To2":   {1:2},  # transfer stress to MuscleContractionSolver gamma
        "connectedSlotsTerm2To1":   None,   # transfer nothing back

        "Term1": { # fibers (FastMonodomainSolver)
          "MultipleInstances": { 
            "ranksAllComputedInstances":    list(range(n_ranks)),
            "nInstances":                   1,

            "instances": [{
              "ranks": [0],

              "StrangSplitting": {
                "timeStepWidth":            variables.dt_splitting,
                "logTimeStepWidthAsKey":    "dt_splitting",
                "durationLogKey":           "duration_splitting",
                "timeStepOutputInterval":   100,
                "connectedSlotsTerm1To2":   None,
                "connectedSlotsTerm2To1":   None,

                "Term1": { # reaction term
                  "MultipleInstances": {
                    "nInstances":   variables.fb_x * variables.fb_y,

                    "instances": [{
                      "ranks": [0],

                      "Heun": {
                        "timeStepWidth":            variables.dt_0D,
                        "logTimeStepWidthAsKey":    "dt_0D",
                        "durationLogKey":           "duration_0D",
                        "timeStepOutputInterval":   100,

                        "initialValues":                [],
                        "dirichletBoundaryConditions":  {},
                        "dirichletOutputFilename":      None,
                        "inputMeshIsGlobal":            True,
                        "checkForNanInf":               False,
                        "nAdditionalFieldVariables":    0,
                        "additionalSlotNames":          [],
                        "OutputWriter":                 [],

                        "CellML": {
                          "modelFilename":          variables.input_dir + "hodgkin_huxley-razumova.cellml",
                          "meshName":               "fiber{}".format(variables.get_fiber_no(fiber_x, fiber_y)), 
                          "stimulationLogFilename": "out/" + scenario_name + "stimulation.log",

                          "statesInitialValues":                        [],
                          "initializeStatesToEquilibrium":              False,
                          "initializeStatesToEquilibriumTimeStepWidth": 1e-4,
                          "optimizationType":                           "vc",
                          "approximateExponentialFunction":             True,
                          "compilerFlags":                              "-fPIC -O3 -march=native -Wno-deprecated_declarations -shared",
                          "maximumNumberOfThreads":                     0,

                          "setSpecificStatesCallEnableBegin":       variables.specific_states_call_enable_begin,
                          "setSpecificStatesCallFrequency":         variables.specific_states_call_frequency,
                          "setSpecificStatesRepeatAfterFirstCall":  0.01,
                          "setSpecificStatesFrequencyJitter":       [0] ,
                          "setSpecificStatesCallInterval":          0,
                          "setSpecificStatesFunction":              None,
                          "additionalArgument":                     None, 

                          "mappings": {
                            ("parameter", 0):               "membrane/i_Stim",
                            ("parameter", 1):               "Razumova/l_hs",
                            ("parameter", 2):               ("constant", "Razumova/rel_velo"),
                            ("connectorSlot", "vm"):        "membrane/V",
                            ("connectorSlot", "stress"):    "Razumova/activestress",
                            ("connectorSlot", "alpha"):     "Razumova/activation",
                            ("connectorSlot", "lambda"):    "Razumova/l_hs",
                            ("connectorSlot", "ldot"):      "Razumova/rel_velo"
                          },
                          "parametersInitialValues": [0.0, 1.0, 0.0],
                        },
                      }
                    } for fiber_x in range(variables.fb_x) for fiber_y in range(variables.fb_y)] 
                  }
                },

                "Term2": { # diffusion term
                  "MultipleInstances": {
                    "nInstances": variables.fb_x * variables.fb_y, 

                    "OutputWriter": [
                      {
                        "format":             "Paraview",
                        "outputInterval":     int(1.0 / variables.dt_3D * variables.output_interval),
                        "filename":           "out/" + scenario_name + "/fibers",
                        "fileNumbering":      "incremental",
                        "binary":             True,
                        "fixedFormat":        False,
                        "onlyNodalValues":    True,
                        "combineFiles":       True
                      }
                    ],

                    "instances": [{
                      "ranks": [0],

                      "ImplicitEuler": {
                        "timeStepWidth":            variables.dt_1D,
                        "logTimeStepWidthAsKey":    "dt_1D",
                        "durationLogKey":           "duration_1D",
                        "timeStepOutputInterval":   100,

                        "nAdditionalFieldVariables":    4,
                        "additionalSlotNames":          ["stress", "alpha", "lambda", "ldot"],

                        "solverName":                       "diffusionSolver",
                        "timeStepWidthRelativeTolerance":   1e-10,

                        "dirichletBoundaryConditions":      {},
                        "dirichletOutputFilename":          None,
                        "inputMeshIsGlobal":                True,
                        "checkForNanInf":                   False,
                        "OutputWriter":                     [],

                        "FiniteElementMethod": {
                          "meshName":           "fiber{}".format(variables.get_fiber_no(fiber_x, fiber_y)),
                          "inputMeshIsGlobal":  True,
                          "solverName":         "diffusionSolver",
                          "prefactor":          variables.diffusion_prefactor,
                          "slotName":           "vm"
                        }
                      }
                    } for fiber_x in range(variables.fb_x) for fiber_y in range(variables.fb_y)]
                  }
                }
              }
            }]
          },

          "fiberDistributionFile":                              variables.fiber_distribution_file,
          "firingTimesFile":                                    variables.firing_times_file,
          "valueForStimulatedPoint":                            20.0,
          "onlyComputeIfHasBeenStimulated":                     True,
          "disableComputationWhenStatesAreCloseToEquilibrium":  True,
          "neuromuscularJunctionRelativeSize":                  0.1,
          "generateGPUSource":                                  True,
          "useSinglePrecision":                                 False
        },

        "Term2": { # solid mechanics (MuscleContractionSolver)
          "MuscleContractionSolver": {
            "Pmax":                         variables.pmax,
            "slotNames":                    ["lambda", "ldot", "gamma", "T"],
            "dynamic":                      True,

            "numberTimeSteps":              1,
            "timeStepOutputInterval":       100,
            "lambdaDotScalingFactor":       1,
            "enableForceLengthRelation":    True,
            "mapGeometryToMeshes":          [],

            "OutputWriter": [
              {
                "format":             "Paraview",
                "outputInterval":     int(1.0 / variables.dt_3D * variables.output_interval),
                "filename":           "out/" + scenario_name + "/mechanics",
                "fileNumbering":      "incremental",
                "binary":             True,
                "fixedFormat":        False,
                "onlyNodalValues":    True,
                "combineFiles":       True
              }
            ],

            "DynamicHyperelasticitySolver": {
              "durationLogKey":         "duration_3D",
              "logTimeStepWidthAsKey":  "dt_3D",
              "numberTimeSteps":        1,
              "materialParameters":     variables.material_parameters,
              "density":                variables.rho,
              "timeStepOutputInterval": 1,

              "meshName":                   "3Dmesh_quadratic",
              "fiberDirectionInElement":    variables.fiber_direction,
              "inputMeshIsGlobal":          True,
              "fiberMeshNames":             [],
              "fiberDirection":             None,

              "solverName":                 "mechanicsSolver",
              "displacementsScalingFactor":  1.0,
              "useAnalyticJacobian":        True,
              "useNumericJacobian":         False,
              "dumpDenseMatlabVariables":   False,
              "loadFactorGiveUpThreshold":  1,
              "loadFactors":                [],
              "scaleInitialGuess":          False,
              "extrapolateInitialGuess":    True,
              "nNonlinearSolveCalls":       1,

              "dirichletBoundaryConditions":                            {}, #elasticity_dirichlet_bc, #variables.dirichlet_bc,
              "neumannBoundaryConditions":                              {}, #elasticity_neumann_bc, #variables.neumann_bc,
              "updateDirichletBoundaryConditionsFunction":              None,
              "updateDirichletBoundaryConditionsFunctionCallInterval":  1,
              "divideNeumannBoundaryConditionValuesByTotalArea":        True,

              #"initialValuesDisplacements": [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
              "initialValuesVelocities":    [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
              "constantBodyForce":          (0, 0, 0),

              "dirichletOutputFilename":    "out/" + scenario_name + "/dirichlet_output",
              "residualNormLogFilename":    "out/" + scenario_name + "/residual_norm_log.txt",
              "totalForceLogFilename":      "out/" + scenario_name + "/total_force_log.txt",

              "OutputWriter": [
                {
                  "format": "PythonCallback",
                  "callback": callback_function_contraction,
                  "outputInterval": 1,
                }
              ],
              "pressure":       { "OutputWriter": [] },
              "dynamic":        { "OutputWriter": [] },
              "LoadIncrements": { "OutputWriter": [] }
            }
          }
        }
      }
    }
  },
}
