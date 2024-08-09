import subprocess
import os
import shlex
import csv

"""
opendihu_case = os.path.join(os.environ["OPENDIHU_HOME"],"examples/electrophysiology/fibers/fibers_contraction/no_precice/cuboid_muscle/build_release")
os.chdir(opendihu_case)
subprocess.run("pwd")
#os.chdir("OpenDiHu/opendihu/examples/electrophysiology/fibers/fibers_contraction/no_precice/cuboid_muscle/build_release")
"""
os.chdir("build_release")

command = shlex.split("./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin 0.0")
#command = shlex.split("./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin 4.5") #5%
#command = shlex.split("./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin 8.2") #10%
#command = shlex.split("./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin 15.93") #25%

#command = shlex.split("./nearly_incompressible_mooney_rivlin ../settings_force.py nearly_incompressible_mooney_rivlin 15.93")
#command = shlex.split("./compressible_mooney_rivlin ../settings_force.py compressible_mooney_rivlin 15.9118")
#command = shlex.split("./muscle_contraction ../settings_muscle_contraction.py variables.py")
subprocess.run(command)

"""
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

f = open("contractions.csv", "a")
f.write(str(contraction))
f.write(",")
f.close()
"""

#~/Desktop/Bachelorthesis/OpenDiHu/opendihu/examples/electrophysiology/fibers/fibers_contraction/no_precice/cuboid_muscle/build_release
