import subprocess
import os
import shlex
import csv
import time

"""
opendihu_case = os.path.join(os.environ["OPENDIHU_HOME"],"examples/electrophysiology/fibers/fibers_contraction/no_precice/cuboid_muscle/build_release")
os.chdir(opendihu_case)
subprocess.run("pwd")
#os.chdir("OpenDiHu/opendihu/examples/electrophysiology/fibers/fibers_contraction/no_precice/cuboid_muscle/build_release")
"""
os.chdir("build_release")

force = 20.0

for i in range(1):
    individuality_parameter = str(int(time.time()))+str(force)
    command = shlex.split(f"./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin {force} {individuality_parameter}")
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

    f = open("contractions_20N_20ms.csv", "a")
    f.write(str(contraction))
    f.write(",")
    f.close()
