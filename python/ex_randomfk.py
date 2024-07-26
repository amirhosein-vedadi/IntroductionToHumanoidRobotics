# main.py

import numpy as np
from link import Link, FindMother, FindChildren, ForwardKinematics

# Setup Biped Robot
uLINK = []

UX = np.array([1, 0, 0])
UY = np.array([0, 1, 0])
UZ = np.array([0, 0, 1])

uLINK.append(Link('BODY', mass=10, sister=-1, child=1, b=np.array([0, 0, 0.7]), a=UZ, q=0))
uLINK.append(Link('RLEG_J0', mass=5, sister=7, child=2, b=np.array([0, -0.1, 0]), a=UZ, q=0))
uLINK.append(Link('RLEG_J1', mass=1, sister=-1, child=3, b=np.array([0, 0, 0]), a=UX, q=0))
uLINK.append(Link('RLEG_J2', mass=5, sister=-1, child=4, b=np.array([0, 0, 0]), a=UY, q=0))
uLINK.append(Link('RLEG_J3', mass=1, sister=-1, child=5, b=np.array([0, 0, -0.3]), a=UY, q=0))
uLINK.append(Link('RLEG_J4', mass=6, sister=-1, child=6, b=np.array([0, 0, -0.3]), a=UY, q=0))
uLINK.append(Link('RLEG_J5', mass=2, sister=-1, child=-1, b=np.array([0, 0, 0]), a=UX, q=0))
uLINK.append(Link('LLEG_J0', mass=5, sister=-1, child=8, b=np.array([0, 0.1, 0]), a=UZ, q=0))
uLINK.append(Link('LLEG_J1', mass=1, sister=-1, child=9, b=np.array([0, 0, 0]), a=UX, q=0))
uLINK.append(Link('LLEG_J2', mass=5, sister=-1, child=10, b=np.array([0, 0, 0]), a=UY, q=0))
uLINK.append(Link('LLEG_J3', mass=1, sister=-1, child=11, b=np.array([0, 0, -0.3]), a=UY, q=0))
uLINK.append(Link('LLEG_J4', mass=6, sister=-1, child=12, b=np.array([0, 0, -0.3]), a=UY, q=0))
uLINK.append(Link('LLEG_J5', mass=2, sister=-1, child=-1, b=np.array([0, 0, 0]), a=UX, q=0))

FindMother(uLINK, 0)
FindChildren(uLINK)

for n, link in enumerate(uLINK):
    globals()[link.name] = n

uLINK[BODY].p = np.array([0.0, 0.0, 0.65])
uLINK[BODY].R = np.eye(3)

np.random.seed(42)

def print_link_positions():
    print("\nLink positions after forward kinematics:")
    for i, link in enumerate(uLINK):
        print(f"{link.name}: {link.p}")

# Generate random joint angles and perform forward kinematics
print("\n--- New Random Configuration ---")

qR1 = 2/3 * np.pi * (np.random.rand(6) - 0.5)
qR1[3] = np.pi * np.random.rand()

qL1 = np.pi * (np.random.rand(6) - 0.5)
qL1[3] = np.pi * np.random.rand()

print("Random joint angles (in radians):")
print("Right leg:", qR1)
print("Left leg:", qL1)

for n in range(6):
    uLINK[RLEG_J0 + n].q = qR1[n]
    uLINK[LLEG_J0 + n].q = qL1[n]

uLINK[BODY].p = np.array([0.0, 0.0, 0.7])
uLINK[BODY].R = np.eye(3)
ForwardKinematics(uLINK, 0)

print_link_positions()