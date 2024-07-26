# main.py

import numpy as np
from link import Link, FindMother, FindChildren, ForwardKinematics, RPY2R, IK_leg

# Setup Biped Robot
uLINK = []

UX = np.array([1, 0, 0])
UY = np.array([0, 1, 0])
UZ = np.array([0, 0, 1])

uLINK.append(Link('BODY', mass=10, sister=-1, child=1, b=np.array([0, 0, 0.6]), a=UZ, q=0))
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

# uLINK[BODY].p = np.array([0.0, 0.0, 0.6])
# uLINK[BODY].R = np.eye(3)

def print_joint_angles():
    print("\nJoint angles (in radians):")
    for i in range(6):
        print(f"RLEG_J{i}: {uLINK[RLEG_J0 + i].q}")
        print(f"LLEG_J{i}: {uLINK[LLEG_J0 + i].q}")

# Raise the right leg by 5 cm
target_foot_pos = uLINK[RLEG_J5].p + np.array([0, -0.1, -0.6])  # Raise by 5 cm
Rfoot = Link('Rfoot', p=target_foot_pos, R=uLINK[RLEG_J5].R)

# Calculate inverse kinematics for the right leg
qR = IK_leg(uLINK[BODY], -0.1, 0.3, 0.3, Rfoot)

# Apply the calculated joint angles
for n in range(6):
    uLINK[RLEG_J0 + n].q = qR[n]

# Perform forward kinematics to update link positions
ForwardKinematics(uLINK, 0)

# Print the resulting joint angles
print_joint_angles()