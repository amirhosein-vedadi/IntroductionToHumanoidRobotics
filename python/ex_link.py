import numpy as np
from link import Link, FindMother, FindChildren, ForwardKinematics

if __name__ == "__main__":
    # Create a list to hold the Link objects
    uLINK = []

    # Constants
    ToDeg = 180 / np.pi
    ToRad = np.pi / 180
    UX = np.array([1, 0, 0])
    UY = np.array([0, 1, 0])
    UZ = np.array([0, 0, 1])

    # Add Link objects to the list
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

    # Find mother link from sister and child data
    FindMother(uLINK, 0)

    # Find children
    FindChildren(uLINK)

    # Substitute the ID to the link name variables. For example, BODY=0.
    for n, link in enumerate(uLINK):
        globals()[link.name] = n

    uLINK[BODY].p = np.array([0.0, 0.0, 0.65])
    uLINK[BODY].R = np.eye(3)

    # Initialize velocities
    uLINK[BODY].v = np.zeros(3)
    uLINK[BODY].w = np.zeros(3)
    for link in uLINK:
        link.dq = 0  # Joint speed [rad/s]

    # Print the link information
    print('[[[[[ uLINK struct was set as following ]]]]]')
    print('-------------------------------------')
    print('ID     name    sister child   mass')
    print('-------------------------------------')
    for i, link in enumerate(uLINK):
        print(f"{i:2d}  {link}")

    # Print the route from BODY to LLEG_J5
    route = FindRoute(uLINK, LLEG_J5)
    print("\nRoute from BODY to LLEG_J5:")
    print([uLINK[i].name for i in route])