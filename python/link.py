import numpy as np

class Link:
    def __init__(self, name, mass=0, sister=-1, child=-1, mother=-1, 
                 p=None, R=None, v=None, w=None, q=0, dq=0, ddq=0, 
                 a=None, b=None, vertex=None, face=None, c=None, I=None):
        self.name = name
        self.mass = mass
        self.sister = sister
        self.child = child
        self.mother = mother
        self.p = p if p is not None else np.zeros(3)
        self.R = R if R is not None else np.eye(3)
        self.v = v if v is not None else np.zeros(3)
        self.w = w if w is not None else np.zeros(3)
        self.q = q
        self.dq = dq
        self.ddq = ddq
        self.a = a if a is not None else np.zeros(3)
        self.b = b if b is not None else np.zeros(3)
        self.vertex = vertex
        self.face = face
        self.c = c if c is not None else np.zeros(3)
        self.I = I if I is not None else np.zeros((3, 3))

    def __str__(self):
        return f"{self.name:8s} {self.sister:2d} {self.child:2d} {self.mass:7.2f}"

    @classmethod
    def print_link_tree(cls, links, j):
        if j != -1:
            print(f"j={j:2d} : {links[j].name}")  # Print current link's name
            cls.print_link_tree(links, links[j].child)  # Print child's name
            cls.print_link_tree(links, links[j].sister)  # Print sister's name

    @classmethod
    def total_mass(cls, links, j):
        if j == -1:
            return 0
        else:
            return (links[j].mass + 
                    cls.total_mass(links, links[j].sister) + 
                    cls.total_mass(links, links[j].child))

def FindMother(links, j):
    if j != -1:
        if j == 0:
            links[j].mother = -1
        if links[j].child != -1:
            links[links[j].child].mother = j
            FindMother(links, links[j].child)
        if links[j].sister != -1:
            links[links[j].sister].mother = links[j].mother
            FindMother(links, links[j].sister)

def FindRoute(links, to):
    # return the list of joint numbers connecting ROOT to 'to'
    i = links[to].mother
    if i == 0:
        return [to]
    else:
        return FindRoute(links, i) + [to]
    
def FindChildren(links):
    num = len(links)

    for n in range(num):
        links[n].sister = -1
        links[n].child = -1

    for n in range(num):
        mom = links[n].mother
        if mom != -1:
            if links[mom].child == -1:
                links[mom].child = n  # I am the first child.
            else:
                eldest_sister = links[mom].child
                SetYoungestSister(links, eldest_sister, n)

def SetYoungestSister(links, eldest, youngest):
    if links[eldest].sister == -1:
        links[eldest].sister = youngest
    else:
        SetYoungestSister(links, links[eldest].sister, youngest)

def Rodrigues(w, dt):
    """
    Compute the rotation matrix from angular velocity vector and time step
    w: angular velocity vector (3D vector)
    dt: time step
    """
    norm_w = np.linalg.norm(w)
    if norm_w < np.finfo(float).eps:
        return np.eye(3)
    else:
        wn = w / norm_w  # rotation axis (unit vector)
        th = norm_w * dt  # amount of rotation (rad)
        w_wedge = np.array([
            [0, -wn[2], wn[1]],
            [wn[2], 0, -wn[0]],
            [-wn[1], wn[0], 0]
        ])
        R = np.eye(3) + w_wedge * np.sin(th) + np.dot(w_wedge, w_wedge) * (1 - np.cos(th))
        return R
    
def ForwardKinematics(links, j):
    if j == -1:
        return
    if j != 0:
        mom = links[j].mother
        links[j].p = links[mom].R @ links[j].b + links[mom].p
        links[j].R = links[mom].R @ Rodrigues(links[j].a, links[j].q)
    
    ForwardKinematics(links, links[j].sister)
    ForwardKinematics(links, links[j].child)

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

    # Perform forward kinematics
    ForwardKinematics(uLINK, 0)

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

    # Print positions after forward kinematics
    print("\nPositions after forward kinematics:")
    for i, link in enumerate(uLINK):
        print(f"{link.name}: {link.p}")