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


def Rroll(angle):
    """Rotation matrix for a roll (rotation around x-axis)."""
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Rpitch(angle):
    """Rotation matrix for a pitch (rotation around y-axis)."""
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def RPY2R(rpy):
    """Convert roll-pitch-yaw angles to a rotation matrix."""
    roll, pitch, yaw = rpy
    R = np.array([
        [np.cos(yaw) * np.cos(pitch), np.cos(yaw) * np.sin(pitch) * np.sin(roll) - np.sin(yaw) * np.cos(roll), np.cos(yaw) * np.sin(pitch) * np.cos(roll) + np.sin(yaw) * np.sin(roll)],
        [np.sin(yaw) * np.cos(pitch), np.sin(yaw) * np.sin(pitch) * np.sin(roll) + np.cos(yaw) * np.cos(roll), np.sin(yaw) * np.sin(pitch) * np.cos(roll) - np.cos(yaw) * np.sin(roll)],
        [-np.sin(pitch), np.cos(pitch) * np.sin(roll), np.cos(pitch) * np.cos(roll)]
    ])
    return R

def IK_leg(Body, D, A, B, Foot):
    r = Foot.R.T @ (Body.p + Body.R @ np.array([0, D, 0]) - Foot.p)  # Relative position
    C = np.linalg.norm(r)
    c5 = (C**2 - A**2 - B**2) / (2.0 * A * B)
    
    if c5 >= 1:
        q5 = 0.0
    elif c5 <= -1:
        q5 = np.pi
    else:
        q5 = np.arccos(c5)  # knee pitch
    
    q6a = np.arcsin((A / C) * np.sin(np.pi - q5))  # ankle pitch sub
    
    q7 = np.arctan2(r[1], r[2])  # ankle roll -pi/2 < q(6) < pi/2
    if q7 > np.pi / 2:
        q7 -= np.pi
    elif q7 < -np.pi / 2:
        q7 += np.pi
    
    q6 = -np.arctan2(r[0], np.sign(r[2]) * np.sqrt(r[1]**2 + r[2]**2)) - q6a  # ankle pitch
    
    R = Body.R.T @ Foot.R @ Rroll(-q7) @ Rpitch(-q6 - q5)  # hipZ*hipX*hipY
    q2 = np.arctan2(-R[0, 1], R[1, 1])  # hip yaw
    cz = np.cos(q2)
    sz = np.sin(q2)
    q3 = np.arctan2(R[2, 1], -R[0, 1] * sz + R[1, 1] * cz)  # hip roll
    q4 = np.arctan2(-R[2, 0], R[2, 2])  # hip pitch
    
    q = np.array([q2, q3, q4, q5, q6, q7])
    return q

def InverseKinematics_LM(to, Target, links):
    idx = FindRoute(links, to)
    wn_pos = 1/0.3
    wn_ang = 1/(2*np.pi)
    We = np.diag([wn_pos, wn_pos, wn_pos, wn_ang, wn_ang, wn_ang])
    Wn = np.eye(len(idx))

    ForwardKinematics(links, 0)
    err = CalcVWerr(Target, links[to])
    Ek = err.T @ We @ err
    # print(err)

    for n in range(10):
        J = CalcJacobian(idx, links)
        Jh = J.T @ We @ J + Wn * (Ek + 0.002)  # Hk + wn
        gerr = J.T @ We @ err  # gk
        dq = np.linalg.solve(Jh, gerr)  # new

        MoveJoints(idx, dq, links)
        ForwardKinematics(links, 0)
        err = CalcVWerr(Target, links[to])
        Ek2 = err.T @ We @ err
        if Ek2 < 1E-12:
            break
        elif Ek2 < Ek:
            Ek = Ek2
        else:
            MoveJoints(idx, -dq, links)  # revert
            ForwardKinematics(links, 0)
            break

    return np.linalg.norm(err)

def CalcVWerr(Cref, Cnow):
    perr = Cref.p - Cnow.p
    Rerr = Cnow.R.T @ Cref.R
    werr = Cnow.R @ rot2omega(Rerr)
    return np.concatenate((perr, werr))

def rot2omega(R):
    el = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    norm_el = np.linalg.norm(el)
    if norm_el > np.finfo(float).eps:
        w = np.arctan2(norm_el, np.trace(R) - 1) / norm_el * el
    elif R[0, 0] > 0 and R[1, 1] > 0 and R[2, 2] > 0:
        w = np.zeros(3)
    else:
        w = np.pi/2 * np.array([R[0, 0] + 1, R[1, 1] + 1, R[2, 2] + 1])
    return w

def CalcJacobian(idx, links):
    jsize = len(idx)
    target = links[idx[-1]].p  # absolute target position
    J = np.zeros((6, jsize))

    for n in range(jsize):
        j = idx[n]
        a = links[j].R @ links[j].a  # joint axis vector in world frame
        J[:, n] = np.concatenate((np.cross(a, target - links[j].p), a))
    return J

def MoveJoints(idx, dq, links):
    for n in range(len(idx)):
        j = idx[n]
        links[j].q += dq[n]

def SetJointAngles(ulink, idx, q):
    for n, j in enumerate(idx):
        ulink[j].q = q[n]
    
    ForwardKinematics(ulink, 0)  # In Python, we typically start indexing from 0
