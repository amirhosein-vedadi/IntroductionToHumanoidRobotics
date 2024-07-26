import numpy as np
import matplotlib.pyplot as plt
from link import Link, FindMother, ForwardKinematics, FindRoute, SetJointAngles, InverseKinematics_LM, IK_leg, CalcJacobian

# Global variables
uLINK = []
ToRad = np.pi / 180
ToDeg = 180 / np.pi

# Constants
BODY = 0
RLEG_J0, RLEG_J1, RLEG_J2, RLEG_J3, RLEG_J4, RLEG_J5 = range(1, 7)
LLEG_J0, LLEG_J1, LLEG_J2, LLEG_J3, LLEG_J4, LLEG_J5 = range(7, 13)

def SetupBipedRobot():
    global uLINK
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
    uLINK[BODY].p = np.array([0.0, 0.0, 0.65])
    uLINK[BODY].R = np.eye(3)
    ForwardKinematics(uLINK, 0)

    for link in uLINK:
        link.v = np.zeros(3)
        link.w = np.zeros(3)
        link.dq = 0

def print_joint_angles():
    print("\nJoint angles (in radians):")
    for i in range(6):
        print(f"RLEG_J{i}: {uLINK[RLEG_J0 + i].q}")
        print(f"LLEG_J{i}: {uLINK[LLEG_J0 + i].q}")

def main():
    global uLINK
    SetupBipedRobot()

    idx = FindRoute(uLINK, RLEG_J5)
    SetJointAngles(uLINK, idx, np.array([0, 0, -25, 50, -25, 0]) * ToRad)
    print_joint_angles()
    Rfoot = uLINK[RLEG_J5]

    Height = np.linalg.norm(uLINK[RLEG_J0].p - uLINK[RLEG_J5].p)
    LegLength = 0.6
    SingularPoint = np.sqrt(0.6**2 - Height**2)

    print('*** Levenberg-Marquardt ***')

    xd_m = np.arange(0, 0.405, 0.005)
    Nstep = len(xd_m)
    q_m = np.zeros((Nstep, 6))
    analy_q_m = np.zeros((Nstep, 6))
    x_m = np.zeros(Nstep)
    manip_m = np.zeros(Nstep)

    for n in range(Nstep):
        # Create a copy of the position of Rfoot
        Rfoot_position = Rfoot.p.copy()
        Rfoot_position[0] = xd_m[n]  # Update only the x position

        # Create a temporary link for the foot with the updated position
        Rfoot_temp = Link(name='Rfoot_temp', mass=Rfoot.mass, sister=Rfoot.sister,
                          child=Rfoot.child, b=Rfoot.b, a=Rfoot.a, q=Rfoot.q)
        Rfoot_temp.p = Rfoot_position
        Rfoot_temp.R = Rfoot.R  # Keep the same orientation

        rerr_norm = InverseKinematics_LM(RLEG_J5, Rfoot_temp, uLINK)

        x_m[n] = uLINK[RLEG_J5].p[0]
        q_m[n, :] = [uLINK[i].q for i in idx]
        analy_q_m[n, :] = IK_leg(uLINK[BODY], -0.1, 0.3, 0.3, Rfoot_temp)

        J = CalcJacobian(idx, uLINK)
        
        manip_m[n] = abs(np.linalg.det(J))

        print(f'Error: {rerr_norm:8.3e}')

    # Plotting
    plt.figure()
    plt.plot(xd_m, ToDeg * q_m[:, [2, 3, 4]], linewidth=2)
    plt.plot(xd_m, ToDeg * analy_q_m[:, [2, 3, 4]], ':')
    plt.axvline(x=SingularPoint, color='m', linestyle='-.')
    plt.legend(['RLEG_J2', 'RLEG_J3', 'RLEG_J4'])
    plt.ylabel('q [deg]')
    plt.xlabel('x [m]')

    plt.figure()
    plt.subplot(211)
    plt.plot(xd_m, x_m, xd_m, xd_m, 'r--')
    plt.axvline(x=SingularPoint, color='m', linestyle='-.')
    plt.ylabel('x [m]')

    plt.subplot(212)
    plt.plot(xd_m, manip_m)
    plt.axvline(x=SingularPoint, color='m', linestyle='-.')
    plt.ylabel('Manipulability')
    plt.xlabel('xd [m]')

    plt.show()

if __name__ == "__main__":
    main()

