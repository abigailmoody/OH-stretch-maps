
a0 = 0.5291772

class Spectroscopic_Map:
    # Uses TIP4P as the default when not defined elsewhere
    def __init__(self):
        return None
    def w_map(self, E):
        w = 3760.2 - (3541.7 * E) - (152677 * (E ** 2))
        return w
    def mu_map(self, E):
        mu = 0.1646 + (11.39 * E) + (63.41 * (E ** 2))
        return mu
    def x_map(self, w):
        x = 0.19285 - (1.7261e-5 * w)
        return x
    def p_map(self, w):
        p = 1.6466 + (5.7692e-4 * w)
        return p
    def intra_map(self, Ei, Ej, xi, xj, pi, pj):
        intra = ((-1361 + (27165 * (Ei + Ej))) * xi * xj) - (1.887 * pi * pj)
        return intra

class TIP3P_Map(Spectroscopic_Map):
    def w_map(self, E):
        w = 3742.81 - (4884.72 * E) - (65278.36 * (E**2))
        return w
    def mu_map(self, E):
        mu = 0.12 + (12.28 * E)
        return mu
    def x_map(self, w):
        x = (0.1019 - (9.0611e-6 * w)) / a0
        return x

class SPCE_Map(Spectroscopic_Map):
    def w_map(self, E):
        w = 3762 - (5060 * E) - (86225 * (E**2))
        return w
    def mu_map(self, E):
        mu = 0.7112 + (75.58 * E)
        return mu
    def x_map(self, w):
        x = 0.1934 - (1.75e-5 * w)
        return x
    def p_map(self, w):
        p = 1.611 + (5.893e-4 * w)
        return p
    def intra_map(self, Ei, Ej, xi, xj, pi, pj):
        intra = ((-1789 + (23852 * (Ei + Ej))) * xi * xj) - (1.966 * pi * pj)
        return intra

class TIP4P_Map(Spectroscopic_Map):
    def w_map(self, E):
        w = 3760.2 - (3541.7 * E) - (152677 * (E ** 2))
        return w
    def mu_map(self, E):
        mu = 0.1646 + (11.39 * E) + (63.41 * (E ** 2))
        return mu
    def x_map(self, w):
        x = 0.19285 - (1.7261e-5 * w)
        return x
    def p_map(self, w):
        p = 1.6466 + (5.7692e-4 * w)
        return p
    def intra_map(self, Ei, Ej, xi, xj, pi, pj):
        intra = ((-1361 + (27165 * (Ei + Ej))) * xi * xj) - (1.887 * pi * pj)
        return intra


class Alcohol_Map(Spectroscopic_Map):
    # Uses aliphatic alcohol map as the default when not defined elsewhere
    def __init__(self):
        return None
    def w_map(self, E):
        w = 3744 - 7239*E - 52826*E**2
        return w
    def x_map(self, w):
        x = 0.1912 - 1.7158e-5*w
        return x
    def mu_map(self, E):
        mu = 0.1053 + 16.02*E
        return mu

def oNP_Map(ground_truth='B3LYP'):
    if ground_truth.lower() == 'aimnet':
        return oNP_Map_AIMNet()
    else:
        return oNP_Map_B3LYP()
            
class oNP_Map_B3LYP(Alcohol_Map):
    def w_map(self, E):
        w = 3315.51 - 3363.21*E - 4037.08*E**2
        return w

class oNP_Map_AIMNet(Alcohol_Map):
    def w_map(self, E):
        w = 3335.86 - 6636.19*E + 50016.47*E**2
        return w

def pNP_Map(ground_truth='B3LYP'):
    if ground_truth.lower() == 'aimnet':
        return pNP_Map_AIMNet()
    else:
        return pNP_Map_B3LYP()
            
class pNP_Map_B3LYP(Alcohol_Map):
    def w_map(self, E):
        w = 3623.3 - 2545.92*E - 105486.8*E**2
        return w

class pNP_Map_AIMNet(Alcohol_Map):
    def w_map(self, E):
        w = 3612.94 - 3823.67*E - 86909.22*E**2
        return w
        
