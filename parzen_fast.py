from __future__ import print_function
import numpy as np
import math as m
import matplotlib.pyplot as plt

def gauss(sample_phi,phi,sample_psi,psi,sigma):
    return np.exp((-1.0*((sample_phi-phi)**2 + (sample_psi-psi)**2))/(2.0*(sigma**2)))
# print gauss(58.7+shiftPHI[l],phi,sample_psi[k]+shiftPSI[l],psi,sigma)

def gauss_dist( std, x ):
    return np.exp( -x*x/2/std/std )

def write_parzen_data( angle_filename, distr_fn, sigma ):
    with open( distr_fn, 'w' ) as dfn:
        dfn.write( "#Loop_length\tposition\tpsi\tphi\tcounts\n" )

    afn = open( angle_filename, 'r' )
    next(afn) # skip column names

    data = list(map( float, next(afn).strip().split() ))

    loop, pos = data[0], data[1]
    s_psi = [data[2]]
    s_phi = [data[3]]

    for line in afn:
        data = list(map( float, line.strip().split() ))
        if loop != data[0] or pos != data[1]:
            print('\t', loop, pos)
            write_parzen_windows( loop, pos, np.array(s_psi), np.array(s_phi), distr_fn, sigma )
            loop, pos = data[0], data[1]
            s_psi = [data[2]]
            s_phi = [data[3]]
        else:
            s_psi.append( data[2] )
            s_phi.append( data[3] )

    write_parzen_windows( loop, pos, np.array(s_psi), np.array(s_phi), distr_fn, sigma )
    afn.close()


def write_parzen_windows( loop, pos, sample_psi, sample_phi, distr_filename, sigma ):
    Xn = bin_sizeX = 72
    Yn = bin_sizeY = 72
    cell_x = np.linspace(-180, 180, bin_sizeX+1)
    cell_y = np.linspace(-180, 180, bin_sizeY+1)

    [X, Y] = np.meshgrid(cell_x, cell_y)

    cell = np.zeros_like(X)

    n = angle_listSize = len( sample_psi )

    min_dists = np.full( (Xn+1, Yn+1, n), np.inf )
    for ii in [-360, 0, 360]:
        for jj in [-360, 0, 360]:
            dists = (cell_x.reshape(Xn+1, 1, 1) - sample_psi.reshape(1, 1, n) + ii)**2 + \
                    (cell_y.reshape(1, Yn+1, 1) - sample_phi.reshape(1, 1, n) + jj)**2
            min_dists = np.minimum( dists, min_dists )

    cell = np.sum(gauss_dist( sigma, np.sqrt( min_dists ) ), axis=2)

    with open( distr_filename, 'a') as dfn:
        for i in range( bin_sizeX):
            for j in range( bin_sizeY ):
                dfn.write("{}\t{}\t{}\t{}\t{}\n".format(int(loop), int(pos),
                                                        Y[i][j], X[i][j],
                                                        cell[i][j]) )

if __name__ == "__main__":

    #for sigma in np.arange(0.1, 10.1, 0.1):
    #    print("Extend w/ sigma {}".format( sigma ))
    #    dist_path = "distributions2/Extend_Phi_Psi-distribution-sigma-{:.2f}_fast.txt".format( sigma )
    #    write_parzen_data( "angles2/Extend_Phi_Psi-angles.txt", dist_path, sigma )

    #for sigma in np.arange(0.1, 10.1, 0.1):
    #    print("Kink w/ sigma {}".format( sigma ))
    sigma = 10.0
    dist_path = "distributions2/Kink_Phi_Psi-distribution-sigma-{:.2f}_yee.txt".format( sigma )
    write_parzen_data( "angles2/Kink_Phi_Psi-angles.txt", dist_path, sigma )

