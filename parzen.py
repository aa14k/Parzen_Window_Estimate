import numpy as np
import math as m
import matplotlib.pyplot as plt

def gauss(sample_phi,phi,sample_psi,psi,sigma):
    return np.exp((-1.0*((sample_phi-phi)**2 + (sample_psi-psi)**2))/(2.0*(sigma**2)))
# print gauss(58.7+shiftPHI[l],phi,sample_psi[k]+shiftPSI[l],psi,sigma)

def write_parzen_data( angle_filename, distr_fn ):
    with open( distr_fn, 'w' ) as dfn:
        dfn.write( "#Loop_length\tposition\tpsi\tphi\tcounts\n" )
    
    afn = open( angle_filename, 'r' )
    next(afn) # skip column names
    
    data = map( float, next(afn).strip().split() )
    loop, pos = data[0], data[1]
    s_psi = [data[2]]
    s_phi = [data[3]]
    
    for line in afn:
        data = map( float, line.strip().split() )
        if loop != data[0] or pos != data[1]:
            print loop, pos
            write_parzen_windows( loop, pos, np.array(s_psi), np.array(s_phi), distr_fn )
            loop, pos = data[0], data[1]
            s_psi = [data[2]]
            s_phi = [data[3]]
        else:
            s_psi.append( data[2] )
            s_phi.append( data[3] )
    
    write_parzen_windows( loop, pos, np.array(s_psi), np.array(s_phi), distr_fn )
    afn.close()
 

def write_parzen_windows( loop, pos, sample_psi, sample_phi, distr_filename ):
    bin_sizeX = 72
    bin_sizeY = 72
    cell_x = np.linspace(-180, 180, bin_sizeX+1)
    cell_y = np.linspace(-180, 180, bin_sizeY+1)
    
    [X, Y] = np.meshgrid(cell_x, cell_y)
    
    cell = np.zeros_like(X)
    sigma = 1.5
    
    angle_listSize = len( sample_psi )
    
    shiftPhi1 = sample_phi + 360.0
    shiftPhi2 = sample_phi + 360.0
    shiftPhi3 = sample_phi + 360.0
    shiftPhi4 = sample_phi 
    shiftPhi5 = sample_phi
    shiftPhi6 = sample_phi 
    shiftPhi7 = sample_phi - 360.0
    shiftPhi8 = sample_phi - 360.0
    shiftPhi9 = sample_phi - 360.0
    
    shiftPsi1 = sample_psi - 360.0
    shiftPsi2 = sample_psi 
    shiftPsi3 = sample_psi + 360.0
    shiftPsi4 = sample_psi - 360.0
    shiftPsi5 = sample_psi 
    shiftPsi6 = sample_psi + 360.0
    shiftPsi7 = sample_psi - 360.0
    shiftPsi8 = sample_psi 
    shiftPsi9 = sample_psi + 360.0
    
    for i in range(bin_sizeX):
        psi = i*360.0/bin_sizeX - 180.0
        for j in range(bin_sizeY):
            phi = j*360.0/bin_sizeY - 180.0
            for k in range(angle_listSize):
                distance1 = gauss(shiftPhi1[k],phi,shiftPsi1[k],psi,sigma)
                distance2 = gauss(shiftPhi2[k],phi,shiftPsi2[k],psi,sigma)
                distance3 = gauss(shiftPhi3[k],phi,shiftPsi3[k],psi,sigma)
                distance4 = gauss(shiftPhi4[k],phi,shiftPsi4[k],psi,sigma)
                distance5 = gauss(shiftPhi5[k],phi,shiftPsi5[k],psi,sigma)
                distance6 = gauss(shiftPhi6[k],phi,shiftPsi6[k],psi,sigma)
                distance7 = gauss(shiftPhi7[k],phi,shiftPsi7[k],psi,sigma)
                distance8 = gauss(shiftPhi8[k],phi,shiftPsi8[k],psi,sigma)
                distance9 = gauss(shiftPhi9[k],phi,shiftPsi9[k],psi,sigma)
                
                cell[i,j] = cell[i,j] + max(distance1, distance2, distance3,
                                            distance4, distance5, distance6,
                                            distance7, distance8, distance9)

                
    with open( distr_filename, 'a') as dfn:
        for i in range( bin_sizeX):
            for j in range( bin_sizeY ):
                dfn.write("{}\t{}\t{}\t{}\t{}\n".format(int(loop), int(pos),
                                                        Y[i][j], X[i][j],
                                                        cell[i][j]) )

if __name__ == "__main__":
    write_parzen_data( "Kink_Phi_Psi-angles.txt", "Kink_Phi_Psi-distribution.txt" )
