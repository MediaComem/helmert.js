from cgi import print_directory
from traceback import print_list
import numpy as np
import copy
import matplotlib.pyplot as plt
import random
np.set_printoptions(threshold=np.inf)


def Cardan2R3D(alpha_deg,beta_deg,gamma_deg):
    """Matrice de rotation a partir des angles de Cardan
    avec la convention X->Y->Z

    Parameters
    ----------
        alpha_deg : float
            angle de rotation autour de l'axe X'
        beta_deg : float
            angle de rotation autour de l'axe Y'
        gamma_deg : float
            angle de rotation autour de l'axe Z'
    
    Returns
    -------
        R : 2D numpy array (3,3)
            matrice de rotation 3D
    """
    
    a = alpha_deg*np.pi/180;
    b = beta_deg*np.pi/180;
    c = gamma_deg*np.pi/180;
    
    Rz = np.array([[np.cos(c),np.sin(c),0],
          [-np.sin(c),np.cos(c),0],
          [0,0,1]])
    
    Ry = np.array([[np.cos(b),0,-np.sin(b)],
          [0,1,0],
          [np.sin(b),0,np.cos(b)]])
 
    Rx = np.array([[1,0,0],
          [0,np.cos(a),np.sin(a)],
          [0,-np.sin(a),np.cos(a)]])      
      
    R = Rz@Ry@Rx
    return R;

def R3D2Cardan(R):
    """Angles de Cardan, avec la convention X->Y->Z,
    d'une matrice de rotation 3D

    Parameters
    ----------
        R : 2D numpy array (3,3)
            matrice de rotation 3D
    
    Returns
    -------
        alpha_deg : float
            angle de rotation autour de l'axe X'
        beta_deg : float
            angle de rotation autour de l'axe Y'
        gamma_deg : float
            angle de rotation autour de l'axe Z'
    """
    
    if R[2][0] == 1 or R[2][0] == -1:
        c = 0
        if R[2][0] == -1:
            b = np.pi/2
            a = c+np.arctan2(R[0][1],R[0][2])
        else:                
            b = -np.pi/2
            a = -c+np.arctan2(-R[0][1],-R[0][2])
            
    else:
        b = -np.arcsin(R[2][0])
        a = np.arctan2(R[2][1]/np.cos(b),R[2][2]/np.cos(b))
        c = np.arctan2(R[1][0]/np.cos(b),R[0][0]/np.cos(b))
        
    alpha_deg = -a*180/np.pi
    beta_deg = -b*180/np.pi
    gamma_deg = -c*180/np.pi
    
    return alpha_deg,beta_deg,gamma_deg


class HelmertTransformation3D:
    
    def __init__(self,id):  
        self.id = id
        self.successful = False
        self.iter = '99 (99)'
        self.JD_UTC = np.nan

        self.points_global = {}
        self.points_local = {}
        self.points_common = {}
        
        self.helmert3DParam = {'tX':[0.0,'unknown'],
                           'tY':[0.0,'unknown'],
                           'tZ':[0.0,'unknown'],
                           'rX':[0.0,'unknown'],
                           'rY':[0.0,'unknown'],
                           'rZ':[0.0,'unknown']
                           }                      


        self.helmert3DParamPrecision = {'s0':[np.nan,'unknown'],
                           'stX':[np.nan,'unknown'],
                           'stY':[np.nan,'unknown'],
                           'stZ':[np.nan,'unknown'],
                           'srX':[np.nan,'unknown'],
                           'srY':[np.nan,'unknown'],
                           'srZ':[np.nan,'unknown']
                           }                              
        
        self.sigma_X = 0.001
        self.sigma_Y = 0.001
        self.sigma_Z = 0.001
        
        self.sigma_x = 0.001
        self.sigma_y = 0.001
        self.sigma_z = 0.001
        
        self.sigma = 0.001

    def setParamTo0(self):
        self.helmert3DParam['tX'][0]=0.0
        self.helmert3DParam['tY'][0]=0.0
        self.helmert3DParam['tZ'][0]=0.0
        self.helmert3DParam['rX'][0]=0.0
        self.helmert3DParam['rY'][0]=0.0
        self.helmert3DParam['rZ'][0]=0.0
        

    #import fichier points global
    def importPointsGlobal(self,path):
        #numero:[X,Y,Z]
        #numero:[0,1,2]        
        self.points_global = {}
        with open(path) as fp:
           line = fp.readline().strip()
           while line:
               pt,X,Y,Z = line.split('\t')
               self.points_global.update({pt:[float(X),float(Y),float(Z)]})
               line = fp.readline().strip()
                   
    def importPointsLocal(self,path):
        #numero:[x,y,z]
        #numero:[0,1,2]        
        self.points_local = {}
        with open(path) as fp:
           line = fp.readline().strip()
           while line:
               pt,x,y,z = line.split('\t')
               self.points_local.update({pt:[float(x),float(y),float(z)]})
               line = fp.readline().strip()
               
               
               
    def obs_eq_X(self,param,x,y,z):        
                    
        tX = param['tX'][0]
        tY = param['tY'][0]
        tZ = param['tZ'][0]       
        rX = param['rX'][0]
        rY = param['rY'][0]
        rZ = param['rZ'][0]       

        R = Cardan2R3D(rX,rY,rZ);

        t_=np.array([[tX],
                    [tY],
                    [tZ]])

        x_=np.array([[x],
                     [y],
                     [z]])
        
        X_ = t_+R@x_
            
        return X_[0][0]               

    def obs_eq_Y(self,param,x,y,z):        
                    
        tX = param['tX'][0]
        tY = param['tY'][0]
        tZ = param['tZ'][0]       
        rX = param['rX'][0]
        rY = param['rY'][0]
        rZ = param['rZ'][0]       

        R = Cardan2R3D(rX,rY,rZ);

        t_=np.array([[tX],
                    [tY],
                    [tZ]])

        x_=np.array([[x],
                     [y],
                     [z]])
        
        X_ = t_+R@x_
            
        return X_[1][0]               


    def obs_eq_Z(self,param,x,y,z):        
                    
        tX = param['tX'][0]
        tY = param['tY'][0]
        tZ = param['tZ'][0]       
        rX = param['rX'][0]
        rY = param['rY'][0]
        rZ = param['rZ'][0]       

        R = Cardan2R3D(rX,rY,rZ);

        t_=np.array([[tX],
                    [tY],
                    [tZ]])

        x_=np.array([[x],
                     [y],
                     [z]])
        
        X_ = t_+R@x_
            
        return X_[2][0]               





    def computeHelmert(self,dict_points, P):
        '''
        Calculate Helmert transformation (3 translations, 3 rotations) with SVD method (least-squares rigid motion using SVD)

        Parameters
        ----------
        dict_points : dict
            Dictionnary of common points.
        P : np.array
            Matrix of weights.

        Returns
        -------
        None.

        '''
        ################################################################
        #Compute local and global centroid    
        ################################################################
 
        global_X = 0
        global_Y = 0
        global_Z = 0
            
        local_x = 0
        local_y = 0
        local_z = 0 

        poids_X = 0
        poids_Y = 0
        poids_Z = 0
            
        poids_x = 0
        poids_y = 0
        poids_z = 0
        
        k = 0
        for key, values in dict_points.items():
            global_X += values[0]*P[k][k]
            global_Y += values[1]*P[k+1][k+1]
            global_Z += values[2]*P[k+2][k+2]
            
            poids_X += P[k][k]
            poids_Y += P[k+1][k+1]
            poids_Z += P[k+2][k+2]
            
            local_x += values[3]*P[k][k]
            local_y += values[4]*P[k+1][k+1]
            local_z += values[5]*P[k+2][k+2]

            poids_x += P[k][k]
            poids_y += P[k+1][k+1]
            poids_z += P[k+2][k+2]
            
            k += 3
        
        #Global centroid
        centroid_X = global_X/poids_X
        centroid_Y = global_Y/poids_Y
        centroid_Z = global_Z/poids_Z

        #Local centroid
        centroid_x = local_x/poids_x
        centroid_y = local_y/poids_y
        centroid_z = local_z/poids_z
        
        global_centroid = np.array([[centroid_X],
                                    [centroid_Y],
                                    [centroid_Z]])
        
        local_centroid = np.array([[centroid_x],
                                    [centroid_y],
                                    [centroid_z]])

        nbr_pts_communs = len(dict_points)
        W_SVD = np.eye(nbr_pts_communs)

        
        ################################################################
        #Compute centered vectors
        ################################################################
        
        global_matrix = np.zeros(shape=(3,nbr_pts_communs))
        local_matrix = np.zeros(shape=(3,nbr_pts_communs))        
        
        k = 0
        for key, values in dict_points.items():
            reduit_X = values[0] - centroid_X
            reduit_Y = values[1] - centroid_Y
            reduit_Z = values[2] - centroid_Z

            reduit_x = values[3] - centroid_x
            reduit_y = values[4] - centroid_y
            reduit_z = values[5] - centroid_z            

        
            global_matrix[0][k] = reduit_X
            global_matrix[1][k] = reduit_Y
            global_matrix[2][k] = reduit_Z

            local_matrix[0][k] = reduit_x
            local_matrix[1][k] = reduit_y
            local_matrix[2][k] = reduit_z
            
            k += 1            
        
        ################################################################    
        #Compute covariance matrix
        ################################################################
        
        S = local_matrix@W_SVD@global_matrix.T
        #Compute singular value decomposition
        U, s, VT = np.linalg.svd(S)
        print(np.linalg.svd(S))
        V = VT.T
        
        #Determinant
        det_VUT = np.linalg.det(V@U.T)
        one_matrix = np.eye(len(V))
        
        one_matrix[-1][-1] = det_VUT
        
        ################################################################
        #Compute rotation
        ################################################################
        R = V@one_matrix@U.T
        
        #Calculate degree angles
        alpha_deg,beta_deg,gamma_deg = R3D2Cardan(R)

        self.helmert3DParam['rX'][0] = alpha_deg
        self.helmert3DParam['rY'][0] = beta_deg
        self.helmert3DParam['rZ'][0] = gamma_deg 
        
        ################################################################
        #Compute translation
        ################################################################
        
        t = global_centroid - R@local_centroid

        t_x = t[0][0]
        t_y = t[1][0]
        t_z = t[2][0]
 
        self.helmert3DParam['tX'][0] = t_x
        self.helmert3DParam['tY'][0] = t_y
        self.helmert3DParam['tZ'][0] = t_z        
         

    def computeResiduals(self, dict_points, param_RANSAC):
        '''
        Compute residuals

        Parameters
        ----------
        dict_points : dict
            Dictionnary of common points.
        param_RANSAC : dict
            Dictionnary of RANSAC parameters.

        Returns
        -------
        v : array
            Array of residuals.

        '''
        nbr_observations = len(dict_points)*3
        v = np.zeros(shape=(nbr_observations,1))
        k = 0
        for key,value in dict_points.items():
            X = value[0]
            Y = value[1]
            Z = value[2]
            x = value[3]
            y = value[4]
            z = value[5]
            Xr = self.obs_eq_X(self.helmert3DParam,x,y,z)
            Yr = self.obs_eq_Y(self.helmert3DParam,x,y,z)
            Zr = self.obs_eq_Z(self.helmert3DParam,x,y,z)
            value[6] = X-Xr
            value[7] = Y-Yr
            value[8] = Z-Zr  
            
            v[k] = value[6]
            v[k+1] = value[7]
            v[k+2] = value[8]
            
            #Check if residual is in tolerance
            if np.fabs(v[k]) <= param_RANSAC['tol'] and np.fabs(v[k+1]) <= param_RANSAC['tol'] and np.fabs(v[k+2]) <= param_RANSAC['tol']:
                value[10] = True
            else:
                value[10] = False
            k += 3  
        return v

        

    def estimateHelmert3DSVD(self, nb_iterations_ransac):
        '''
        Estimate Helmert 3D with singular values decomposition method

        Returns
        -------
        None.

        '''
      
        print('\nCOMPUTE HELMERT3D WITH SVD')
        #create common points dictionnary
        self.points_common = {}
        for key,value in self.points_global.items():
            if key in self.points_local.keys():
                X = value[0]
                Y = value[1]
                Z = value[2]
                x = self.points_local[key][0]
                y = self.points_local[key][1]
                z = self.points_local[key][2]
                self.points_common.update({key:[X,Y,Z,x,y,z,np.nan,np.nan,np.nan,np.nan, False]}) 
                                              

        if len(self.points_common) < 3:
            print('\nNOT ENOUGH POINTS FOR HELMERT3D\n')
            self.successful = False
            return

        #Compute weight matrix
        P = np.eye(len(self.points_common))

        #RANSAC parameters
        param_RANSAC = {}
        param_RANSAC.update({'nbr_iter_RANSAC':nb_iterations_ransac})
        param_RANSAC.update({'tol':self.sigma*3})
        
        #RANSAC results
        res_RANSAC = {}
        res_RANSAC.update({'max_nbr_inliers':0})
        res_RANSAC.update({'tX_best':0.0})
        res_RANSAC.update({'tY_best':0.0})
        res_RANSAC.update({'tZ_best':0.0})     
        res_RANSAC.update({'rX_best':0.0})
        res_RANSAC.update({'rY_best':0.0})
        res_RANSAC.update({'rZ_best':0.0})
        
        for iter_RANSAC in range(0,param_RANSAC['nbr_iter_RANSAC']):            
            
            dict_points_ransac = {}
            
            print("RANSAC:{:d}, #inliers={:d}".format(iter_RANSAC,res_RANSAC['max_nbr_inliers']))
                          
            #Choose 3 random points
            no1 = random.choice(list(self.points_common.keys()))
            no2 = random.choice(list(self.points_common.keys()))
            while no1==no2:
                no2 = random.choice(list(self.points_common.keys()))    
            no3 = random.choice(list(self.points_common.keys()))
            while no1==no3 or no2==no3:
                no3 = random.choice(list(self.points_common.keys())) 

            dict_points_ransac.update({no1:self.points_common[no1]})
            dict_points_ransac.update({no2:self.points_common[no2]})
            dict_points_ransac.update({no3:self.points_common[no3]})
            
            #Compute Helmert adjustment with 3 random points
            self.computeHelmert(dict_points_ransac, P)
            
            #compute residuals
            v = self.computeResiduals(self.points_common, param_RANSAC)
            
            #Count inliers
            nbr_inliers = 0
            for i in v:
                if np.fabs(i) <= param_RANSAC['tol']:
                    nbr_inliers += 1
            
            #Update Helmert parameters
            if  nbr_inliers >  res_RANSAC['max_nbr_inliers']:
                res_RANSAC['max_nbr_inliers'] = nbr_inliers
                res_RANSAC.update({'pt1_best':no1})
                res_RANSAC.update({'pt2_best':no2})
                res_RANSAC.update({'pt3_best':no3})
                res_RANSAC.update({'tX_best':self.helmert3DParam['tX'][0]})
                res_RANSAC.update({'tY_best':self.helmert3DParam['tY'][0]})
                res_RANSAC.update({'tZ_best':self.helmert3DParam['tZ'][0]})    
                res_RANSAC.update({'rX_best':self.helmert3DParam['rX'][0]})
                res_RANSAC.update({'rY_best':self.helmert3DParam['rY'][0]})
                res_RANSAC.update({'rZ_best':self.helmert3DParam['rZ'][0]})

        #Refaire une classification inliers/outliers avec les 3 points qui donnaient le plus d'inliers
        dict_points_best = {}
        dict_points_best.update({res_RANSAC['pt1_best']:self.points_common[res_RANSAC['pt1_best']]})
        dict_points_best.update({res_RANSAC['pt2_best']:self.points_common[res_RANSAC['pt2_best']]})
        dict_points_best.update({res_RANSAC['pt3_best']:self.points_common[res_RANSAC['pt3_best']]})            
        #Compute Helmert adjustment with 3 best points
        self.computeHelmert(dict_points_best, P)        
        #compute residuals
        v = self.computeResiduals(self.points_common, param_RANSAC)        
        
        
        #Selection of common points
        dict_final_inliers_points = {}
        
        for key,value in self.points_common.items():
            if value[10] == True:
                dict_final_inliers_points.update({key:value})
        
        #Compute weight matrix
        P = np.eye(len(dict_final_inliers_points)*3)
        self.computeHelmert(dict_final_inliers_points, P)


    def estimateHelmert3DSVD_minimum(self):
        '''
        Estimate Helmert 3D with singular values decomposition method

        Returns
        -------
        None.

        '''
        
        print('\nCOMPUTE HELMERT3D WITH SVD')
        #create common points dictionnary
        self.points_common = {}
        for key,value in self.points_global.items():
            if key in self.points_local.keys():
                X = value[0]
                Y = value[1]
                Z = value[2]
                x = self.points_local[key][0]
                y = self.points_local[key][1]
                z = self.points_local[key][2]
                self.points_common.update({key:[X,Y,Z,x,y,z,np.nan,np.nan,np.nan,np.nan, False]}) 
                                              

        if len(self.points_common) < 3:
            print('\nNOT ENOUGH POINTS FOR HELMERT3D\n')
            self.successful = False
            return

        #Compute weight matrix
        P = np.eye(len(self.points_common)*3)

        #Compute Helmert adjustment with 3 best points
        self.computeHelmert(self.points_common, P)        

    

    def plot2DMap(self,typePoints='local'):
        if typePoints == 'local':
            plt.figure()
            plt.title('Points local')
            
            for key, data in self.points_local.items():

                plt.plot(data[0], data[1], 'ob')
                plt.text(data[0], data[1], key)
        
                plt.xlabel('x local [m]')
                plt.ylabel('y local [m]')

        if typePoints == 'global':
            plt.figure()
            plt.title('Points global')  
            plt.xlabel('x global [m]')
            plt.ylabel('y global [m]')
            
            for key, data in self.points_global.items():
                plt.plot(data[0], data[1], 'ob')
                plt.text(data[0], data[1], key)
                


    def exportAllLocalPointsINGlobalCoord(self,path):
        file_out = open(path,'w')
        print('\nLOCAL POINTS in GLOBAL COORDINATES')
        for key,value in self.points_local.items():
            x = value[0]
            y = value[1]
            z = value[2]
            
            X = self.obs_eq_X(self.helmert3DParam,x,y,z)
            Y = self.obs_eq_Y(self.helmert3DParam,x,y,z)
            Z = self.obs_eq_Z(self.helmert3DParam,x,y,z)
            
            print('{}\t{:0.6f}\t{:0.6f}\t{:0.6f}'.format(key,X,Y,Z))
            file_out.write('{}\t{:0.6f}\t{:0.6f}\t{:0.6f}\n'.format(key,X,Y,Z))
            
        file_out.close()
        
    def exportParameters(self,path):
        file_out = open(path,'w')
        tX = self.helmert3DParam['tX'][0]
        tY = self.helmert3DParam['tY'][0]
        tZ = self.helmert3DParam['tZ'][0]
        rX = self.helmert3DParam['rX'][0]
        rY = self.helmert3DParam['rY'][0]
        rZ = self.helmert3DParam['rZ'][0]
        file_out.write('tX\t{:0.6f}\n'.format(tX))
        file_out.write('tY\t{:0.6f}\n'.format(tY))
        file_out.write('tZ\t{:0.6f}\n'.format(tZ))
        file_out.write('rX\t{:0.6f}\n'.format(rX))
        file_out.write('rY\t{:0.6f}\n'.format(rY))
        file_out.write('rZ\t{:0.6f}\n'.format(rZ))
            
        file_out.close()
            
            
    def printParameters(self):
        print('\nPARAMETERS')
        print('tX = {:+0.6f} [m] +/- {:+0.6f} [m]'.format(self.helmert3DParam['tX'][0],self.helmert3DParamPrecision['stX'][0]))
        print('tY = {:+0.6f} [m] +/- {:+0.6f} [m]'.format(self.helmert3DParam['tY'][0],self.helmert3DParamPrecision['stY'][0]))
        print('tZ = {:+0.6f} [m] +/- {:+0.6f} [m]'.format(self.helmert3DParam['tZ'][0],self.helmert3DParamPrecision['stZ'][0]))
        print('rX = {:+0.6f} [deg] +/- {:+0.6f} [deg]'.format(self.helmert3DParam['rX'][0],self.helmert3DParamPrecision['srX'][0]))
        print('rY = {:+0.6f} [deg] +/- {:+0.6f} [deg]'.format(self.helmert3DParam['rY'][0],self.helmert3DParamPrecision['srY'][0]))
        print('rZ = {:+0.6f} [deg] +/- {:+0.6f} [deg]'.format(self.helmert3DParam['rZ'][0],self.helmert3DParamPrecision['srZ'][0]))


    def printResiduals(self):
        print('\nRESIDUALS')
        for key,value in self.points_common.items():
            X = value[0]
            Y = value[1]
            Z = value[2]
            vX = value[6]
            vY = value[7]
            vZ = value[8]
            print('{:5s}: X:{:+9.4f} Y:{:+9.4f} Z:{:+9.4f} vX:{:+6.4f} vY:{:+6.4f} vZ:{:+6.4f}'.format(key,X,Y,Z,vX,vY,vZ))



def local2global(helm, dict_points_local):
    points_global = {}

    tX = helm.helmert3DParam['tX'][0]
    tY = helm.helmert3DParam['tY'][0]
    tZ = helm.helmert3DParam['tZ'][0]       
    rX = helm.helmert3DParam['rX'][0]
    rY = helm.helmert3DParam['rY'][0]
    rZ = helm.helmert3DParam['rZ'][0]       

    R = Cardan2R3D(rX,rY,rZ);

    t_=np.array([[tX],
                [tY],
                [tZ]])

    for key, data in dict_points_local.items():
        no_point = key
        x = data[0]
        y = data[1]
        z = data[2]
        
        x_=np.array([[x],
                     [y],
                     [z]])
    
        X_ = t_+R@x_
        
        X = X_[0][0]
        Y = X_[1][0]
        Z = X_[2][0]
    
        points_global.update({no_point:[X,Y,Z]})
    return points_global            


def global2local(helm, dict_points_global):
    points_local = {}

    tX = helm.helmert3DParam['tX'][0]
    tY = helm.helmert3DParam['tY'][0]
    tZ = helm.helmert3DParam['tZ'][0]       
    rX = helm.helmert3DParam['rX'][0]
    rY = helm.helmert3DParam['rY'][0]
    rZ = helm.helmert3DParam['rZ'][0]       

    R = Cardan2R3D(rX,rY,rZ);

    t_=np.array([[tX],
                [tY],
                [tZ]])

    for key, data in dict_points_global.items():
        no_point = key
        X = data[0]
        Y = data[1]
        Z = data[2]
        
        X_=np.array([[X],
                     [Y],
                     [Z]])
    
        x_ = R@(X_-t_)
        
        x = x_[0][0]
        y = x_[1][0]
        z = x_[2][0]
    
        points_local.update({no_point:[x,y,z]})
    return points_local 
