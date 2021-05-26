import numpy as np

# ========================================================= #
# ===  make__referenceTrajectory                        === #
# ========================================================= #

def make__referenceTrajectory():

    # ------------------------------------------------- #
    # --- [1] make coordinate                       --- #
    # ------------------------------------------------- #
    r0 = 1.05
    sp = np.linspace( 0.0, 1.0, 101 )
    th = sp * 0.2*np.pi - np.pi/3.0
    xp = r0*np.cos( th )
    yp = r0*np.sin( th )
    zp = sp*0.0

    Data = np.concatenate( [xp[:,None],yp[:,None],zp[:,None]], axis=1 )

    # ------------------------------------------------- #
    # --- [2] save coordinate                       --- #
    # ------------------------------------------------- #
    outFile   = "dat/trajectory_reference.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    return()



# ========================================================= #
# ===  make__referenceBField                            === #
# ========================================================= #

def make__referenceBField():

    x_,y_,z_ = 0,1,2
    eps      = 1.e-4

    # ------------------------------------------------- #
    # --- [1] make coordinates                      --- #
    # ------------------------------------------------- #

    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum   = [ -1.5, 1.5, 101 ]
    x2MinMaxNum   = [ -1.5, 1.5, 101 ]
    x3MinMaxNum   = [ -1.0, 1.0,  11 ]
    ret           = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                       x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    radii         = np.sqrt( ret[...,x_]**2 + ret[...,y_]**2 + ret[...,z_]**2 )
    radii[np.where( radii < eps )] = eps
    Data          = np.zeros( (ret.shape[0],ret.shape[1],ret.shape[2],6) )
    Data[...,0:3] = ret[...,0:3]
    Data[...,3]   = 0.1 * Data[...,x_] / radii
    Data[...,4]   = 0.1 * Data[...,y_] / radii
    Data[...,5]   = +1.0 + 0.2 * np.exp( - 2.0 * radii**2 )
    # r_cbinv     = 1.0 / radii
    # vect        = ret * np.repeat( r_cbinv[:,:,:,None], 3, axis=3 )
    # Data        = np.concatenate( [ret,vect], axis=3 )


    outFile   = "dat/bfield_reference.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    return()


# ========================================================= #
# ===  make__inputTrajectory                            === #
# ========================================================= #

def make__inputTrajectory():

    # ------------------------------------------------- #
    # --- [1] make coordinate                       --- #
    # ------------------------------------------------- #
    r0 = 1.0
    sp = np.linspace( 0.0, 1.0, 101 )
    th = sp * 0.2*np.pi - np.pi/3.0
    xp = r0*np.cos( th ) + 0.05
    yp = r0*np.sin( th ) + 0.05
    zp = sp*0.0

    Data = np.concatenate( [xp[:,None],yp[:,None],zp[:,None]], axis=1 )

    # ------------------------------------------------- #
    # --- [2] save coordinate                       --- #
    # ------------------------------------------------- #
    outFile   = "dat/trajectory_input.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    return()


# ========================================================= #
# ===  make__inputGrid                                  === #
# ========================================================= #

def make__inputGrid():

    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [  0.6,  0.8,  21 ]
    x2MinMaxNum = [ -0.8, -0.6,  21 ]
    x3MinMaxNum = [ -0.01, 0.01, 11 ]
    grid        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    outFile   = "dat/coord_xy.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=grid )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__inputTrajectory()
    make__inputGrid()
    make__referenceTrajectory()
    make__referenceBField()
