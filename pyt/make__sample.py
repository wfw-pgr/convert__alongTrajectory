import numpy as np

# ========================================================= #
# ===  make__sampleTrajectory                           === #
# ========================================================= #

def make__sampleTrajectory():

    # ------------------------------------------------- #
    # --- [1] make coordinate                       --- #
    # ------------------------------------------------- #
    r0 = 1.0
    sp = np.linspace( 0.0, 1.0, 101 )
    th = sp * 0.2*np.pi - np.pi/3.0
    xp = r0*np.cos( th )
    yp = r0*np.sin( th )
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
# ===  make__sampleBField                               === #
# ========================================================= #

def make__sampleBField():

    x_,y_,z_ = 0,1,2
    eps      = 1.e-4

    # ------------------------------------------------- #
    # --- [1] make coordinates                      --- #
    # ------------------------------------------------- #

    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.0, 1.0, 101 ]
    x2MinMaxNum = [ -1.0, 1.0, 101 ]
    x3MinMaxNum = [ -1.0, 1.0,  11 ]
    ret         = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    radii       = np.sqrt( ret[...,x_]**2 + ret[...,y_]**2 + ret[...,z_]**2 )
    radii[np.where( radii < eps )] = eps
    r_cbinv     = 1.0 / radii
    vect        = ret * np.repeat( r_cbinv[:,:,:,None], 3, axis=3 )
    Data        = np.concatenate( [ret,vect], axis=3 )

    outFile   = "dat/bfield_input.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    return()
    

# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__sampleTrajectory()
    make__sampleBField()
