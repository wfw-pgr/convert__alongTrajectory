import numpy as np

# ========================================================= #
# ===  show__tangentNormal.py                           === #
# ========================================================= #

def show__tangentNormal():

    x_,y_,z_ = 0, 1, 2

    # ------------------------------------------------- #
    # --- [1] load reference trajectory             --- #
    # ------------------------------------------------- #
    inpFile  = "dat/trajectory_reference.dat"
    import nkUtilities.load__pointFile as lpf
    traj_ref = lpf.load__pointFile( inpFile=inpFile, returnType="point" )


    # ------------------------------------------------- #
    # --- [2] tangent & normal vector               --- #
    # ------------------------------------------------- #
    #  -- [2-1] tangent vector  :: s                --  #
    ds           = np.roll( traj_ref, -1, axis=0 ) - traj_ref
    ds[-1,:]     = ds[-2,:]
    slen         = np.cumsum( np.sqrt( np.sum( ds**2, axis=1 ) ) )
    svec         = ds / np.repeat( ( np.sqrt( np.sum( ds**2, axis=1 ) ) )[:,None], 3, axis=1 )

    #  -- [2-2] normal vector   :: n                --  #
    nvec         = np.zeros( (svec.shape[0],3) )
    nvec[:,0]    = + np.copy( svec[:,1] )
    nvec[:,1]    = - np.copy( svec[:,0] )
    nvec[:,2]    =   0.0
    
    # ------------------------------------------------- #
    # --- [2] plot trajectory                       --- #
    # ------------------------------------------------- #
    import nkUtilities.plot1D       as pl1
    import nkUtilities.load__config as lcf
    pngFile = "png/tangent_normal.png"
    config  = lcf.load__config()
    fig     = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=traj_ref[:,x_], yAxis=traj_ref[:,y_], \
                   label="trajectory-Ref", linestyle="--" )
    fig.add__arrow( xAxis=traj_ref[:,x_], yAxis=traj_ref[:,y_], \
                    uvec =svec    [:,x_], vvec =svec    [:,y_],
                    color="blue"   , width=0.005, scale=0.02, nvec=20 )
    fig.add__arrow( xAxis=traj_ref[:,x_], yAxis=traj_ref[:,y_], \
                    uvec =nvec    [:,x_], vvec =nvec    [:,y_],
                    color="magenta", width=0.005, scale=0.02, nvec=20 )
    fig.set__axis()
    fig.save__figure()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    show__tangentNormal()
