import numpy as np

# ========================================================= #
# ===  convert__alongTrajectory.py                      === #
# ========================================================= #

def convert__alongTrajectory():
    
    # ------------------------------------------------- #
    # --- [1] load bfield                           --- #
    # ------------------------------------------------- #

    inpFile = "dat/bfield_input.dat"
    import nkUtilities.load__pointFile as lpf
    bfield  = lpf.load__pointFile( inpFile=inpFile, returnType="structured" )

    inpFile = "dat/trajectory_input.dat"
    import nkUtilities.load__pointFile as lpf
    trajec  = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    inpFile = "dat/trajectory_ref.dat"
    import nkUtilities.load__pointFile as lpf
    traj_ref= lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # ------------------------------------------------- #
    # --- [2] coordinate to be interpolated         --- #
    # ------------------------------------------------- #
    
    ds           = np.roll( trajec, -1, axis=0 ) - trajec
    ds[-1,:]     = ds[-2,:]
    slen         = np.cumsum( np.sqrt( np.sum( ds**2, axis=1 ) ) )
    svec         = ds / np.repeat( ( np.sqrt( np.sum( ds**2, axis=1 ) ) )[:,None], 3, axis=1 )

    nvec         = np.zeros( (svec.shape[0],3) )
    nvec[:,0]    = + np.copy( svec[:,1] )
    nvec[:,1]    = - np.copy( svec[:,0] )
    nvec[:,2]    =   0.0

    nc_MinMaxNum = [ -0.2, +0.2, 21 ]
    ncoord       = np.linspace( nc_MinMaxNum[0], nc_MinMaxNum[1], nc_MinMaxNum[2] )

    Ls           = trajec.shape[0]
    Ln           = ncoord.shape[0]
    traj_inP_xy  = np.zeros( (Ls,Ln,3) )
    for inc,Anc in enumerate(ncoord):
        traj_inP_xy[:,inc,:] = trajec + Anc*nvec
    sco, nco     = np.meshgrid( slen, ncoord, indexing="ij" )
    zco          = np.zeros_like(sco)
    traj_inP_sn  = np.concatenate( [sco[:,:,None],nco[:,:,None],zco[:,:,None]], axis=2 )
    
    # ------------------------------------------------- #
    # --- [3] interpolation                         --- #
    # ------------------------------------------------- #
    import nkInterpolator.interpolate__tricubic as itc
    grid    = bfield[:,:,:,0:3]
    points  = np.reshape( traj_inP_xy , (-1,3) )
    bf      = np.zeros( (points.shape[0],3) )
    for ik in range(3):
        bt        = bfield[:,:,:,ik+3]
        gridData  = np.concatenate( [grid,bt[:,:,:,None]],axis=3 )
        pointData = np.concatenate( [points,np.zeros( (points.shape[0],1))], axis=1 )
        ret       = itc.interpolate__tricubic( gridData=gridData, pointData=pointData )
        bf[:,ik]  = np.copy( ret[:,3] )
    b_along  = np.concatenate( [points,bf], axis=1 )
    b_along  = np.reshape( b_along, (Ls,Ln,6) )
    b_along  = np.concatenate( [traj_inP_sn,b_along],axis=2 )


    # ------------------------------------------------- #
    # --- [4] desired coordinate in xy              --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.0, 1.0, 21 ]
    x2MinMaxNum = [ -1.0, 1.0, 21 ]
    x3MinMaxNum = [ -1.0, 1.0, 21 ]
    ret         = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )

    
    
    
    return()

# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    convert__alongTrajectory()
